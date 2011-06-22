#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "hstio.h"

#include "acs.h"
#include "acsinfo.h"
#include "acserr.h"

/* number of bias cols */
#define NBIAS_COLS 24

/* prototypes for functions in this file */
static int destripe(ACSInfo * acs);
static int calc_mean_std(const int len, const double array[], const double sig,
                         double *mean, double *std);
static int bias_col_mean_std(const int arr_rows, const int arr_cols, const double array[],
                             double col_means[], double col_stds[]);
static int sub_bias_col_means(const int arr_rows, const int arr_cols, const int bias_cols,
                              const double bias_means[], double array[]);
static int find_good_rows(const int arr_rows, const int arr_cols, const double * array,
                          char * good_rows, int * num_good_rows);
static int calc_bias_mean_std(const int arr_rows, const int arr_cols,
                              const double * array, const char * good_rows,
                              double * mean, double * std, int * pix_used);
static int remove_stripes(const int arr_rows, const int arr_cols, 
                          char * good_rows[NAMPS], double * ampdata[NAMPS]);
static int make_amp_array(const int arr_rows, const int arr_cols, SingleGroup * im,
                          int amp, double * array);
static int unmake_amp_array(const int arr_rows, const int arr_cols, SingleGroup * im,
                            int amp, double * array);

int doDestripe(ACSInfo * acs) {
  extern int status;
  
  Hdr phdr;		/* primary header for input image */
  
  /* functions from lib */
  int LoadHdr (char *, Hdr *);
  void PrSwitch (char *, int);
  void TimeStamp (char *, char *);
  
  PrSwitch("blevcorr", PERFORM);
  trlmessage("Performing stripe removal and bias level subtraction.");
  
  if (destripe(acs)) {
    return status;
  }
  
  PrSwitch("blevcorr", COMPLETE);
  
  if (acs->printtime) {
    TimeStamp("BLEVCORR complete", acs->rootname);
  }
  
  return status;
}

static int destripe(ACSInfo * acs) {
  extern int status;
  
  /* iteration variables */
  int i,j,k;
  
  /* structures to hold data for both chips */
  SingleGroup chip2, chip1;
  
  /* amp array size variables */
  int arr_rows, arr_cols;
  
  /* array of arrays for each amp's data in order of AMPSORDER */
  double * ampdata[NAMPS];
  
  /* arrays of bias column means and standard deviation */
  double bias_col_means[NBIAS_COLS];
  double bias_col_stds[NBIAS_COLS];
  
  /* array of arrays designating whether a row is usable or not. 0: bad, 1:good */
  char * good_rows[NAMPS];
  
  /* array of number of good rows for each amp */
  int num_good_rows[NAMPS];
  
  /* bias pixel mean, standard deviation, and number of good pixels */
  double bias_mean, bias_std;
  int good_bias_pix;
  
  /* holder of bias means from each amp, saved here so I can put the
   * MEANBLEV keyword in the science extension headers */
  double bias_mean_arr[NAMPS];
  
  /* output filename. will become the new acs->input for the next steps */
  char outname[ACS_LINE+1];
  
  int PutKeyFlt(Hdr *, char *, float, char *);
  int blevHistory(ACSInfo *, Hdr *, int, int);
  
  /* get data and put it into individual amp arrays, with the amp in the
   * lower left corner */
  initSingleGroup(&chip2);
  initSingleGroup(&chip1);
  
  getSingleGroup(acs->input, 1, &chip2);
  if (hstio_err()) {
    return (status = OPEN_FAILED);
  }
  
  getSingleGroup(acs->input, 2, &chip1);
  if (hstio_err()) {
    freeSingleGroup(&chip2);
    return (status = OPEN_FAILED);
  }
  
  /* figure out the size of individual amp arrays 
   * should be 2068 rows by 2072 columns */
  arr_rows = chip2.sci.data.ny;
  arr_cols = chip2.sci.data.nx/2;
  
  /* allocate space for the amp arrays */
  for (i = 0; i < NAMPS; i++) {
    ampdata[i] = malloc(arr_rows * arr_cols * sizeof(double));
    good_rows[i] = malloc(arr_rows * sizeof(char));
  }
  
  /* copy data from SingleGroup structs to amp arrays */
  for (i = 0; i < NAMPS; i++) {
    if (i < 2) {
      make_amp_array(arr_rows, arr_cols, &chip1, i, ampdata[i]);
    } else {
      make_amp_array(arr_rows, arr_cols, &chip2, i, ampdata[i]);
    }
  }
  
  /* subtract each column's mean as computed after removing the mean of each
   * bias row, ignoring bias rows near saturated pixels, and doing 
   * sigma rejection of outlying bias pixels. */
  for (i = 0; i < NAMPS; i++) {
    if (bias_col_mean_std(arr_rows, arr_cols, ampdata[i], 
                          bias_col_means, bias_col_stds)) {
      return status;
    }
    if (sub_bias_col_means(arr_rows, arr_cols, NBIAS_COLS, bias_col_means,
                           ampdata[i])) {
      return status;
    }
  }
  
  /* for each amp figure out which rows to use and how many good rows there are */
  for (i = 0; i < NAMPS; i++) {
    find_good_rows(arr_rows, arr_cols, ampdata[i], good_rows[i], &num_good_rows[i]);
  }
  
  /* for each amp, figure out the mean of the good bias pixels with "sigma"
   * clipping, then subtract that mean from all that amp's data. */
  for (i = 0; i < NAMPS; i++) {
    if (calc_bias_mean_std(arr_rows, arr_cols, ampdata[i], good_rows[i],
                           &bias_mean, &bias_std, &good_bias_pix)) {
      return status;
    }
    
    /* subtract the mean from all the pixels in the image. */
    for (j = 0; j < arr_rows; j++) {
      for (k = 0; k < arr_cols; k++) {
        ampdata[i][arr_cols*j + k] -= bias_mean;
      }
    }
    
    /* report bias level subtracted to user */
    sprintf(MsgText, "     bias level of %.6g was subtracted for AMP %c.", bias_mean, AMPSORDER[i]);
    trlmessage(MsgText); 
    
    bias_mean_arr[i] = bias_mean;
  }
  
  /* add MEANBLEV keyword to science extension headers */
  if (PutKeyFlt (&chip1.sci.hdr, "MEANBLEV", (bias_mean_arr[0] + bias_mean_arr[1])/2.,
                 "mean of bias levels subtracted")) {
    return (status);
  }
  if (PutKeyFlt (&chip2.sci.hdr, "MEANBLEV", (bias_mean_arr[2] + bias_mean_arr[3])/2.,
                 "mean of bias levels subtracted")) {
    return (status);
  }
  
  /* remove stripes */
  if (remove_stripes(arr_rows,arr_cols,good_rows,ampdata)) {
    return status;
  }
  
  /* copy modified data back to SingleGroup structs */
  for (i = 0; i < NAMPS; i++) {
    if (i < 2) {
      unmake_amp_array(arr_rows, arr_cols, &chip1, i, ampdata[i]);
    } else {
      unmake_amp_array(arr_rows, arr_cols, &chip2, i, ampdata[i]);
    }
  }
  
  /* save new file */
  /* make output file name */
  if (MkName(acs->input, "_raw", "_strp_tmp", "", outname, ACS_LINE)) {
    return status;
  }
  
  /* copy outname to acs->input */
  strcpy(acs->input, outname);
  
  /* update header history and keywords */
  if (blevHistory(acs, chip2.globalhdr, YES, NO)) {
    return status;
  }
  
  /* output destriped data to temp file */
  putSingleGroup(outname, 1, &chip2, 0);
  putSingleGroup(outname, 2, &chip1, 0);
  
  /* free allocated arrays */
  for (i = 0; i < NAMPS; i++) {
    free(ampdata[i]);
    free(good_rows[i]);
  }
  
  freeSingleGroup(&chip2);
  freeSingleGroup(&chip1);
  
  return status;
}

/* calculate stripe corrections and remove stripes */
static int remove_stripes(const int arr_rows, const int arr_cols, 
                          char * good_rows[NAMPS], double * ampdata[NAMPS]) {
  extern int status;
  
  /* iteration variables */
  int i,j,k;
  
  /* number of rows skipped and fixed */
  int num_skipped = 0;
  int num_fixed = 0;
  
  /* summation holders */
  double sum_mean;
  double sum_std;
  int sum_num;
  
  /* arrays of means/std devs, one element for each amp */
  double ampmeans[NAMPS];
  double ampstds[NAMPS];
  
  /* arrays of row means/stds */
  double * rowmeans;
  
  /* interpolation points */
  double interp_pt;
  int interp_lower;
  double interp_dist;
  double interp_mean;
  
  /* allocate arrays as necessary */
  rowmeans = malloc(arr_rows * sizeof(double));
  
  for (i = 0; i < arr_rows; i++) {
    /* calculate mean of this row for each amp */
    for (k = 0; k < NAMPS; k++) {
      sum_mean = 0.0;
      
      for (j = 0; j < NBIAS_COLS; j++) {
        sum_mean += ampdata[k][arr_cols*i + j];
      }
      ampmeans[k] = sum_mean / (double) NBIAS_COLS;
    }
    
    /* calculate the std dev of this row for each amp */
    for (k = 0; k < NAMPS; k++) {
      sum_std = 0.0;
      
      for (j = 0; j < NBIAS_COLS; j++) {
        sum_std += pow(ampdata[k][arr_cols*i + j] - ampmeans[k], 2);
      }
      ampstds[k] = sqrt(sum_std / (double) (NBIAS_COLS - 1));
    }
    
    /* check whether this row meets usability standards */
    for (k = 0; k < NAMPS; k++) {
      if (ampstds[k] >= 7.5) {
        good_rows[k][i] = 0;
      }
    }
    
    /* calculate number of amp rows we like */
    sum_mean = 0.0;
    sum_num = 0;
    for (k = 0; k < NAMPS; k++) {
      if (good_rows[k][i] == 1) {
        sum_mean += good_rows[k][i]*ampmeans[k];
        sum_num++;
      }
    }
    
    /* calculate mean of amp row means */
    rowmeans[i] = 0.0;
    
    if (sum_num >= 2) {
      rowmeans[i] = sum_mean / (double) sum_num;
    
      num_fixed++;
    } else {
      num_skipped++;
    }
  }
  
  /* perform actual destriping correction */
  for (j = 0; j < arr_cols; j++) {
    for (i = 0; i < arr_rows; i++) {
      interp_pt = ((double) i) + (((double) j) - 11.5) / ((double) (arr_cols + 146));
      interp_lower = (int) floor(interp_pt);
      
      /* check for out of bounds interp_lower */
      if (interp_lower < 0) {
        interp_lower = 0;
      } else if (interp_lower > arr_rows - 2) {
        interp_lower = arr_rows - 2;
      }
      
      interp_dist = interp_pt - (double) interp_lower;
      
      interp_mean = rowmeans[interp_lower] + 
                    interp_dist * (rowmeans[interp_lower+1] - rowmeans[interp_lower]);
      
      for (k = 0; k < NAMPS; k++) {
        ampdata[k][arr_cols*i + j] -= interp_mean;
      }
    }
  }
  
  free(rowmeans);
  
  return status;
}

/* calculate the mean and standard deviation of all the bias pixels with
 * sigma clipping. clipping is done according to absolute deviation instead
 * of absolute deviation, but the returned std is standard deviation. */
static int calc_bias_mean_std(const int arr_rows, const int arr_cols,
                              const double * array, const char * good_rows,
                              double * mean, double * std, int * pix_used) {
  extern int status;
  
  /* iteration variables */
  int i,j;
  
  /* sigma clipping level */
  double sig_lev = 6.0;
  
  /* summation holders */
  double sum_mean, sum_std;
  
  /* holders for mean, std results */
  double temp_mean, temp_std;
  
  /* holders of total number of pixels used */
  int ntot1, ntot2;
  
  /* array of pixel good/bad */
  char good_pix[arr_rows * NBIAS_COLS];
  
  /* number of sigma clipping iterations performed */
  int nits = 0;
  
  /* fill good_pix, and add up good pixels */
  ntot1 = 0;
  for (j = 0; j < NBIAS_COLS; j++) {
    for (i = 0; i < arr_rows; i++) {
      good_pix[NBIAS_COLS*i + j] = good_rows[i];
      
      if (good_rows[i] == 1) {
        ntot1++;
      }
    }
  }
  
  /* iterate, performing "sigma" clipping at least 5 times but a max of 10 */
  do {
    ntot2 = ntot1;
    ntot1 = 0;
    nits++;
    
    /* calculate absolute deviation and current number of good pixels*/
    sum_mean = 0.0;
    for (j = 0; j < NBIAS_COLS; j++) {
      for (i = 0; i < arr_rows; i++) {
        if (good_pix[NBIAS_COLS*i +j] == 1) {
          sum_mean += array[arr_cols*i + j];
          ntot1++;
        }
      }
    }
    
    temp_mean = sum_mean / (double) ntot1;

    /* calculate absolute deviation */
    sum_std = 0.0;
    for (j = 0; j < NBIAS_COLS; j++) {
      for (i = 0; i < arr_rows; i++) {
        if (good_pix[NBIAS_COLS*i +j] == 1) {
          sum_std += fabs(array[arr_cols*i + j] - temp_mean);
        }
      }
    }

    temp_std = sum_std/(double) ntot1;
    
    ntot1 = 0;
    
    for (j = 0; j < NBIAS_COLS; j++) {
      for (i = 0; i < arr_rows; i++) {
        if ((good_pix[NBIAS_COLS*i + j] == 1) &&
            (fabs(array[arr_cols*i+j] - temp_mean) > sig_lev*temp_std)) {
          good_pix[NBIAS_COLS*i + j] = 0;
        } else if (good_pix[NBIAS_COLS*i + j] == 1) {
          ntot1++;
        }
      }
    }
  } while ((nits <= 5) || (nits <= 10 && ntot1 < ntot2));
  
  /* calculate final mean and standard deviation of good pixels */
  /* calculate absolute deviation and current number of good pixels*/
  sum_mean = 0.0;
  for (j = 0; j < NBIAS_COLS; j++) {
    for (i = 0; i < arr_rows; i++) {
      if (good_pix[NBIAS_COLS*i +j] == 1) {
        sum_mean += array[arr_cols*i + j];
      }
    }
  }
  
  temp_mean = sum_mean / (double) ntot1;
  
  sum_std = 0.0;
  for (j = 0; j < NBIAS_COLS; j++) {
    for (i = 0; i < arr_rows; i++) {
      if (good_pix[NBIAS_COLS*i +j] == 1) {
        sum_std += pow(array[arr_cols*i + j] - temp_mean,2);
      }
    }
  }
  
  temp_std = sqrt(sum_std / (double) ntot1);
  
  *mean = temp_mean;
  *std = temp_std;
  *pix_used = ntot1;
  
  return status;
}

/* for the given array of amp data, figure out which rows are good and how
 * many good rows there are. 
 */
static int find_good_rows(const int arr_rows, const int arr_cols, const double * array,
                          char * good_rows, int * num_good_rows) {
  extern int status;
  
  /* iteration variables */
  int i,j;
  
  /* means of data rows - row_ref for that row */
  double row_means[arr_rows];
  
  /* some kind of weird mean of the bias rows */
  double row_ref[arr_rows];
  
  double sum;
  
  /* never use the first row */
  good_rows[0] = 0;
  
  *num_good_rows = 0;
  
  /* calculate "row_ref" for the first row */
  sum = 0.0;
  for (j = 0; j < NBIAS_COLS; j++) {
    sum += array[j];
  }
  
  row_ref[0]= sum / (double) NBIAS_COLS;
  
  row_means[0] = 0.0;
  
  for (i = 1; i < arr_rows; i++) {
    sum = 0.0;
    
    for (j = 0; j < NBIAS_COLS; j++) {
      sum += array[arr_cols*i + j] + array[arr_cols*(i-1) + j];
    }
    
    row_ref[i] = sum / (double) (2*NBIAS_COLS);
    
    sum = 0.0;
    
    for (j = NBIAS_COLS; j < arr_cols; j++) {
      sum += array[arr_cols*i + j] - row_ref[i];
    }
    
    row_means[i] = sum / (double) (arr_cols - NBIAS_COLS);
    
    if (row_means[i] > 50) {
      good_rows[i] = 0;
    } else if (array[arr_cols*i + 25] > 25000) {
      good_rows[i] = 0;
    } else if (array[arr_cols*i + 26] > 25000) {
      good_rows[i] = 0;
    } else if (array[arr_cols*i + 27] > 25000) {
      good_rows[i] = 0;
    } else if (array[arr_cols*i + 28] > 25000) {
      good_rows[i] = 0;
    } else if (array[arr_cols*i + 29] > 25000) {
      good_rows[i] = 0;
    } else if (array[arr_cols*i + 30] > 25000) {
      good_rows[i] = 0;
    } else {
      good_rows[i] = 1;
      (*num_good_rows)++;
    }
  }
  
  return status;
}

/* calculate the mean and standard deviation of each bias column. these are
 * calculated after removing the mean of each bias row, removing rows adjacent
 * to saturated pixels, and doing 4 iterations of 5-sigma clipping.
 */
static int bias_col_mean_std(const int arr_rows, const int arr_cols, const double array[],
                             double col_means[], double col_stds[]) {
  extern int status;
  
  /* iteration variables */
  int i,j;
  
  /* sigma level at which to do clipping */
  double sig_lev = 5.0;
  
  double sum;
  double mean;
  
  /* array of bias columns with row by row mean removed */
  double minus_row_mean[arr_rows*NBIAS_COLS];
  
  /* number of rows not skipped due to saturated pixels */
  int not_skipped;
  
  /* a single bias column */
  double bias_col[arr_rows];
  
  /* values returned by calc_mean_std */
  double temp_mean, temp_std;
  
  /* subtract the average from each row (of the bias columns) where the average
   * is calculated using columns 13 - 24 */
  for (i = 0; i < arr_rows; i++) {
    /* loop over columns 13-24 for this row and calculate mean */
    sum = 0.0;
    for (j = 12; j < NBIAS_COLS; j++) {
      sum += array[arr_cols*i + j];
    }
    mean = sum / 12.0;
    
    /* subtract 13-24 mean from all columns in this row */
    for (j = 0; j < NBIAS_COLS; j++) {
      minus_row_mean[NBIAS_COLS*i + j] = array[arr_cols*i + j] - mean;
    }
  }
  
  /* loop over columns finding the 5-sigma clipped mean and standard deviation
   * of each while ignoring rows affected by nearby saturated pixels */
  for (j = 0; j < NBIAS_COLS; j++) {
    not_skipped = 0;
    
    for (i = 0; i < arr_rows; i++) {
      if ((array[arr_cols*i + 24] < 25000) &&
          (array[arr_cols*i + 25] < 25000) &&
          (array[arr_cols*i + 26] < 25000) &&
          (array[arr_cols*i + 27] < 25000) &&
          (array[arr_cols*i + 28] < 25000) &&
          (array[arr_cols*i + 29] < 25000)) {
        bias_col[not_skipped] = minus_row_mean[NBIAS_COLS*i + j];
        not_skipped ++;
      }
    }
    
    calc_mean_std(not_skipped, bias_col, sig_lev, &temp_mean, &temp_std);
    
    col_means[j] = temp_mean;
    col_stds[j] = temp_std;
  }
  
  return status;
}

/* calculate the mean and standard deviation of a 1D array with sigma clipping.
 * len specifies the length of the array and sig the level of clipping. 
 */
static int calc_mean_std(const int len, const double array[], const double sig,
                         double *mean, double *std) {
  extern int status;

  /* iteration variables */
  int i,j;
  int clip_it = 4;
  
  double sum_mean;
  double sum_std;
  double temp_mean;
  double temp_std;
  int nused;          /* number of pixels used for average (len - skipped)
  
  /* calculate the initial mean */
  sum_mean = 0.0;
  for (i=0; i < len; i++) {
    sum_mean += array[i];
  }
  
  temp_mean = sum_mean / (double) len;
  
  /* calculate initial standard deviation */
  sum_std = 0.0;
  for (i=0; i < len; i++) {
    sum_std += pow(array[i] - temp_mean,2);
  }
  
  temp_std = sqrt(sum_std / (double) (len - 1));
  
  /* do sigma clipping iterations */
  for (j = 0; j < clip_it; j++) {
    nused = 0;
    sum_mean = 0.0;
    sum_std = 0.0;
    
    for (i = 0; i < len; i++) {
      if (fabs(array[i] - temp_mean) <= sig*temp_std) {
        sum_mean += array[i];
        sum_std += pow(array[i] - temp_mean, 2);
        nused++;
      }
    }
    
    temp_mean = sum_mean / (double) nused;
    temp_std = sqrt(sum_std / (double) (nused - 1));
  }
  
  *mean = temp_mean;
  *std = temp_std;
  
  return status;
}

/* subtract the reference column means from the reference columns */
static int sub_bias_col_means(const int arr_rows, const int arr_cols, const int bias_cols,
                              const double bias_means[], double array[]) {
  extern int status;
  
  /* iteration variables */
  int i,j;
  
  for (j = 0; j < bias_cols; j++) {
    for (i = 0; i < arr_rows; i++) {
      array[arr_cols*i + j] = array[arr_cols*i + j] - bias_means[j];
    }
  }
  
  for (j = bias_cols; j < arr_cols; j++) {
    for (i = 0; i < arr_rows; i++) {
      array[arr_cols*i + j] = array[arr_cols*i + j];
    }
  }
  
  return status;
}

/* 
 * make_amp_array returns an array view of the data readout through the
 * specified amp in which the amp is at the lower left hand corner.
 * based on make_amp_array from dopcte.c.
 */
static int make_amp_array(const int arr_rows, const int arr_cols, SingleGroup *im,
                          int amp, double * array) {
  
  extern int status;
  
  /* iteration variables */
  int i, j;
  
  /* variables for the image row/column we want */
  int r,c;
  
  for (i = 0; i < arr_rows; i++) {
    for (j = 0; j < arr_cols; j++) {
      if (amp == AMP_A) {
        r = im->sci.data.ny - i - 1;
        c = j;
      } else if (amp == AMP_B) {
        r = im->sci.data.ny - i - 1;
        c = im->sci.data.nx - j - 1;
      } else if (amp == AMP_C) {
        r = i;
        c = j;
      } else if (amp == AMP_D) {
        r = i;
        c = im->sci.data.nx - j -1;
      } else {
        trlerror("Amp number not recognized, must be 0-3.");
        status = ERROR_RETURN;
        return status;
      }
      
      array[i*arr_cols + j] = Pix(im->sci.data, c, r);
    }
  }
  
  return status;
}

/*
 * unmake_amp_array does the opposite of make_amp_array, it takes amp array
 * views and puts them back into the single group in the right order.
 */
static int unmake_amp_array(const int arr_rows, const int arr_cols, SingleGroup *im,
                            int amp, double * array) {
  
  extern int status;
  
  /* iteration variables */
  int i, j;
  
  /* variables for the image row/column we want */
  int r,c;
  
  for (i = 0; i < arr_rows; i++) {
    for (j = 0; j < arr_cols; j++) {
      if (amp == AMP_A) {
        r = im->sci.data.ny - i - 1;
        c = j;
      } else if (amp == AMP_B) {
        r = im->sci.data.ny - i - 1;
        c = im->sci.data.nx - j - 1;
      } else if (amp == AMP_C) {
        r = i;
        c = j;
      } else if (amp == AMP_D) {
        r = i;
        c = im->sci.data.nx - j -1;
      } else {
        trlerror("Amp number not recognized, must be 0-3.");
        status = ERROR_RETURN;
        return status;
      }
      
      Pix(im->sci.data, c, r) = (float) array[i*arr_cols + j];
    }
  }
  
  return status;
}
