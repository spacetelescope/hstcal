#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "hstio.h"

#include "acs.h"
#include "acsinfo.h"
#include "err.h"

/* number of bias cols */
#define NBIAS_COLS 24

/* number of overscan rows */
#define NOSCN_ROWS 20

/* prototypes for functions in this file */
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
                          char * good_rows[NAMPS], double * ampdata[NAMPS],
                          int * num_fixed, int * num_skipped);
static int make_amp_array(const int arr_rows, const int arr_cols, SingleGroup * im,
                          int amp, double * array);
static int unmake_amp_array(const int arr_rows, const int arr_cols, SingleGroup * im,
                            int amp, double * array);


/* remove signal dependent bias shift from post-SM4 full frame WFC data.
 * based on ISR http://www.stsci.edu/hst/acs/documents/isrs/isr1202.pdf
 * chip2 is amps C & D, chip1 is amps A & B. */
int bias_shift_corr(ACSInfo *acs, SingleGroup *chip2, SingleGroup *chip1) {
  extern int status;
  
  int i, j, k;     /* iteration variables */
  
  /* array for amp data and amp data + gap pixels */
  double * ampdata, * ampdata_gap;
  
  const double serial_freq = 1000./22.;    /* serial pixel frequency */
  const double parallel_shift = 3.212;     /* parallel shift time */
  
  /* number of virtual pixels at end of each row */
  const int ngap_pix = (int) (serial_freq * parallel_shift + 0.5);
  
  const int arr_rows = chip2->sci.data.ny;
  const int arr_cols = chip2->sci.data.nx/2;
  
  /* total number of real and gap pixels per quadrant */
  const int nquad_pix = (arr_cols + ngap_pix) * arr_rows;
  
  /* array of true DC bias levels */
  double * dc_bias_levels;
  
  /* arrays below are in order of amp A, B, C, D */
  
  /* time constant of AC high pass filter in external pre-amp. */
  const double time_const[NAMPS] = {37.290251, 36.180001, 37.867770, 42.461249};

  /* ratio of DC offset shift and pixel signal */
  const double dc_ratio[NAMPS] = {0.3, 0.3, 0.3, 0.3};
  
  /* DSI sensitivity */
  const double dsi_sens[NAMPS] = {2.9188919e-3, 10.805754e-3, 12.432145e-3, 3.8596253e-3};
    
  /* factor combining time constant and clocking frequency */
  double factor;
  
  /* summation variables */
  double sum;
  int num;
  double magic_square_mean;
  
  /* allocate space for data arrays */
  ampdata = malloc(arr_rows * arr_cols * sizeof(double));
  ampdata_gap = malloc(nquad_pix * sizeof(double));
  
  /* allocate space for true DC bias levels array */
  dc_bias_levels = malloc((nquad_pix + 1) * sizeof(double));
  
  for (i = 0; i < NAMPS; i++) {
    /* put a single amp's data into ampdata */
    if (i < 2) {
      make_amp_array(arr_rows, arr_cols, chip1, i, ampdata);
    } else {
      make_amp_array(arr_rows, arr_cols, chip2, i, ampdata);
    }
       
    /* calculate "magic square mean" */
    sum = 0.0;
    num = 0;
    
    for (j = 2057; j <= 2066; j++) {
      for (k = 13; k <= 22; k++) {
        sum += ampdata[arr_cols*j + k];
        num++;
      }
    }
    
    magic_square_mean = sum / (double) num;
    
    /* make amp + gap array */
    for (j = 0; j < arr_rows; j++) {
      for (k = 0; k < (arr_cols + ngap_pix); k++) {
        if (k < arr_cols) {
          ampdata_gap[(arr_cols + ngap_pix)*j + k] = ampdata[arr_cols*j + k];
        } else {
          ampdata_gap[(arr_cols + ngap_pix)*j + k] = magic_square_mean;
        }
      }
    }

    /* calculate true DC bias levels */
    factor = 1.0 - exp(-1.0 / (time_const[i] * serial_freq));
    
    dc_bias_levels[0] = magic_square_mean * dc_ratio[i];
    
    for (j = 1; j < nquad_pix + 1; j++) {
      dc_bias_levels[j] = ampdata_gap[j-1] * factor * dc_ratio[i] + 
                            (1.0 - factor) * dc_bias_levels[j-1];
    }
    
    /* calculate correction to data */
    for (j = 0; j < nquad_pix; j++) {
      ampdata_gap[j] = (ampdata_gap[j] - dsi_sens[i] * dc_bias_levels[j+1]) -
                          (10./22.) * (dc_bias_levels[j+1] - dc_bias_levels[j]);
    }
    
    /* copy corrected data back to ampdata */
    for (j = 0; j < arr_rows; j++) {
      for (k = 0; k < arr_cols; k++) {
        ampdata[arr_cols*j + k] = ampdata_gap[(arr_cols + ngap_pix)*j + k];
      }
    }
    
    /* re-calculate "magic square mean" */
    sum = 0.0;
    num = 0;
    
    for (j = 2057; j <= 2066; j++) {
      for (k = 13; k <= 22; k++) {
        sum += ampdata[arr_cols*j + k];
        num++;
      }
    }
    
    magic_square_mean = sum / (double) num;
    /* Keep track of the total bias correction applied to each AMP region */
    acs->blev[i] += magic_square_mean;
    /* Report to the user the contribution to the bias level correction 
       made by this processing step.
    */
    sprintf(MsgText, "Bias shift correcting for bias level in Amp %c of %0.4f electrons.",AMPSORDER[i], magic_square_mean);
    trlmessage(MsgText);
    
    /* subtract "magic square mean" from data*/
    for (j = 0; j < arr_rows; j++) {
      for (k = 0; k < arr_cols; k++) {
        ampdata[arr_cols*j + k] -= magic_square_mean;
      }
    }
    
    /* copy modified data back to SingleGroup structs */
    if (i < 2) {
      unmake_amp_array(arr_rows, arr_cols, chip1, i, ampdata);
    } else {
      unmake_amp_array(arr_rows, arr_cols, chip2, i, ampdata);
    }
  }
  
  free(ampdata);
  free(ampdata_gap);
  free(dc_bias_levels);
  
  return status;
}


/* remove amplifier cross-talk */
void cross_talk_corr(ACSInfo *acs, SingleGroup *im) {
  /* iteration variables */
  int i, j;
  
  /* cross talk scaling constant */
  double cross_scale = 9.1e-5;
  
  double temp;
  
  const int arr_rows = im->sci.data.ny;
  const int arr_cols = im->sci.data.nx;

  for (i = 0; i < arr_rows; i++) {
    for (j = 0; j < arr_cols; j++) {
      temp = Pix(im->sci.data, arr_cols-j-1, i) * cross_scale;
      
      Pix(im->sci.data, j, i) += (float) temp;
    }
  }
}


/* remove stripes from post-SM4 full frame WFC data using information in
 * the prescan regions. chip2 is amps C & D, chip1 is amps A & B. */
int doDestripe(ACSInfo *acs, SingleGroup *chip2, SingleGroup *chip1) {
  extern int status;

  /* iteration variables */
  int i, j, k;

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
  
  /* number of rows ultimately worked on and fixed */
  int rows_fixed;
  int rows_skipped;
  
  /* character array for holding history messages */
  char history[ACS_LINE];

  /* bias pixel mean, standard deviation, and number of good pixels */
  double bias_mean, bias_std;
  int good_bias_pix;

  /* holder of bias means from each amp, saved here so I can put the
   * MEANBLEV keyword in the science extension headers */
  double bias_mean_arr[NAMPS];
  
  int PutKeyFlt(Hdr *, char *, float, char *);
  int blevHistory(ACSInfo *, Hdr *, int, int);
  int MkName (char *, char *, char *, char *, char *, int);

  /* figure out the size of individual amp arrays
   * should be 2068 rows by 2072 columns */
  arr_rows = chip2->sci.data.ny;
  arr_cols = chip2->sci.data.nx/2;

  /* allocate space for the amp arrays */
  for (i = 0; i < NAMPS; i++) {
    ampdata[i] = malloc(arr_rows * arr_cols * sizeof(double));
    good_rows[i] = malloc(arr_rows * sizeof(char));
  }

  /* copy data from SingleGroup structs to amp arrays */
  for (i = 0; i < NAMPS; i++) {
    if (i < 2) {
      make_amp_array(arr_rows, arr_cols, chip1, i, ampdata[i]);
    } else {
      make_amp_array(arr_rows, arr_cols, chip2, i, ampdata[i]);
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
    sprintf(MsgText, "     bias level of %.6g electrons was subtracted for AMP %c.",
            bias_mean, AMPSORDER[i]);
    trlmessage(MsgText);

    acs->blev[i] += bias_mean;
    bias_mean_arr[i] = bias_mean;
  }

  /* add MEANBLEV keyword to science extension headers */
  if (PutKeyFlt (&chip1->sci.hdr, "MEANBLEV", (bias_mean_arr[0] + bias_mean_arr[1])/2.,
                 "mean of bias levels subtracted in electrons")) {
    return (status);
  }
  if (PutKeyFlt (&chip2->sci.hdr, "MEANBLEV", (bias_mean_arr[2] + bias_mean_arr[3])/2.,
                 "mean of bias levels subtracted in electrons")) {
    return (status);
  }

  /* remove stripes */
  if (remove_stripes(arr_rows, arr_cols, good_rows, ampdata, &rows_fixed, &rows_skipped)) {
    return status;
  }
  
  /* add history keywords about rows fixed and rows skipped */
  sprintf(history, "DESTRIPE: number of rows fixed per amp: %i", rows_fixed);
  addHistoryKw(chip2->globalhdr, history);
  if (hstio_err()) {
    return (status = HEADER_PROBLEM);
  }
  
  sprintf(history, "DESTRIPE: number of rows skipped per amp: %i", rows_skipped);
  addHistoryKw(chip2->globalhdr, history);
  if (hstio_err()) {
    return (status = HEADER_PROBLEM);
  }

  /* copy modified data back to SingleGroup structs */
  for (i = 0; i < NAMPS; i++) {
    if (i < 2) {
      unmake_amp_array(arr_rows, arr_cols, chip1, i, ampdata[i]);
    } else {
      unmake_amp_array(arr_rows, arr_cols, chip2, i, ampdata[i]);
    }
  }

  /* free allocated arrays */
  for (i = 0; i < NAMPS; i++) {
    free(ampdata[i]);
    free(good_rows[i]);
  }

  return status;
}

/* calculate stripe corrections and remove stripes */
static int remove_stripes(const int arr_rows, const int arr_cols,
                          char * good_rows[NAMPS], double * ampdata[NAMPS],
                          int * num_fixed, int * num_skipped) {
  extern int status;

  /* iteration variables */
  int i,j,k;

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
  
  *num_fixed = 0;
  *num_skipped = 0;

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
      if (ampstds[k] >= 15) {
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

      (*num_fixed)++;
    } else {
      (*num_skipped)++;
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
 * of standard deviation, but the returned std is standard deviation. */
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

    if (ntot1 != 0) {
      temp_mean = sum_mean / (double) ntot1;
    } else {
      temp_mean = 0;
    }

    /* calculate absolute deviation */
    sum_std = 0.0;
    for (j = 0; j < NBIAS_COLS; j++) {
      for (i = 0; i < arr_rows; i++) {
        if (good_pix[NBIAS_COLS*i +j] == 1) {
          sum_std += fabs(array[arr_cols*i + j] - temp_mean);
        }
      }
    }

    if (ntot1 != 0) {
      temp_std = sum_std/(double) ntot1;
    } else {
      temp_std = 0;
    }

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

  if (ntot1 != 0) {
    temp_mean = sum_mean / (double) ntot1;
  } else {
    temp_mean = 0;
  }

  sum_std = 0.0;
  for (j = 0; j < NBIAS_COLS; j++) {
    for (i = 0; i < arr_rows; i++) {
      if (good_pix[NBIAS_COLS*i +j] == 1) {
        sum_std += pow(array[arr_cols*i + j] - temp_mean,2);
      }
    }
  }

  if (ntot1 != 0) {
    temp_std = sqrt(sum_std/(double) (ntot1 - 1));
  } else {
    temp_std = 0;
  }

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
  
  /* mean of the bias region */
  double bias_mean;
  
  /* overall mean of science area of amp */
  double amp_mean;

  /* means of data rows - row_ref for that row */
  double row_means[arr_rows];

  /* some kind of weird mean of the bias rows */
  double row_ref[arr_rows];

  double sum;
  
  /* calculate the bias mean */
  sum = 0.0;
  
  for (i = 0; i < arr_rows-NOSCN_ROWS; i++) {
    for (j = 0; j < NBIAS_COLS; j++) {
      sum += array[arr_cols*i + j];
    }
  }
  
  bias_mean = sum / (double) ((arr_rows - NOSCN_ROWS) * NBIAS_COLS);
  
  /* calculate the science area mean */
  sum = 0.0;
  
  for (i = 0; i < arr_rows-NOSCN_ROWS; i++) {
    for (j = NBIAS_COLS; j < arr_cols; j++) {
      sum += (array[arr_cols*i + j] - bias_mean);
    }
  }
  
  amp_mean = sum / (double) ((arr_cols - NBIAS_COLS) * (arr_rows - NOSCN_ROWS));

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

    if (abs(row_means[i] - amp_mean) > 100) {
      good_rows[i] = 0;
    } else if (array[arr_cols*i + 24] > 70000) {
      good_rows[i] = 0;
    } else if (array[arr_cols*i + 25] > 70000) {
      good_rows[i] = 0;
    } else if (array[arr_cols*i + 26] > 70000) {
      good_rows[i] = 0;
    } else if (array[arr_cols*i + 27] > 70000) {
      good_rows[i] = 0;
    } else if (array[arr_cols*i + 28] > 70000) {
      good_rows[i] = 0;
    } else if (array[arr_cols*i + 29] > 70000) {
      good_rows[i] = 0;
    } else {
      good_rows[i] = 1;
      (*num_good_rows)++;
    }
  }
  
  //printf("%i   %i   %8.2f   %8.2f\n",*num_good_rows,2068-*num_good_rows,bias_mean,amp_mean);
  
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
      if ((array[arr_cols*i + 24] < 50000) &&
          (array[arr_cols*i + 25] < 50000) &&
          (array[arr_cols*i + 26] < 50000) &&
          (array[arr_cols*i + 27] < 50000) &&
          (array[arr_cols*i + 28] < 50000) &&
          (array[arr_cols*i + 29] < 50000)) {
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
  int nused;          /* number of pixels used for average (len - skipped) */

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
