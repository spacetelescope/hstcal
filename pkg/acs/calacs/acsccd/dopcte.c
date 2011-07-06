#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "hstio.h"

#include "acs.h"
#include "acsinfo.h"
#include "acserr.h"

#include "pcte.h"

static int get_amp_array_size(const ACSInfo *acs, SingleGroup *x,
                              const int amp, char *amploc, char *ccdamp,
                              int *xsize, int *ysize, int *xbeg, 
                              int *xend, int *ybeg, int *yend);
static int make_amp_array(const ACSInfo *acs, const SingleGroup *im, const int amp, 
                          const int arr1, const int arr2,
                          const int xbeg, const int ybeg,
                          double amp_sci_array[arr1*arr2], 
                          double amp_err_array[arr1*arr2]);
static int unmake_amp_array(const ACSInfo *acs, const SingleGroup *im, const int amp, 
                            const int arr1, const int arr2,
                            const int xbeg, const int ybeg,
                            double amp_sci_array[arr1*arr2], 
                            double amp_err_array[arr1*arr2]);

/*
 * Perform a pixel based CTE correction on the SCI data extension of ACS CCD
 * data. Parameters of the CTE characterization are read from the PCTE reference
 * file.
 *
 * Originally adapted from code written by Jay Anderson and rewritten by
 * Pey-Lian Lim as a standalone application that operated on FLT files.
 *
 * MRD 10 Mar 2011
 */

int doPCTE (ACSInfo *acs, SingleGroup *x) {
  
  /* arguments:
   ACSInfo *acs     i: calibration switches, etc
   SingleGroup *x    io: image to be calibrated; written to in-place
   */
  
	extern int status;
  
  /* cte scaling factor */
  double cte_frac;
  
  /* interpolated cte profile shape parameters */
  double chg_leak[MAX_TAIL_LEN*NUM_LOGQ];
  double chg_open[MAX_TAIL_LEN*NUM_LOGQ];
  
  /* interpolated profile charge parameters */
  double dtde_q[MAX_PHI];
  
  /* arrays interpolated at each parameterized charge level */
  double chg_leak_lt[MAX_TAIL_LEN*NUM_LEV];
  double chg_open_lt[MAX_TAIL_LEN*NUM_LEV];
  double dpde_l[NUM_LEV];
  int tail_len[NUM_LEV];
  
  /* structure to hold CTE parameters from file */
  CTEParams pars;
  
  /* temporary variable used during final error calculation */
  double temp_err;
  
  /* iteration variable */
  int i, j, k, m;
  
  char ccdamp[strlen(AMPSTR1)+1]; /* string to hold amps on current chip */
  int numamps;                    /* number of amps on chip */
  int amp;                        /* index amp A:0, B:1, etc. */
  char * amploc;                  /* pointer to amp character in AMPSORDER */
  int amp_arr1, amp_arr2;         /* integers to hold C array dimension sizes */
  int amp_xsize, amp_ysize;       /* integers to hold amp array size (x/y in CCD coords) */
  int amp_xbeg, amp_xend;         /* integers for beginning and end of amp arrays on chip*/
  int amp_ybeg, amp_yend;
  
  /* make arrays to hold data amp by amp.
   * these can be large arrays so it's best to declare them as pointers and
   * get the space for the ararys using malloc */
  double * amp_sci_arr; /* original sci data */
  double * amp_err_arr; /* original err data */
  double * amp_sig_arr; /* decomposed signal */
  double * amp_nse_arr; /* decomposed readout error */
  double * amp_cor_arr; /* cte corrected data */
  
  /* functions from calacs/lib */
  void parseWFCamps (char *acsamps, int chip, char *ccdamp);
  void TimeStamp (char *message, char *rootname);
  int PutKeyDbl (Hdr *hd, char *keyword, double value, char *comment);
  
  /* test whether this we've been given ACS WFC data, since the CTE algorithm
   * is currently only valid for WFC data. */
  if (acs->detector != WFC_CCD_DETECTOR) {
    trlerror("(pctecorr) only valid for WFC CCD data, PCTECORR should be OMIT for all others.");
    return (status = ERROR_RETURN);
  }
  
  if (acs->printtime) {
    TimeStamp("Starting CTE correction...","");
  }
  
  /****************** read and calculate parameters of CTE model **************/
  if (PixCteParams(acs->pcte.name, &pars)) {
    return (status);
  }
  
  if (CompareCteParams(x, &pars)) {
    return (status);
  }
  sprintf(MsgText,"(pctecorr) Max read noise amplitude PCTERNCL: %f",pars.rn_clip);
  trlmessage(MsgText);
  sprintf(MsgText,"(pctecorr) Readout simulation iterations PCTESMIT: %i",pars.sim_nit);
  trlmessage(MsgText);
  sprintf(MsgText,"(pctecorr) Number of readout shifts PCTESHFT: %i",pars.shft_nit);
  trlmessage(MsgText);
  
  cte_frac = CalcCteFrac(acs->expstart, acs->detector);
  sprintf(MsgText,"(pctecorr) CTE_FRAC: %f",cte_frac);
  trlmessage(MsgText);
  
  /* also add cte_frac as header keyword */
  if (acs->chip == 2) {
    if (PutKeyDbl(x->globalhdr,"PCTEFRAC",cte_frac,"CTE scaling factor")) {
      trlerror("(pctecorr) Error writing PCTEFRAC to image header");
      return (status = HEADER_PROBLEM);
    }
  }
  
  if (InterpolatePsi(pars.chg_leak, pars.psi_node, chg_leak, chg_open)) {
    return (status);
  }
  
  if (InterpolatePhi(pars.dtde_l, pars.q_dtde, pars.shft_nit, dtde_q)) {
    return (status);
  }
  
  if (FillLevelArrays(chg_leak, chg_open, dtde_q, pars.levels, chg_leak_lt,
                      chg_open_lt, dpde_l, tail_len)) {
    return (status);
  }
  /*********** done reading and calculating CTE model parameters **************/
  
  /* need to figure out which amps are on this chip */
  ccdamp[0] = '\0'; /* "reset" the string for reuse */
  parseWFCamps(acs->ccdamp, acs->chip, ccdamp);
  
  /* loop over amps on this chip and do CTE correction */
  numamps = strlen(ccdamp);
  for (i = 0; i < numamps; i++) {
    sprintf(MsgText,"(pctecorr) Performing CTE correction for amp %c",ccdamp[i]);
    trlmessage(MsgText);
    
    /* get the amp letter and number where A:0, B:1, etc. */
    amploc = strchr(AMPSORDER, ccdamp[i]);
    amp = *amploc - AMPSORDER[0];
    
    /* get amp array size */
    if (get_amp_array_size(acs, x, amp, amploc, ccdamp,
                           &amp_xsize, &amp_ysize, &amp_xbeg, 
                           &amp_xend, &amp_ybeg, &amp_yend)) {
      return (status);
    }
    
    amp_arr1 = amp_ysize;
    amp_arr2 = amp_xsize;
    
    /* allocate space to hold this amp's data in its various forms */
    amp_sci_arr = (double *) malloc(amp_arr1 * amp_arr2 * sizeof(double));
    amp_err_arr = (double *) malloc(amp_arr1 * amp_arr2 * sizeof(double));
    amp_sig_arr = (double *) malloc(amp_arr1 * amp_arr2 * sizeof(double));
    amp_nse_arr = (double *) malloc(amp_arr1 * amp_arr2 * sizeof(double));
    amp_cor_arr = (double *) malloc(amp_arr1 * amp_arr2 * sizeof(double));
    
    /* read data from the SingleGroup into an array containing data from
     * just one amp */
    if (make_amp_array(acs, x, amp, amp_arr1, amp_arr2, amp_xbeg, amp_ybeg, 
                       amp_sci_arr, amp_err_arr)) {
      return (status);
    }
    
    /* convert data from DN to electrons */
    for (k = 0; k < amp_arr1; k++) {
      for (m = 0; m < amp_arr2; m++) {
        amp_sci_arr[k*amp_arr2 + m] *= (double) acs->atodgain[amp];
      }
    }
    
    /* do some smoothing on the data so we don't amplify the read noise. 
     * data should be in electrons. */
    if (DecomposeRN(amp_arr1, amp_arr2, amp_sci_arr,
                    pars.rn_clip, amp_sig_arr, amp_nse_arr)) {
      return (status);
    }
    
    /* perform CTE correction */
    if (FixYCte(amp_arr1, amp_arr2, amp_sig_arr, amp_cor_arr, cte_frac,
                pars.sim_nit, pars.shft_nit, pars.levels, dpde_l, tail_len,
                chg_leak_lt, chg_open_lt)) {
      return (status);
    }
    
    /* add readout noise back and convert corrected data back to DN.
     * add 10% correction to error in quadrature. */
    for (k = 0; k < amp_arr1; k++) {
      for (m = 0; m < amp_arr2; m++) {
        amp_cor_arr[k*amp_arr2 + m] = (amp_cor_arr[k*amp_arr2 + m] 
        + amp_nse_arr[k*amp_arr2 + m]) / acs->atodgain[amp];
        
        temp_err = 0.1 * fabs(amp_cor_arr[k*amp_arr2 + m] - (amp_sci_arr[k*amp_arr2 + m]/acs->atodgain[amp]));
        amp_err_arr[k*amp_arr2 + m] = sqrt(pow(amp_err_arr[k*amp_arr2 + m],2) + pow(temp_err,2));
      }
    }
    
    /* put the CTE corrected data back into the SingleGroup structure */
    if (unmake_amp_array(acs, x, amp, amp_arr1, amp_arr2, amp_xbeg, amp_ybeg, 
                         amp_cor_arr, amp_err_arr)) {
      return (status);
    }
    
    /* free space used by our amp arrays */
    free(amp_sci_arr);
    free(amp_err_arr);
    free(amp_sig_arr);
    free(amp_nse_arr);
    free(amp_cor_arr);
  }
  
  if (acs->printtime) {
    TimeStamp("CTE corrections complete...","");
  }
   
  return (status);
}
  
/* 
 * returns the x/y dimensions for an array that holds data readout through a
 * single amp. currently only works for ACS WFC data.
 *
 * the standalone version has the array size hard wired since _flt files will
 * always have 2048 x 2048 amp regions starting at pixel 0. Here we want to be
 * a bit more careful because the overscan regions are still part of the data.
 *
 * the logic for figuring out the amp regions has been copied from doblev.
 * - MRD 14 Mar 2011
 */
static int get_amp_array_size(const ACSInfo *acs, SingleGroup *x,
                              const int amp, char *amploc, char *ccdamp,
                              int *xsize, int *ysize, int *xbeg, int *xend, 
                              int *ybeg, int *yend) {
  extern int status;
  
  int bias_loc;
  int bias_ampx, bias_ampy;
  int bias_orderx[4] = {0,1,0,1};
  int bias_ordery[4] = {0,0,1,1};
  
  int trimx1, trimx2, trimy1, trimy2;
  
  if (acs->detector == WFC_CCD_DETECTOR) {
    /* Copy out overscan info for ease of reference in this function*/
    trimx1 = acs->trimx[0];
    trimx2 = acs->trimx[1];
    trimy1 = acs->trimy[0];
    trimy2 = acs->trimy[1];
    
    bias_loc = *amploc - ccdamp[0];
    bias_ampx = bias_orderx[bias_loc];
    bias_ampy = bias_ordery[bias_loc];
    
    /* Compute range of pixels affected by each amp */	
    *xbeg = trimx1 + acs->ampx * bias_ampx;
    *xend = (bias_ampx == 0 && acs->ampx != 0) ? acs->ampx + trimx1 : x->sci.data.nx - trimx2;
    *ybeg = trimy1 + acs->ampy* bias_ampy;
    *yend = (bias_ampy == 0 && acs->ampy != 0) ? acs->ampy +trimy1 : x->sci.data.ny - trimy2; 
    
    /* Make sure that xend and yend do not extend beyond the bounds of the
     image... WJH 8 Sept 2000
     */
    if (*xend > x->sci.data.nx) *xend = x->sci.data.nx;
    if (*yend > x->sci.data.ny) *yend = x->sci.data.ny;
    *xsize = *xend - *xbeg;
    *ysize = *yend - *ybeg;	
    
  } else {
    sprintf(MsgText,"(pctecorr) Detector not supported: %i",acs->detector);
    trlerror(MsgText);
    status = ERROR_RETURN;
    return status;
  }
  
  return status;
}

/* 
 * make_amp_array returns an array view of the data readout through the
 * specified amp in which the amp is at the lower left hand corner.
 */
static int make_amp_array(const ACSInfo *acs, const SingleGroup *im, const int amp, 
                          const int arr1, const int arr2,
                          const int xbeg, const int ybeg,
                          double amp_sci_array[arr1*arr2], 
                          double amp_err_array[arr1*arr2]) {
  
  extern int status;
  
  /* iteration variables */
  int i, j;
  
  /* variables for the image row/column we want */
  int r,c;
  
  if (acs->detector == WFC_CCD_DETECTOR) {
    for (i = 0; i < arr1; i++) {
      for (j = 0; j < arr2; j++) {
        if (amp == AMP_A) {
          r = ybeg + arr1 - i - 1;
          c = xbeg + j;
        } else if (amp == AMP_B) {
          r = ybeg + arr1 - i - 1;
          c = xbeg + arr2 - j - 1;
        } else if (amp == AMP_C) {
          r = ybeg + i;
          c = xbeg + j;
        } else if (amp == AMP_D) {
          r = ybeg + i;
          c = xbeg + arr2 - j -1;
        } else {
          trlerror("Amp number not recognized, must be 0-3.");
          status = ERROR_RETURN;
          return status;
        }
        
        amp_sci_array[i*arr2 + j] = Pix(im->sci.data, c, r);
        amp_err_array[i*arr2 + j] = Pix(im->err.data, c, r);
      }
    }
  } else {
    sprintf(MsgText,"(pctecorr) Detector not supported: %i",acs->detector);
    trlerror(MsgText);
    status = ERROR_RETURN;
    return status;
  }
  
  return status;
}

/*
 * unmake_amp_array does the opposite of make_amp_array, it takes amp array
 * views and puts them back into the single group in the right order.
 */
static int unmake_amp_array(const ACSInfo *acs, const SingleGroup *im, const int amp, 
                            const int arr1, const int arr2,
                            const int xbeg, const int ybeg,
                            double amp_sci_array[arr1*arr2], 
                            double amp_err_array[arr1*arr2]) {
  
  extern int status;
  
  /* iteration variables */
  int i, j;
  
  /* variables for the image row/column we want */
  int r,c;
  
  if (acs->detector == WFC_CCD_DETECTOR) {
    for (i = 0; i < arr1; i++) {
      for (j = 0; j < arr2; j++) {
        if (amp == AMP_A) {
          r = ybeg + arr1 - i - 1;
          c = xbeg + j;
        } else if (amp == AMP_B) {
          r = ybeg + arr1 - i - 1;
          c = xbeg + arr2 - j - 1;
        } else if (amp == AMP_C) {
          r = ybeg + i;
          c = xbeg + j;
        } else if (amp == AMP_D) {
          r = ybeg + i;
          c = xbeg + arr2 - j -1;
        } else {
          trlerror("Amp number not recognized, must be 0-3.");
          status = ERROR_RETURN;
          return status;
        }
        
        Pix(im->sci.data, c, r) = (float) amp_sci_array[i*arr2 + j];
        Pix(im->err.data, c, r) = (float) amp_err_array[i*arr2 + j];
      }
    }
  } else {
    sprintf(MsgText,"(pctecorr) Detector not supported: %i",acs->detector);
    trlerror(MsgText);
    status = ERROR_RETURN;
    return status;
  }
  
  return status;
}

