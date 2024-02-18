#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "hstcal.h"
#include "hstio.h"

#include "acs.h"
#include "acsinfo.h"
#include "hstcalerr.h"
#include "trlbuf.h"

#include "pcte.h"


int get_amp_array_size_acs_cte(const ACSInfo *acs, SingleGroup *x,
                              char *amploc, char *ccdamp,
                              int *xsize, int *ysize, int *xbeg,
                              int *xend, int *ybeg, int *yend);
static int make_amp_array(const ACSInfo *acs, const SingleGroup *im,
                          const int amp,
                          const int arr1, const int arr2,
                          const int xbeg, const int ybeg,
                          double amp_sci_array[arr1*arr2],
                          double amp_err_array[arr1*arr2]);
static int unmake_amp_array(const ACSInfo *acs, const SingleGroup *im,
                            const int amp,
                            const int arr1, const int arr2,
                            const int xbeg, const int ybeg,
                            double amp_sci_array[arr1*arr2],
                            double amp_err_array[arr1*arr2]);


/* Returns the x/y dimensions for an array that holds data readout through a
   single amp. currently only works for ACS WFC data.

   the standalone version has the array size hard wired since _flt files will
   always have 2048 x 2048 amp regions starting at pixel 0. Here we want to be
   a bit more careful because the overscan regions are still part of the data.

   the logic for figuring out the amp regions has been copied from doblev.
   - MRD 14 Mar 2011
*/
int get_amp_array_size_acs_cte(const ACSInfo *acs, SingleGroup *x,
                              char *amploc, char *ccdamp,
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
        *xbeg = (trimx1 + acs->ampx) * bias_ampx;
        *xend = (bias_ampx == 0 && acs->ampx != 0) ? acs->ampx + trimx1 : x->sci.data.nx;
        *ybeg = (trimy1 + acs->ampy) * bias_ampy;
        *yend = (bias_ampy == 0 && acs->ampy != 0) ? acs->ampy + trimy1 : x->sci.data.ny;
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


/* Make_amp_array returns an array view of the data readout through the
   specified amp in which the amp is at the lower left hand corner.
*/
static int make_amp_array(const ACSInfo *acs, const SingleGroup *im,
                          const int amp,
                          const int arr1, const int arr2,
                          const int xbeg, const int ybeg,
                          double amp_sci_array[arr1*arr2],
                          double amp_err_array[arr1*arr2]) {

    extern int status;

    /* iteration variables */
    int i, j;

    /* variables for the image row/column we want */
    int r, c;

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


/* unmake_amp_array does the opposite of make_amp_array, it takes amp array
   views and puts them back into the single group in the right order.
*/
static int unmake_amp_array(const ACSInfo *acs, const SingleGroup *im,
                            const int amp,
                            const int arr1, const int arr2,
                            const int xbeg, const int ybeg,
                            double amp_sci_array[arr1*arr2],
                            double amp_err_array[arr1*arr2]) {

    extern int status;

    /* iteration variables */
    int i, j;

    /* variables for the image row/column we want */
    int r, c;

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
