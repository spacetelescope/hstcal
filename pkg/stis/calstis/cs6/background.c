/*
    This file contains routines that handle background processing:
        CalcBack
        X1DBack
	SmoothBack



   Revision history:
   ----------------
   20 Feb 97  -  Implemented (I.Busko)
   11 Apr 97  -  OPR 33790: flag background values that had more than 30%
                 of their pixels affected by serious DQ flags (IB).
   14 Apr 97  -  OPR 33792: compute background errors using fit scatter
                 instead of input pixel errors (IB).
   15 Apr 97  -  OPR 33789: perform sigma-clipping in background fit (IB).
   28 Apr 97  -  Fixed error in background box coordinate computation (IB).
   28 Jul 97  -  Print 1-indexed debug coordinates (IB).
   10 Apr 98  -  Remove debug file (IB).
   14 Oct 98  -  Add X1D_BAD_BACKGROUND macro (IB)
   16 Dec 98  -  Output background error: rms around fit / sqrt(n)
   08 Sep 99  -  Fixed sigma-clip test (IB)
   13 Jan 00  -  Scattered light correction (IB)
   07 Mar 01  -  Wrong polynomial order defaults to zero order (IB)
   16 Oct 01  -  Out-of-aperture pixels were being counted in (IB)
   05 Aug 02  -  Background smoothing (IB)
   12 Jan 05  -  Turn off background smoothing for first-order FUV if
                 there are too few points to fit a polynomial (presumably
                 due to Lyman Alpha avoidance with a large aperture) (PEH)
   07 Feb 06  -  Add an argument to Interp2D (PEH)
   28 Aug 13  -  Use integer loop index and limits in X1DBack (PEH)
*/

# include <stdio.h>
# include <math.h>
# include <float.h>

# include "xtables.h"
# include "hstio.h"

# include "stis.h"
# include "stisdef.h"
# include "err.h"
# include "stisdq.h"
# include "calstis6.h"

# define MAX_CLIP	5	/* max. number of sigma-clip iterations */

static void X1DBack (StisInfo6 *, SingleGroup *, FloatHdrData *,
                     FloatHdrData *, double, double, double, double *,
                     double *, double *, double *, double *, double *,
                     double *, int *, int *, int *, double, double *,
                     double *, double *, double *, int *, int);
static int PolynomialFilter (float *, int, int, int, int, int, int);
static void BoxcarFilter (float *, int, int, int);
static float median (float *, int);
static int isValid (int, int, int, int, int, int);
static void copyArray (float *, float *, int);


/* Compute background coefficients. Zero, first and nth order polynomials
   are supported. The zero-order background is computed as the unweighted
   average of all eligible pixels in both background extraction boxes. The
   first order polynomial is fitted to eligible pixels in each background
   extraction box. The independent variable is the pixel index in the A2
   direction (in physical pixel array coordinates, not reference
   coordinates) and the dependent variable is the pixel content.

   7th order is used exclusively by the scattered light correction
   algorithm. It was implemented later and the code had to be tweaked
   in order to support it. Note that sigma clipping is effectively
   turned off when scattered light correction is taking place.

   The background DQ flag is set according to the total fraction of flagged
   pixels in the background extraction boxes.

*/

int CalcBack (StisInfo6 *sts, XtractInfo *xtr, SingleGroup *in,
              FloatHdrData *ssgx, FloatHdrData *ssgy, int ipix,
              double ysp, int debug) {

/* arguments:
StisInfo6 *sts         i: calibration switches and info
XtractInfo *xtr        i: extraction parameters
SingleGroup *in	       i: input image
FloatHdrData ssgx;     i: small-scale distortion in X (not used)
FloatHdrData ssgy;     i: small-scale distortion in Y (not used)
int ipix;              i: index of image pixel in the A1 direction
double ysp;            i: center of spectrum extraction box
int debug;             i: debug control
*/
	double rpix;            /* index of reference pixel in A1 direction */
	double xcent[2];	/* center of background boxes */
	double ycent[2];
	double bksize[2];       /* sizes of background boxes */
	double bkoffset[2];     /* offsets of background boxes */
	double npts;		/* accumulators used for zero and first */
	double sumx;		/* order only */
	double sumy;
	double sumx2;
	double sumy2;
	double sumxy;
	double sumvar;
	double delta;
	double sigma;		/* scatter around current fit */
	double new_sigma;	/* updated scatter around current fit */
	int    nbck;		/* total number of pixels in backgr. boxes */
	int    nfbck;		/* discarded pixels due to flagging */
	int    nsbck;		/* discarded pixels due to sigma-clipping */
	double  *yval;		/* data for polynomial fit */
	double  *xval;
	double  *wval;
	int    ndata;		/* # of data points in polynomial fit */
	int    backord, iclip, maxclip, i;

	int FitPoly (double *, double *, double *, int, int, double *);

	/* See if there are command-line overrides. Overrides are disabled
           in case the scatterd light correction algorithm is active.
        */
	if (!(sts->scatter)) {
	    for (i = 0; i < 2; i++) {
	        if (sts->bksize[i] == NO_SIZE)
	            bksize[i] = xtr->bksize[i];
	        else
	            bksize[i] = sts->bksize[i];
	        if (sts->bkoffset[i] == NO_SIZE)
	            bkoffset[i] = xtr->bkoffset[i];
	        else
	            bkoffset[i] = sts->bkoffset[i];
	    }
	    if (sts->bkord == NO_ORDER)
	        backord = xtr->backord;
	    else
	        backord = sts->bkord;
	} else {
	    /* Size and offset values are supposed to be re-defined by
               previous call to DefineBackRegions function.
            */
	    for (i = 0; i < 2; i++) {
	        bksize[i]   = sts->bksize[i];
	        bkoffset[i] = sts->bkoffset[i];
	    }
	    backord = BACKP;
	}

	/* Translate image pixel index into reference pixel index. */
	rpix = (ipix - sts->ltv[0]) / sts->ltm[0];

	/* Compute centers of background boxes. */
	xcent[0] = rpix + bkoffset[0] * sts->sin_bktilt;
	ycent[0] = ysp  + bkoffset[0] * sts->cos_bktilt;
	xcent[1] = rpix + bkoffset[1] * sts->sin_bktilt;
	ycent[1] = ysp  + bkoffset[1] * sts->cos_bktilt;

	sigma = DBL_MAX / 3.1; /*make sure first 3*sigma computation succeeds*/
	sts->ebck   = 0.0;
	sts->bck[0] = 0.0;
	sts->bck[1] = 0.0;;

	/* Main sigma-clip loop. */
	maxclip = (sts->scatter) ? 1 : MAX_CLIP;

	for (iclip = 1; iclip <= maxclip; iclip++) {

	    /* Clear all accumulators. */
	    npts   = 0.0;
	    sumx   = 0.0;
	    sumy   = 0.0;
	    sumx2  = 0.0;
	    sumy2  = 0.0;
	    sumxy  = 0.0;
	    sumvar = 0.0;
	    new_sigma = 0.0;
	    nbck  = 0;
	    nfbck = 0;
	    nsbck = 0;

	    /* Alloc memory used in 7th order background fit. This must be
               done even when the scattered light correction algorithm is
               not active, the reason being to avoid a possible addressing
               error when the calling sequence to X1DBack gets executed.
	       Note that in this case the sigma-clip loop is executed
               only once.
            */
	    if ((xval = (double *) calloc (200, sizeof (double))) == NULL)
	        return (OUT_OF_MEMORY);
	    if ((yval = (double *) calloc (200, sizeof (double))) == NULL)
	        return (OUT_OF_MEMORY);
	    if ((wval = (double *) calloc (200, sizeof (double))) == NULL)
	        return (OUT_OF_MEMORY);
	    ndata = 0;

	    /* Extract background in each box. */

	    X1DBack (sts, in, ssgx, ssgy, xcent[0], ycent[0], bksize[0],
                     &npts, &sumx, &sumy, &sumx2, &sumy2, &sumxy, &sumvar,
                     &nbck, &nfbck, &nsbck, sigma, &new_sigma,
                     xval, yval, wval, &ndata, debug);

	    X1DBack (sts, in, ssgx, ssgy, xcent[1], ycent[1], bksize[1],
                     &npts, &sumx, &sumy, &sumx2, &sumy2, &sumxy, &sumvar,
                     &nbck, &nfbck, &nsbck, sigma, &new_sigma,
                     xval, yval, wval, &ndata, debug);

	    /* If more than 30% of the background pixels were rejected
               by sigma-clip, or flagged, set background's 11th flag bit.
            */
	    if ((float)(nfbck+nsbck) / (float)nbck > 0.33333F)
	        sts->dqbck = X1D_BAD_BACKGROUND;
	    else
	        sts->dqbck = 0;

	    /* Compute background coefficients. Only unweighted fits are
               implemented in this version, so sumvar is never used.
            */
	    switch (backord) {
	    case 0:
	    default:
	        if (npts > 0.0) {
	            sts->bck[0]  = sumy / npts;
	            if (npts > 1.0)
	                sts->vbck[0] = (sumy2 - sumy*sumy/npts) /
                                       (npts - 1.0);
	            else
	                sts->vbck[0] = 0.0;
	        } else {
	            sts->bck[0]  = 0.0;
	            sts->vbck[0] = 0.0;
	        }
	        sts->bck[1]  = 0.0;
	        sts->vbck[1] = 0.0;
	        break;
	    case 1:
	        delta = npts * sumx2 - sumx * sumx;
	        if (fabs(delta) > 0.0) {
	            sts->bck[0] = (sumx2 * sumy - sumx * sumxy) / delta;
	            sts->bck[1] = (sumxy * npts - sumx * sumy)  / delta;
	            sts->vbck[0] = sumx2 / delta;
	            sts->vbck[1] = npts  / delta;
	        } else {
	            sts->bck[0]  = 0.0;
	            sts->vbck[0] = 0.0;
	            sts->bck[1]  = 0.0;
	            sts->vbck[1] = 0.0;
	        }
	        break;
	    case BACKP:
	        /* This code does NOT compute the coefficient variances !
                */
	        if (ndata > BACKP) {
	            if (!FitPoly (xval, yval, wval, ndata, BACKP, sts->bck)) {
	                for (i = 0; i < BACKP+3; i++)
	                    sts->vbck[i] = 0.0;
	            } else {
	                for (i = 0; i < BACKP+3; i++) {
	                    sts->bck[i]  = 0.0;
	                    sts->vbck[i] = 0.0;
	                }
	            }
	        } else {
	            for (i = 0; i < BACKP+3; i++) {
	                sts->bck[i]  = 0.0;
	                sts->vbck[i] = 0.0;
	            }
	        }
	        break;
/*
	    default:
	        sts->bck[0]  = 0.0;
	        sts->vbck[0] = 0.0;
	        sts->bck[1]  = 0.0;
	        sts->vbck[1] = 0.0;
	        break;
*/
	    }

	    /* Free memory used in 7th order background fit. */
	    free (wval);
	    free (yval);
	    free (xval);

	    /* Update sigma for next iteration. */
	    if (npts > 1.0) {
	        sigma = sqrt (new_sigma / (npts - 1.0));
	        sts->ebck = sigma;
	    }
	    else
	        break;

	    /* Exit the loop if no pixels where discarded by sigma-clip. */
	    if (nsbck == 0 && iclip > 1)
	        break;
	}

	return (0);
}




/*  Extract background pixels and update fitting accumulators in a given
    background extraction box. The interpolation routine does *not* take
    into account that interpolated pixels may have a tilt respect to the
    input pixel array.
*/

static void X1DBack (StisInfo6 *sts, SingleGroup *in, FloatHdrData *ssgx,
                     FloatHdrData *ssgy, double xcenter, double ycenter,
                     double size, double *npts, double *sumx, double *sumy,
                     double *sumx2, double *sumy2, double *sumxy,
                     double *sumvar, int *nbck, int *nfbck, int *nsbck,
                     double sigma, double *new_sigma, double *xval,
                     double *yval, double *wval, int *ndata, int debug) {

/* arguments:
StisInfo6 *sts         i:  calibration switches and info
SingleGroup *in	       i:  input image
FloatHdrData ssgx;     i:  small-scale distortion in X (not used)
FloatHdrData ssgy;     i:  small-scale distortion in Y (not used)
double xcenter;        i:  center of spectrum extraction box in A1 direction
double ycenter;        i:  center of spectrum extraction box in A2 direction
double size;           i:  size of box
double *npts...        o:  accumulators
int *nbck;             o:  total number of accumulated pixels
int *nfbck;            o:  number of flagged pixels
int *nsbck;            o:  number of discarded pixels in this s-clip iteration
double sigma;          i:  scatter around current fit
double *new_sigma;     io: updated scatter
double  *yval;         io: data for polynomial fit
double  *xval;
double  *wval;
int    ndata;          io: # of data points in polynomial fit
int debug;             i: debug control
*/
	double x1, y1, y2;	/* end points of background extraction box */
	double xx, yy;		/* coordinates of current pixel */
	float  oSci, oErr;	/* interpolated values from image array */
	short  oDQ;
	double residual;	/* residual background */
	double tcent;		/* center of trace in image coordinates */
	double hold;
        int k, isize;           /* loop index and limit */

	/* Compute endpoints of background extraction box. Notice that
           the 0.5 pixel offsets below, used to make the endpoint coordinates
           point to the pixel center, only make sense for small tilt angles.
           In this case most of the 0.5 pixel correction goes into the A2
           direction. A more sophisticated correction is needed for large
           angles.
        */
	x1 = xcenter - (size / 2.0 * sts->sin_bktilt);
	y1 = ycenter - (size / 2.0 * sts->cos_bktilt) + 0.5;
	y2 = ycenter + (size / 2.0 * sts->cos_bktilt) - 0.5;

	/* Translate reference pixel coordinates into physical indices. */
	x1    = x1      * sts->ltm[0] + sts->ltv[0];
	y1    = y1      * sts->ltm[1] + sts->ltv[1];
	y2    = y2      * sts->ltm[1] + sts->ltv[1];
	tcent = ycenter * sts->ltm[1] + sts->ltv[1];

	/* Scan background box, extract pixels and update accumulators. */
        xx = x1;
        yy = y1;
        isize = NINT(size);
	for (k = 0;  k < isize;  k++) {

	    /* Interpolate, checking for out of bounds. */
	    Interp2D (in, sts->sdqflags, xx, yy, 1.0, WGT_VARIANCE,
			&oSci, &oErr, &oDQ);

	    /* Update accumulators, discarding flagged data and outliers.
               Notice that background is computed as a function of physical
               pixel coordinates. Two separate sets of tests are used because
               each background processing algorithm has a different set of
               requirements.
            */
	    (*nbck)++;
	    if (!(sts->scatter)) {

	        /* for zero-th and 1st order */

                /* Fixed 16-Oct-01: out-of-aperture dq bit must be
                   explicitly ORed in.
                */
	        if (!(oDQ & (sts->sdqflags | DETECTORPROB))) {
	            residual = fabs (oSci - (sts->bck[0] + sts->bck[1] * yy));
	            if (residual <= (3.0 * sigma)) {
	                *npts   += 1.0;
	                *sumvar += oErr * oErr;
	                *sumx   += yy;
	                *sumy   += oSci;
	                *sumx2  += yy * yy;
	                *sumxy  += yy * oSci;
	                *sumy2  += oSci * oSci;
	                *new_sigma += residual * residual;
	            } else
	                (*nsbck)++;
	        } else
	            (*nfbck)++;
	    } else {

	        /* for 7th order */
	        if (!(oDQ & DETECTORPROB & sts->sdqflags)) {
	            xval[*ndata] = yy;
	            yval[*ndata] = (double)oSci;
	            hold = yy - tcent;
	            if (hold < 0)
	                hold = -hold;
	            if (hold > 0)
	                hold = 1.0 / sqrt (hold);
	            else
	                hold = 1.0;
	            wval[*ndata] = hold;
	            (*ndata)++;
	        } else
	            (*nfbck)++;
	    }
            xx += sts->sin_bktilt;
            yy += sts->cos_bktilt;
	}
}


/*  Smooths the background.

    The algorithm actually used to smooth the background is chosen based on
    the specific detector/grating combination:

    - CCD:  median smooth with a running window 9 pixels wide, then fit a
            n-th degree polynomial. Optionally replace the median by a simple
            average.

    - MAMA 1st order: fit a n-th degree polynomial.

    - MAMA echelle: average with a running window 31 pixels wide.

    This is documented in OPR #46279 and private e-mail with Ralph Bohlin.

    Note that the Lya avoidance applies only to first order MAMA.

*/

int SmoothBack (StisInfo6 *sts, RowContents *row_contents) {

	int i, status;

	if (sts->bks_mode == BKS_OFF) {
	    return (status = 0);
	}

	if (sts->detector == CCD_DETECTOR) {

	    BoxcarFilter (row_contents->back, row_contents->npts, 9,
                          sts->bks_mode);

	    if ((status = PolynomialFilter (row_contents->back,
                                           row_contents->npts,
                                           sts->bks_order,
                                            0, 0, 0, 0))) {
	        return (status);
	    }

	} else if (((sts->detector == NUV_MAMA_DETECTOR) ||
                    (sts->detector == FUV_MAMA_DETECTOR))
                     && !sts->echelle) {

	    status = PolynomialFilter (row_contents->back,
	                               row_contents->npts,
	                               sts->bks_order,
	                               sts->avoid1a, sts->avoid2a,
	                               sts->avoid1b, sts->avoid2b);
	    if (status < 0 && sts->detector == FUV_MAMA_DETECTOR) {
		printf ("Warning  There was not enough background data" \
			" to fit a polynomial, probably due\n");
		printf ("Warning  to the use of a large aperture;" \
			" background smoothing will not be done.\n");
		sts->bks_mode = BKS_OFF;
		return 0;
	    } else if (status != 0) {
		return (status);
	    }

        } else if (((sts->detector == NUV_MAMA_DETECTOR) ||
                    (sts->detector == FUV_MAMA_DETECTOR))
                     && sts->echelle) {

	    BoxcarFilter (row_contents->back, row_contents->npts, 31,
                          BKS_AVERAGE);
	    BoxcarFilter (row_contents->back, row_contents->npts, 31,
                          BKS_AVERAGE);
	}

	/* Re-compute net array. */

	for (i = 0; i < row_contents->npts; i++) {
	    row_contents->net[i] = row_contents->gross[i] -
                                   row_contents->back[i];
	}

	return (0);
}


static int PolynomialFilter (float *array, int size, int order,
                             int avoid1a, int avoid2a, int avoid1b,
                             int avoid2b) {

	double *ax, *ay, *aw, *coeff;
	int i, k, status;

	int FitPoly (double *, double *, double *, int, int, double *);
	void ComputePoly (double *, int, double *, int, double *);

	/* Adjust values found by the Lya-finding function
           as per requirements.
        */

	if (avoid1a != 0) {
	    avoid1a += 5;
	    avoid2a -= 5;
	    avoid1b += 5;
	    avoid2b -= 5;
	}

	/* The polynomial functions require double arrays. */

	if ((ax = (double *)calloc(size, sizeof (double))) == NULL)
	    return (OUT_OF_MEMORY);
	if ((ay = (double *)calloc(size, sizeof (double))) == NULL)
	    return (OUT_OF_MEMORY);
	if ((aw = (double *)calloc(size, sizeof (double))) == NULL)
	    return (OUT_OF_MEMORY);

	if ((coeff = (double *)calloc(order+3, sizeof (double))) == NULL)
	    return (OUT_OF_MEMORY);

	/* Here we build arrays with the subset of data points to
           be included in the fit. We cannot use the weight array
           in the fit functions for that.
	*/
	k = 0;
	for (i = 0; i < size; i++) {
	    if (isValid (i, size, avoid1a, avoid2a, avoid1b, avoid2b)) {
	        ax[k] = (double)i;
	        ay[k] = (double)array[i];
	        aw[k] = 1.0;
	        k++;
	    }
	}

	/* Fit and compute polynomial over subset. */

	if ((status = FitPoly (ax, ay, aw, k, order, coeff)) != 0) {
	    free (aw);
	    free (ay);
	    free (ax);
	    return (status);
	}
	ComputePoly (ax, k, coeff, order, ay);

	/* Re-fill the output array with smoothed data. The voids
           are left with the original input data.
        */
	k = 0;
	for (i = 0; i < size; i++) {
	    if (isValid (i, size, avoid1a, avoid2a, avoid1b, avoid2b))
	        array[i] = (float)ay[k++];
	}

	free (aw);
	free (ay);
	free (ax);

	return (0);
}


static void BoxcarFilter (float *input_array, int size, int wsize, int mode) {

	float *array;
	double sum;
	int i, j, w1, w2, vw2, npts;

	array = (float *) malloc (size * sizeof (float));
	copyArray (input_array, array, size);

	for (i = 0; i < size; i++) {
	    switch (mode) {
	        case BKS_AVERAGE:

	            if (i == 0) {
	                 sum  = 0.0;
	                 npts = 0;
	                 w1   = 0;
	                 w2   = wsize / 2;
	                 for (j = 0; j < w2; j++) {
	                     sum += array[j];
	                     npts++;
	                 }
	            } else {
	                w2++;
	                if (w2 < size) {
	                    sum += array[w2];
	                    npts++;
	                }

	                w1 = w2 - wsize;

	                if (w1 > 0) {
	                    sum -= array[w1-1];
	                    npts--;
	                }
	            }

	            input_array[i] = sum / npts;

	            break;

	        case BKS_MEDIAN:

	            if (i == 0) {
	                 w1  = 0;
	                 w2  = wsize / 2;
	                 vw2 = w2;
	            } else {
	                w2++;
	                vw2++;
	                if (w2 >= size) {
	                    w2--;
	                }

	                w1 = vw2 - wsize;

	                if (w1 < 0) {
	                    w1 = 0;
	                }
	            }

	            input_array[i] = median (array + w1, w2 - w1 + 1);

	            break;
	    }
	}

	free (array);
}


static void copyArray (float *input_array, float *array, int size) {

	int i;

	for (i = 0; i < size; i++) {
	    array[i] = input_array[i];
	}
}


static float median (float *window, int size) {

	float MedianFloat (float *, int, int);

	return (MedianFloat (window, size, 0));
}


static int isValid (int index, int size, int avoid1a, int avoid2a,
                                         int avoid1b, int avoid2b) {
	int include;

	/* The logic here is: if near array endpoints, always avoid.
           Then check for the avoidance region; if it is not present,
           include always. Only then check for the avoidance limits.
        */

	if (index < 10 || index > size-10)
	    return (0);

	if (avoid1a == 0)
	    return (1);

	include = 1;

	if (index > avoid1a && index < avoid2a)
	    include = 0;
	if (index > avoid1b && index < avoid2b)
	    include = 0;

	return (include);
}
