/*  This file contains:
	BlevDrift	fit to virtual overscan
	DriftEval	evaluate the fit for a given line number
	DriftMean	return the mean value for the drift

	and the following local functions:
	  VMedianY	take median in a column, to reject outliers
	  DriftInit	zero sums for a new fit
	  DriftAccum	increment sums for a new (line, value) point
	  DriftFit	compute the coefficients of the fit
	  DriftSet	set coefficients of fit to a given value
*/

# include <stdio.h>
# include <stdlib.h>
#include "hstcal.h"
# include "hstio.h"
# include <math.h>
# include "wf3.h"
# include "wf3info.h"
# include "hstcalerr.h"
# include "wf3dq.h"		/* for GOODPIXEL */

# define NELEM_SUMS  5		/* size of sums array */

static int VMedianY (SingleGroup *, int, int *, short, double *, double *);
static void cleanDriftFit (double *, int *, int *, float);
static void DriftInit (double, int);
static void DriftAccum (int, double);
static int DriftFit (void);
static void DriftSet (double);

/* These are used by DriftInit, ... DriftEval. */
static double sums[NELEM_SUMS];
static double slope;		/* intercept = 0. */
static double zero_column;
static int middle_col;

/* This routine fits a line to the virtual overscan values as a function
   of column number.

   When evaluating the fit using DriftEval, the independent variable
   should be X pixel number minus the middle of biassect.  We achieve this
   by ignoring the intercept of the fit and subtracting the middle of
   biassect from the X pixel number in DriftEval.

   Warren Hack, 1998 June 2:
	Revised to use ACS header files... Otherwise basically unchanged...
   Howard Bushouse, 2000 Aug 29:
	Revised for WFC3 use.
   H.Bushouse, 2002 March 1:
	Corrected computation of zerocol to be units of input, untrimmed
	image pixels, rather than output, trimmed image pixels.
   H.Bushouse, 2008 Sep 17:
	Upgraded to reject outliers from parallel overscan array before
	fitting (as is already done for serial overscan fit). Added new
	routine cleanDriftFit to reject outliers (equivalent to serial
	cleanBiasFit). Added readnoise as input argument, which is used in 
	cleanDriftFit.
   H.Bushouse, 2009 Jan 16:
	Upgraded the methods used in cleanDriftFit to compute the mean and
	standard deviation of unrejected values so that they use only the
	good values returned by VMedianY. Also added checks for potential 
	divide-by-zero conditions.
*/

int BlevDrift (SingleGroup *in, int *vx, int *vy, int trimx1, int *biassect,
	       int *driftcorr, short sdqflags, float rn) {

/* arguments:
SingleGroup *in      i: image to be calibrated
int vx[2], vy[2]     i: range of pixel numbers for virtual overscan region
int trimx1           i: width to trim off beginning of line
int biassect[4]      i: section to use for finding bias level
int *driftcorr       i: true if correction can be applied
float rn             i: readnoise (units of DN)
*/

	extern int status;
	double *scratch;
	double value;		/* median of values in column */
	double zerocol;		/* zero point for fit */
	int midcol;		/* middle column of overscan region */
	int i;			/* X pixel in input image */
	double *biasvals;	/* intermediate array for bias level values */
	int    *biasmask;	/* mask array for biasvals:0 means don't use */
	char nodriftcorr[] =
		{"No correction for slope will be applied."};

	*driftcorr = 0;		/* initial value */

	if (vx[1] <= vx[0] || vy[1] <= vy[0]) {
	    trlmessage ("(blevcorr) No virtual overscan region;");
	    trlmessage (nodriftcorr);
	    DriftSet (0.);
	    return (status);
	}

	/* Allocate space for temporary arrays */
	if ((scratch = malloc ((vy[1]-vy[0]+1) * sizeof (double))) == NULL) {
	    trlwarn ("Out of memory in BlevDrift.");
	    return (status = OUT_OF_MEMORY);
	}
	biasvals = (double *) calloc (vx[1]+1, sizeof(double));
	biasmask = (int *) calloc (vx[1]+1, sizeof(int));

	/* Initialize for fitting.  The first argument is the middle of
	   the section used for determining the bias level, in units of
	   pixel number in the input, untrimmed, image.
	*/
	
	/* HAB 1-Mar-2002: Modified zerocol computation so that it is
	** pixel units of the input, untrimmed image, rather than the
	** STIS approach of pixel units in the output, trimmed image. */
	/*zerocol = (double) (biassect[0] + biassect[1]) / 2. - trimx1;*/
	zerocol = (double) (biassect[0] + biassect[1]) / 2.;
	midcol = (vx[0] + vx[1]) / 2;
	DriftInit (zerocol, midcol);

	/* For each column, determine the value of the virtual overscan,
	** and accumulate sums. */
	for (i = vx[0];  i <= vx[1];  i++) {

	    /* Compute median of column */
	    if (VMedianY (in, i, vy, sdqflags, &value, scratch))
		continue;
	    biasvals[i] = value;
	    biasmask[i] = 1;
	}

	/* Analyze biasvals for outliers and mask them by setting the
	** corresponding value in the biasmask array to zero. */
	cleanDriftFit (biasvals, biasmask, vx, rn);

	/* Now that we have read in all bias level values from the overscan
	** regions and thrown out any outliers (cosmic-ray hits), we can
	** accumulate the fitting statistics now. */
	for (i = vx[0];  i <= vx[1];  i++) {
	     if (biasmask[i] == 1)
		 DriftAccum (i, biasvals[i]);
	}

	/* Fit a curve to the values found. */
	if (DriftFit()) {
	    trlwarn ("(blevcorr) Singular fit to virtual overscan; ");
	    trlwarn (nodriftcorr);
	    DriftSet (0.);
	} else {
	    *driftcorr = 1;
	}

	/* Free up local memory */
	free (scratch);
	free (biasvals);
	free (biasmask);

	return (status);
}

/* This routine takes the median of the SCI data values at x = i
   over y = vy[0] to vy[1] inclusive.  Note that the independent variable
   i is the column number in the input image.
   Using a median instead of a mean is to reject outliers.

   This function does not set status.  The function value will be zero
   if there is at least one data value not flagged as bad, but if there
   is no data, the function value will be NO_GOOD_DATA.


   H. Bushouse, 2006 June 20: Fixed bug in logic for rejecting flagged
   pixels using the sdqflags variable.

*/

static int VMedianY (SingleGroup *in, int i, int *vy, short sdqflags, 
		     double *median, double *scratch) {

/* arguments:
SingleGroup *in    i: image to be calibrated
int i              i: X pixel location in input image
int vy[2]          i: range of pixel numbers in Y direction
double *median     o: median of values in column
double scratch[]  io: scratch space
*/

	int j;		/* loop indices */
	int nvals;		/* number of good pixels */
	int inplace = 1;	/* OK to sort in-place */

	double MedianDouble (double *, int, int);

	nvals = 0;
	for (j = vy[0];  j < vy[1];  j++) {
	    if (DQPix (in->dq.data, i, j) == GOODPIXEL ||
		!(DQPix (in->dq.data, i, j) & sdqflags)) {
		scratch[nvals] = Pix (in->sci.data, i, j);
		nvals++;
	    }
	}

	if (nvals < 1)
	    return (NO_GOOD_DATA);

	*median = MedianDouble (scratch, nvals, inplace);

	return (0);
}

/* This routine uses iterative sigma clipping to reject outliers from
   the array of bias values.
*/

void cleanDriftFit (double *barray, int *bmask, int *vx, float rn) {

        int j;
        double bsum, bmean, sdev, svar;
        double s;
        int nsum;
        float clip;
        int nrej=0;

        bsum = bmean = sdev = svar = 0.0;
	nsum = 0;

        for (j=vx[0]; j <= vx[1]; j++) {
	     if (bmask[j] == 1) {
		 bsum += barray[j];
		 nsum++;
	     }
	}
	if (nsum == 0) return;

        bmean = bsum / nsum;

        for (j=vx[0]; j <= vx[1]; j++) {
	     if (bmask[j] == 1) {
		 s = barray[j] - bmean;
		 svar += (s*s);
	     }
        }
        sdev = sqrt (svar/nsum);

        /* Reset stddev from mean to stddev of poisson dist centered
        ** on mean. This will keep cosmic-ray hits or bleeding from
        ** bright sources into the overscan columns from biasing
        ** the fit. HAB 19-Feb-2004
        */
        if (sdev > sqrt(bmean)) sdev = sqrt(bmean);
	clip = 3.5;

        /* With statistics in hand, ID and flag outliers */
        for (j=vx[0]; j <= vx[1]; j++) {
             if (barray[j] > clip*sdev+bmean) {
                 bmask[j] = 0;
                 nrej++;
             }
        }

        /* Recompute the mean based on clipped values. */
        bsum = bmean = 0.0;
        nsum = 0;
        for (j=vx[0]; j <= vx[1]; j++) {
             /* if value has not already been thrown out, use it
             ** to compute new mean ... */
             if (bmask[j] != 0) {
                 bsum += barray[j];
                 nsum++;
             }
        }
	if (nsum == 0) return;
        bmean = bsum / nsum;
	clip = 2.0;

        /* With statistics in hand, ID and flag outliers based on
        ** readnoise as sigma to further refine the value ... */
        for (j=vx[0]; j <= vx[1]; j++) {
             if (barray[j] > clip*rn+bmean) {
		 if (bmask[j] != 0) {
                     bmask[j] = 0;
                     nrej++;
		 }
             }
        }

        sprintf (MsgText,
		"(blevcorr) Rejected %d bias values from parallel fit.",nrej);
        trlmessage (MsgText);
}

/* This routine resets the sums and coefficients to zero, and it copies
   zerocol to zero_column and midcol to middle_col.  Note that zerocol
   should be an X pixel number in the input, untrimmed, image.
*/

static void DriftInit (double zerocol, int midcol) {

/* argument:
double zerocol     i: fit will be forced to zero at this column number
int midcol         i: middle column of image
*/

	int i;

	for (i = 0;  i < NELEM_SUMS;  i++)
	    sums[i] = 0.;
	slope = 0.;
	zero_column = zerocol;
	middle_col = midcol;
}

/* This routine increments the sums for an unweighted linear fit of
   virtual overscan as a function of i, the column number.  We subtract
   middle_col from i before accumulating the sums to reduce numerical
   roundoff problems, but we will ignore middle_col when evaluating
   the fit because we're ignoring the intercept of the fit and using
   a different location as the zero point.
*/

static void DriftAccum (int i, double value) {

/* arguments:
int i          i: independent variable, X pixel number
double value   i: dependent variable, the virtual overscan at this column
*/

	double x;

	/* Subtract middle column number to reduce numerical problems. */
	x = (double) (i - middle_col);

	sums[0]++;
	sums[1] += x;
	sums[2] += value;
	sums[3] += (x * value);
	sums[4] += (x * x);
}

/* This routine computes the coefficients of fit (slope only).
   If there is no data, the function value will be -1, if the fit is
   singular, the function value will be -2; otherwise, the function will
   return zero.
*/

static int DriftFit (void) {

	double d;
	double xmean, ymean;	/* mean values of column and value */

	if (sums[0] < 1.)
	    return (-1);

	d = sums[4] - sums[1] * sums[1] / sums[0];
	if (d == 0.)
	    return (-2);

	xmean = sums[1] / sums[0];
	ymean = sums[2] / sums[0];

	slope = (sums[3] - xmean * ymean * sums[0]) / d;

        sprintf (MsgText, "Computed a parallel fit with slope of %g", slope);
        trlmessage (MsgText);

	return (0);
}

/* This routine sets the coefficients of fit to specific values so that
   DriftEval will return "value" regardless of the column number.
*/

static void DriftSet (double value) {

	slope = value;
}

/* This routine evaluates the fit at a particular column number i,
   the zero-indexed X pixel coordinate AFTER overscan subtraction,
   i.e. the output column number.  (This is the only routine in this
   file where it matters whether the X pixel number is in the input
   or output file, because we're ignoring the intercept.)
   Note that we only include the slope of the fit, not the intercept;
   this is to ignore the overall bias level.
   Note also that subtract zero_column here, so it will be the origin.
*/

double DriftEval (double i) {

	double xi;		/* i - zero_column */

	xi = i - zero_column;

	return (xi * slope);
}

/* This routine returns the value of the drift averaged over an output
   line.  Since we're currently using just a linear fit, this will be
   the value of the drift at the column in the middle of the output image.
   The argument nx is the size of the first axis of the output image.
*/

double DriftMean (double nx) {

	double DriftEval (double);

	return (DriftEval (nx / 2.) );
}
