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
# include <math.h>
# include "hstio.h"
# include "stis.h"
# include "calstis1.h"
# include "hstcalerr.h"
# include "stisdq.h"

# define NELEM_SUMS  5		/* size of sums array */

static int rejectHighValues (SingleGroup *, int *, int *, float);
static int VMedianY (SingleGroup *, short, int, int *, double *, double *);
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

   Extreme outliers in the virtual overscan region will be replaced,
   but since the input image was opened read-only, the input file itself
   will not be modified.

   Phil Hodge, 1997 Nov 4:
	Function created.

   Phil Hodge, 1998 Sept 28:
	In VMedianY, delete k, and use nvals instead of k.

   Phil Hodge, 1998 Oct 16:
	Include sdqflags in the calling sequences of BlevDrift and VMedianY;
	use sdqflags in the latter, instead of requiring that there be no
	data quality flag set.

   Phil Hodge, 2004 July 28:
	Add function rejectHighValues, to filter out hot columns from the
	virtual overscan region.

   Phil Hodge, 2005 May 12:
	In function rejectHighValues, change the test for rejecting outliers
	from:
	    if (scratch[k] < 4.5 * mad) {
	to:
	    if (scratch[k] <= 4.5 * mad) {

   Phil Hodge, 2007 May 23:
	In BlevDrift, change:
	zerocol = (double) (biassect[0] + biassect[2]) / 2. - trimx1;
	to:
	zerocol = (double) (biassect[0] + biassect[1]) / 2. - trimx1;
*/

int BlevDrift (SingleGroup *in, short sdqflags, int *vx, int *vy,
		int trimx1, int *biassect, float blev_clip, int *driftcorr) {

/* arguments:
SingleGroup *in     io: image to be calibrated
                        (outliers in virtual overscan will be replaced)
short sdqflags       i: "serious" data quality flags
int vx[2], vy[2]     i: range of pixel numbers for virtual overscan region
int trimx1           i: width to trim off beginning of line
int biassect[2]      i: section to use for finding bias level
float blev_clip      i: criterion for clipping in virtual overscan region
int *driftcorr       o: true if correction can be applied
*/

	double *scratch;
	double value;		/* median of values in column */
	double zerocol;		/* zero point for fit */
	int midcol;		/* middle column of overscan region */
	int i;			/* X pixel in input image */
	char nodriftcorr[] =
		{"Warning  no correction for slope will be applied.\n"};

	*driftcorr = 0;		/* initial value */

	if (vx[1] <= vx[0] || vy[1] <= vy[0]) {
	    printf ("Warning  (blevcorr) No virtual overscan region; \\\n");
	    printf ("%s", nodriftcorr);
	    DriftSet (0.);
	    return (0);
	}

	if ((scratch = malloc ((vy[1]-vy[0]+1) * sizeof (double))) == NULL) {
	    printf ("ERROR    Out of memory in BlevDrift.\n");
	    return (OUT_OF_MEMORY);
	}

	/* Initialize for fitting.  The first argument is the middle of
	   the section used for determining the bias level, but shifted
	   to pixel number in the output, trimmed, image.
	*/
	zerocol = (double) (biassect[0] + biassect[1]) / 2. - trimx1;
	midcol = (vx[0] + vx[1]) / 2;
	DriftInit (zerocol, midcol);

	/* Replace extreme outliers in the virtual overscan region. */
	if (rejectHighValues (in, vx, vy, blev_clip))
	    return (OUT_OF_MEMORY);

	/* For each column, determine the value of the virtual overscan,
	   and accumulate sums.
	*/
	for (i = vx[0];  i <= vx[1];  i++) {
	    /* median of column */
	    if (VMedianY (in, sdqflags, i, vy, &value, scratch))
		continue;
	    DriftAccum (i, value);			/* increment sums */
	}

	/* Fit a curve to the values found. */
	if (DriftFit()) {
	    printf (
	"Warning  (blevcorr) Singular fit to virtual overscan; \\\n");
	    printf ("%s", nodriftcorr);
	    DriftSet (0.);
	} else {
	    *driftcorr = 1;
	}

	free (scratch);

	return (0);
}

/* Look for extreme outliers in the virtual overscan region, and
   replace them with neighboring values.  This is needed because some
   data have hot columns extending into the virtual overscan region,
   but VMedianY only looks at values within a single column, and if
   all the values in the column are high, they might not be rejected.

   The reason for actually replacing the values rather than just
   flagging them in the data quality array is that sdqflags might
   be set to zero, in which case the flagged data would be included
   by VMedianY even though the data were flagged.

   The criterion for rejecting an outlier is that the value be greater
   than (median + blev_clip * sigma), where 'median' is the median of all
   the values in the virtual overscan region (excluding the parallel
   overscan regions on the ends), and 'sigma' is the RMS of the values
   in the virtual overscan region (after rejecting 4.5*MAD outliers,
   where MAD is the median of the absolute value of deviations of the
   values from 'median').
*/

static int rejectHighValues (SingleGroup *in, int *vx, int *vy,
		float blev_clip) {

/* arguments:
SingleGroup *in   io: image to be calibrated (virtual overscan can be modified)
int vx[2], vy[2]   i: range of pixel numbers for virtual overscan region
float blev_clip    i: criterion for clipping in virtual overscan region
*/

	double *scratch;	/* for a copy of the virtual overscan data */
	int nscratch;		/* length of scratch array */
	double median;		/* median of the entire virtual overscan */
	double mad;		/* median of absolute values of deviations */
	float replacement;	/* replacement value for outliers */
	int i, j;		/* loop indices for column and row */
	int k;			/* index in scratch */
	int inplace = 1;	/* sort in-place */
	int i0, i1;		/* for finding replacement */
	int ii;
	double sumsq, n;	/* for finding sigma (rms) */
	double sigma;

	nscratch = (vx[1]-vx[0]+1) * (vy[1]-vy[0]+1);
	scratch = (double *)malloc (nscratch * sizeof (double));
	if (scratch == NULL) {
	    printf ("ERROR    Out of memory in BlevDrift.\n");
	    return (OUT_OF_MEMORY);
	}
	k = 0;
	for (j = vy[0];  j <= vy[1];  j++) {
	    for (i = vx[0];  i <= vx[1];  i++) {
		scratch[k] = Pix (in->sci.data, i, j);
		k++;
	    }
	}
	/* Compute the median and MAD. */
	median = MedianDouble (scratch, nscratch, inplace);
	for (k = 0;  k < nscratch;  k++)
	    scratch[k] = fabs (scratch[k] - median);
	mad = MedianDouble (scratch, k, inplace);

	/* Find the standard deviation.  At this point scratch contains
	   the absolute values of deviations from the median, which is
	   OK since we'll be squaring the elements of scratch.
	*/
	n = 0.;
	sumsq = 0.;
	for (k = 0;  k < nscratch;  k++) {
	    if (scratch[k] <= 4.5 * mad) {		/* approx 3 sigma */
		sumsq += scratch[k] * scratch[k];
		n += 1.;
	    }
	}
	sigma = sqrt (sumsq / n);

	/* Reject outliers.  */
	for (j = vy[0];  j <= vy[1];  j++) {
	    for (i = vx[0];  i <= vx[1];  i++) {
		if (Pix (in->sci.data, i, j) > median + blev_clip * sigma) {
		    i0 = i - 10;
		    i0 = i0 < vx[0] ? vx[0] : i0;
		    i1 = i0 + 20;
		    if (i1 > vx[1]) {
			i1 = vx[1];
			i0 = i1 - 20;
		    }
		    for (ii = i0, k = 0;  ii <= i1;  ii++, k++)
			scratch[k] = Pix (in->sci.data, ii, j);
		    replacement = MedianDouble (scratch, k, inplace);
		    Pix (in->sci.data, i, j) = replacement;
		}
	    }
	}

	free (scratch);
	return (0);
}

/* This routine takes the median of the SCI data values at x = i
   over y = vy[0] to vy[1] inclusive.  Note that the independent variable
   i is the column number in the input image.
   Using a median instead of a mean is to reject outliers.

   This function does not set status.  The function value will be zero
   if there is at least one data value not flagged as bad, but if there
   is no data, the function value will be NO_GOOD_DATA.
*/

static int VMedianY (SingleGroup *in, short sdqflags, int i, int *vy,
		double *median, double *scratch) {

/* arguments:
SingleGroup *in    i: image to be calibrated
short sdqflags     i: "serious" data quality flags
int i              i: X pixel location in input image
int vy[2]          i: range of pixel numbers in Y direction
double *median     o: median of values in column
double scratch[]  io: scratch space
*/

	int j;			/* loop index */
	int nvals;		/* number of good pixels */
	int inplace = 1;	/* OK to sort in-place */

	nvals = 0;
	for (j = vy[0];  j <= vy[1];  j++) {
	    if (!(sdqflags & DQPix (in->dq.data, i, j))) {
		scratch[nvals] = Pix (in->sci.data, i, j);
		nvals++;
	    }
	}

	if (nvals < 1)
	    return (NO_GOOD_DATA);

	*median = MedianDouble (scratch, nvals, inplace);

	return (0);
}

/* This routine resets the sums and coefficients to zero, and it copies
   zerocol to zero_column and midcol to middle_col.  Note that zerocol
   should be an X pixel number in the output, trimmed, image.
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

double DriftEval (int i) {

	double xi;		/* i - zero_column */

	xi = (double)i - zero_column;

	return (xi * slope);
}

/* This routine returns the value of the drift averaged over an output
   line.  Since we're currently using just a linear fit, this will be
   the value of the drift at the column in the middle of the output image.
   The argument nx is the size of the first axis of the output image.
*/

double DriftMean (int nx) {

	double DriftEval (int);

	return (DriftEval (nx / 2));
}
