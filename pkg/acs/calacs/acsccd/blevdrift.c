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
# include "hstio.h"
# include "acs.h"
# include "acsinfo.h"
# include "err.h"
# include "acsdq.h"		/* for GOODPIXEL */

# define NELEM_SUMS  5		/* size of sums array */

static int VMedianY (SingleGroup *, int, int *, short, double *, double *);
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
*/

int BlevDrift (SingleGroup *in, int *vx, int *vy,
		int trimx1, int *biassect, int *driftcorr, short sdqflags) {

/* arguments:
SingleGroup *in      i: image to be calibrated
int vx[2], vy[2]     i: range of pixel numbers for virtual overscan region
int trimx1           i: width to trim off beginning of line
int biassect[2]      i: section to use for finding bias level
int *driftcorr       i: true if correction can be applied
*/

	extern int status;
	double *scratch;
	double value;		/* median of values in column */
	double zerocol;		/* zero point for fit */
	int midcol;		/* middle column of overscan region */
	int i;			/* X pixel in input image */
	char nodriftcorr[] =
		{"No correction for slope will be applied."};

	*driftcorr = 0;		/* initial value */

	if (vx[1] <= vx[0] || vy[1] <= vy[0]) {
	    trlmessage ("(blevcorr) No virtual overscan region;");
		trlmessage (nodriftcorr);
	    DriftSet (0.);
	    return (status);
	}

	if ((scratch = malloc ((vy[1]-vy[0]+1) * sizeof (double))) == NULL) {
	    trlwarn ("Out of memory in BlevDrift.");
	    return (status = OUT_OF_MEMORY);
	}

	/* Initialize for fitting.  The first argument is the middle of
	   the section used for determining the bias level, but shifted
	   to pixel number in the output, trimmed, image.
	*/
	
	zerocol = (double) (biassect[0] + biassect[1]) / 2.;
	midcol = (vx[0] + vx[1]) / 2;
	DriftInit (zerocol, midcol);

	/* For each column, determine the value of the virtual overscan,
	   and accumulate sums.
	*/
	for (i = vx[0];  i <= vx[1];  i++) {
	    if (VMedianY (in, i, vy, sdqflags, &value, scratch))	/* median of column */
		continue;
	    DriftAccum (i, value);			/* increment sums */
	}

	/* Fit a curve to the values found. */
	if (DriftFit()) {
	    trlwarn ("(blevcorr) Singular fit to virtual overscan; ");
	    trlwarn (nodriftcorr);
	    DriftSet (0.);
	} else {
	    *driftcorr = 1;
	}

	free (scratch);

	return (status);
}

/* This routine takes the median of the SCI data values at x = i
   over y = vy[0] to vy[1] inclusive.  Note that the independent variable
   i is the column number in the input image.
   Using a median instead of a mean is to reject outliers.

   This function does not set status.  The function value will be zero
   if there is at least one data value not flagged as bad, but if there
   is no data, the function value will be NO_GOOD_DATA.
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
	for (j = vy[0];  j <= vy[1];  j++) {
	    if (DQPix (in->dq.data, i, j) == GOODPIXEL || !(DQPix (in->dq.data, i, j) & sdqflags)) {
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
