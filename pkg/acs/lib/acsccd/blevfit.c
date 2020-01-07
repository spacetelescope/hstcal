/*  This file contains:
	BlevInit        zero sums for a new fit
	BlevAccum       increment sums for a new (line, bias) point
	BlevFit         compute the coefficients of the fit
	BlevSet         set coefficients of fit to a given value
	BlevEval        evaluate the fit for a given line number
    BlevResults     returns values for slope and intercept from fit
*/
# include <float.h>
#include "hstcal.h"
# include "acs.h"        /* for MsgText */

# define NELEM_SUMS  5		/* size of sums array */

/* This sets epsilon based on compilers 
	minimum defined float value which will
	be compared against a double. 
*/
# define EPSILON	FLT_MIN 

static double sums[NELEM_SUMS];
static double slope = 0., intercept = 0.;
static int middle_line;

/* This routine resets the sums and coefficients to zero, and it copies
   midline to middle_line.

   Phil Hodge, 1997 Nov 4:
	Change length of sums array from 6 to 5.
	
   Warren Hack, 1998 July 28:
   	Unchanged STIS code used for ACS...
*/

void BlevInit (int midline) {

/* argument:
int midline        i: middle line of image, will be subtracted from line number
*/

	int i;

	for (i = 0;  i < NELEM_SUMS;  i++)
	    sums[i] = 0.;
	slope = 0.;
	intercept = 0.;
	middle_line = midline;
}

/* This routine increments the sums for an unweighted linear fit of
   biaslevel as a function of line number j.  middle_line will be
   subtracted from j before incrementing the sums to reduce numerical
   problems.  The fit is linear.
*/

void BlevAccum (int j, double biaslevel) {

/* arguments:
int j              i: independent variable, line number in the image
double biaslevel   i: dependent variable, the bias level for this line
*/

	double xj;		/* = j - middle_line */

	xj = (double) (j - middle_line);

	sums[0]++;
	sums[1] += xj;
	sums[2] += biaslevel;
	sums[3] += (xj * biaslevel);
	sums[4] += (xj * xj);
}

/* This routine computes the coefficients of fit (slope and intercept).
   If there is no data, the function value will be -1, if the fit is
   singular, the function value will be -2; otherwise, the function will
   return zero.
*/

int BlevFit (void) {

	double d;
	double xmean, ymean;	/* mean values of line and biaslevel */

	if (sums[0] < 1.)
	    return (-1);

	d = sums[4] - sums[1] * sums[1] / sums[0];
	if (d < EPSILON && d > -EPSILON)
	    return (-2);

	xmean = sums[1] / sums[0];
	ymean = sums[2] / sums[0];

	slope = (sums[3] - xmean * ymean * sums[0]) / d;
	intercept = ymean - slope * xmean;

    sprintf (MsgText, "Computed a fit with slope of %g and intercept of %g", slope, intercept);
    trlmessage (MsgText);
    
	return (0);
}

void BlevResults(double *s, double *i) {

    *s = slope;
    *i = intercept;
}

/* This routine sets the coefficients of fit to specific values so that
   BlevEval will return "value" regardless of the line number.
*/

void BlevSet (double value) {

	slope = 0.;
	intercept = value;
}

/* This routine evaluates the fit at a particular line number j. */

double BlevEval (double j) {

	double xj;		/* = j - middle_line */

	xj = (j - middle_line);

	return (xj * slope + intercept);
}
