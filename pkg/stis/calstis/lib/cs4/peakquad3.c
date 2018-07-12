# include <stdio.h>
# include "hstcalerr.h"

/* This routine finds the peak of the quadratic that passes through
   three uniformly spaced points.  The peak can be either a maximum
   or a minimum.
*/

int PeakQuad3 (double *y, double *x) {

/* arguments:
double *y      i: array of three points
double *x      o: location of axis relative to middle point,
		in units of the spacing of the independent variable
*/

	double denominator;

	denominator = y[0] - 2. * y[1] + y[2];
	if (denominator == 0.)
	    return (NO_GOOD_DATA);

	*x = (y[0] - y[2]) / (2. * denominator);

	return (0);
}
