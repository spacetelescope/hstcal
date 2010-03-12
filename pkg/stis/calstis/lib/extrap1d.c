# include <stdio.h>
# include "stis.h"

static double extrapolate (double, double, double, double, double);

/* This is a replacement for interp1d that uses linear extrapolation if
   the argument falls outside the range of the tabulated values.


   Ivo Busko, 2002 Apr 20
	Adapted from interp1d.c
*/

double extrap1d (double x, double wl[], double f[], int n, int *starti) {

/* arguments:
double x        i: the value at which the function is to be evaluated
double wl[]     i: the array of independent variable values
double f[]      i: the array of dependent variable values
int n           i: size of wl and f arrays
int *starti     io: begin here for finding nearest wl to x
*/
	double wl1, wl2, f1, f2;
	double y;		/* return value */

        double interp1d (double, double *, double *, int, int *);

	if (n == 1 || x <= wl[0]) {
	    wl1 = wl[0];
	    wl2 = wl[1];
	    f1  = f[0];
	    f2  = f[1];
	    y = extrapolate (x, wl1, wl2, f1, f2);
	} else if (x >= wl[n-1]) {
	    wl1 = wl[n-2];
	    wl2 = wl[n-1];
	    f1  = f[n-2];
	    f2  = f[n-1];
	    y = extrapolate (x, wl1, wl2, f1, f2);
	} else {
	    y = interp1d (x, wl, f, n, starti);
	}

	return (y);
}

static double extrapolate (double x, double x1, double x2, 
                           double y1, double y2) {

	double slope, intercept;

	slope     = (y2 - y1) / (x2 - x1);
	intercept = y1 - slope * x1;

        return (slope * x + intercept);
}

