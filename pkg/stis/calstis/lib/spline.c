# include <stdlib.h>

# include "stis.h"

/* This file contains splint_nr and spline_nr.
   These routines were copied from Numerical Recipes (in Fortran,
   first edition, pp 88, 89) and were translated to C.
*/

static void spline_nr (double *, double *, double *, double *, int);

/* ya[i] is the function value at xa[i], for i = 1..n.  We will compute
   y[i], the estimate of the function at x[i], for i = 1..nelem, using
   a natural spline interpolation on the xa,ya.
   Note that for values of x[i] outside the range (x[0], x[n-1]),
   the values of y[i] will be extrapolated and can be quite wild.

   Phil Hodge, 1997 Oct 22:
	Functions created.
*/

int splint_nr (double *xa, double *ya, int n,
		double *x, double *y, int nelem) {

	double *ua;		/* scratch space */
	double *y2a;		/* array of second derivatives at xa */
	double a, b, h;
	int i;
	int k, klo, khi;

	ua = calloc (n, sizeof (double));
	y2a = calloc (n, sizeof (double));
	if (ua == NULL || y2a == NULL)
	    return (111);

	/* Compute the spline fit for xa vs ya. */
	spline_nr (xa, ya, y2a, ua, n);

	/* Evaluate the fit at each element of x. */
	for (i = 0;  i < nelem;  i++) {

	    /* Find the segment containing x[i]. */
	    klo = 0;
	    khi = n - 1;
	    while (khi - klo > 1) {
		k = (klo + khi) / 2;
		if (xa[k] > x[i])
		    khi = k;
		else
		    klo = k;
	    }

	    /* Evaluate. */

	    h = xa[khi] - xa[klo];
	    if (h == 0.)
		return (113);

	    a = (xa[khi] - x[i]) / h;
	    b = (x[i] - xa[klo]) / h;

	    y[i] = a * ya[klo] + b * ya[khi] +
		((a*a*a - a) * y2a[klo] + (b*b*b - b) * y2a[khi]) * h * h / 6.;
	}

	free (ua);
	free (y2a);

	return (0);
}

static void spline_nr (double *xa, double *ya, double *y2a, double *ua, int n) {

/* arguments:
double xa[n]     i: array of X coordinates
double ya[n]     i: Y values at corresponding xa
double y2a[n]    o: second derivatives needed by splint
double ua[n]     o: scratch space
int n            i: size of arrays
*/

	double sig, p;
	int i;

	/* "natural" spline */
	y2a[0] = 0.;
	y2a[n-1] = 0.;

	/* This is the decomposition loop of the tridiagonal algorithm.
	   y2a and ua are used for temporary storage of the decomposed factors.
	*/
	for (i = 1;  i < n-1;  i++) {
	    sig = (xa[i] - xa[i-1]) / (xa[i+1] - xa[i-1]);
	    p = sig * y2a[i-1] + 2.;
	    y2a[i] = (sig - 1.) / p;
	    ua[i] = (6. *
		((ya[i+1] - ya[i]) / (xa[i+1] - xa[i]) -
		 (ya[i] - ya[i-1]) / (xa[i] - xa[i-1])) /
		(xa[i+1] - xa[i-1]) - sig * ua[i-1]) / p;
	}

	/* This is the backsubstitution loop of the tridiagonal algorithm. */
	for (i = n-2;  i >= 0;  i--)
	    y2a[i] = y2a[i] * y2a[i+1] + ua[i];
}
