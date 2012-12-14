# include <stdio.h>

/* DET2x2 is the determinant of a 2 x 2 matrix:
	| a  b |
	| c  d |
*/
# define DET2x2(a,b,c,d) ((a) * (d) - (b) * (c))

/* This routine uses Aitken's method for cubic interpolation.
   (See F. S. Acton, 1970, Numerical Methods that Work,
    or Abramowitz and Stegun, Handbook of Mathematical Functions, AMS 55.)
   The function value is the interpolated value.

   The independent variable array (wl) is assumed to be monotonically
   increasing.  If x is outside the range of wl values, the nearest
   element (i.e. the first or last) in the dependent variable array (f)
   will be returned.

   Before the first call to this routine, starti can be initialized to
   some value such as one.  This is the starting index for searching for
   the minimum difference between x and elements of wl; starti will be
   updated by this function.

   Phil Hodge, 1999 Mar 25:
	Rewrite section to search for new starti; add BinarySearch.

   Phil Hodge, 1999 Sept 3
	Add another check for index i out of range.
	Convert from quadratic to cubic interpolation.
*/

static int BinarySearch (double, double [], int);

double interp1d (double x, double wl[], double f[], int n, int *starti) {

/* arguments:
double x        i: the value at which the function is to be evaluated
double wl[]     i: the array of independent variable values
double f[]      i: the array of dependent variable values
int n           i: size of wl and f arrays
int *starti     io: begin here for finding nearest wl to x
*/

	double p0, p1, p2, p3;	/* differences between x and wl elements */
	double f0, f1, f2, f3;	/* nearest elements in f array */
	double wl0, wl1, wl2, wl3;	/* nearest elements in wl array */
	double temp1, temp2, temp3;
	double temp12, temp13;
	double y;		/* interpolated value */
	int i;			/* index for which wl[i] is closest to x */

	if (n == 1 || x <= wl[0]) {

	    y = f[0];

	} else if (x >= wl[n-1]) {

	    y = f[n-1];

	} else if (n == 2) {

	    p0 = (wl[1] - x) / (wl[1] - wl[0]);
	    p1 = 1. - x;
	    y = f[0] * p0 + f[1] * p1;

	} else if (n == 3) {

	    p0 = wl[0] - x;
	    p1 = wl[1] - x;
	    p2 = wl[2] - x;

	    y = ((f[0] * p1 - f[1] * p0) * p2 / (wl[1] - wl[0]) -
		 (f[0] * p2 - f[2] * p0) * p1 / (wl[2] - wl[0])) /
			(wl[2] - wl[1]);

	} else {

	    /* Find the index i such that x is between wl[i] and wl[i+1].
		We want i and i+1 to be the middle elements of a set of four,
		so we limit i to the range from one through n-3 inclusive.
	    */
	    i = *starti;

	    if (i < 1 || i > n-3 ||		/* i is out of range? */
		x < wl[i-1] || x > wl[i+2]) {	/* i is not close? */

		i = BinarySearch (x, wl, n);

	    } else {

		/* starting value is close, so check nearby points */
		if (wl[i-1] <= x && x < wl[i]) {
		    if (i > 1) {
			i--;
		    }
		} else if (wl[i+1] <= x && x <= wl[i+2]) {
		    if (i < n - 3) {
			i++;
		    }
		}
	    }

	    *starti = i;		/* update the starting location */

	    wl0 = wl[i-1];
	    wl1 = wl[i];
	    wl2 = wl[i+1];
	    wl3 = wl[i+2];

	    f0 = f[i-1];
	    f1 = f[i];
	    f2 = f[i+1];
	    f3 = f[i+2];

	    p0 = wl0 - x;
	    p1 = wl1 - x;
	    p2 = wl2 - x;
	    p3 = wl3 - x;

	    temp1 = DET2x2 (f0, p0,
	                    f1, p1) / (wl1 - wl0);

	    temp2 = DET2x2 (f0, p0,
	                    f2, p2) / (wl2 - wl0);

	    temp3 = DET2x2 (f0, p0,
	                    f3, p3) / (wl3 - wl0);

	    temp12 = DET2x2 (temp1, p1,
	                     temp2, p2) / (wl2 - wl1);

	    temp13 = DET2x2 (temp1, p1,
	                     temp3, p3) / (wl3 - wl1);

	    y = DET2x2 (temp12, p2,
	                temp13, p3) / (wl3 - wl2);
	}

	return (y);
}

/* This function does a binary search and returns the index i such that
   x is between wl[i] and wl[i+1], except that i is restricted to the
   range from 1 to n-3 inclusive.
*/

static int BinarySearch (double x, double wl[], int n) {

	int low, high;		/* range of elements to consider */
	int k;			/* middle element between low and high */

	low = 1;
	high = n - 3;

	while (high - low > 1) {

	    k = (low + high) / 2;
	    if (x < wl[k]) {
		high = k;
	    } else if (x > wl[k+1]) {
		low = k;
	    } else {
		high = k;
		low = k;
	    }
	}

	if (low < high) {
	    if (wl[high] <= x)
		return (high);
	    else
		return (low);
	} else {
	    return (low);
	}
}

