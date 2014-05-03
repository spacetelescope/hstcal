# include <stdio.h>
# include "stis.h"

/* This routine uses linear interpolation.
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

   Phil Hodge, 2000 Nov 8
	Convert from cubic to linear interpolation.
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

	double p0, p1;		/* differences between x and wl elements */
	double y;		/* interpolated value */
	int i;			/* index for which wl[i] is closest to x */

	if (n == 1 || x <= wl[0]) {

	    y = f[0];

	} else if (x >= wl[n-1]) {

	    y = f[n-1];

	} else {

	    /* Find the index i such that x is between wl[i] and wl[i+1];
		i can have any value from zero to n-2, inclusive.
	    */
	    i = *starti;

	    if (i < 0 || i > n-2) {		/* i is out of range? */

		i = BinarySearch (x, wl, n);

	    } else if (i > 0 && wl[i-1] <= x && x < wl[i]) {

		i--;			/* x is in the previous interval */

	    } else if (i < n - 2 && wl[i+1] <= x && x <= wl[i+2]) {

		i++;			/* x is in the next interval */

	    } else if (x < wl[i] || x > wl[i+1]) {

		i = BinarySearch (x, wl, n);

	    }

	    /* Update the starting index, in case we call this function
		again.
	    */
	    *starti = i;

	    /* linear interpolation */
	    p0 = (wl[i+1] - x) / (wl[i+1] - wl[i]);
	    p1 = 1. - p0;
	    y = f[i] * p0 + f[i+1] * p1;
	}

	return (y);
}

/* This function does a binary search and returns the index i such that
   x is between wl[i] and wl[i+1], except that i is restricted to the
   range from 0 to n-2 inclusive.
*/

static int BinarySearch (double x, double wl[], int n) {

	int low, high;		/* range of elements to consider */
	int k;			/* middle element between low and high */

	low = 0;
	high = n - 2;

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
