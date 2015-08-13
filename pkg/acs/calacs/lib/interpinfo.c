/* This file contains the following:
	InterpInfo
	InterpDQInfo
*/

# include <math.h>		/* sqrt */

/* This routine determines the appropriate array index i and weights
   p and q for linear interpolation.  If the array is called a, and ai
   is the independent variable (in units of the array index), then
   the interpolated value may be computed as:  p * a[i] + q * a[i+1].
*/

void InterpInfo (float ai, int npts, int *i, float *p, float *q) {

/* arguments:
float ai        i: independent variable in same units as i
int npts        i: size of array within which *i is an index
int *i          o: array index close to ai
float *p, *q    o: weights for linear interpolation
*/

	*i = (int) ai;
	*i = (*i < 0) ? 0 : *i;
	*i = (*i >= npts - 1) ? (npts - 2) : *i;
	*q = ai - *i;
	*p = 1.0F - *q;
}

/* This routine determines which array indexes i1 and i2 to use for
   ORing the data quality information.
*/

void InterpDQInfo (float ai, int npts, int *i1, int *i2, int *num_i) {

/* arguments:
float ai        i: independent variable in units of array index
int npts        i: size of array within which *i1 and *i2 are indexes
int *i1, *i2    o: array indexes close to ai
int *num_i      o: 1 or 2, i.e. use just *i1 or both *i1 and *i2
*/

	*i1 = (int) ai;
	*i2 = *i1 + 1;

	if (ai == (float)(*i1))
	    *num_i = 1;			/* ai is an integer, so just use i1 */
	else
	    *num_i = 2;

	if (*i1 <= 0) {
	    *i1 = 0;
	    *num_i = 1;			/* off the low end; just use i1 */
	}

	if (*i1 >= npts - 1) {
	    *i1 = npts - 1;
	    *num_i = 1;			/* off the high end; just use i1 */
	}
}
