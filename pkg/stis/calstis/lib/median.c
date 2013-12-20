/* This file contains the following global symbols:
	MedianDouble
	MedianFloat
*/

# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include "stis.h"

static int CompareDouble (const void *vp, const void *vq);
static int CompareFloat (const void *vp, const void *vq);

/* These two routines (one for double precision, one for single) return
   the median of a floating-point array.

   If inplace = 0, the input array is copied to a scratch array, and the
   scratch array is sorted; otherwise, the input array is sorted in-place.
*/

double MedianDouble (double *v, int n, int inplace) {

/* arguments:
double v[n]    i: input array (io if inplace=1)
int n          i: size of array
int inplace    i: sort the input array in-place?
*/

	double *vt;		/* scratch for a copy of v */
	double median;
	int i;

	if (n < 1) {
	    printf ("Warning  (MedianDouble) No data.\n");
	    return (0.);
	} else if (n == 1) {
	    return (v[0]);
	} else if (n == 2) {
	    return ((v[0] + v[1]) / 2.);
	}

	if (inplace) {
	    vt = v;
	} else {
	    if ((vt = malloc (n * sizeof (double))) == NULL) {
		printf ("Warning  (MedianDouble) Out of memory.\n");
		return (v[0]);
	    }
	    memcpy (vt, v, n * sizeof (double));
	}

	qsort (vt, n, sizeof (double), CompareDouble);

	i = n / 2;
	if (n == i * 2)
	    median = (vt[i-1] + vt[i]) / 2.;
	else
	    median = vt[i];

	if (!inplace)
	    free (vt);

	return (median);
}

float MedianFloat (float *v, int n, int inplace) {

/* arguments:
double v[n]    i: input array (io if inplace=1)
int n          i: size of array
int inplace    i: sort the input array in-place?
*/

	float *vt;
	float median;
	int i;

	if (n < 1) {
	    printf ("Warning  (MedianFloat) No data.\n");
	    return (0.F);
	} else if (n == 1) {
	    return (v[0]);
	} else if (n == 2) {
	    return ((v[0] + v[1]) / 2.F);
	}

	if (inplace) {
	    vt = v;
	} else {
	    if ((vt = malloc (n * sizeof (float))) == NULL) {
		printf ("Warning  (MedianFloat) Out of memory.\n");
		return (v[0]);
	    }
	    memcpy (vt, v, n * sizeof (float));
	}

	qsort (vt, n, sizeof (float), CompareFloat);

	i = n / 2;
	if (n == i * 2)
	    median = (vt[i-1] + vt[i]) / 2.F;
	else
	    median = vt[i];

	if (!inplace)
	    free (vt);

	return (median);
}

/* These routines are used to sort the array. */

static int CompareDouble (const void *vp, const void *vq) {

	const double *p = vp;
	const double *q = vq;

	if (*p > *q)
	    return (1);
	else if (*p < *q)
	    return (-1);
	else
	    return (0);
}

static int CompareFloat (const void *vp, const void *vq) {

	const float *p = vp;
	const float *q = vq;

	if (*p > *q)
	    return (1);
	else if (*p < *q)
	    return (-1);
	else
	    return (0);
}
