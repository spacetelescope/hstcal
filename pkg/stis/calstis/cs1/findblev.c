# include <stdio.h>
# include <stdlib.h>		/* for calloc and qsort */

# include "hstio.h"
# include "stis.h"
# include "err.h"	/* for NO_GOOD_DATA */
# include "stisdq.h"

/* This routine determines the bias level for one line of an image
   by taking the median of the values in the overscan region.

   Phil Hodge, 1998 Oct 16:
	Use sdqflags instead of requiring that there be no dq flag set.
*/

int FindBlev (SingleGroup *x, short sdqflags, int j, int *biassect,
		double *biaslevel, int *npix) {

/* arguments:
SingleGroup *x      i: needed for science data and data quality
short sdqflags      i: "serious" data quality flags
int j               i: line of x data to use to get overscan
int biassect[2]     i: beginning and end of region to use for overscan
double *biaslevel   o: median bias level for current (j) line
int *npix           o: number of pixels used to compute bias level
*/

	double *over;	/* values extracted from overscan region */
	int nvals;	/* number of good pixels extracted from overscan */
	int i;
	int inplace = 1;	/* sort the array in-place */

	/* Allocate space for the overscan, and copy out good data. */

	over = calloc (biassect[1]-biassect[0]+1, sizeof (double));

	nvals = 0;
	for (i = biassect[0];  i <= biassect[1];  i++) {
	    if (!(sdqflags & DQPix (x->dq.data, i, j))) {
		over[nvals] = Pix (x->sci.data, i, j);
		nvals++;
	    }
	}

	*npix = nvals;

	if (nvals < 1) {
	    free (over);
	    return (NO_GOOD_DATA);
	}

	/* Find the median. */
	*biaslevel = MedianDouble (over, nvals, inplace);

	free (over);

	return (0);
}
