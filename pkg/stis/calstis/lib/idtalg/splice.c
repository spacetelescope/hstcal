# include <stdio.h>
# include <stdlib.h>

# include "xtables.h"

# include "stis.h"
# include "hstcalerr.h"
# include "calstis6.h"
# include "idtalg.h"

static int FindIndex (double *, int, double);


/*
   splice --  Merges echelle orders.

   The input data lives in an array of RowContents elements. The output
   wavelength and flux arrays are allocated by this function but must
   be fred by the caller (their exact size is only known after the
   entire spectrum is spliced together). The NET column in the RowContents
   structure is used for filling up the flux array.

   This function also takes care of the conversion between counts/s
   to counts units.

   This function was written with the assumption that all spectral
   orders carry the same number of pixels, which is a reasonable
   assumption as long as calstis6std doesn't change its behavior.


   Revision history:
   ----------------
   13 Mar 00  -  Implemented (I.Busko)
*/

int Splice (RowContents **x1d, int nrows, double exptime, Spliced *merge) {
/* arguments
RowContents **x1d;		i: extracted 1-D data array
int nrows;			i: number of rows in x1d
double exptime;			i: exposure time
Spliced *merge;			o: output spliced spectrum
*/
	double **twave, **tflux, *wmerge, *fmerge;
	double midw;
	int i, i1, i2, i3, j, j1, j2, jstep, k, npt;

	double **Alloc2DArrayD (int, int);
	void Free2DArrayD (double **, int);

	npt = x1d[0]->npts;

	/* Check for reversed wavelenghts across orders. */

	if (x1d[nrows-1]->wave[0] < x1d[0]->wave[0]) {
	    j1 = nrows-1;
	    j2 = 0;
	    jstep = -1;
	} else {
	    j1 = 0;
	    j2 = nrows-1;
	    jstep = 1;
	}

	/* Alloc working memory. We need two buffers because the final
           size of the spliced spectrum will only be know after the
           entire splicing operation is complete.
        */

	twave = Alloc2DArrayD (npt, nrows);
	tflux = Alloc2DArrayD (npt, nrows);
	wmerge = (double *) malloc (npt * nrows * sizeof (double));
	fmerge = (double *) malloc (npt * nrows * sizeof (double));
	if (twave == NULL || tflux == NULL || wmerge == NULL || fmerge == NULL)
	    return (OUT_OF_MEMORY);

	/* Flip/copy spectral orders into working memory. Also translate
           from counts/s to counts.
        */

	for (j = j1, k = 0; k < nrows; j += jstep, k++) {
	    if (x1d[j]->wave[1] < x1d[j]->wave[0]) {
	        for (i1 = 0, i2 = npt-1; i1 < npt; i1++, i2--) {
	            twave[k][i1] = x1d[j]->wave[i2];
	            tflux[k][i1] = x1d[j]->net[i2] * exptime;
	        }
	    } else {
	        for (i1 = 0; i1 < npt; i1++) {
	            twave[k][i1] = x1d[j]->wave[i1];
	            tflux[k][i1] = x1d[j]->net[i1] * exptime;
	        }
	    }
	}

	/* Splice orders together. */

	i1 = 0;
	k  = 0;
	for (j = 0; j < nrows-1; j++) {

	    /* Find where wavelegth in current order gets larger than
               1st wavelength in next order.
            */
	    i2 = FindIndex (twave[j], npt, twave[j+1][0]);

	    /* Find where wavelegth in next order gets larger than
               last wavelength in current order.
            */
	    i3 = FindIndex (twave[j+1], npt, twave[j][npt-1]);

	    /* i2 and i3 define overlap region. Now compute midpoint. */

	    midw = (twave[j+1][i3] - twave[j][i2]) / 2. + twave[j][i2];

	    /* Re-compute indices to point to midpoint wavelength
               in both orders.
            */
	    i2 = FindIndex (twave[j],   npt, midw);
	    i3 = FindIndex (twave[j+1], npt, midw);

	    /* Ignore first 4 and last 3 pixels of each order,
               which may be bad.
            */
	    i2 = (i2 < (npt-3)) ? i2 : npt-4;
	    i3 = (i3 > 3)       ? i3 : 4;

	    /* Copy current order to output arrays. */
	    for (i = i1; i <= i2; i++, k++) {
	        wmerge[k] = twave[j][i];
	        fmerge[k] = tflux[j][i];
	    }

	    /* Set initial index for next iteration. */
	    i1 = i3;
	}

	/* Copy last row. */
	for (i = i1; i < npt; i++, k++) {
	    wmerge[k] = twave[nrows-1][i];
	    fmerge[k] = tflux[nrows-1][i];
	}

	/* Free first buffer. */

	Free2DArrayD (tflux, nrows);
	Free2DArrayD (twave, nrows);

	/* Alloc output array with exact size and copy data. */

	merge->npts = k;
	merge->wmerge = (double *) malloc (merge->npts * sizeof (double));
	merge->fmerge = (double *) malloc (merge->npts * sizeof (double));

	for (i = 0; i < k; i++) {
	    merge->wmerge[i] = wmerge[i];
	    merge->fmerge[i] = fmerge[i];
	}

	/* Free second buffer. */

	free (wmerge);
	free (fmerge);

	return (STIS_OK);
}


/*
   Finds the index of a given value in array. Array is assumed to be
   ordered in increasing order. Zero is returned if value is not in
   the range of values encompassed by the array elements.
*/

static int FindIndex (double *array, int n, double value) {

/* arguments:
double *array;		i: the array
int n;			i: the array size
double value;		i: the value to look for

returns: the index of the array element immediatley larger than 'value'
*/
	int i;

	for (i = 0; i < n; i++) {
	    if (array[i] > value)
	        return (i);
	}
	return (0);
}
