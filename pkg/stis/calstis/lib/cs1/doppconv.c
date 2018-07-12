# include <stdio.h>
# include <stdlib.h>		/* calloc */
# include <math.h>		/* sqrt */
# include "hstio.h"

# include "stisdq.h"		/* for GOODPIXEL and DETECTORPROB */
# include "hstcalerr.h"

/* This routine convolves the science, error, and data quality arrays
   with a Doppler smearing function in the first axis direction, leaving
   the result in the input SingleGroup.

   Phil Hodge, 1998 Mar 13:
	Rewrite; change calling sequence.

   Phil Hodge, 1998 Apr 2:
	Use DETECTORPROB instead of zero to fill edges of DQ array.

   Phil Hodge, 1998 Aug 6:
	Add border to calling sequence, for correcting the flat field;
	remove the section for flagging left and right edges as off
	the detector.
*/

int DoppConv (SingleGroup *a, int border, float *ds, int nds, int d0) {

/* arguments:
SingleGroup *a     io: science data array
int border         i: magnitude of Doppler shift (high-res pixels)
float ds[]         i: Doppler smearing array
int nds            i: size of Doppler smearing array
int d0             i: index in ds corresponding to zero Doppler shift;
			d0 = (nds - 1) / 2
*/

	int nx;			/* length of first axis */
	float sum;		/* for convolving */
	float *x;		/* scratch space */
	short sum_dq;		/* for convolving data quality */
	short *xdq;		/* scratch space */
	int nscr;		/* size of scratch space */
	int i, j, k;
	int kmin, kmax;		/* range of indexes in ds */

	nx = a->sci.data.nx;
	nscr = nx + nds - 1;

	/* Find the range of non-zero elements in ds, to reduce the range
	   over which we need to do the convolution.
	*/
	kmin = nds - 1;		/* initial values */
	kmax = 0;
	for (k = 0;  k < nds;  k++) {
	    if (ds[k] > 0.) {		/* there will be no negative values */
		if (k < kmin)
		    kmin = k;
		if (k > kmax)
		    kmax = k;
	    }
	}

	x = calloc (nscr, sizeof(float));
	if (x == NULL) {
	    printf ("ERROR    (DoppConv) can't allocate memory.\n");
	    return (OUT_OF_MEMORY);
	}

	/* science data */
	for (j = 0;  j < a->sci.data.ny;  j++) {	/* each row */

	    /* Extract the row into the scratch buffer.  The border region
		included here is for flat fielding, to account for the loss
		of counts due to the flight software's effective subarray.
	    */
	    for (i = border;  i < nx-border;  i++)
		x[i+d0] = Pix (a->sci.data, i, j);

	    /* Convolve. */
	    for (i = 0;  i < nx;  i++) {
		sum = 0.;
		for (k = kmin;  k <= kmax;  k++)
		    sum += x[i+2*d0-k] * ds[k];
		Pix (a->sci.data, i, j) = sum;
	    }
	}

	/* error array */
	for (j = 0;  j < a->err.data.ny;  j++) {

	    for (i = 0;  i < nx;  i++)
		x[i+d0] = Pix (a->err.data, i, j);

	    for (i = 0;  i < nx;  i++) {
		/* sum of squares */
		sum = 0.;
		for (k = kmin;  k <= kmax;  k++)
		    sum += x[i+2*d0-k] * x[i+2*d0-k] * ds[k] * ds[k];
		Pix (a->err.data, i, j) = sqrt (sum);
	    }
	}

	free (x);

	/* data quality array */

	xdq = calloc (nscr, sizeof(short));
	if (xdq == NULL) {
	    printf ("ERROR    (DoppConv) can't allocate memory.\n");
	    return (OUT_OF_MEMORY);
	}

	for (j = 0;  j < a->dq.data.ny;  j++) {

	    for (i = 0;  i < nx;  i++)
		xdq[i+d0] = DQPix (a->dq.data, i, j);

	    /* OR the data quality values of all pixels touched
		due to the Doppler shift. */
	    for (i = 0;  i < nx;  i++) {
		sum_dq = GOODPIXEL;
		for (k = kmin;  k <= kmax;  k++)
		    sum_dq |= xdq[i+2*d0-k];
		DQSetPix (a->dq.data, i, j, sum_dq);
	    }
	}

	free (xdq);

	return (0);
}
