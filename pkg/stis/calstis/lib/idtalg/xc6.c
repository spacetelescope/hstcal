/* This file contains functions that handle complex arrays and
   Fast Fourier Transform operations:

   Float2Cmplx     - copies **float array into real part of cplx array
   Cmplx2Float     - copies real part of cplx array into **float array
   FFTConvolve     - convolve two cplx arrays with real data
   FFTConvolve2    - convolve two cplx arrays with real and cplx data


   Revision history:
   ----------------
   03 Mar 00  -  Borrowed from cs4 (I.Busko)
   03 Mar 00  -  Replaced SingleGroup structure by Image structure (IB)
   03 Mar 00  -  Add Phil's FFTConvolve function (IB)
   03 Mar 00  -  Replace malloc by calloc (IB)
   04 Mar 00  -  Add FFTShift (IB)
   23 Mar 00  -  Add Float2Cmplx, Cmplx2Float, FFTConvolve2, InitCmplxArray (IB)
   25 Jul 00  -  Remove functions that were moved to the calstis library (IB)
*/

# include <stdio.h>
# include <stdlib.h>
# include "hstio.h"

# include "stis.h"
# include "hstcalerr.h"
# include "idtalg.h"

/*static int CDebug (char *, CmplxArray *);*/ /* Not used */


/* Copies the array data into the real part of z->data, centering it.
   Eventual areas around edges are zeroed. Also sets the entire imaginary
   part of z->data to zero. The output complex array must be allocated by
   the caller, and must be of same size or larger than the input array.
*/

void Float2Cmplx (float **in, int inx, int iny, CmplxArray *z) {

/* arguments:
float **in        i: input data
int inx           i: X axis dimension of input array
int iny           i: Y axis dimension of input array
CmplxArray *z     o: complex data, copied from input
*/
	int i, j, i1, j1, ik, jk;

	i1 = (z->nx - inx) / 2;
	j1 = (z->ny - iny) / 2;

	for (j = 0;  j < z->ny;  j++) {
	    for (i = 0;  i < z->nx;  i++) {
		RPIX2D (z, i, j) = 0.0;
		IPIX2D (z, i, j) = 0.0;
	    }
	}

	for (j = j1, jk = 0;  jk < iny;  j++, jk++) {
	    for (i = i1, ik = 0;  ik < inx;  i++, ik++)
		RPIX2D (z, i, j) = in[jk][ik];
	}
}



/* Copies the real part of z->data into the output array. If the input
   array is larger than the output array, only its central subarray is
   copied. The output array must be allocated by the caller, and must
   be of same size or smaller than the input complex array.
*/

void Cmplx2Float (CmplxArray *z, float **out, int onx, int ony) {

/* arguments:
CmplxArray *z     i: input complex data
float **out       o: output data
int onx           i: X axis dimension of output array
int ony           i: Y axis dimension of output array
*/
	int i, j, i1, j1, ik, jk;

	i1 = (z->nx - onx) / 2;
	j1 = (z->ny - ony) / 2;

	for (j = j1, jk = 0;  jk < ony;  j++, jk++) {
	    for (i = i1, ik = 0;  ik < onx;  i++, ik++)
	        out[jk][ik] = RPIX2D (z, i, j);
	}
}



int FFTConvolve (CmplxArray *z1, CmplxArray *z2) {

/* arguments:
CmplxArray *z1    i: input data 1 copied into real part of array
CmplxArray *z2    i: input data 2 copied into real part of array
CmplxArray *z1    o: convolution of 1 and 2 in real part of array
CmplxArray *z2    o: FT of input z2
*/

	float temp;
	int i, j;		/* loop indexes */
	int nx, ny;		/* image size */
	int status;

	nx = z1->nx;
	ny = z1->ny;

	/* Take the forward Fourier transform of each input array, in-place. */
	if ((status = fft2d (z1)))
	    return (status);
	if ((status = fft2d (z2)))
	    return (status);

	/* Multiply the two complex arrays, leaving the product in z1. */
	for (j = 0;  j < ny;  j++) {
	    for (i = 0;  i < nx;  i++) {
		temp =	RPIX2D (z1, i, j) * RPIX2D (z2, i, j) -
			IPIX2D (z1, i, j) * IPIX2D (z2, i, j);
		IPIX2D (z1, i, j) =
			RPIX2D (z1, i, j) * IPIX2D (z2, i, j) +
			IPIX2D (z1, i, j) * RPIX2D (z2, i, j);
		RPIX2D (z1, i, j) = temp;
	    }
	}

	/* Shift to center. */
	FFTShift (z1);

	/* Take the inverse Fourier transform of z1, in-place. */
	if ((status = ifft2d (z1)))
	    return (status);

	return (0);
}


int FFTConvolve2 (CmplxArray *z1, CmplxArray *z2) {

/* arguments:
CmplxArray *z1    i: input data 1 copied into real part of array
CmplxArray *z2    i: FT of input data 2
CmplxArray *z1    o: convolution of 1 and 2 in real part of array
*/
	float temp;
	int i, j;		/* loop indexes */
	int nx, ny;		/* image size */
	int status;

	nx = z1->nx;
	ny = z1->ny;

	/* Take the forward Fourier transform of 1st input array, in-place. */

	if ((status = fft2d (z1)))
	    return (status);

	/* Multiply the two complex arrays, leaving the product in z1. */

	for (j = 0;  j < ny;  j++) {
	    for (i = 0;  i < nx;  i++) {
		temp =	RPIX2D (z1, i, j) * RPIX2D (z2, i, j) -
			IPIX2D (z1, i, j) * IPIX2D (z2, i, j);
		IPIX2D (z1, i, j) =
			RPIX2D (z1, i, j) * IPIX2D (z2, i, j) +
			IPIX2D (z1, i, j) * RPIX2D (z2, i, j);
		RPIX2D (z1, i, j) = temp;
	    }
	}

	/* Shift to center. */

	FFTShift (z1);

	/* Take the inverse Fourier transform of z1, in-place. */

	if ((status = ifft2d (z1)))
	    return (status);

	return (0);
}


/**************************************************************************/

/* Not used */
/*
static int CDebug (char *name, CmplxArray *z) {

	SingleGroup out;
	int i, j;

	initSingleGroup (&out);
	if (allocSingleGroup (&out, z->nx, z->ny) == -1)
	    return (OUT_OF_MEMORY);
	for (j = 0; j < z->ny; j++) {
	    for (i = 0; i < z->nx; i++) {
	        Pix (out.sci.data, i, j) = RPIX2D (z, i, j);
	        Pix (out.err.data, i, j) = IPIX2D (z, i, j);
	    }
	}

	printf ("Writing %s image.\n", name);

	if (putSingleGroup (name, 1, &out, 0))
	    return (OPEN_FAILED);
	freeSingleGroup (&out);

	printf ("Done writing.\n");

	return (STIS_OK);
}
*/
