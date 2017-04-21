/* This file contains functions that handle complex arrays and
   Fast Fourier Transform operations:

   FFTShift        - can be used to re-center cross correlation result
   AllocCmplxArray - complex array memory management
   FreeCmplxArray
   InitCmplxArray
   fft2d           - forward 2-D FFT
   ifft2d          - inverse 2-D FFT
   CpyCmplx        - copy a cmplx array


   Revision history:
   ----------------
   21 Jul 2000  -  Extracted from cx4/ and cs6/idtalg/

*/

# include <stdio.h>
# include <stdlib.h>

# include "stis.h"
# include "err.h"


/*
   To be used before transforming FT back to data space, to shift image
   to center of field. Works for arrays with even dimensions.
*/

void FFTShift (CmplxArray *z) {

/* arguments:
CmplxArray *z       io: FT array to be shifted
*/
	int i, j;

        for (j = 1; j < z->ny; j+=2) {
            for (i = 0; i < z->nx-1; i+=2) {
	        RPIX2D (z, i, j) = - RPIX2D (z, i, j);
	        IPIX2D (z, i, j) = - IPIX2D (z, i, j);
	    }
        }
        for (j = 0; j < z->ny-1; j+=2) {
            for (i = 1; i < z->nx; i+= 2) {
	        RPIX2D (z, i, j) = - RPIX2D (z, i, j);
	        IPIX2D (z, i, j) = - IPIX2D (z, i, j);
	    }
        }

}


/* AllocCmplxArray allocates memory for an array of complex numbers.
   Memory is also allocated for the work space used by the NCAR FFT
   routines written by Paul Swarztrauber.
*/

int AllocCmplxArray (CmplxArray *z, int nx, int ny) {

	if (z->allocated) {
	    printf ("Error:  complex array already allocated.\n");
	    return (ERROR_RETURN);
	}

	if (ny < 1)
	    ny = 1;

	z->data = calloc (nx * ny * 2, sizeof(float));
	z->workx = calloc ((15 + 4 * nx), sizeof (float));
	z->worky = calloc ((15 + 4 * ny), sizeof (float));
	z->nx = nx;
	z->ny = ny;

	if (z->data == NULL || z->workx == NULL || z->worky == NULL)
	    return (OUT_OF_MEMORY);

	z->allocated = 1;

	/* Compute trig functions for use by the NCAR FFT routines. */
	CFFTI (&nx, z->workx);
	if (ny > 1)
	    CFFTI (&ny, z->worky);

	return (0);
}

/* This routine frees memory allocated by AllocCmplxArray. */

void FreeCmplxArray (CmplxArray *z) {

	if (z->allocated) {
	    free (z->data);
	    free (z->workx);
	    free (z->worky);
	    z->nx = 0;
	    z->ny = 0;
	    z->allocated = 0;
	}
}


/* Initializes a CmplxArray structure. */

void InitCmplxArray (CmplxArray *z) {

	z->allocated = 0;
	z->data = NULL;
	z->workx = NULL;
	z->worky = NULL;
	z->nx = 0;
	z->ny = 0;
}


/* Forward Fourier transform of 2-D data.
   On input, z contains an array of complex data.
   On output, z will contain the fourier transform of the input data.
*/

int fft2d (CmplxArray *z) {

/* argument:
CmplxArray *z       io: complex array
*/

	CmplxArray scr;		/* scratch space for FT of columns */
	int i, j;		/* loop indexes */
	int nx, ny;		/* image size */
	int status;

	nx = z->nx;
	ny = z->ny;

	/* Allocate a 1-D array for a column.  Note that scr.nx = z->ny. */
	InitCmplxArray (&scr);
	if ((status = AllocCmplxArray (&scr, ny, 1)))
	    return (status);

	for (j = 0;  j < ny;  j++)	/* transform each line */
	    CFFTF (&nx, &(RPIX2D(z,0,j)), z->workx);

	for (i = 0;  i < nx;  i++) {	/* transform each column */

	    for (j = 0;  j < ny;  j++) {
		RPIX1D (&scr, j) = RPIX2D (z, i, j);
		IPIX1D (&scr, j) = IPIX2D (z, i, j);
	    }
	    CFFTF (&scr.nx, scr.data, scr.workx);
	    for (j = 0;  j < ny;  j++) {
		RPIX2D (z, i, j) = RPIX1D (&scr, j);
		IPIX2D (z, i, j) = IPIX1D (&scr, j);
	    }
	}

	FreeCmplxArray (&scr);
	return (0);
}

/* Inverse Fourier transform of 2-D data.
   The normalization (dividing by nx * ny) is done in this routine.
   Note that z is both input and output.
*/

int ifft2d (CmplxArray *z) {

/* argument:
CmplxArray *z       io: complex array
*/

	CmplxArray scr;		/* scratch space for FT of columns */
	int i, j;		/* loop indexes */
	int nx, ny;		/* image size */
	int status;

	nx = z->nx;
	ny = z->ny;

	/* Allocate a 1-D array for a column. */
	InitCmplxArray (&scr);
	if ((status = AllocCmplxArray (&scr, ny, 1)))
	    return (status);

	for (j = 0;  j < ny;  j++)	/* transform each line */
	    CFFTB (&nx, &(RPIX2D(z,0,j)), z->workx);

	for (i = 0;  i < nx;  i++) {	/* transform each column */

	    for (j = 0;  j < ny;  j++) {
		RPIX1D (&scr, j) = RPIX2D (z, i, j);
		IPIX1D (&scr, j) = IPIX2D (z, i, j);
	    }
	    CFFTB (&scr.nx, scr.data, scr.workx);
	    for (j = 0;  j < ny;  j++) {
		RPIX2D (z, i, j) = RPIX1D (&scr, j) / (nx * ny);
		IPIX2D (z, i, j) = IPIX1D (&scr, j) / (nx * ny);
	    }
	}

	FreeCmplxArray (&scr);
	return (0);
}


void CpyCmplx (CmplxArray *zin, CmplxArray *zout) {

	int i,j;

	for (j = 0;  j < zout->ny;  j++) {
	    for (i = 0; i < zout->nx;  i++) {
		RPIX2D (zout, i, j) = RPIX2D (zin, i, j);
		IPIX2D (zout, i, j) = IPIX2D (zin, i, j);
	    }
	}
}
