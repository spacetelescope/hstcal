# include <stdio.h>
# include <stdlib.h>
# include "stis.h"
# include "err.h"
# include "idtalg.h"


/*
   Memory management functions.



   Revision history:
   ----------------
   22 Feb 00  -  Implemented (I.Busko)

*/

float **Alloc2DArrayF (int x, int y) {

	float **array;
	int i;

	array = (float **) malloc (y * sizeof (float *));
	if (array == NULL) {
            printf ("Not enough memory to allocate array.\n");
	    return (NULL);
	}
	for (i = 0; i < y; i++) {
	    array[i] = (float *) calloc (x, sizeof (float));
	    if (array[i] == NULL) {
                printf ("Not enough memory to allocate array.\n");
	        return (NULL);
	    }
	}
	return (array);
}


void Free2DArrayF (float **array, int y) {

	int i;

	for (i = 0; i < y; free (array[i++]));
	free (array);
}


double **Alloc2DArrayD (int x, int y) {

	double **array;
	int i;

	array = (double **) malloc (y * sizeof (double *));
	if (array == NULL) {
            printf ("Not enough memory to allocate array.\n");
	    return (NULL);
	}
	for (i = 0; i < y; i++) {
	    array[i] = (double *) calloc (x, sizeof (double));
	    if (array[i] == NULL) {
                printf ("Not enough memory to allocate array.\n");
	        return (NULL);
	    }
	}
	return (array);
}


void Free2DArrayD (double **array, int y) {

	int i;

	for (i = 0; i < y; free (array[i++]));
	free (array);
}


void InitImage (Image *image) {

	image->im    = NULL;
	image->pix   = NULL;
	image->naxis = 0;
	image->nx    = 0;
	image->ny    = 0;
	image->npix  = 0;
	image->imio  = 0;
}

int Alloc1DImage (Image *image, int n) {

	image->im    = NULL;
	image->naxis = 1;
	image->nx    = n;
	image->ny    = 0;
	image->npix  = n;
	image->imio  = 0;
	image->pix = (float *) calloc (image->npix, sizeof (float));
	if (image->pix == NULL)
	    return (OUT_OF_MEMORY);
	return (STIS_OK);
}


int Alloc2DImage (Image *image, int nx, int ny) {

	image->im    = NULL;
	image->naxis = 2;
	image->nx    = nx;
	image->ny    = ny;
	image->npix  = nx * ny;
	image->imio  = 0;
	image->pix = (float *) calloc (image->npix, sizeof (float));
	if (image->pix == NULL)
	    return (OUT_OF_MEMORY);
	return (STIS_OK);
}


void FreeImage (Image *image) {

	if (image->im != NULL) {
	    c_imtclose (image->im);
	    image->im = NULL;
	}
	image->naxis = 0;
	image->nx    = 0;
	image->ny    = 0;
	image->npix  = 0;
	if (image->pix != NULL && !(image->imio))
	    free (image->pix);
	image->pix = NULL;
	image->imio = 0;
}


void FreeScatter (ScatterFunctions *scf) {

	int i;

	for (i = 0; i < scf->nsc; free (scf->scfunc[i++].values));
	free (scf->scfunc);
	for (i = 0; i < scf->nrp; free (scf->rpfunc[i++].wavelengths));
	for (i = 0; i < scf->nrp; free (scf->rpfunc[i++].values));
	free (scf->rpfunc);
	if (scf->spsf != NULL)
	    free (scf->spsf);
	if (scf->xdisp != NULL)
	    free (scf->xdisp);
}
