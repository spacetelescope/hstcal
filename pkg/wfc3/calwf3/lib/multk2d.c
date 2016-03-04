# include "hstio.h"

/* The science data and the error array values are multiplied by k;
   the data quality array is not modified.
*/

int multk2d (SingleGroup *a, float k) {

/* arguments:
SingleGroup *a   io: input data; output product
float k          i: multiply a by this constant
*/

	extern int status;

	int i, j;
	int dimx, dimy;

	if (k == 1.)
	    return (status);

	dimx = a->sci.data.nx;
	dimy = a->sci.data.ny;
    
	for (j = 0;  j < dimy;  j++) {
	    for (i = 0;  i < dimx;  i++) {
		Pix (a->sci.data, i, j) = k * Pix (a->sci.data, i, j);
		Pix (a->err.data, i, j) = k * Pix (a->err.data, i, j);
	    }
	}

	return (status);
}

