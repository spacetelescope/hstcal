# include "hstio.h"
# include "hstcal.h"

/* The science data and the error array values are multiplied by k;
   the data quality array is not modified.

   M.D. De La Pena: 05 June 2018
   Created AvgSciVal from AvgSciValLine (multk1d.c) and put the new routine 
   here to complement the multk2d routine.
*/

int multk2d (SingleGroup *a, float k) {

/* arguments:
SingleGroup *a   io: input data; output product
float k          i: multiply a by this constant
*/

	extern int status;

	if (k == 1.)
	    return (status);

	/* science data */
    int dimx = a->sci.data.nx;
    int dimy = a->sci.data.ny;
    
    {unsigned int i, j;
    for (j = 0;  j < dimy;  j++) {
        for (i = 0;  i < dimx;  i++) {

            /* science array */
            Pix (a->sci.data, i, j) = k * Pix (a->sci.data, i, j);

            /* error array */
            Pix (a->err.data, i, j) = k * Pix (a->err.data, i, j);
        }
    }}

	return (status);
}

/* Compute the average of all the good pixels in the image, as
   well as a weight value.
*/ 
void AvgSciVal (SingleGroup *y, short sdqflags, double *mean, double *weight) {

    double sum  = 0.0;
    int numgood = 0;     /* number of good pixels */
    short flagval;       /* data quality flag value */

    int dimx = y->sci.data.nx;
    int dimy = y->sci.data.ny;

    {unsigned int i, j;
    for (j = 0; j < dimy; j++) {
        for (i = 0;  i < dimx;  i++) {
            flagval = DQPix (y->dq.data, i, j);

            /* no serious flag bit set */
            if ( ! (sdqflags & flagval) ) {
               sum += Pix (y->sci.data, i, j);
               numgood++;
            }
        }
    }}

    *mean   = 0.0;
    *weight = 0.0;
    if (numgood > 0) {
        *mean   = sum / (double) numgood;
        *weight = (double) numgood / (double)(dimx * dimy);
    } 
}
