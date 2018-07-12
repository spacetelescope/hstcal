# include "hstio.h"

/* The constant k is added to the science data;
   the error and data quality arrays are not modified.
*/
int addk1d (SingleGroupLine *a, float k) {

/* arguments:
SingleGroupLine *a	io: input data; output sum
float k          	i: constant to be added to a
*/

	int i;
	int dimx;

	if (k == 0.)
	    return (0);

	/* Use this to remove pointer dereferencing for every pixel */
	dimx = a->sci.tot_nx;
    
	/* Only the science data are modified. */
	for (i = 0;  i < dimx;  i++)
	     a->sci.line[i] = k + a->sci.line[i];

	return (0);
}

int addkline (float *a, float k, int dimx) {

/* arguments:
float *a	       io: input data; output sum
float k			i: constant to be added to a
int dimx		i: size of array
*/

	int i;

	if (k == 0.)
	    return (0);

	/* Only the science data are modified. */
	for (i = 0;  i < dimx;  i++)
	     a[i] = k + a[i];

	return (0);
}

