# include "c_iraf.h"
# include "hstio.h"

# include "stis.h"
# include "stisdef.h"

/* The science data and the error array values are multiplied by k;
   the data quality array is not modified.
*/

int multk2d (SingleGroup *a, float k) {

/* arguments:
SingleGroup *a   io: input data; output product
float k          i: multiply a by this constant
*/

	int i, j;

	if (k == 1.)
	    return (0);

	/* science data */
	for (j = 0;  j < a->sci.data.ny;  j++)
	    for (i = 0;  i < a->sci.data.nx;  i++)
		Pix (a->sci.data, i, j) = k * Pix (a->sci.data, i, j);

	/* error array */
	for (j = 0;  j < a->err.data.ny;  j++)
	    for (i = 0;  i < a->err.data.nx;  i++)
		Pix (a->err.data, i, j) = k * Pix (a->err.data, i, j);

	return (0);
}
