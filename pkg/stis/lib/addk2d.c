# include "c_iraf.h"
# include "hstio.h"

# include "stis.h"
# include "stisdef.h"

/* The constant k is added to the science data;
   the error and data quality arrays are not modified.
*/

int addk2d (SingleGroup *a, float k) {

/* arguments:
SingleGroup *a   io: input data; output sum
float k          i: constant to be added to a
*/

	int i, j;

	if (k == 0.)
	    return (0);

	/* Only the science data are modified. */
	for (j = 0;  j < a->sci.data.ny;  j++)
	    for (i = 0;  i < a->sci.data.nx;  i++)
		Pix (a->sci.data, i, j) = k + Pix (a->sci.data, i, j);

	return (0);
}
