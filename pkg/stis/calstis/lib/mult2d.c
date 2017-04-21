# include <math.h>
# include "c_iraf.h"
# include "hstio.h"

# include "stis.h"
# include "err.h"	/* SIZE_MISMATCH */
# include "stisdef.h"

/* Multiply two SingleGroup triplets, leaving the result in the first.

   (*a) *= (*b)

   The science data arrays are multiplied together; the error arrays are
   combined; the data quality arrays are combined.
*/

int mult2d (SingleGroup *a, SingleGroup *b) {

/* arguments:
SingleGroup *a   io: input data; output product
SingleGroup *b   i: second input data
*/

	int i, j;
	float a_sci, b_sci;	/* science data values from a and b */
	float a_db;		/* a value * b error */
	float b_da;		/* b value * a error */
	short dqa, dqb, dqab;	/* data quality for a, b, combined */

	if (a->sci.data.nx != b->sci.data.nx ||
	    a->sci.data.ny != b->sci.data.ny)
	    return (SIZE_MISMATCH);

	/* error array and science data */
	for (j = 0;  j < a->err.data.ny;  j++) {
	    for (i = 0;  i < a->err.data.nx;  i++) {

		a_sci = Pix (a->sci.data, i, j);
		b_sci = Pix (b->sci.data, i, j);
		Pix (a->sci.data, i, j) = a_sci * b_sci;

		a_db = a_sci * Pix (b->err.data, i, j);
		b_da = b_sci * Pix (a->err.data, i, j);
		Pix (a->err.data, i, j) = sqrt (a_db * a_db + b_da * b_da);
	    }
	}

	/* data quality */
	for (j = 0;  j < a->dq.data.ny;  j++) {
	    for (i = 0;  i < a->dq.data.nx;  i++) {
		dqa = DQPix (a->dq.data, i, j);
		dqb = DQPix (b->dq.data, i, j);
		dqab = dqa | dqb;
		DQSetPix (a->dq.data, i, j, dqab);
	    }
	}

	return (0);
}
