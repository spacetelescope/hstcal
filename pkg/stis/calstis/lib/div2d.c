# include <math.h>
# include "c_iraf.h"
# include "hstio.h"

# include "stis.h"
# include "stisdq.h"	/* CALIBDEFECT = bad pixel value for divide by zero */
# include "err.h"	/* SIZE_MISMATCH */
# include "stisdef.h"

/* Divide the first SingleGroup triplet by the second, leaving the result
   in the first.

   (*a) /= (*b)

   The science data arrays are divided; the error arrays are combined;
   the data quality arrays are combined.
*/

int div2d (SingleGroup *a, SingleGroup *b) {

/* arguments:
SingleGroup *a   io: input data; output quotient
SingleGroup *b   i: second input data
*/

	int i, j;
	float a_sci, b_sci;	/* science data values from a and b */
	float a_expr, b_expr;	/* expressions involving errors in a and b */
	short dqa, dqb, dqab;	/* data quality for a, b, combined */

	if (a->sci.data.nx != b->sci.data.nx ||
	    a->sci.data.ny != b->sci.data.ny)
	    return (SIZE_MISMATCH);

	/* error array and science data */
	for (j = 0;  j < a->sci.data.ny;  j++) {
	    for (i = 0;  i < a->sci.data.nx;  i++) {

		a_sci = Pix (a->sci.data, i, j);
		b_sci = Pix (b->sci.data, i, j);

		if (b_sci == 0.) {
		    /* Flag divide by zero as lost data. */
		    dqa = CALIBDEFECT | DQPix (a->dq.data, i, j);
		    DQSetPix (a->dq.data, i, j, dqa);
		} else {
		    Pix (a->sci.data, i, j) = a_sci / b_sci;
		    a_expr = Pix (a->err.data, i, j) / b_sci;
		    b_expr = Pix (b->err.data, i, j)  * a_sci / b_sci / b_sci;
		    Pix (a->err.data, i, j) =
				sqrt (a_expr * a_expr + b_expr * b_expr);
		}
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
