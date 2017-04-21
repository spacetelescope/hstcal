# include <math.h>
# include "hstio.h"

# include "acsdq.h"	/* CALIBDEFECT = bad pixel value for divide by zero */
# include "err.h"	/* SIZE_MISMATCH */

/* Divide the first SingleGroup triplet by the second, leaving the result
   in the first.

   (*a) /= (*b)

   The science data arrays are divided; the error arrays are combined;
   the data quality arrays are combined.
*/

int div1d (SingleGroup *a, int line, SingleGroupLine *b) {

/* arguments:
SingleGroup *a   	 io: input data; output quotient
int line		  	  i: row number from input data to divide
SingleGroupLine *b    i: second input data
*/

	extern int status;

	int i;
	float a_sci, b_sci;	/* science data values from a and b */
	float a_expr, b_expr;	/* expressions involving errors in a and b */
	short dqa, dqb, dqab;	/* data quality for a, b, combined */
    int sci_dimx, dq_dimx ;

	if (a->sci.data.nx != b->sci.tot_nx )
	    return (status = SIZE_MISMATCH);

    sci_dimx = a->sci.data.nx;
	/* error array and science data */
	for (i = 0;  i < sci_dimx;  i++) {

		a_sci = Pix (a->sci.data, i, line);
		b_sci = b->sci.line[i];

		if (b_sci == 0.) {
		    /* Flag divide by zero as lost data. */
		    dqa = CALIBDEFECT | DQPix (a->dq.data, i, line);
		    DQSetPix (a->dq.data, i, line, dqa);
		} else {
		    Pix (a->sci.data, i, line) = a_sci / b_sci;
		    a_expr = Pix (a->err.data, i, line) / b_sci;
		    b_expr = b->err.line[i]  * a_sci / b_sci / b_sci;
		    Pix (a->err.data, i, line) =
				sqrt (a_expr * a_expr + b_expr * b_expr);
		}
	}


	/* data quality */
    dq_dimx = a->dq.data.nx;
	for (i = 0;  i < dq_dimx;  i++) {
	
		dqa = DQPix (a->dq.data, i, line);
		dqb = b->dq.line[i];
		dqab = dqa | dqb;
		DQSetPix (a->dq.data, i, line, dqab);
	}

	return (status);
}
