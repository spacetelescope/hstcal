# include <math.h>
# include "hstio.h"
# include "err.h"	/* SIZE_MISMATCH */

/* Multiply SingleGroup triplet by SingleGroupLine triplet, leaving 
	the result in the first.

   (*a) *= (*b)

   The science data arrays are multiplied together; the error arrays are
   combined; the data quality arrays are combined.
*/

int mult1d (SingleGroup *a, int line, SingleGroupLine *b) {

/* arguments:
SingleGroup *a   	io: input data; output product
int line			 i: line in input data to operate on
SingleGroupLine *b   i: second input data
*/

	extern int status;

	int i;
	float a_sci, b_sci;	/* science data values from a and b */
	float a_db;		/* a value * b error */
	float b_da;		/* b value * a error */
	short dqa, dqb, dqab;	/* data quality for a, b, combined */
    int     err_dimx, dq_dimx;

	if (a->sci.data.nx != b->sci.tot_nx)
	    return (status = SIZE_MISMATCH);

	/* error array and science data */
    err_dimx = a->err.data.nx;
	for (i = 0;  i < err_dimx;  i++) {

		a_sci = Pix (a->sci.data, i, line);
		b_sci = b->sci.line[i];
		Pix (a->sci.data, i, line) = a_sci * b_sci;

		a_db = a_sci * b->err.line[i];
		b_da = b_sci * Pix (a->err.data, i, line);
		Pix (a->err.data, i, line) = sqrt (a_db * a_db + b_da * b_da);
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

int multlines (SingleGroupLine *a, SingleGroupLine *b) {

/* arguments:
SingleGroupLine *a   	io: input data; output product
SingleGroupLine *b   	 i: second input data
*/

	extern int status;

	int i;
	float a_sci, b_sci;	/* science data values from a and b */
	float a_db;		/* a value * b error */
	float b_da;		/* b value * a error */
	short dqa, dqb, dqab;	/* data quality for a, b, combined */
    int     err_dimx, dq_dimx;

	if (a->sci.tot_nx != b->sci.tot_nx)
	    return (status = SIZE_MISMATCH);

	/* error array and science data */
    err_dimx = a->err.tot_nx;
	for (i = 0;  i < err_dimx;  i++) {

		a_sci = a->sci.line[i];
		b_sci = b->sci.line[i];
		a->sci.line[i] = a_sci * b_sci;

		a_db = a_sci * b->err.line[i];
		b_da = b_sci * a->err.line[i];
		a->err.line[i] = sqrt (a_db * a_db + b_da * b_da);
	}

	/* data quality */
    dq_dimx = a->dq.tot_nx;
    for (i = 0;  i < dq_dimx;  i++) {
		dqa = a->dq.line[i];
		dqb = b->dq.line[i];
		dqab = dqa | dqb;
		a->dq.line[i] = dqab;
    }


	return (status);
}
