# include <math.h>
# include "hstio.h"
# include "hstcalerr.h"	/* SIZE_MISMATCH */

/* Subtract the second SingleGroupLine from the first, leaving the
   result in the first.

   (*a) -= (*b)

   The science data arrays are subtracted; the error arrays are combined;
   the data quality arrays are combined.
	
*/

int sub1d (SingleGroup *a, int line, SingleGroupLine *b) {

/* arguments:
SingleGroup *a   	io: input data; output difference
int line			 i: line of input data to subtract 1-d data from
SingleGroupLine *b   i: second input data
*/

	extern int status;

	int i;
	float da, db;		/* errors for a and b */
	short dqa, dqb, dqab;	/* data quality for a, b, combined */
    int dimx;

	if (a->sci.data.nx != b->sci.tot_nx)
	    return (status = SIZE_MISMATCH);

	/* science, error, and DQ data */
    dimx = a->sci.data.nx;
	for (i = 0;  i < dimx;  i++) {
    	/* science array */
		Pix(a->sci.data, i, line) =
			Pix (a->sci.data, i, line) - b->sci.line[i];

    	/* error array */
		da = Pix (a->err.data, i, line);
		db = b->err.line[i];
		Pix (a->err.data, i, line) = sqrt (da * da + db * db);

    	/* data quality */
		dqa = DQPix (a->dq.data, i, line);
		dqb = b->dq.line[i];
		dqab = dqa | dqb;
		DQSetPix (a->dq.data, i, line, dqab);
	}

	return (status);
}
