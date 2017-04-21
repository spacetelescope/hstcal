# include <math.h>
# include "c_iraf.h"
# include "hstio.h"

# include "stis.h"
# include "err.h"	/* SIZE_MISMATCH */
# include "stisdef.h"

/* Subtract the second SingleGroup triplet from the first, leaving the
   result in the first.

   (*a) -= (*b)

   The science data arrays are subtracted; the error arrays are combined;
   the data quality arrays are combined.
*/

int sub2d (SingleGroup *a, SingleGroup *b) {

/* arguments:
SingleGroup *a   io: input data; output difference
SingleGroup *b   i: second input data
*/

	int i, j;
	float da, db;		/* errors for a and b */
	short dqa, dqb, dqab;	/* data quality for a, b, combined */

	if (a->sci.data.nx != b->sci.data.nx ||
	    a->sci.data.ny != b->sci.data.ny)
	    return (SIZE_MISMATCH);

	/* science data */
	for (j = 0;  j < a->sci.data.ny;  j++) {
	    for (i = 0;  i < a->sci.data.nx;  i++) {
		Pix (a->sci.data, i, j) =
			Pix (a->sci.data, i, j) - Pix (b->sci.data, i, j);
	    }
	}

	/* error array */
	for (j = 0;  j < a->err.data.ny;  j++) {
	    for (i = 0;  i < a->err.data.nx;  i++) {
		da = Pix (a->err.data, i, j);
		db = Pix (b->err.data, i, j);
		Pix (a->err.data, i, j) = sqrt (da * da + db * db);
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
