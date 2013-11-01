# include <math.h>

# include "hstio.h"	/* Defines HST I/O functions */
# include "wf3.h"
# include "wf3info.h"

/* MATH: Contains general utilities for performing arithmetic
** on "NicmosGroup" two-dimensional data groups. Errors and DQ are
** updated.  SAMP and TIME arrays are unchanged.
**
** Revision history:
** H.Bushouse	Oct. 2000	Initial CALNICA to CALWF3 port.
** H.Bushouse	10-Apr-2002	Added *_noref routines to enable skipping
**				of IR array reference pixels.
** H.Bushouse	01-Aug-2006	Updated adiv and adiv_noref routines to
**				avoid division by zero.
*/

void aadd (SingleNicmosGroup *a, SingleNicmosGroup *b) {

/* Add two SingleNicmosGroups, leaving the result in the first.

   (*a) += (*b)

   The science data arrays are added together; the error arrays are
   combined; the data quality arrays are "or'ed".
*/

	int i, j;	/* array indexes */
	float aerr;	/* a error */
	float berr;	/* b error */

	for (j=0; j < a->sci.data.ny; j++) {
	     for (i=0; i < a->sci.data.nx; i++) {

		  /* error data */
		  aerr = Pix(a->err.data,i,j);
		  berr = Pix(b->err.data,i,j);
		  Pix(a->err.data,i,j) = sqrt(aerr*aerr + berr*berr);

		  /* science data */
		  Pix(a->sci.data,i,j) += Pix(b->sci.data,i,j);

		  /* data quality */
		  DQSetPix(a->dq.data,i,j,
		  DQPix(a->dq.data,i,j) | DQPix(b->dq.data,i,j) );

	     }
	}
}

void asub (SingleNicmosGroup *a, SingleNicmosGroup *b) {

/* Subtract two SingleNicmosGroups, leaving the result in the first.

   (*a) -= (*b)

   The science data arrays are subtracted; the error arrays are
   combined; the data quality arrays are "or'ed".
*/

	int i, j;	/* array indexes */
	float aerr;	/* a error */
	float berr;	/* b error */

	for (j=0; j < a->sci.data.ny; j++) {
	     for (i=0; i < a->sci.data.nx; i++) {

		  /* error data */
		  aerr = Pix(a->err.data,i,j);
		  berr = Pix(b->err.data,i,j);
		  Pix(a->err.data,i,j) = sqrt(aerr*aerr + berr*berr);

		  /* science data */
		  Pix(a->sci.data,i,j) -= Pix(b->sci.data,i,j);

		  /* data quality */
		  DQSetPix(a->dq.data,i,j,
		  DQPix(a->dq.data,i,j) | DQPix(b->dq.data,i,j) );

	     }
	}
}

void amul (SingleNicmosGroup *a, SingleNicmosGroup *b) {

/* Multiply two SingleNicmosGroups, leaving the result in the first.

   (*a) *= (*b)

   The science data arrays are multiplied together; the error arrays are
   combined; the data quality arrays are "or'ed".
*/

	int i, j;	/* array indexes */
	float a_db;	/* a value * b error */
	float b_da;	/* b value * a error */

	for (j=0; j < a->sci.data.ny; j++) {
	     for (i=0; i < a->sci.data.nx; i++) {

		  /* error data */
		  a_db = Pix(a->sci.data,i,j) * Pix(b->err.data,i,j);
		  b_da = Pix(b->sci.data,i,j) * Pix(a->err.data,i,j);
		  Pix(a->err.data,i,j) = sqrt (a_db*a_db + b_da*b_da);

		  /* science data */
		  Pix(a->sci.data,i,j) *= Pix(b->sci.data,i,j);

		  /* data quality */
		  DQSetPix(a->dq.data,i,j,
		  DQPix(a->dq.data,i,j) | DQPix(b->dq.data,i,j) );

	     }
	}
}

void adiv (SingleNicmosGroup *a, SingleNicmosGroup *b) {

/* Divide two SingleNicmosGroups, leaving the result in the first.

   (*a) /= (*b)

   The science data arrays are divided; the error arrays are
   combined; the data quality arrays are "or'ed".
*/

	int i, j;	/* array indexes */
	float asci;	/* a science */
	float bsci;	/* b science */
	float bsci2;	/* b science squared */
	float aerr;	/* a error */
	float berr;	/* b error */

	for (j=0; j < a->sci.data.ny; j++) {
	     for (i=0; i < a->sci.data.nx; i++) {

		  /* error and science data */
		  asci = Pix(a->sci.data,i,j);
		  bsci = Pix(b->sci.data,i,j);
		  bsci2 = bsci*bsci;
		  aerr = Pix(a->err.data,i,j);
		  berr = Pix(b->err.data,i,j);

		  if (bsci != 0) {
		      Pix(a->sci.data,i,j) = asci / bsci;
		      Pix(a->err.data,i,j) =
		      sqrt(aerr*aerr/bsci2 + asci*asci*berr*berr/(bsci2*bsci2));
		  } else {
		      Pix(a->sci.data,i,j) = 0.0;
		      Pix(a->err.data,i,j) = 0.0;
		  }

		  /* data quality */
		  DQSetPix(a->dq.data,i,j,
			DQPix(a->dq.data,i,j) | DQPix(b->dq.data,i,j) );

	     }
	}
}

void asubk (SingleNicmosGroup *a, float b) {

/* Subtract in-place a constant from a SingleNicmosGroup.

   (*a) -= b

   The constant is subtracted from the science data array only;
   the error and data quality arrays are unchanged.
*/

	int i, j;	/* array indexes */

	for (j=0; j < a->sci.data.ny; j++) {
	     for (i=0; i < a->sci.data.nx; i++) {

		  /* science data */
		  Pix(a->sci.data,i,j) -= b;

	     }
	}
}

void amulk (SingleNicmosGroup *a, float b) {

/* Multiply in-place a SingleNicmosGroup by a constant.

   (*a) *= b

   The science and error data arrays are multiplied by the constant;
   the data quality array is unchanged.
*/

	int i, j;	/* array indexes */

	for (j=0; j < a->sci.data.ny; j++) {
	     for (i=0; i < a->sci.data.nx; i++) {

		  /* science data */
		  Pix(a->sci.data,i,j) *= b;

		  /* error data */
		  Pix(a->err.data,i,j) *= b;

	     }
	}
}

void aor (SingleNicmosGroup *a, SingleNicmosGroup *b) {

/* Take the logical OR of two SingleNicmosGroup DQ arrays, leaving the result 
** in the first.

   (*a.dq) = (*a.dq) | (*b.dq)

   The science and error data arrays are unchanged.
*/

	int i, j;	/* array indexes */

	/* data quality */
	for (j=0; j < a->dq.data.ny; j++) {
	     for (i=0; i < a->dq.data.nx; i++) {

		  DQSetPix (a->dq.data,i,j,
			DQPix (a->dq.data,i,j) | DQPix (b->dq.data,i,j) );
	     }
	}
}

void aadd_noref (WF3Info *wf3, SingleNicmosGroup *a, SingleNicmosGroup *b) {

/* Add two SingleNicmosGroups, leaving the result in the first.

   (*a) += (*b)

   The science data arrays are added together; the error arrays are
   combined; the data quality arrays are "or'ed".
*/

	int i, j;	/* array indexes */
	int ibeg, iend;	/* loop limits */
	int jbeg, jend;	/* loop limits */
	float aerr;	/* a error */
	float berr;	/* b error */

	ibeg = wf3->trimx[0]; iend = a->sci.data.nx - wf3->trimx[1];
	jbeg = wf3->trimy[0]; jend = a->sci.data.ny - wf3->trimy[1];

	for (j = jbeg; j < jend; j++) {
	     for (i = ibeg; i < iend; i++) {

		  /* error data */
		  aerr = Pix(a->err.data,i,j);
		  berr = Pix(b->err.data,i,j);
		  Pix(a->err.data,i,j) = sqrt(aerr*aerr + berr*berr);

		  /* science data */
		  Pix(a->sci.data,i,j) += Pix(b->sci.data,i,j);

		  /* data quality */
		  DQSetPix(a->dq.data,i,j,
		  DQPix(a->dq.data,i,j) | DQPix(b->dq.data,i,j) );

	     }
	}
}

void asub_noref (WF3Info *wf3, SingleNicmosGroup *a, SingleNicmosGroup *b) {

/* Subtract two SingleNicmosGroups, leaving the result in the first.

   (*a) -= (*b)

   The science data arrays are subtracted; the error arrays are
   combined; the data quality arrays are "or'ed".
*/

	int i, j;	/* array indexes */
	int ibeg, iend;	/* loop limits */
	int jbeg, jend;	/* loop limits */
	float aerr;	/* a error */
	float berr;	/* b error */

	ibeg = wf3->trimx[0]; iend = a->sci.data.nx - wf3->trimx[1];
	jbeg = wf3->trimy[0]; jend = a->sci.data.ny - wf3->trimy[1];

	for (j = jbeg; j < jend; j++) {
	     for (i = ibeg; i < iend; i++) {

		  /* error data */
		  aerr = Pix(a->err.data,i,j);
		  berr = Pix(b->err.data,i,j);
		  Pix(a->err.data,i,j) = sqrt(aerr*aerr + berr*berr);

		  /* science data */
		  Pix(a->sci.data,i,j) -= Pix(b->sci.data,i,j);

		  /* data quality */
		  DQSetPix(a->dq.data,i,j,
		  DQPix(a->dq.data,i,j) | DQPix(b->dq.data,i,j) );

	     }
	}
}

void amul_noref (WF3Info *wf3, SingleNicmosGroup *a, SingleNicmosGroup *b) {

/* Multiply two SingleNicmosGroups, leaving the result in the first.

   (*a) *= (*b)

   The science data arrays are multiplied together; the error arrays are
   combined; the data quality arrays are "or'ed".
*/

	int i, j;	/* array indexes */
	int ibeg, iend;	/* loop limits */
	int jbeg, jend;	/* loop limits */
	float a_db;	/* a value * b error */
	float b_da;	/* b value * a error */

	ibeg = wf3->trimx[0]; iend = a->sci.data.nx - wf3->trimx[1];
	jbeg = wf3->trimy[0]; jend = a->sci.data.ny - wf3->trimy[1];

	for (j = jbeg; j < jend; j++) {
	     for (i = ibeg; i < iend; i++) {

		  /* error data */
		  a_db = Pix(a->sci.data,i,j) * Pix(b->err.data,i,j);
		  b_da = Pix(b->sci.data,i,j) * Pix(a->err.data,i,j);
		  Pix(a->err.data,i,j) = sqrt (a_db*a_db + b_da*b_da);

		  /* science data */
		  Pix(a->sci.data,i,j) *= Pix(b->sci.data,i,j);

		  /* data quality */
		  DQSetPix(a->dq.data,i,j,
		  DQPix(a->dq.data,i,j) | DQPix(b->dq.data,i,j) );

	     }
	}
}

void adiv_noref (WF3Info *wf3, SingleNicmosGroup *a, SingleNicmosGroup *b) {

/* Divide two SingleNicmosGroups, leaving the result in the first.

   (*a) /= (*b)

   The science data arrays are divided; the error arrays are
   combined; the data quality arrays are "or'ed".
*/

	int i, j;	/* array indexes */
	int ibeg, iend;	/* loop limits */
	int jbeg, jend;	/* loop limits */
	float asci;	/* a science */
	float bsci;	/* b science */
	float bsci2;	/* b science squared */
	float aerr;	/* a error */
	float berr;	/* b error */

	ibeg = wf3->trimx[0]; iend = a->sci.data.nx - wf3->trimx[1];
	jbeg = wf3->trimy[0]; jend = a->sci.data.ny - wf3->trimy[1];

	for (j = jbeg; j < jend; j++) {
	     for (i = ibeg; i < iend; i++) {

		  /* error and science data */
		  asci = Pix(a->sci.data,i,j);
		  bsci = Pix(b->sci.data,i,j);
		  bsci2 = bsci*bsci;
		  aerr = Pix(a->err.data,i,j);
		  berr = Pix(b->err.data,i,j);

		  if (bsci != 0) {
		      Pix(a->sci.data,i,j) = asci / bsci;
		      Pix(a->err.data,i,j) =
		      sqrt(aerr*aerr/bsci2 + asci*asci*berr*berr/(bsci2*bsci2));
		  } else {
		      Pix(a->sci.data,i,j) = 0.0;
		      Pix(a->err.data,i,j) = 0.0;
		  }

		  /* data quality */
		  DQSetPix(a->dq.data,i,j,
			DQPix(a->dq.data,i,j) | DQPix(b->dq.data,i,j) );

	     }
	}
}

void asubk_noref (WF3Info *wf3, SingleNicmosGroup *a, float b) {

/* Subtract in-place a constant from a SingleNicmosGroup.

   (*a) -= b

   The constant is subtracted from the science data array only;
   the error and data quality arrays are unchanged.
*/

	int i, j;	/* array indexes */
	int ibeg, iend;	/* loop limits */
	int jbeg, jend;	/* loop limits */

	ibeg = wf3->trimx[0]; iend = a->sci.data.nx - wf3->trimx[1];
	jbeg = wf3->trimy[0]; jend = a->sci.data.ny - wf3->trimy[1];

	for (j = jbeg; j < jend; j++) {
	     for (i = ibeg; i < iend; i++) {

		  /* science data */
		  Pix(a->sci.data,i,j) -= b;

	     }
	}
}

void amulk_noref (WF3Info *wf3, SingleNicmosGroup *a, float b) {

/* Multiply in-place a SingleNicmosGroup by a constant.

   (*a) *= b

   The science and error data arrays are multiplied by the constant;
   the data quality array is unchanged.
*/

	int i, j;	/* array indexes */
	int ibeg, iend;	/* loop limits */
	int jbeg, jend;	/* loop limits */

	ibeg = wf3->trimx[0]; iend = a->sci.data.nx - wf3->trimx[1];
	jbeg = wf3->trimy[0]; jend = a->sci.data.ny - wf3->trimy[1];

	for (j = jbeg; j < jend; j++) {
	     for (i = ibeg; i < iend; i++) {

		  /* science data */
		  Pix(a->sci.data,i,j) *= b;

		  /* error data */
		  Pix(a->err.data,i,j) *= b;

	     }
	}
}

void aor_noref (WF3Info *wf3, SingleNicmosGroup *a, SingleNicmosGroup *b) {

/* Take the logical OR of two SingleNicmosGroup DQ arrays, leaving the result 
** in the first.

   (*a.dq) = (*a.dq) | (*b.dq)

   The science and error data arrays are unchanged.
*/

	int i, j;	/* array indexes */
	int ibeg, iend;	/* loop limits */
	int jbeg, jend;	/* loop limits */

	ibeg = wf3->trimx[0]; iend = a->sci.data.nx - wf3->trimx[1];
	jbeg = wf3->trimy[0]; jend = a->sci.data.ny - wf3->trimy[1];

	for (j = jbeg; j < jend; j++) {
	     for (i = ibeg; i < iend; i++) {
		  DQSetPix (a->dq.data,i,j,
			DQPix (a->dq.data,i,j) | DQPix (b->dq.data,i,j) );
	     }
	}
}

