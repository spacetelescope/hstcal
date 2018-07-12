/* This file contains the following:
	unbin2d
	InterpInfo
	InterpDQInfo
*/

# include <stdio.h>
# include <stdlib.h>		/* calloc */
# include <string.h>		/* strncmp */
# include <math.h>		/* sqrt */

# include "c_iraf.h"
# include "hstio.h"

# include "stis.h"
# include "hstcalerr.h"
# include "stisdef.h"

static void InterpInfo (float, int, int *, float *, float *);
static void InterpDQInfo (float, int, int *, int *, int *);

/* This routine takes an input data array and expands it by linear
   interpolation, writing to the output array.  The calling routine
   must allocate the output SingleGroup (setting its size) and free
   it when done.

   The coordinate keywords in the output extension headers will be updated.

   Note that, in contrast to bin2d, this routine does not include
   any offset (i.e. xcorner, ycorner) in either input or output.

   Note that the errors are *not* increased by the factor
   sqrt (binx * biny), which would be necessary in order to make
   this routine the true inverse of bin2d.

   Phil Hodge, 1998 Oct 5:
	Change status values to GENERIC_ERROR_CODE or HEADER_PROBLEM.
*/

/* The computation of errors is not what one would normally do for
   linear interpolation, but it's reasonable in this context, which is
   that unbin2d should be the inverse of bin2d (except for the factor
   sqrt (binx*biny) mentioned above).
*/

int unbin2d (SingleGroup *a, SingleGroup *b) {

/* arguments:
SingleGroup *a        i: input data
SingleGroup *b        o: output data
*/

	int status;

	double block[2];	/* number of input pixels for one output */
	double offset[2] = {0., 0.};	/* offset of binned image */
	float p, q, r, s;	/* for interpolating */
	float xoffset, yoffset;	/* for getting location of output in input */
	float ai, aj;		/* location in input corresponding to m,n */
	float value;		/* interpolated value */
	float e1, e2, e3, e4;	/* errors at four corners */
	int inx, iny;		/* size of input array */
	int onx, ony;		/* size of output array */
	int binx, biny;		/* number of output pixels per input pixel */
	int m, n;		/* pixel index in output array */
	int i, j;		/* pixel index in input array */
	int i1, i2, num_i;	/* for ORing data quality */
	int j1, j2, num_j;	/* for ORing data quality */
	short dq;		/* OR of the data quality arrays */

	inx = a->sci.data.nx;
	iny = a->sci.data.ny;
	onx = b->sci.data.nx;
	ony = b->sci.data.ny;

	binx = onx / inx;
	biny = ony / iny;
	if (binx * inx != onx || biny * iny != ony) {
	    printf ("ERROR    (unbin2d) bin ratio is not an integer.\n");
	    return (GENERIC_ERROR_CODE);
	}

	xoffset = (float)(binx - 1) / 2.0F;
	yoffset = (float)(biny - 1) / 2.0F;

	if (binx == 1 && biny == 1) {

	    /* Same size, so just copy. */

	    /* Copy the science data. */
	    for (n = 0;  n < ony;  n++)
		for (m = 0;  m < onx;  m++)
		    Pix (b->sci.data, m, n) = Pix (a->sci.data, m, n);

	    /* Copy the error array. */
	    for (n = 0;  n < ony;  n++)
		for (m = 0;  m < onx;  m++)
		    Pix (b->err.data, m, n) = Pix (a->err.data, m, n);

	    /* Copy the data quality array. */
	    for (n = 0;  n < ony;  n++)
		for (m = 0;  m < onx;  m++)
		    DQSetPix (b->dq.data, m, n, DQPix(a->dq.data,m,n));

	} else if (binx == 1) {

	    /* Interpolate in Y. */

	    /* Science data array. */
	    for (n = 0;  n < ony;  n++) {
		aj = ((float)n - yoffset) / (float)biny;
		InterpInfo (aj, iny, &j, &r, &s);
		for (m = 0;  m < onx;  m++) {
		    value = r * Pix (a->sci.data, m, j) +
			    s * Pix (a->sci.data, m, j+1);
		    Pix (b->sci.data, m, n) = value;
		}
	    }

	    /* Error array. */
	    for (n = 0;  n < ony;  n++) {
		aj = ((float)n - yoffset) / (float)biny;
		InterpInfo (aj, iny, &j, &r, &s);
		for (m = 0;  m < onx;  m++) {
		    e1 = Pix (a->err.data, m, j);
		    e2 = Pix (a->err.data, m, j+1);
		    value = r * e1*e1 + s * e2*e2;
		    Pix (b->err.data, m, n) = sqrt (value);
		}
	    }

	    /* Data quality array. */
	    for (n = 0;  n < ony;  n++) {
		aj = ((float)n - yoffset) / (float)biny;
		InterpDQInfo (aj, iny, &j1, &j2, &num_j);
		for (m = 0;  m < onx;  m++) {
		    if (num_j == 1)
			dq = DQPix (a->dq.data, m, j1);
		    else
			dq = DQPix (a->dq.data, m, j1) |
			     DQPix (a->dq.data, m, j2);
		    DQSetPix (b->dq.data, m, n, dq);
		}
	    }

	} else if (biny == 1) {

	    /* Interpolate in X. */

	    /* Science data array. */
	    for (n = 0;  n < ony;  n++) {
		for (m = 0;  m < onx;  m++) {
		    ai = ((float)m - xoffset) / (float)binx;
		    InterpInfo (ai, inx, &i, &p, &q);
		    value = p * Pix (a->sci.data, i, n) +
			    q * Pix (a->sci.data, i+1, n);
		    Pix (b->sci.data, m, n) = value;
		}
	    }

	    /* Error array. */
	    for (n = 0;  n < ony;  n++) {
		for (m = 0;  m < onx;  m++) {
		    ai = ((float)m - xoffset) / (float)binx;
		    InterpInfo (ai, inx, &i, &p, &q);
		    e1 = Pix (a->err.data, i, n);
		    e2 = Pix (a->err.data, i+1, n);
		    value = p * e1*e1 + q * e2*e2;
		    Pix (b->err.data, m, n) = sqrt (value);
		}
	    }

	    /* Data quality array. */
	    for (n = 0;  n < ony;  n++) {
		for (m = 0;  m < onx;  m++) {
		    ai = ((float)m - xoffset) / (float)binx;
		    InterpDQInfo (ai, inx, &i1, &i2, &num_i);
		    if (num_i == 1)
			dq = DQPix (a->dq.data, i1, n);
		    else
			dq = DQPix (a->dq.data, i1, n) |
			     DQPix (a->dq.data, i2, n);
		    DQSetPix (b->dq.data, m, n, dq);
		}
	    }

	} else {

	    /* Science data array. */
	    for (n = 0;  n < ony;  n++) {
		aj = ((float)n - yoffset) / (float)biny;
		InterpInfo (aj, iny, &j, &r, &s);
		for (m = 0;  m < onx;  m++) {
		    ai = ((float)m - xoffset) / (float)binx;
		    InterpInfo (ai, inx, &i, &p, &q);
		    value = p * r * Pix (a->sci.data, i,   j) +
			    q * r * Pix (a->sci.data, i+1, j) +
			    p * s * Pix (a->sci.data, i,   j+1) +
			    q * s * Pix (a->sci.data, i+1, j+1);
		    Pix (b->sci.data, m, n) = value;
		}
	    }

	    /* Error array. */
	    for (n = 0;  n < ony;  n++) {
		aj = ((float)n - yoffset) / (float)biny;
		InterpInfo (aj, iny, &j, &r, &s);
		for (m = 0;  m < onx;  m++) {
		    ai = ((float)m - xoffset) / (float)binx;
		    InterpInfo (ai, inx, &i, &p, &q);
		    e1 = Pix (a->err.data, i,   j);
		    e2 = Pix (a->err.data, i+1, j);
		    e3 = Pix (a->err.data, i,   j+1);
		    e4 = Pix (a->err.data, i+1, j+1);
		    value = p * r * e1*e1 + q * r * e2*e2 +
			    p * s * e3*e3 + q * s * e4*e4;
		    Pix (b->err.data, m, n) = sqrt (value);
		}
	    }

	    /* Data quality array. */
	    for (n = 0;  n < ony;  n++) {
		InterpDQInfo (aj, iny, &j1, &j2, &num_j);
		for (m = 0;  m < onx;  m++) {
		    ai = ((float)m - xoffset) / (float)binx;
		    InterpDQInfo (ai, inx, &i1, &i2, &num_i);
		    if (num_i == 1 && num_j == 1)
			dq = DQPix (a->dq.data, i1, j1);
		    else if (num_i == 1)
			dq = DQPix (a->dq.data, i1, j1) |
			     DQPix (a->dq.data, i1, j2);
		    else if (num_j == 1)
			dq = DQPix (a->dq.data, i1, j1) |
			     DQPix (a->dq.data, i2, j1);
		    else
			dq = DQPix (a->dq.data, i1, j1) |
			     DQPix (a->dq.data, i2, j1) |
			     DQPix (a->dq.data, i1, j2) |
			     DQPix (a->dq.data, i2, j2);
		    DQSetPix (b->dq.data, m, n, dq);
		}
	    }
	}

	/* Copy the headers. */
	copyHdr (b->globalhdr, a->globalhdr);
	if (hstio_err())
	    return (HEADER_PROBLEM);
	copyHdr (&b->sci.hdr, &a->sci.hdr);
	if (hstio_err())
	    return (HEADER_PROBLEM);
	copyHdr (&b->err.hdr, &a->err.hdr);
	if (hstio_err())
	    return (HEADER_PROBLEM);
	copyHdr (&b->dq.hdr, &a->dq.hdr);
	if (hstio_err())
	    return (HEADER_PROBLEM);

	/* Update the coordinate parameters that depend on the binning. */
	block[0] = 1. / (double)binx;
	block[1] = 1. / (double)biny;
	if ((status = BinCoords (&a->sci.hdr, block, offset,
                                 &b->sci.hdr, &b->err.hdr, &b->dq.hdr)))
	    return (status);

	return (0);
}

/* This routine determines the appropriate array index i and weights
   p and q for linear interpolation.  If the array is called a, and ai
   is the independent variable (in units of the array index), then
   the interpolated value may be computed as:  p * a[i] + q * a[i+1].
*/

static void InterpInfo (float ai, int npts, int *i, float *p, float *q) {

/* arguments:
float ai        i: independent variable in same units as i
int npts        i: size of array within which *i is an index
int *i          o: array index close to ai
float *p, *q    o: weights for linear interpolation
*/

	*i = (int) ai;
	*i = (*i < 0) ? 0 : *i;
	*i = (*i >= npts - 1) ? (npts - 2) : *i;
	*q = ai - *i;
	*p = 1.0F - *q;
}

/* This routine determines which array indexes i1 and i2 to use for
   ORing the data quality information.
*/

static void InterpDQInfo (float ai, int npts, int *i1, int *i2, int *num_i) {

/* arguments:
float ai        i: independent variable in units of array index
int npts        i: size of array within which *i1 and *i2 are indexes
int *i1, *i2    o: array indexes close to ai
int *num_i      o: 1 or 2, i.e. use just *i1 or both *i1 and *i2
*/

	*i1 = (int) ai;
	*i2 = *i1 + 1;

	if (ai == (float)(*i1))
	    *num_i = 1;			/* ai is an integer, so just use i1 */
	else
	    *num_i = 2;

	if (*i1 <= 0) {
	    *i1 = 0;
	    *num_i = 1;			/* off the low end; just use i1 */
	}

	if (*i1 >= npts - 1) {
	    *i1 = npts - 1;
	    *num_i = 1;			/* off the high end; just use i1 */
	}
}
