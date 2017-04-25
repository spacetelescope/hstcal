/* This file contains the following:
	unbinsect
*/

# include <stdlib.h>		/* calloc */
# include <string.h>		/* strncmp */
# include <math.h>		/* sqrt */
# include "hstio.h"  /* for SingleGroupLine definitions */

# include "wf3.h"
# include "wf3info.h"
# include "hstcalerr.h"

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

   The computation of errors is not what one would normally do for
   linear interpolation, but it's reasonable in this context, which is
   that unbin2d should be the inverse of bin2d (except for the factor
   sqrt (binx*biny) mentioned above).

   Revision history:
   H. Bushouse	27-Jan-2009	Added check to computations of err value to make
				sure argument to sqrt() is positive.
*/

int unbinsect (WF3sect *a, int update, WF3sect *b ) {

/* arguments:
WF3sect *a        i: input binned data
WF3sect *b        o: output expanded data
*/

	extern int status;

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

	int BinCoords (Hdr *, double *, double *, Hdr *, Hdr *, Hdr *);
	void InterpInfo (float, int, int *, float *, float *);
	void InterpDQInfo (float, int, int *, int *, int *);

	inx = a->npix;
	iny = a->nlines;
	onx = b->npix;
	ony = b->nlines;

	binx = onx / inx;
	biny = ony / iny;
	if (binx * inx != onx || biny * iny != ony) {
	    trlerror ("ERROR    (unbinsect) bin ratio is not an integer.");
	    return (status = INVALID_VALUE);
	}

	xoffset = (float)(binx - 1) / 2.0F;
	yoffset = (float)(biny - 1) / 2.0F;

	if (binx == 1 && biny == 1) {

	    /* Same size, so just copy. */

	    /* Copy the science data. */
	    /* Copy the error array. */
	    /* Copy the data quality array. */
	    for (n = 0;  n < ony;  n++) {
		    for (m = 0;  m < onx;  m++) {
		        b->sci[n].line[m] = a->sci[n].line[m];
		        b->err[n].line[m] = a->err[n].line[m];
		        b->dq[n].line[m]  = a->dq[n].line[m];
            	    }
            }

	} else if (binx == 1) {

	    /* Interpolate in Y. */
	    for (n = 0;  n < ony;  n++) {
		 aj = ((float)n - yoffset) / (float)biny;
		 InterpInfo (aj, iny, &j, &r, &s);
		 InterpDQInfo (aj, iny, &j1, &j2, &num_j);

		 for (m = 0;  m < onx;  m++) {
        	      /* Science data array. */
		      value = r * a->sci[j].line[m] + s * a->sci[j+1].line[m];
		      b->sci[n].line[m] = value;

        	      /* Error array. */
		      e1 = a->err[j].line[m];
		      e2 = a->err[j+1].line[m];
		      value = r * e1*e1 + s * e2*e2;
		      b->err[n].line[m] = (value > 0.0) ? sqrt (value) : 0.0;

        	      /* Data quality array. */
		      if (num_j == 1)
			  dq = a->dq[j1].line[m];
		      else
	    		  dq = a->dq[j1].line[m] | a->dq[j2].line[m];
		      b->dq[n].line[m] = dq;
		 }
	    }

	} else if (biny == 1) {

	    /* Interpolate in X. */
	    for (n = 0;  n < ony;  n++) {
		 for (m = 0;  m < onx;  m++) {
		      ai = ((float)m - xoffset) / (float)binx;
		      InterpInfo (ai, inx, &i, &p, &q);

		      /* Science data array. */
		      value = p * a->sci[n].line[i] + q * a->sci[n].line[i+1];
		      b->sci[n].line[m] = value;

		      /* Error array. */
		      e1 = a->err[n].line[i];
		      e2 = a->err[n].line[i+1];
		      value = p * e1*e1 + q * e2*e2;
		      b->err[n].line[m] = (value > 0.0) ? sqrt (value) : 0.0;

		      /* Data quality array. */
		      InterpDQInfo (ai, inx, &i1, &i2, &num_i);
		      if (num_i == 1)
			  dq = a->dq[n].line[i1];
		      else
			  dq = a->dq[n].line[i1] | a->dq[n].line[i2];
		      b->dq[n].line[m] = dq;
		 }
	    }

	} else {

	    for (n = 0;  n < ony;  n++) {
		 aj = ((float)n - yoffset) / (float)biny;
		 InterpInfo (aj, iny, &j, &r, &s);
		 InterpDQInfo (aj, iny, &j1, &j2, &num_j);
		 for (m = 0;  m < onx;  m++) {
		      ai = ((float)m - xoffset) / (float)binx;
		      InterpInfo (ai, inx, &i, &p, &q);
                
		      /* Science data array. */
		      value = p * r * a->sci[j].line[i] +
			      q * r * a->sci[j].line[i+1] +
			      p * s * a->sci[j+1].line[i] +
			      q * s * a->sci[j+1].line[i+1];
		      b->sci[n].line[m] = value;

        	      /* Error array. */
		      e1 = a->err[j].line[i];
		      e2 = a->err[j].line[i+1];
		      e3 = a->err[j+1].line[i];
		      e4 = a->err[j+1].line[i+1];
		      value = p * r * e1*e1 + q * r * e2*e2 + p * s * e3*e3 +
			      q * s * e4*e4;
		      b->err[n].line[m] = (value > 0.0) ? sqrt (value) : 0.0;

        	      /* Data quality array. */
		      InterpDQInfo (ai, inx, &i1, &i2, &num_i);
		      if (num_i == 1 && num_j == 1)
			  dq = a->dq[j1].line[i1];
		      else if (num_i == 1)
			  dq = a->dq[j1].line[i1] | a->dq[j2].line[i1];
		      else if (num_j == 1)
			  dq = a->dq[j1].line[i1] | a->dq[j1].line[i2];
		      else
			  dq = a->dq[j1].line[i1] | a->dq[j1].line[i2] |
			       a->dq[j2].line[i1] | a->dq[j2].line[i2];
		      b->dq[n].line[m] = dq;
		 }
	    }
	}

	if (update == YES) {
	    /* Copy the headers. */
	    copyHdr (b->globalhdr, a->globalhdr);
	    if (hstio_err())
		return (status = 1001);
	    copyHdr (&b->sci->hdr, &a->sci->hdr);
	    if (hstio_err())
		return (status = 1001);
	    copyHdr (&b->err->hdr, &a->err->hdr);
	    if (hstio_err())
		return (status = 1001);
	    copyHdr (&b->dq->hdr, &a->dq->hdr);
	    if (hstio_err())
		return (status = 1001);

	    /* Update the coordinate parameters that depend on the binning. */
	    block[0] = 1. / (double)binx;
	    block[1] = 1. / (double)biny;
	    if (BinCoords (&a->sci->hdr, block, offset, &b->sci->hdr,
			   &b->err->hdr, &b->dq->hdr))
		return (status);
	}
	return (status);
}
