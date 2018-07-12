# include <stdio.h>
# include <math.h>		/* for fabs */
# include "hstio.h"		/* for Hdr */

/* This routine gets from the extension header the keywords
   LTV1, LTV2, LTM1_1, and LTM2_2.  These are converted to
   bin and corner (zero indexed integers).

check:
   NOTE:  LTV1 for the CCD uses the beginning of the illuminated
   portion as the origin, not the beginning of the overscan region.
   Thus the xcorner that we compute has the same origin as ltv1,
   which is what we want, but it differs from the CENTERA1 header
   keyword, which has the beginning of the overscan region as origin.
*/

int GetCorner (Hdr *hdr, int rsize, int *bin, int *corner) {

/* arguments:
Hdr *hdr         i: extension header
int rsize        i: size of reference pixel in units of high-res pixels
int bin[2]       o: pixel size in X and Y
int corner[2]    o: corner of subarray in X and Y
*/

	extern int status;

	double ltm[2];		/* diagonal elements of MWCS matrix */
	double ltv[2];		/* MWCS linear transformation vector */

	int GetLT (Hdr *, double *, double *);
	int FromLT (int, double *, double *, int *, int *);

	/* Get the ltm and ltv values from the header. */
	if (GetLT (hdr, ltm, ltv))
	    return (status);

	/* Compute the corner location and bin size from ltm & ltv. */
	if (FromLT (rsize, ltm, ltv, bin, corner))
	    return (status);

	return (status);
}
