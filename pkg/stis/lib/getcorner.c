# include <stdio.h>
# include <math.h>		/* for fabs */
# include "c_iraf.h"
# include "hstio.h"		/* for Hdr */

# include "stis.h"
# include "stisdef.h"

static void FromLT (int, double *, double *, int *, int *);

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

	int status;

	double ltm[2];		/* diagonal elements of MWCS matrix */
	double ltv[2];		/* MWCS linear transformation vector */

	/* Get the ltm and ltv values from the header. */
	if ((status = GetLT (hdr, ltm, ltv)))
	    return (status);

	/* Compute the corner location and bin size from ltm & ltv. */
	FromLT (rsize, ltm, ltv, bin, corner);

	return (0);
}

/* This routine uses the LTV and LTM keyword values to compute the corner
   location and pixel size.  For the MAMA detectors, the corner location
   and pixel size are in units of high-res pixels, in contrast to the
   header keywords CENTERA1, SIZAXIS1, etc.

   ltm[0] and ltm[1] are assumed to be greater than zero.
*/

static void FromLT (int rsize, double *ltm, double *ltv,
	int *bin, int *corner) {

/* arguments:
int rsize        i: reference pixel size, 1 or 2
double ltm[2]    i: diagonal elements of MWCS matrix
double ltv[2]    i: MWCS linear transformation vector
int bin[2]       o: pixel size in X and Y
int corner[2]    o: corner of subarray in X and Y
*/

	double dbinx, dbiny, dxcorner, dycorner;

	dbinx = (double)rsize / ltm[0];
	dbiny = (double)rsize / ltm[1];

	dxcorner = (dbinx - rsize) / 2. - dbinx * ltv[0];
	dycorner = (dbiny - rsize) / 2. - dbiny * ltv[1];

	/* Round off to the nearest integer. */
	corner[0] = NINT (dxcorner);
	corner[1] = NINT (dycorner);
	bin[0] = NINT (dbinx);
	bin[1] = NINT (dbiny);
}
