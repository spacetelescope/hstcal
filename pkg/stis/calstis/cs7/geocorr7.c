# include <stdio.h>
# include <math.h>		/* for sqrt in Interp7 */

# include "c_iraf.h"
# include "hstio.h"

# include "stis.h"
# include "calstis7.h"
# include "hstcalerr.h"
# include "stisdq.h"
# include "stisdef.h"

static int GeoOnePoint (StisInfo7 *, double, double,
	DistInfo *, FloatHdrData *, FloatHdrData *,
	double *, double *);
static void GeoEval (DistInfo *, double, double, double *, double *);
static int InverseMap (StisInfo7 *, double, double,
	DistInfo *, FloatHdrData *, FloatHdrData *,
	double *, double *);
static int UpdateCrpix (SingleGroup *, double []);
static void ApplyLT (double, double, double [], double [],
			double *, double *);
static void InverseLT (double, double, double [], double [],
			double *, double *);

/* This routine corrects for geometric distortion of an image.
   Note that this is only used for OBSTYPE = IMAGING.

   The mapping from distorted to undistorted coordinates is done as follows:

	output image size is CXSIZE, CYSIZE from IDC table, but scaled
		according to binning of input image

	for each pixel in undistorted image:
	    convert from image pixels to reference pixels
	    subtract the reference point
		(dist->cxref, dist->cyref = CXREF, CYREF from IDC table)
	    apply the mapping (DIRECTION="INVERSE", CXij, CYij)
	    add the reference point
		(dist->xref, dist->yref = XREF, YREF from IDC table)
	    convert from reference pixels to image pixels

   Phil Hodge, 1998 Dec 17:
	Set sdqflags to zero locally before calling Interp2D.

   Phil Hodge, 1999 Sept 8:
	Add total_offset to crpix, not to input pixel location [ix,iy].

   Phil Hodge, 2000 Jan 12:
	Remove coord_o and slit from the calling sequence, and don't call
	DataMasked.

   Phil Hodge, 2000 Dec 5:
	Evaluate the coefficients directly (new function GeoEval), rather
	than using gsurfit.

   Phil Hodge, 2006 Feb 7:
	Add err_algorithm to the calling sequence of Interp2D.
*/

int GeoCorr7 (StisInfo7 *sts, DistInfo *dist,
	FloatHdrData *ssgx, FloatHdrData *ssgy,
	SingleGroup *in, SingleGroup *out) {

/* arguments:
StisInfo7 *sts             i: calibration switches and info
DistInfo *dist             i: distortion coefficients
FloatHdrData *ssgx, *ssgy  i: small-scale distortions in X and Y
SingleGroup *in            i: input data
SingleGroup *out           o: output data
*/

	int status;

	double ix, iy;		/* pixel coordinates in input image */
	double ox, oy;		/* pixel coordinates in output image */
	double crpix[2];	/* corrected for distortion */

	int i, j;		/* pixel coordinates in output array */
	double jacobian;	/* for conserving flux */
	short dq;		/* an interpolated data quality value */
	short sdqflags;		/* local value, set to zero */

	jacobian = 1.;		/* do not conserve flux */

	/* Fill each pixel of the output array. */
	for (j = 0;  j < out->sci.data.ny;  j++) {
	    for (i = 0;  i < out->sci.data.nx;  i++) {

		ox = (double) i;
		oy = (double) j;

		/* Find the point [ix,iy] corresponding to [ox,oy]. */
		if ((status = GeoOnePoint (sts, ox, oy, dist,
                                           ssgx, ssgy, &ix, &iy)))
		    return (status);

		/* interpolate */
		sdqflags = 0;		/* instead of sts->sdqflags */
		Interp2D (in, sdqflags, ix, iy, jacobian,
				sts->err_algorithm,
				&Pix(out->sci.data,i,j),
				&Pix(out->err.data,i,j), &dq);
		DQSetPix (out->dq.data, i, j, dq);
	    }
	}

	/* Apply the mapping to CRPIXi, and update in output.
	   (total_offset is currently zero)
	*/
	if ((status = InverseMap (sts,
		sts->crpix[0] + sts->total_offset[0],
		sts->crpix[1] + sts->total_offset[1],
                dist, ssgx, ssgy, &crpix[0], &crpix[1])))
	    return (status);
	if ((status = UpdateCrpix (out, crpix)))
	    return (status);

	return (0);
}

/* This routine evaluates the mapping from a single point in the output
   image to the corresponding point in the input image.
*/

static int GeoOnePoint (StisInfo7 *sts, double ox, double oy,
	DistInfo *dist, FloatHdrData *ssgx, FloatHdrData *ssgy,
	double *ix, double *iy) {

/* arguments:
StisInfo7 *sts             i: calibration switches and info
double ox, oy              i: pixel location in output
DistInfo *dist             i: distortion coefficients
FloatHdrData *ssgx, *ssgy  i: small-scale distortions in X and Y
double *ix, *iy            o: location in input corresponding to [ox,oy]
*/

	double ix_r, iy_r;	/* coords in input, but in ref coords */
	double ox_r, oy_r;	/* coords in output, but in ref coords */

	int x, y;		/* a pixel location in input array */

	/* Convert to reference coordinates. */
	InverseLT (ox, oy, sts->ltm, sts->ltv, &ox_r, &oy_r);

	/* Find corresponding pixel location in input.  Note
	   that we have reference coordinates at this point,
	   not necessarily coords in input image.
	*/
	GeoEval (dist, ox_r, oy_r, &ix_r, &iy_r);

	/* Add offsets due to small-scale geometric distortion. */
	if (sts->sgeocorr == PERFORM) {
	    x = NINT(ix_r);  y = NINT(iy_r);
	    /* within the image? */
	    if (x < 0)
		x = 0;
	    if (x > ssgx->data.nx - 1)
		x = ssgx->data.nx - 1;
	    if (y < 0)
		y = 0;
	    if (y > ssgx->data.ny - 1)
		y = ssgx->data.ny - 1;

	    ix_r += Pix (ssgx->data, x, y);
	    iy_r += Pix (ssgy->data, x, y);
	}

	/* Convert from reference coords to input pixel coords. */
	ApplyLT (ix_r, iy_r, sts->ltm, sts->ltv, ix, iy);

	return (0);
}


/* This routine evaluates the polynomial coefficients at one point in
   the undistorted output array to give x and y in the distorted input
   array.  Both input and output are reference pixel coordinates,
   not necessarily image pixels.

   Note that the way the coefficients are interpreted, the reference
   point gets mapped to the same pixel position.  The reference point
   is (xref,yref), which is not necessarily the same as the reference
   pixel (CRPIX1,CRPIX2).
*/

static void GeoEval (DistInfo *dist,
		double ox, double oy, double *ix, double *iy) {

/* arguments:
DistInfo *dist    i: distortion coefficients
double ox, oy     i: pixel coordinates (reference pixels), undistorted
double *ix, *iy   o: pixel coordinates (reference pixels), distorted
*/

	double x, y;		/* difference from reference pixel */
	double xpow[MAX_ORDER+1], ypow[MAX_ORDER+1];
	int i, j, k;

	/* Subtract the reference point, and scale to arcseconds. */
	x = (ox - dist->cxref) * dist->scale;
	y = (oy - dist->cyref) * dist->scale;

	/* set up arrays of powers of x and y */
	xpow[0] = 1.;
	ypow[0] = 1.;
	for (i = 1;  i <= dist->norder;  i++) {
	    xpow[i] = xpow[i-1] * x;
	    ypow[i] = ypow[i-1] * y;
	}

	/* initial values;
	   this assumes xcoeff[0] and ycoeff[0] are always zero,
	   but we're going to add those coefficients anyway.
	*/
	*ix = dist->xref;
	*iy = dist->yref;

	for (i = 0;  i <= dist->norder;  i++) {

	    for (j = 0;  j <= i;  j++) {

		k = WHICH_COEFF (i, j);

		*ix += dist->xcoeff[k] * xpow[j] * ypow[i-j];
		*iy += dist->ycoeff[k] * xpow[j] * ypow[i-j];
	    }
	}
}

/* This routine evaluates the inverse mapping from a single point in
   the input image to the corresponding point in the output image.
*/

static int InverseMap (StisInfo7 *sts, double ix, double iy,
	DistInfo *dist, FloatHdrData *ssgx, FloatHdrData *ssgy,
	double *ox, double *oy) {

	int status;

	double tempx, tempy;
	double x, y;

	/* first approximation */
	if ((status = GeoOnePoint (sts, ix, iy, dist, ssgx, ssgy,
                                   &tempx, &tempy)))
	    return (status);
	x = ix + (ix - tempx);
	y = iy + (iy - tempy);

	/* next iteration */
	if ((status = GeoOnePoint (sts, x, y, dist, ssgx, ssgy,
                                   &tempx, &tempy)))
	    return (status);
	x = ix + (x - tempx);
	y = iy + (y - tempy);

	/* one more */
	if ((status = GeoOnePoint (sts, x, y, dist, ssgx, ssgy,
                                   &tempx, &tempy)))
	    return (status);
	*ox = ix + (x - tempx);
	*oy = iy + (y - tempy);

	return (0);
}

static int UpdateCrpix (SingleGroup *out, double crpix[]) {

	int status;

	/* Note:  add one to crpix to convert to one indexing. */

	if ((status = Put_KeyD (&out->sci.hdr, "CRPIX1", crpix[0]+1., "")))
	    return (status);
	if ((status = Put_KeyD (&out->sci.hdr, "CRPIX2", crpix[1]+1., "")))
	    return (status);

	if ((status = Put_KeyD (&out->err.hdr, "CRPIX1", crpix[0]+1., "")))
	    return (status);
	if ((status = Put_KeyD (&out->err.hdr, "CRPIX2", crpix[1]+1., "")))
	    return (status);

	if ((status = Put_KeyD (&out->dq.hdr, "CRPIX1", crpix[0]+1., "")))
	    return (status);
	if ((status = Put_KeyD (&out->dq.hdr, "CRPIX2", crpix[1]+1., "")))
	    return (status);

	return (0);
}

/* This routine transforms reference pixel coordinates to image pixel
   coordinates.  The values of ix_r and iy_r input to this routine
   are in the reference coordinate system (unbinned CCD or low-res MAMA,
   and full detector); the values of ix and iy output from this routine
   are pixel coordinates in the current input image.
*/

static void ApplyLT (double ix_r, double iy_r,
		double ltm[], double ltv[],
		double *ix, double *iy) {

/* arguments:
double ix_r, iy_r  i: pixel coordinates in input image, in reference coords
double ltm[2]      i: linear transformation, diagonal of matrix part
double ltv[2]      i: linear transformation, vector part
double *ix, *iy    o: pixel coordinates in input image
*/

	*ix = ix_r * ltm[0] + ltv[0];
	*iy = iy_r * ltm[1] + ltv[1];
}

/* This is the inverse of ApplyLT; it transforms the image coordinates to
   reference pixel coordinates.
*/

static void InverseLT (double ox, double oy,
		double ltm[], double ltv[],
		double *ox_r, double *oy_r) {

/* arguments:
double ox, oy        i: pixel coordinates in output image
double ltm[2]        i: linear transformation, diagonal of matrix part
double ltv[2]        i: linear transformation, vector part
double *ox_r, *oy_r  o: pixel coordinates in output image, reference coords
*/

	*ox_r = (ox - ltv[0]) / ltm[0];
	*oy_r = (oy - ltv[1]) / ltm[1];
}
