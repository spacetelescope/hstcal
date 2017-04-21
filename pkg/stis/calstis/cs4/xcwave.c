/* This file contains:
	XCWave
	SumSpec
	ChopTemplate
	PixToWl
	GetMaxPixel
*/

# include <stdio.h>
# include <stdlib.h>	/* calloc */
# include <string.h>	/* memcpy */
# include <math.h>

# include "stis.h"
# include "calstis4.h"
# include "err.h"

static int SumSpec (double [], double [], int, int, double [],
		double, double, double, double, double,
		double [], int);
static void ChopTemplate (double [], short [], int, short);
static double PixToWl (int, double [],
		double, double, double, double, double, double);
static void GetMaxPixel (int, double [], int, double, double);

/* This is reassigned by calling GetMaxPixel in XCWave, and it is used
   by PixToWl.  For image pixel coordinates (dispersion direction)
   greater than this, PixToWl will return the maximum wavelength for
   the prism, rather than evaluating the dispersion relation.  This is
   needed because the dispersion relation turns over for large pixel
   coordinates.
*/
static double max_pixel = 2048.;	/* image pixel units */

/* The size of the cross correlation array must be at least three pixels
   because we will fit a quadratic to the three values nearest the peak.
*/
# define MIN_RANGE       3

/* This is the number of extra pixels that we chop off the ends of the
   template spectrum, beyond the point where the data begin to be flagged
   as good.
*/
# define CHOP_EXTRA     10

/* This routine finds the shift of a spectrum in the dispersion direction
   by cross correlating the data with a reference spectrum, i.e. the
   calibration lamp spectrum read from a reference file.

   If the data are at a larger pixel number than the calibration spectrum,
   then the shift will be positive.

   Phil Hodge, 1998 Dec 11:
	Add sts to calling sequence, for wl_range and dbg; write
	wavelength, convolved template, observed spectrum to debug file.

   Phil Hodge, 2000 Jan 5:
	Move ConvSlit and FindWL to separate files.

   Phil Hodge, 2000 May 2:
	Add ChopTemplate, to set the template spectrum to zero at the
	ends of the array.
	Print the DQ value to the debug file for every pixel, not just
	when DQ is non-zero.

   Phil Hodge, 2001 Feb 28:
	Add disp to calling sequence.  In SumSpec, call PixToWl to
	convert pixel coordinate to a wavelength.  Both these changes
	were added to support non-linear wavelengths (e.g. prism data).
	Remove normalization of template spectrum in SumSpec.

   Phil Hodge, 2010 July 29:
	In PixToWl, call prismDisp instead of explicitly evaluating the
	prism dispersion relation.
*/

int XCWave (StisInfo4 *sts, double wl[], double flux[], int nelem,
	double v[], short qv[], int nwl, short sdqflags,
	double slitwidth,
	DispRelation *disp, double crpix, double crval, double cdelt,
	double *shift) {

/* arguments:
StisInfo4 *sts      i: calibration switches and info
double wl[]         i: array of wavelengths for the reference spectrum
double flux[]       i: array of fluxes, the reference spectrum
int nelem           i: size of wl and flux arrays
double v[]          i: array containing observed spectrum
short qv[]          i: data quality flags corresponding to v
int nwl             i: size of v and qv arrays
short sdqflags      i: "serious" data quality flags
double slitwidth    i: width of slit (dispersion direction) in pixels
DispRelation *disp  i: dispersion relation (currently only used for prism)
double crpix        i: reference pixel number in v array (zero indexed)
double crval        i: wavelength (Angstroms) at the reference pixel
double cdelt        i: Angstroms per pixel
double *shift       o: the shift, in pixels
*/

	int status;

	int range;		/* range for cross correlation */
	double *tspec;		/* resampled reference spectrum */
	double c7_shift;	/* shift in units of calstis7 pixels */
	int i;			/* loop index for debug file */

	int ConvSlit (double, double [], int);
	int XCPeak (StisInfo4 *, double *, short *, double *,
		int, int, short, double *);

	range = sts->wl_range;
	if (range > nwl / 2)
	    range = nwl / 2;
	if (range < MIN_RANGE)
	    range = MIN_RANGE;
	if (range / 2 * 2 == range)
	    range++;			/* must be odd */

	tspec = (double *) calloc (nwl, sizeof (double));
	if (tspec == NULL)
	    return (OUT_OF_MEMORY);

	/* Set the value of max_pixel, for use by PixToWl. */
	GetMaxPixel (sts->disp_type, disp->coeff, nwl,
			sts->ltm[0], sts->ltv[0]);

	/* Resample reference spectrum to match observed data. */
	if ((status = SumSpec (wl, flux, nelem,
			sts->disp_type, disp->coeff, crpix, crval, cdelt,
			sts->ltm[0], sts->ltv[0],
                        tspec, nwl))) {
	    free (tspec);
	    return (status);
	}

	/* Convolve the reference spectrum with the slit width. */
	if ((status = ConvSlit (slitwidth, tspec, nwl))) {
	    free (tspec);
	    return (status);
	}

	/* Chop off (by setting to zero) the template spectrum at the
	   ends of the array where the data are flagged as bad.
	*/
	ChopTemplate (tspec, qv, nwl, sdqflags);

	/* Write info to debug file. */
	if (sts->dbg != NULL) {
	    double wavelength;
	    fprintf (sts->dbg,
"# (XCWave) pixel, wavelength, convolved template, observed spectrum, DQ:\n");
	    for (i = 0;  i < nwl;  i++) {
		wavelength = PixToWl (sts->disp_type, disp->coeff, crpix,
			crval, cdelt, sts->ltm[0], sts->ltv[0], (double)i);
		fprintf (sts->dbg, "%d %.4f %.6g %.6g %d\n",
			i + 1,			/* one indexed */
			wavelength, tspec[i], v[i], qv[i]);
	    }
	}

	/* Do the cross correlation and find the shift. */
	if ((status = XCPeak (sts, v, qv, tspec, nwl, range, sdqflags,
                              &c7_shift))) {
	    free (tspec);
	    return (status);
	}

	*shift = c7_shift;

	free (tspec);

	return (0);
}

/* Integrate the reference (template) spectrum to match the pixel spacing
   of the observed data.
*/

static int SumSpec (double wl[], double flux[], int nelem,
		int disp_type, double coeff[],
		double crpix, double crval, double cdelt,
		double ltm, double ltv,
		double tspec[], int nwl) {

/* arguments:
double wl[nelem+1]  i: wavelengths for the template spectrum (edges of pixels)
double flux[nelem]  i: template spectrum flux
int nelem           i: size of flux array
int disp_type       i: distinguishes grating from prism
double coeff[]      i: dispersion coefficients (currently used only for prism)
double crpix        i: reference pixel number in v array (zero indexed)
double crval        i: wavelength (Angstroms) at the reference pixel
double cdelt        i: Angstroms per pixel
double tspec[nwl]   o: template spectrum, integrated over image pixels
int nwl             i: size of v and qv arrays
*/

	/* wavelengths at edges of a pixel in the observed spectrum */
	double wl_left, wl_right;

	/* indices in input template spectrum corresponding to
	   wl_left, wl_right */
	int jl, jr;

	int i, j;
	void FindWL (double, double *, int, int *);

	/* -0.5 and +0.5 are the pixel coordinates at the left and right
	   edges of the first pixel.  (1 is spectral order)
	*/
	wl_left = PixToWl (disp_type, coeff, crpix, crval, cdelt,
			ltm, ltv, -0.5);
	wl_right = PixToWl (disp_type, coeff, crpix, crval, cdelt,
			ltm, ltv, 0.5);

	/* Find jl, the element in template spectrum containing wl_left. */
	jl = -1;				/* initial value */
	FindWL (wl_left, wl, nelem, &jl);
	jr = jl;				/* updated in loop */

	for (i = 0;  i < nwl;  i++) {

	    /* integrate template spectrum over the range wl_left to wl_right */

	    FindWL (wl_right, wl, nelem, &jr);		/* update jr */

	    if (jr >= nelem)		/* past the end of the template? */
		break;

	    if (jl >= 0) {
		if (jl == jr) {
		    tspec[i] = flux[jl] * (wl_right - wl_left);
		} else {
		    tspec[i] = flux[jl] * (wl[jl+1] - wl_left);
		    for (j = jl+1;  j < jr;  j++)
			tspec[i] += flux[j] * (wl[j+1] - wl[j]);
		    tspec[i] += flux[jr] * (wl_right - wl[jr]);
		}
	    }

	    jl = jr;			/* next pixel of template spectrum */
	    wl_left = wl_right;
	    wl_right = PixToWl (disp_type, coeff, crpix, crval, cdelt,
			ltm, ltv, i+1.5);
	}

	return (0);
}

/* This routine sets the template spectrum to zero at the endpoints where
   the observed wavecal spectrum is flagged as bad.  Starting at the left
   endpoint, each element of tspec is set to zero until an element is
   encountered that is not flagged as bad.  Then the same is done starting
   at the right endpoint.  An additional few pixels (CHOP_EXTRA) of tspec
   will also be set to zero at each end, to allow for MSM slop.

   If all of tspec would be set to zero, due to the range of pixels flagged
   by qv, this routine returns without actually changing tspec.
*/

static void ChopTemplate (double tspec[], short qv[], int nwl,
		short sdqflags) {

/* arguments:
double tspec[]     io: template spectrum
short qv[]         i: data quality flags for observed wavecal
int nwl            i: size of arrays
short sdqflags     i: "serious" data quality flags
*/

	int i;
	int first_good, last_good;

	first_good = 0;				/* initial value */
	for (i = 0;  i < nwl;  i++) {
	    if (!(qv[i] & sdqflags)) {		/* not flagged as bad? */
		first_good = i;
		break;
	    }
	}

	last_good = nwl - 1;
	for (i = nwl - 1;  i >= 0;  i--) {
	    if (!(qv[i] & sdqflags)) {
		last_good = i;
		break;
	    }
	}

	first_good += CHOP_EXTRA;
	if (first_good >= nwl)
	    return;
	last_good -= CHOP_EXTRA;
	if (last_good < 0)
	    return;

	for (i = 0;  i < first_good;  i++)
	    tspec[i] = 0.;
	for (i = last_good + 1;  i < nwl;  i++)
	    tspec[i] = 0.;
}

/* This function returns the wavelength corresponding to an input image
   pixel coordinate.  This can be used for 2-D rectified or prism data,
   but not for non-2-D rectified echelle data (see maketemplate.c).
*/

static double PixToWl (int disp_type, double coeff[],
		double crpix, double crval, double cdelt,
		double ltm, double ltv, double pixel) {

/* arguments:
int disp_type       i: 2-D rectified (i.e. first order), prism, or echelle
double coeff[]      i: coefficients for dispersion relation (for PRISM only)
double crpix        i: these three coordinate parameters are already set
double crval        i:    to be used for zero-indexed image pixels
double cdelt        i:
double ltm, ltv     i: for converting from image to reference pixels
double pixel        i: X coordinate, zero-indexed image pixel units
*/

	double x_ref0;	/* pixel converted to reference coordinates */
	double wl;	/* value to be returned */

	if (disp_type == RECTIFIED) {

	    wl = crval + (pixel - crpix) * cdelt;

	} else if (disp_type == PRISM_DISP) {

	    if (pixel > max_pixel) {

		/* dispersion relation is unreliable at this point */
		wl = MAX_PRISM_WAVELENGTH;

	    } else {

		/* Convert to reference coords, but leave as zero-indexed. */
		x_ref0 = (pixel - ltv) / ltm;

		wl = prismDisp (coeff, x_ref0);

		if (wl > MAX_PRISM_WAVELENGTH)
		    wl = MAX_PRISM_WAVELENGTH;
	    }

	} else {

	    wl = 0.;		/* shouldn't happen */
	}

	return (wl);
}

/* This routine sets the value of max_pixel, which is global to this file.
   max_pixel is used by PixToWl.
*/

static void GetMaxPixel (int disp_type, double coeff[],
		int nwl, double ltm, double ltv) {

	int i;
	double wl;
	double prev_wl;

	if (disp_type != PRISM_DISP)
	    return;

	/* initial value, pixel 0 */
	prev_wl = PixToWl (disp_type, coeff, 0., 0., 0., ltm, ltv, 0.);

	for (i = 1;  i < nwl;  i++) {

	    wl = PixToWl (disp_type, coeff, 0., 0., 0., ltm, ltv, (double)i);

	    if (wl >= MAX_PRISM_WAVELENGTH || wl < prev_wl) {
		max_pixel = (double)i;
		return;
	    }
	}
}
