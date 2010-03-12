# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <math.h>

# include <hstio.h>	/* this is only used for writing debug image */
# include "../stis.h"
# include "calstis4.h"
# include "../stiserr.h"

static int AddTrace (StisInfo4 *, LampInfo *, DispRelation *,
		SpTrace *, double [], CmplxArray *);
static void ReadWidth (StisInfo4 *, double []);
static void debugimg (char *, CmplxArray *);
static int ESumSpec (double [], double [], int,
		DispRelation *, int, double, double,
		double [], int);
static double PixToWl (DispRelation *, int, double, double, double);

/* This routine reads the 1-D trace table (spectral trace, 1dt) and
   uses it to fill in the data for a template lamp image.

   Phil Hodge, 2000 Jan 5

   Phil Hodge, 2000 July 21:
	cmplx.h is now in ../
   Phil Hodge, 2004 July 23:
	Add to template in-place, rather than assigning values, because
	the latter does not work properly for slits that are not small
	compared with the spacing between orders.  Use slit_angle.
*/

int MakeTemplate (StisInfo4 *sts, LampInfo *lamp, DispRelation *disp,
		SpTrace *trace, CmplxArray *clamp) {

/* arguments:
StisInfo4 *sts      i: info
LampInfo *lamp      i: template lamp spectrum
DispRelation *disp  i: dispersion relation
SpTrace *trace      i: list of 1-D spectral traces
CmplxArray *clamp   o: 2-D complex array, values will be assigned
*/

	SpTrace *trace_o;	/* spectral trace for current spectral order */
	double slitwidth[2];	/* aperture size in image pixels */
	int i, j;
	int status;

	/* Read the aperture size and convert to pixels. */
	ReadWidth (sts, slitwidth);

	/* Initialize the template image to zero. */
	for (j = 0;  j < clamp->ny;  j++) {
	    for (i = 0;  i < clamp->nx;  i++) {
		RPIX2D (clamp, i, j) = 0.F;
		IPIX2D (clamp, i, j) = 0.F;
	    }
	}

	trace_o = trace;
	while (trace_o != NULL) {

	    if (status = AddTrace (sts, lamp, disp, trace_o, slitwidth, clamp))
		return (status);

	    trace_o = trace_o->next;
	}

	if (sts->dbgfile[0] != '\0')
	    debugimg (sts->dbgfile, clamp);	/* copy clamp to an image */

	return (0);
}

static int AddTrace (StisInfo4 *sts, LampInfo *lamp, DispRelation *disp,
		SpTrace *trace, double slitwidth[], CmplxArray *clamp) {

/* arguments:
StisInfo4 *sts      i: info
LampInfo *lamp      i: template lamp spectrum
DispRelation *disp  i: dispersion relation
SpTrace *trace      i: 1-D spectral trace
double slitwidth[]  i: slit width and height, in image pixels
CmplxArray *clamp   o: 2-D complex array, values will be assigned
*/

	/* reference pixel coordinates */
	int ix_ref;
	double x_ref, y_ref;

	/* image pixel coordinates */
	double y_im;
	int y_low, y_high;	/* range of pixels corresp. to slit height */

	/* scratch for template spectrum, integrated over pixels */
	double *tspec;
	int i, j;		/* loop indexes */
	int status;
	int ConvSlit (double, double [], int);

	double dtilt;	/* for taking into account the slit tilt */
	int i_tilt;	/* i - (nearest integer to dtilt) */

	/* If the middle of the trace is off the image, skip this order. */
	i = trace->a1center / 2;
	y_ref = trace->a2displ[i]  + trace->a2center;
	j = (int)(sts->ltm[1] * y_ref + sts->ltv[1]);
	if (j < 0 || j >= sts->ny)
	    return (0);

	if ((tspec = calloc (sts->nx, sizeof (double))) == NULL)
	    return (OUT_OF_MEMORY);

	/* Integrate the template lamp spectrum over the pixels of
	   the input image.
	*/
	if (status = ESumSpec (lamp->wl, lamp->flux, lamp->nelem,
		disp, trace->sporder,
		sts->ltm[0], sts->ltv[0], tspec, sts->nx))
	    return (status);

	/* Convolve the integrated lamp spectrum with the slit width. */
	if (status = ConvSlit (slitwidth[0], tspec, sts->nx))
	    return (status);

	/* Add the convolved, integrated template spectrum to the
	   template image.
	*/
	for (i = 0;  i < sts->nx;  i++) {

	    x_ref = (i - sts->ltv[0]) / sts->ltm[0];
	    ix_ref = NINT (x_ref);
	    if (ix_ref < 0 || ix_ref >= trace->nelem) {
		continue;
	    } else {
		y_ref = trace->a2center + trace->a2displ[ix_ref];
		y_im = sts->ltm[1] * y_ref + sts->ltv[1];
	    }

	    y_low  = NINT (y_im - slitwidth[1] / 2.);
	    y_high = NINT (y_im + slitwidth[1] / 2.);
	    if (y_high < 0)
		continue;
	    if (y_low >= sts->ny)
		continue;
	    if (y_low < 0)
		y_low = 0;
	    if (y_high >= sts->ny)
		y_high = sts->ny - 1;

	    for (j = y_low;  j <= y_high;  j++) {
		if (sts->slit_angle == 0.) {
		    RPIX2D (clamp, i, j) += tspec[i];
		} else {
		    dtilt = sts->slit_angle * (j - y_im);
		    i_tilt = i - NINT (dtilt);
		    if (i_tilt >= 0 && i_tilt < sts->nx)
			RPIX2D (clamp, i, j) += tspec[i_tilt];
		}
	    }
	}

	free (tspec);

	return (0);
}

/* This routine reads the aperture size in each axis from the APER_FOV
   keyword value.  The size is then converted to image pixels using
   CDELT2 (CD2_2) and LTMi_i from the input extension header.
*/

static void ReadWidth (StisInfo4 *sts, double slitwidth[]) {

/* arguments:
StisInfo4 *sts      i: info
double slitwidth[]  o: slit width and height, in image pixels
*/

	double xwidth, ywidth;		/* aperture size in arcseconds */
	char *wx, *dummy;

	/* default values */
	slitwidth[0] = 1.;
	slitwidth[1] = 1.;

	ywidth = strtod (sts->aper_fov, &wx);
	if (ywidth == 0.) {
	    printf ("Warning  Can't interpret APER_FOV = `%s'\n",
		sts->aper_fov);
	    return;
	}
	wx++;				/* skip over the "X" */
	xwidth = strtod (wx, &dummy);
	if (xwidth == 0.) {
	    printf ("Warning  Can't interpret APER_FOV = `%s'\n",
		sts->aper_fov);
	    return;
	}

	/* 3600 converts cdelt from degrees per pixel to arcsec per pixel. */
	slitwidth[0] = xwidth / (3600. * sts->cdelt[1]) * sts->ltm[0];
	slitwidth[1] = ywidth / (3600. * sts->cdelt[1]) * sts->ltm[1];
}

static void debugimg (char *dbgfile, CmplxArray *z) {

	SciHdrData x;
	int i, j;

	initFloatHdrData (&x);
	allocFloatHdrData (&x, z->nx, z->ny);

	for (j = 0;  j < z->ny;  j++) {
	    for (i = 0;  i < z->nx;  i++)
		Pix (x.data, i, j) = RPIX2D (z, i, j);
	}

	/* Write the data to dbgfile, extension 1, using option = 0. */
	putSci (dbgfile, 1, &x, 0);
	if (hstio_err()) {
	    printf ("Warning  Can't create template image:  %s\n",
			hstio_errmsg());
	    clear_hstioerr();
	}

	freeFloatHdrData (&x);
}

/* Integrate the reference (template) spectrum to match the pixel spacing
   of the observed data.
*/

static int ESumSpec (double wl[], double flux[], int nelem,
		DispRelation *disp, int sporder, double ltm, double ltv,
		double tspec[], int nwl) {

/* arguments:
double wl[nelem+1]  i: wavelengths for the template spectrum (edges of pixels)
double flux[nelem]  i: template spectrum flux
int nelem           i: size of flux array
DispRelation *disp  i: dispersion relation
int sporder         i: spectral order (needed for applying dispersion rel)
double ltm, ltv     i: transformation from reference to image pixel coords
double tspec[nwl]   o: template spectrum, integrated over image pixels
int nwl             i: size of tspec array
*/

	/* wavelengths at edges of a pixel in the integrated template */
	double wl_left, wl_right;

	/* indices in input template spectrum corresponding to
	   wl_left, wl_right */
	int jl, jr;

	int i, j;
	void FindWL (double, double [], int, int *);

	if (disp->ncoeff > 7) {
	    printf ("ERROR   Too many dispersion coefficients\n");
	    return (INTERNAL_ERROR);
	}

	/* -0.5 is the pixel coordinate at the left edge of the first pixel */
	wl_left = PixToWl (disp, sporder, -0.5, ltm, ltv);

	/* Find jl, the element in template spectrum containing wl_left. */
	jl = -1;				/* initial value */
	FindWL (wl_left, wl, nelem, &jl);
	jr = jl;				/* updated in loop */

	for (i = 0;  i < nwl;  i++) {

	    wl_right = PixToWl (disp, sporder, i + 0.5, ltm, ltv);

	    FindWL (wl_right, wl, nelem, &jr);		/* update jr */

	    if (jr >= nelem)		/* past the end of the template? */
		break;

	    /* integrate template spectrum over the range wl_left to wl_right */

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

	    jl = jr;		/* next pixel of template spectrum */
	    wl_left = wl_right;	/* next pixel of integrated spectrum */
	}

	return (0);
}

/* This function inverts the dispersion relation, returning the
   wavelength corresponding to an input image pixel coordinate.
*/

static double PixToWl (DispRelation *disp, int sporder, double pixel,
		double ltm, double ltv) {

/* arguments:
DispRelation *disp  i: dispersion relation
int sporder         i: spectral order
double pixel        i: X coordinate, zero-indexed image pixel units
double ltm, ltv     i: for converting from image to reference pixels
*/

	double m;	/* =sporder */
	double x_ref;	/* pixel converted to reference coordinates */
	double a, b, c;	/* coefficients:  a*x**2 + b*x + c = 0 */
	double discr;	/* discriminant of quadratic; scratch */
	double wl1, wl2;	/* two solutions of quadratic */

	m = (double)sporder;

	/* Convert to reference coords and add one for one-indexing. */
	x_ref = (pixel - ltv) / ltm + 1.;

	a = disp->coeff[2] * m * m + disp->coeff[6] * m;
	b = disp->coeff[1] * m + disp->coeff[4] + disp->coeff[5] * m * m;
	c = disp->coeff[0] + disp->coeff[3] * m - x_ref;

	discr = b * b - 4. * a * c;
	if (discr < 0.)
	    return (0.);
	discr = sqrt (discr);
	wl1 = (-b + discr) / (2. * a);
	wl2 = (-b - discr) / (2. * a);

	if (wl1 > 0. && wl2 > 0.) {

	    if (wl1 > wl2)		/* this is an arbitrary choice */
		return (wl1);
	    else
		return (wl2);

	} else if (wl1 <= 0. && wl2 <= 0.) {

	    return (0.);

	} else if (wl1 > 0.) {

	    return (wl1);

	} else {

	    return (wl2);
	}
}
