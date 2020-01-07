/* This file includes:
	X2DCorr7
	CheckBoundary
	PixelPhys
	ApplyLT
	EvalTrace
	FreeInang
*/

# include <math.h>
# include <stdio.h>
# include <string.h>

# include <c_iraf.h>
# include <hstio.h>

# include "stis.h"
# include "calstis7.h"
# include "hstcalerr.h"
# include "stisdq.h"
# include "stisdef.h"

/* Don't interpolate if the sum of weights is smaller than this value. */
# define TINY_WEIGHT  (0.1)
/* Criterion for being very close to the center of a pixel. */
# define NEAR_CENTER  (0.1)

static void ApplyLT (double, double, double *, double *, double *, double *);
static int CheckBoundary (StisInfo7 *, CoordInfo *, SpTrace **, int, int);
static void EvalTrace (SpTrace *, double, double *, double *);
static void FreeInang (InangInfo *);
static void PixelPhys (CoordInfo *, double, double, double,
		double *, double *, double *, double *);
static void SgeoMsg (StisInfo7 *, int);

/* These are flags to indicate whether we have already printed a warning
   about inconsistency between the number of dispersion coefficients and
   the number of incidence-angle coefficients.
*/
static int warn1, warn2;

/* This routine performs 2-D rectification of a spectrum.
   Note that this is only used for OBSTYPE = SPECTROSCOPIC.

   Note:  The function value will be -1 if the current order is off the
   image or a table row that we need has dummy pedigree.  This may not be
   a fatal error.

   Phil Hodge, 1997 Dec 9:
	Don't call CheckBoundary if first order grating was used.

   Phil Hodge, 1998 Dec 17:
	Set sdqflags to zero locally before calling Interp2D.

   Phil Hodge, 2000 Mar 10:
	Include option to evaluate prism dispersion relation, using
	functions SolveDisp and PrismDisp.

   Phil Hodge, 2000 May 31:
	In EvalDisp, add a seventh term (cubic) in the dispersion relation.

   Phil Hodge, 2000 Aug 9:
	Remove all references to MAMA offset table or coefficients.
	Remove ssgx and ssgy from the calling sequence.

   Phil Hodge, 2001 Feb 22:
	Move prism macros to ../stis.h.

   Phil Hodge, 2001 May 2:
	Move the call to GetDisp from here to Do2Dx, and add disp to
	calling sequence.

   Ivo Busko, 2002 Mar 5:
   	Move EvalDisp and associated functions to its own file.

   Paul Barrett, 2004 Feb 11:
        Correct long slits for tilt angle.

   Paul Barrett and Phil Hodge, 2005 Jan 19:
	Remove scaling factor for binning from correction for tilt angle,
	and include binning in the cross-dispersion direction when applying
	del_tan (the difference in tangents of the angles).

   Phil Hodge, 2006 Feb 7:
	Add err_algorithm to the calling sequence of Interp2D.

   Phil Hodge, 2010 July 30:
	Change the name of EvalDisp to EvalDisp7.
*/

int X2DCorr7 (StisInfo7 *sts,
	CoordInfo *coord_o, DispRelation *disp, ApInfo *slit, int o_extver,
	SingleGroup *in, SingleGroup *out) {

/* arguments:
StisInfo7 *sts             i: calibration switches and info
CoordInfo *coord_o         i: coordinate info
DispRelation *disp         i: dispersion coefficients and other info
ApInfo *slit               i: description of slit
int o_extver               i: zero on first call (for printing info in SgeoMsg)
SingleGroup *in            i: input data
SingleGroup *out           o: output data
*/

	int status;

	double ix, iy;		/* pixel coordinates in input image */
	double ox, oy;		/* pixel coordinates in output image */
	double ix_r, iy_r;	/* coords in input, but in ref coords */
	double dix_r, diy_r;	/* increments in ix_r, iy_r */
	double wl, s;		/* wavelength and position along slit */
	double dwl, ds;		/* increments in wl and s */
	double ydispl;		/* spectrum trace at current pixel */
	double jacobian;	/* ratio of input to output pixel areas */
        double dtor;            /* conversion of degrees to radians */
        double del_tan;         /* difference of tangent angles */
	int disp_type;		/* disperser is grating or prism */
	int i, j;		/* pixel number in output image */
	short dq;		/* data quality at a pixel */
	short sdqflags;		/* local value, set to zero */

	/* approx line number in input image, but in reference coords */
	double line0;

	/* X offset (arcsec) from slit used to measure dispersion coeff */
	double delta;
	int GetApOffset (StisInfo7 *, ApInfo *, char *, double *);

	InangInfo iac;		/* incidence-angle correction coeff */
	int GetInang (StisInfo7 *, RefTab *, int, InangInfo *);

	DispRelation *disp_y;	/* dispersion relation interpolated at y */
	int InterpDisp (DispRelation **, double, DispRelation **);
	void AdjustDisp (StisInfo7 *, DispRelation *,double, InangInfo *,
	                 int *, int *);
	void FreeDisp (DispRelation **);

	SpTrace *trace;		/* list of spectrum traces */
	SpTrace *trace_y;	/* spectrum trace interpolated at y */
	int GetTrace (StisInfo7 *, int, SpTrace **);
	int InterpTrace (SpTrace **, double, SpTrace **);
	void FreeTrace (SpTrace **);

	int DataMasked (StisInfo7 *, CoordInfo *, ApInfo *, SingleGroup *);

	void EvalDisp7 (DispRelation *, int, double, int, double *, double *);

	disp_y = NULL;

	trace = NULL;
	trace_y = NULL;

	iac.allocated = 0;

	sts->x2dcorr_o = PERFORM;	/* initial value for current order */

	if (strcmp (sts->opt_elem, "PRISM") == 0)
	    disp_type = PRISM_DISP;
	else
	    disp_type = GRATING_DISP;

	/* Get offset of slit from that used to measure the dispersion
	   relation, from the APD again.
	*/
	if ((status = GetApOffset (sts, slit, disp->ref_aper, &delta)))
	    return (status);
	if (sts->verbose)
	    printf ("         Delta = %.6g arcsec.\n", delta);

        /* Calculate slit angle correction.
           This is the difference in the tangents of the slit angle
           and the nominal slit angle (= 0.315 degress).
        */
        dtor = acos(0)/90.;
        del_tan = tan(slit->angle*dtor) - tan(REF_ANGLE*dtor);

	/* Get spectrum trace. */
	if ((status = GetTrace (sts, coord_o->sporder, &trace)))
	    return (status);

	if (sts->verbose && sts->trace_rotation != 0.)
           printf ("         trace was rotated by = %.6g degree.\n",
                       sts->trace_rotation);

	/* Get incidence-angle correction coefficients. */
	if ((status = GetInang (sts, &sts->inangtab, coord_o->sporder, &iac)))
	    return (status);

	/* Was the switch reset for the current spectral order due to
	   pedigree = DUMMY in a table row or because no matching row
	   was found?
	*/
	if (sts->x2dcorr_o != PERFORM)
	    return (-1);

	/* Check whether the nominal location (crpix2) of the spectrum in the
	   output image maps back to a point that's outside the input image.
	   Skip this test, however, if a first-order grating was used.
	*/
	if (!sts->first_order) {
	    if ((status = CheckBoundary (sts, coord_o, &trace,
                                         in->sci.data.nx, in->sci.data.ny))) {
		if (status > 0)			/* a real error */
		    return (status);
		printf ("Info     Spectral order %d is outside the image.\n",
			coord_o->sporder);
		return (-1);
	    }
	}

	/* Initialize flags to say that warnings have not been printed. */
	warn1 = 0;  warn2 = 0;

	/* Fill each pixel of the output image. */

	for (j = 0;  j < out->sci.data.ny;  j++) {

	    /* approximate; in reference coordinates */
	    line0 = (j - coord_o->crpix[1]) / sts->ltm[1] +
			coord_o->a2center + sts->total_offset[1];

	    /* Interpolate to get the dispersion coefficients at line0. */
	    if ((status = InterpDisp (&disp, line0, &disp_y)))
		return (status);

            /* Correct each row for slit tilt */
            disp_y->coeff[0] += (j - coord_o->crpix[1])/sts->ltm[1] * del_tan;

	    /* Modify the dispersion coefficients using IAC. */
	    AdjustDisp (sts, disp_y, delta, &iac, &warn1, &warn2);

	    /* Interpolate to get the spectrum trace at line0. */
	    if ((status = InterpTrace (&trace, line0, &trace_y)))
		return (status);

	    for (i = 0;  i < out->sci.data.nx;  i++) {

		ox = (double) i;
		oy = (double) j;

		/* Convert to physical coordinates. */
		PixelPhys (coord_o, ox, oy, sts->hfactor, &wl, &dwl, &s, &ds);

		/* Get corresponding point (ix_r, iy_r) in reference
		   coordinates by evaluating the dispersion coefficients
		   and the spectrum trace; also get the derivatives.
		*/
		EvalDisp7 (disp_y, coord_o->sporder, wl, disp_type,
			&ix_r, &dix_r);

		/* ix_r --> ydispl and derivative */
		EvalTrace (trace_y, ix_r, &ydispl, &diy_r);
		iy_r = trace_y->a2center + ydispl;

/* add small-scale geo */
		if (sts->sgeocorr == PERFORM) {
		    SgeoMsg (sts, o_extver);
		    printf ("Warning  SGEOCORR not implemented yet.\n");
		    PrSwitch ("sgeocorr", SKIPPED);
		}

		/* Convert pixel location from reference coordinates to
		   input pixel units.  (can be out of bounds)
		*/
		ApplyLT (ix_r, iy_r, sts->ltm, sts->ltv, &ix, &iy);

/* compute jacobian */
		jacobian = 1.;	/* stub */

		/* Interpolate, checking for out of bounds. */
		sdqflags = 0;		/* instead of sts->sdqflags */
		Interp2D (in, sdqflags, ix, iy, jacobian,
				sts->err_algorithm,
				&Pix(out->sci.data,i,j),
				&Pix(out->err.data,i,j), &dq);
		DQSetPix (out->dq.data, i, j, dq);
	    }
	}

	/* Flag areas that are beyond the slit or behind an occulting bar. */
	if ((status = DataMasked (sts, coord_o, slit, out)))
	    return (status);

	FreeInang (&iac);
	FreeTrace (&trace_y);
	FreeTrace (&trace);
	FreeDisp (&disp_y);

	return (0);
}

/* This routine maps the nominal location of the spectrum in the output
   image back into the input image.  (The nominal location is crpix2.)
   The spectrum trace at that line is evaluated at the
   middle (in the first axis) of the input image or subarray to get the
   expected location of the spectrum in the input.  If that location is
   outside the boundaries of the input image, this function returns -1;
   if the location is inside the boundaries, zero is returned.

   The only true error condition is when InterpTrace returns an error.
*/

static int CheckBoundary (StisInfo7 *sts, CoordInfo *coord_o, SpTrace **trace,
		int nx, int ny) {

/* arguments:
StisInfo7 *sts      i: calibration switches and info
CoordInfo *coord_o  i: coordinate info
SpTrace **trace     i: list of spectrum traces
int nx, ny          i: size of input image
*/

	int status;

	double line0;		/* approx line number in input in ref coords */
	double ix, iy;		/* pixel coordinates in input image */
	double ix_r, iy_r;	/* coords in input, but in ref coords */
	double ydispl;		/* spectrum trace at current pixel */
	double diy_r;		/* returned by EvalTrace and ignored */

	SpTrace *trace_y;	/* spectrum trace interpolated at y */
	int InterpTrace (SpTrace **, double, SpTrace **);
	void FreeTrace (SpTrace **);

	/* approximate; in reference coordinates */
	line0 = coord_o->a2center + sts->total_offset[1];

	/* Interpolate to get the spectrum trace at line0. */
	trace_y = NULL;
	if ((status = InterpTrace (trace, line0, &trace_y)))
	    return (status);

	/* Get ydispl at the middle of the input image.
	   Note that it could be a subarray in X.
	*/
	ix = nx / 2.;
	ix_r = (ix - sts->ltv[0]) / sts->ltm[0];	/* ref coords */
	EvalTrace (trace_y, ix_r, &ydispl, &diy_r);
	iy_r = trace_y->a2center + ydispl;
	FreeTrace (&trace_y);

	/* Convert from reference coordinates to input pixel units. */
	ApplyLT (ix_r, iy_r, sts->ltm, sts->ltv, &ix, &iy);

	/* Out of bounds in cross dispersion direction? */
	if (iy < 0. || iy >= ny)
	    return (-1);
	else
	    return (0);
}

/* Convert pixel number in output to wavelength and position on slit;
   also compute increments (for dispaxis = 1).
   Note that we divide by the heliocentric correction factor to convert
   the wavelength from heliocentric to observed.
*/

static void PixelPhys (CoordInfo *coord_o,
		double ox, double oy, double hfactor,
		double *wl, double *dwl, double *s, double *ds) {

	*wl = ((ox - coord_o->crpix[0]) * coord_o->cdelt[0] +
			coord_o->crval[0]) / hfactor;

	*s = ((oy - coord_o->crpix[1]) * coord_o->cdelt[1] +
			coord_o->crval[1]);

	*dwl = coord_o->cdelt[0] / hfactor;
	*ds = coord_o->cdelt[1];
}

static void ApplyLT (double ix_r, double iy_r,
		double *ltm, double *ltv,
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


static void EvalTrace (SpTrace *trace_y, double ix_r,
		double *ydispl, double *diy_r) {

	int i;

	/* The nearest integer should be OK since the trace is small and
	   varies slowly.
	*/
	i = NINT (ix_r);

	if (i < 0 || i >= trace_y->nelem) {
	    *ydispl = 0.;
	    *diy_r = 0.;
	} else {
	    *ydispl = trace_y->a2displ[i];
	    *diy_r = 1.;		/* stub */
	}
}


/* This routine frees memory for the incidence-angle coefficients. */

static void FreeInang (InangInfo *iac) {

	if (iac->allocated) {
	    if (iac->ncoeff1 > 0)
		free (iac->coeff1);
	    if (iac->ncoeff2 > 0)
		free (iac->coeff2);
	    iac->allocated = 0;
	}
}

static void SgeoMsg (StisInfo7 *sts, int o_extver) {

	printf ("\n");
	PrSwitch ("sgeocorr", PERFORM);

	if (o_extver == 0) {
	    PrRefInfo ("sdstfile", sts->sdstfile.name,
			sts->sdstfile.pedigree, sts->sdstfile.descrip, "");
	}
}
