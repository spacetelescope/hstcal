# include <stdio.h>

# include "hstio.h"
# include "stis.h"
# include "calstis4.h"

static void RestoreDispCoeff (DispRelation *);
static int A4corrAdjustDisp (StisInfo4 *, SpTrace *, DispRelation *, double);

/* This routine finds the shifts for the current image set.

   The template image clamp is allocated and values are assigned when
   processing the first image set of the input wavecal file.  This
   image (actually its Fourier transform) is then reused for subsequent
   imsets, and it should be deallocated by the calling routine after
   all imsets have been processed.

   Phil Hodge, 2000 Jan 14
   Phil Hodge, 2000 July 21:
	Initialize cwave and clamp; cmplx.h is now in ../
   Phil Hodge, 2007 April 25:
	Do the initialization steps (allocating memory for clamp, etc) if
	clamp.allocated is not true, rather than testing whether extver == 1.
	Remove extver from the calling sequence.
   Phil Hodge, 2011 Jan 5:
	Remove clamp from the calling sequence of EchShift; it is now a
	local variable in this function.
	Move MOCAdjustDisp from wavecal.c to this file, and rename it to
	A4corrAdjustDisp.  Change this function so it does not call GetMOC
	and does not modify the coefficients in-place.
	Add 0 or 1 as a new argument (writedebug) to MakeTemplate.
*/

int EchShift (StisInfo4 *sts, SingleGroup *in,
		LampInfo *lamp, DispRelation *disp, SpTrace *trace,
		double *shiftx, double *shifty) {

/* arguments:
StisInfo4 *sts           i: keyword values, etc.
SingleGroup *in          i: input data
LampInfo *lamp           i: template lamp spectrum
DispRelation *disp       i: dispersion relation
SpTrace *trace           i: list of spectral traces
double *shiftx, *shifty  o: shift from clamp to cwave, in each axis
*/

	CmplxArray clamp;	/* template lamp image (complex array) */
	CmplxArray cwave;	/* input SCI extension data (wavecal) */
	CmplxArray crosscorr;	/* cross correlation of cwave and clamp */
	int status;
	int MakeTemplate (StisInfo4 *, LampInfo *, DispRelation *,
		SpTrace *, CmplxArray *, int);
	int ReadRealImage (SingleGroup *, int [], int [], CmplxArray *);
	int XCWavecal (CmplxArray *, CmplxArray *, CmplxArray *,
		double *, double *);

	InitCmplxArray (&clamp);	/* this sets clamp.allocated = 0 */
	InitCmplxArray (&cwave);
	InitCmplxArray (&crosscorr);

	/* Copy the input image (a subset of the SCI extension) into cwave. */
	if ((status = ReadRealImage (in, sts->wl_sect1, sts->wl_sect2, &cwave)))
	    return (status);

	/* Take the forward Fourier transform of the input wavecal image. */
	if ((status = fft2d (&cwave)))
	    return (status);

	if ((status = AllocCmplxArray (&crosscorr, cwave.nx, cwave.ny)))
	    return (status);
	if ((status = AllocCmplxArray (&clamp, cwave.nx, cwave.ny)))
	    return (status);

	/* Copy coeff_save to coeff, in case the coefficients were modified
	   when processing a previous imset.
	*/
	RestoreDispCoeff (disp);

	/* Create the template image and take the Fourier transform. */
	if ((status = MakeTemplate (sts, lamp, disp, trace, &clamp, 0)))
	    return (status);
	if ((status = fft2d (&clamp)))
	    return (status);

	/* Take the cross correlation of clamp with cwave. */
	if ((status = XCWavecal (&clamp, &cwave, &crosscorr, shiftx, shifty)))
	    return (status);

	/* Modify the coefficients to account for the skew due to shifty. */
	A4corrAdjustDisp (sts, trace, disp, *shifty);

	/* Now do the previous steps again with the modified dispersion
	   coefficients, to improve the shifts.
	   The final argument (writedebug = 1) means that a debug image
	   will be written if one was specified.
	*/
	if ((status = MakeTemplate (sts, lamp, disp, trace, &clamp, 1)))
	    return (status);
	if ((status = fft2d (&clamp)))
	    return (status);
	if ((status = XCWavecal (&clamp, &cwave, &crosscorr, shiftx, shifty)))
	    return (status);

	/* Free memory. */
	FreeCmplxArray (&clamp);
	FreeCmplxArray (&cwave);
	FreeCmplxArray (&crosscorr);

	return (0);
}

/* This copies the dispersion coefficients from coeff_save to coeff. */

static void RestoreDispCoeff (DispRelation *disp) {

	int i;

	for (i = 0;  i < MAX_DISP_COEFF;  i++)
	    disp->coeff[i] = disp->coeff_save[i];
}

/* This routine modifies the dispersion coefficients to account for the
   displacement of the image on the detector from the location that was used
   for measuring the coefficients.
   This only needs to be applied for echelle data.

   Note:  This uses sts->ltm[1], which is in an extension header, so this
   should be called within the loop over imsets.
*/

static int A4corrAdjustDisp (StisInfo4 *sts, SpTrace *trace,
		DispRelation *disp, double shifta2) {

	double ydiff;		/* Y position offset */
	SpTrace *trace_o;	/* spectral trace for current spectral order */
	double a2center;	/* from the trace for order mref */
	double r_shifta2;	/* shifta2 converted to reference pixel size */
	int foundit;		/* true if we have found mref in trace */
	/*int status;*/

	if (disp->a4corr == 0.)
	    return (0);		/* nothing further to do */

	/* Find the trace for spectral order mref, and its a2center. */
	trace_o = trace;
	foundit = 0;
	while (trace_o != NULL) {
	    if (trace_o->sporder == disp->mref) {
		foundit = 1;
		a2center = trace_o->a2center;
	    }
	    trace_o = trace_o->next;
	}

	if (!foundit) {
	    printf (
"Warning  Order %d not found in list of spectral traces; \\\n", disp->mref);
	    printf ("Warning  no A4CORR correction will be applied.\n");
	    return (0);
	}

	/* shifta2 should be accurate, since it's the shift in the dispersion
	   direction that is affected by the a4corr correction.

	   Here's what we're doing:
		ydiff = ypos - yref
		shifta2 = ypos - a2center
	   so:
		ydiff = shifta2 + a2center - yref
	*/
	r_shifta2 = shifta2 / sts->ltm[1];	/* reference pixels */
	ydiff = r_shifta2 + a2center - disp->yref;

	/* Adjust dispersion coefficients using a4corr. */
	disp->coeff[0] = disp->coeff_save[0] -
		(ydiff * sts->cenwave * disp->a4corr);
	disp->coeff[4] = disp->coeff_save[4] + ydiff * disp->a4corr;

	return (0);
}
