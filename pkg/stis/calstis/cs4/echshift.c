# include <stdio.h>

# include <hstio.h>
# include "../stis.h"
# include "calstis4.h"

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
*/

int EchShift (StisInfo4 *sts, SingleGroup *in,
		LampInfo *lamp, DispRelation *disp, SpTrace *trace,
		CmplxArray *clamp, double *shiftx, double *shifty) {

/* arguments:
StisInfo4 *sts           i: keyword values, etc.
SingleGroup *in          i: input data
LampInfo *lamp           i: template lamp spectrum
DispRelation *disp       i: dispersion relation
SpTrace *trace           i: list of spectral traces
CmplxArray *clamp       io: FT of template wavecal, created from lamp spectrum
double *shiftx, *shifty  o: shift from clamp to cwave, in each axis
*/

	CmplxArray cwave;	/* input SCI extension data (wavecal) */
	int status;
	int MakeTemplate (StisInfo4 *, LampInfo *, DispRelation *,
		SpTrace *, CmplxArray *);
	int ReadRealImage (SingleGroup *, int [], int [], CmplxArray *);
	int XCWavecal (CmplxArray *, CmplxArray *, double *, double *);

	InitCmplxArray (&cwave);

	/* Copy the input image (a subset of the SCI extension) into cwave. */
	if (status = ReadRealImage (in, sts->wl_sect1, sts->wl_sect2, &cwave))
	    return (status);

	/* Allocate memory, create the template image (stored in a
	   complex array), and take the forward Fourier transform.
	*/
	if (!clamp->allocated) {
	    InitCmplxArray (clamp);
	    if (status = AllocCmplxArray (clamp, cwave.nx, cwave.ny))
		return (status);
	    if (status = MakeTemplate (sts, lamp, disp, trace, clamp))
		return (status);
	    if (status = fft2d (clamp))
		return (status);
	}

	/* Take the forward Fourier transform of the input wavecal image. */
	if (status = fft2d (&cwave))
	    return (status);

	/* XCWavecal overwrites cwave, but it doesn't modify clamp. */
	if (status = XCWavecal (clamp, &cwave, shiftx, shifty))
	    return (status);

	/* Free memory for the wavecal. */
	FreeCmplxArray (&cwave);

	return (0);
}
