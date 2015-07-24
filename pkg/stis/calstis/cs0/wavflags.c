# include <stdio.h>
# include "hstio.h"
# include "stis.h"
# include "calstis0.h"

static void CCDSanity (int, char *);
static void MAMASanity (int, char *);

/* This routine gets the names of all the reference files (images and
   tables) that are needed for calibrating the science file, based on
   the switches.  The files are not checked at this time to see whether
   they actually exist.

   Some sanity checking is done, and warning messages are printed if
   switch values don't make sense for the current detector.

   Note that in this routine, in contrast to SciFlags, we do not
   distinguish between wav_basic_2d and a 2d_a for calstis1a.

   Phil Hodge, 1997 Dec 10:
	Replace missing with wavref in calling sequence;
	call GetNewRef instead of CheckRef.

   Phil Hodge, 1998 Jan 13:
	Check for APERTAB for PHOTCORR and PCTAB for X2DCORR.

   Phil Hodge, 1998 Mar 9:
	Remove a space from several warning messages.

   Phil Hodge, 1998 Mar 25:
	Don't set sts->wav_basic_2d just because DOPPCORR is PERFORM.

   Phil Hodge, 1998 July 14:
	For x2dcorr and MAMA detector, check MOFFTAB.

   Phil Hodge, 1998 Nov 20:
	Check WCPTAB for wavecal processing.

   Phil Hodge, 2000 Jan 12:
	If the data were taken with an echelle grating, we don't check
	for the existence of SDCTAB and MOFFTAB, since they are only needed
	for the x2dcorr step for first-order data.

   Phil Hodge, 2000 Aug 9:
	Remove test for MOFFTAB.

   Phil Hodge, 2011 July 19:
	Replace PHOTTAB with IMPHTTAB (but irrelevant for a wavecal).
*/

int WavFlags (StisInfo *sts, CalSwitch *wav_sw, Hdr *phdr,
		RefFileInfo *wavref) {

/* arguments:
StisInfo *sts        i: calibration flags and other info
CalSwitch *wav_sw    i: all calibration switches (0 or 1) for wavecal
Hdr *phdr            i: primary header of wavecal
RefFileInfo wavref  io: list of keyword,filename pairs
*/

	int status;

	int refimage_used = 0;	/* = 1 if do bias, dark, flat, or shadcorr */

	int GetNewRef (Hdr *, char *, RefFileInfo *);

	if (sts->detector == CCD_DETECTOR) {
	    if ((status = GetNewRef (phdr, "CCDTAB", wavref)))
		return (status);
	}

	if (wav_sw->doppcorr == PERFORM) {
	    MAMASanity (sts->detector, "DOPPCORR");
	}

	if (wav_sw->lorscorr == PERFORM) {
	    sts->wav_basic_2d = PERFORM;
	    MAMASanity (sts->detector, "LORSCORR");
	}

	if (wav_sw->dqicorr == PERFORM) {
	    sts->wav_basic_2d = PERFORM;
	    if ((status = GetNewRef (phdr, "BPIXTAB", wavref)))
		return (status);
	}

	if (wav_sw->glincorr == PERFORM) {
	    sts->wav_basic_2d = PERFORM;
	    MAMASanity (sts->detector, "GLINCORR");
	}
	if (wav_sw->lflgcorr == PERFORM) {
	    sts->wav_basic_2d = PERFORM;
	    MAMASanity (sts->detector, "LFLGCORR");
	}
	if (wav_sw->glincorr == PERFORM || wav_sw->lflgcorr == PERFORM) {
	    if ((status = GetNewRef (phdr, "MLINTAB", wavref)))
		return (status);
	}

	if (wav_sw->atodcorr == PERFORM) {
	    sts->wav_basic_2d = PERFORM;
	    CCDSanity (sts->detector, "ATODCORR");
	    if ((status = GetNewRef (phdr, "ATODTAB", wavref)))
		return (status);
	}

	if (wav_sw->biascorr == PERFORM) {
	    sts->wav_basic_2d = PERFORM;   refimage_used = 1;
	    CCDSanity (sts->detector, "BIASCORR");
	    if ((status = GetNewRef (phdr, "BIASFILE", wavref)))
		return (status);
	}

	if (wav_sw->darkcorr == PERFORM) {
	    sts->wav_basic_2d = PERFORM;   refimage_used = 1;
	    if ((status = GetNewRef (phdr, "DARKFILE", wavref)))
		return (status);
	}

	if (wav_sw->flatcorr == PERFORM) {
	    sts->wav_basic_2d = PERFORM;   refimage_used = 1;
	    if ((status = GetNewRef (phdr, "PFLTFILE", wavref)))
		return (status);
	    if ((status = GetNewRef (phdr, "DFLTFILE", wavref)))
		return (status);
	    if ((status = GetNewRef (phdr, "LFLTFILE", wavref)))
		return (status);
	}

	if (wav_sw->shadcorr == PERFORM) {
	    sts->wav_basic_2d = PERFORM;   refimage_used = 1;
	    CCDSanity (sts->detector, "SHADCORR");
	    if ((status = GetNewRef (phdr, "SHADFILE", wavref)))
		return (status);
	}

	if (wav_sw->blevcorr == PERFORM) {
	    sts->wav_basic_2d = PERFORM;
	    CCDSanity (sts->detector, "BLEVCORR");
	} else if (wav_sw->blevcorr == OMIT && refimage_used &&
			sts->detector == CCD_DETECTOR) {
	    printf (
	"Warning  For wavecal, should do BLEVCORR to remove overscan \\\n");
	    printf (
	"         before doing other steps that use reference images.\n");
	}

	if (wav_sw->photcorr == PERFORM) {
	    sts->wav_basic_2d = PERFORM;
	    if (sts->obstype != IMAGING_TYPE)
		printf (
"Warning  PHOTCORR = PERFORM in wavecal, but OBSTYPE is not IMAGING.\n");
	    if ((status = GetNewRef (phdr, "IMPHTTAB", wavref)))
		return (status);
	    if ((status = GetNewRef (phdr, "APERTAB", wavref)))
		return (status);
	}

	/* We don't need to check the wavecorr switch!  These are the
	   tables we need for wavecal processing.
	*/
	if ((status = GetNewRef (phdr, "WCPTAB", wavref)))
	    return (status);
	if ((status = GetNewRef (phdr, "LAMPTAB", wavref)))
	    return (status);
	if ((status = GetNewRef (phdr, "APDESTAB", wavref)))
	    return (status);

	/* For echelle data, we need these for the wavecorr step;
	   for first-order data, we need them for the x2dcorr step.
	*/
	if ((status = GetNewRef (phdr, "DISPTAB", wavref)))
	    return (status);
	if ((status = GetNewRef (phdr, "INANGTAB", wavref)))
	    return (status);
	if ((status = GetNewRef (phdr, "SPTRCTAB", wavref)))
	    return (status);

	/* For first-order data, we will do x2dcorr regardless of the
	   wav_sw->x2dcorr value.
	*/
	if (!sts->echelle) {
	    if ((status = GetNewRef (phdr, "SDCTAB", wavref)))
		return (status);
	}

	if (wav_sw->sgeocorr == PERFORM) {
	    if ((status = GetNewRef (phdr, "SDSTFILE", wavref)))
		return (status);
	}
	if (wav_sw->fluxcorr == PERFORM) {	/* really shouldn't be set */
	    if ((status = GetNewRef (phdr, "PHOTTAB", wavref)))
		return (status);
	    if ((status = GetNewRef (phdr, "APERTAB", wavref)))
		return (status);
	    if ((status = GetNewRef (phdr, "PCTAB", wavref)))
		return (status);
	}

	return (0);
}

static void MAMASanity (int detector, char *calswitch) {

	if (detector == CCD_DETECTOR)
	    printf (
	"Warning  %s = PERFORM in wavecal, but detector is CCD.\n",
			calswitch);
}

static void CCDSanity (int detector, char *calswitch) {

	if (detector != CCD_DETECTOR)
	    printf (
	"Warning  %s = PERFORM in wavecal, but detector is MAMA.\n",
			calswitch);
}
