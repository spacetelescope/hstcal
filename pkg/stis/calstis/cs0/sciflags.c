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

   Various flags are set (e.g. sts->basic_2d), depending on which calstis
   module will perform the calibration.  Basic 2-D reduction, calstis1, is
   divided into two parts; dqicorr, atodcorr, and blevcorr constitute the
   calstis1a part, and these are included with the sci_basic_2d_a flag
   rather than with sci_basic_2d which is used for the other calstis1 steps.

   Phil Hodge, 1997 Dec 10:
	Replace missing with sciref in calling sequence;
	call GetNewRef instead of CheckRef.

   Phil Hodge, 1998 Jan 13:
	Check for APERTAB for PHOTCORR and PCTAB for X2DCORR.

   Phil Hodge, 1998 Mar 25:
	Don't set sts->sci_basic_2d just because DOPPCORR is PERFORM.

   Phil Hodge, 1998 June 3:
	Set sci_basic_2d_a instead of sci_basic_2d for BIASCORR.

   Phil Hodge, 1998 July 14:
	For x1dcorr or x2dcorr and MAMA detector, check MOFFTAB.

   Phil Hodge, 1998 Nov 20:
	Check WCPTAB and APDESTAB for wavecal processing.

   Phil Hodge, 2000 Jan 14:
	Check SDCTAB for x1dcorr.

   Phil Hodge, 2000 Aug 9:
	Remove test for MOFFTAB for x2dcorr.

   Phil Hodge, 2000 Oct 5:
	Add reference files for sc2dcorr.

   Phil Hodge, 2000 Nov 9:
	Remove test for MOFFTAB for x1dcorr.

   Phil Hodge, 2003 Jun 20:
	Remove the call to GetNewRef for wbiafile.

   Phil Hodge, 2004 Mar 1:
	Silently reset ctecorr if the detector is not the CCD.

   Phil Hodge, 2011 July 19:
	Replace PHOTTAB with IMPHTTAB.
*/

int SciFlags (StisInfo *sts, CalSwitch *sci_sw, Hdr *phdr,
		RefFileInfo *sciref) {

/* arguments:
StisInfo *sts         i: calibration flags and other info
CalSwitch *sci_sw     i: all calibration switches for science file
Hdr *phdr             i: primary header of science file
RefFileInfo *sciref  io: list of keyword,filename pairs
*/

	int status;

	int refimage_used = 0;	/* = 1 if do bias, dark, flat, or shadcorr */

	int GetNewRef (Hdr *, char *, RefFileInfo *);

	if (sts->detector == CCD_DETECTOR) {
	    if ((status = GetNewRef (phdr, "CCDTAB", sciref)))
		return (status);
	}

	/* Charge transfer inefficiency only affects the CCD. */
	if (sci_sw->ctecorr == PERFORM && sts->detector != CCD_DETECTOR)
	    sci_sw->ctecorr = OMIT;

	if (sci_sw->doppcorr == PERFORM) {
	    MAMASanity (sts->detector, "DOPPCORR");
	}

	if (sci_sw->lorscorr == PERFORM) {
	    sts->sci_basic_2d = PERFORM;
	    MAMASanity (sts->detector, "LORSCORR");
	}

	/* Note that we set sci_basic_2d_a rather than sci_basic_2d. */
	if (sci_sw->dqicorr == PERFORM) {
	    sts->sci_basic_2d_a = PERFORM;
	    if ((status = GetNewRef (phdr, "BPIXTAB", sciref)))
		return (status);
	}

	if (sci_sw->glincorr == PERFORM) {
	    sts->sci_basic_2d = PERFORM;
	    MAMASanity (sts->detector, "GLINCORR");
	}
	if (sci_sw->lflgcorr == PERFORM) {
	    sts->sci_basic_2d = PERFORM;
	    MAMASanity (sts->detector, "LFLGCORR");
	}
	if (sci_sw->glincorr == PERFORM || sci_sw->lflgcorr == PERFORM) {
	    if ((status = GetNewRef (phdr, "MLINTAB", sciref)))
		return (status);
	}

	if (sci_sw->crcorr == PERFORM) {
	    sts->sci_crcorr = PERFORM;
	    if (sts->detector != CCD_DETECTOR) {
		printf (
	"Warning  CRCORR will be omitted because detector is MAMA.\n");
		sts->sci_crcorr = OMIT;
	    }
	    if (sts->nimages < 2) {
		printf (
	"Warning  CRCORR will be omitted because there's only one imset.\n");
		sts->sci_crcorr = OMIT;
	    }
	    /* reference table is checked below, after getting rptcorr */
	}

	if (sci_sw->rptcorr == PERFORM) {
	    if (sts->nimages < 2) {
		printf (
	"Warning  RPTCORR will be omitted because there's only one imset.\n");
	    } else if (sts->detector == CCD_DETECTOR) {
		/* For the CCD, interpret rptcorr to mean crcorr. */
		sts->sci_crcorr = PERFORM;
		sts->sci_rptcorr = OMIT;
	    } else {
		sts->sci_rptcorr = PERFORM;
	    }
	}

	if (sts->sci_crcorr == PERFORM) {
	    if ((status = GetNewRef (phdr, "CRREJTAB", sciref)))
		return (status);
	}

	if (sci_sw->expscorr == PERFORM) {
	    sts->sci_expscorr = PERFORM;
	    if (sts->detector != CCD_DETECTOR) {
		printf (
	"Warning  EXPSCORR will be omitted because detector is MAMA.\n");
		sts->sci_expscorr = OMIT;
	    }
	}

	/* Note that we set sci_basic_2d_a rather than sci_basic_2d. */
	if (sci_sw->atodcorr == PERFORM) {
	    sts->sci_basic_2d_a = PERFORM;
	    CCDSanity (sts->detector, "ATODCORR");
	    if ((status = GetNewRef (phdr, "ATODTAB", sciref)))
		return (status);
	}

	/* Note that we set sci_basic_2d_a rather than sci_basic_2d. */
	if (sci_sw->biascorr == PERFORM) {
	    sts->sci_basic_2d_a = PERFORM;   refimage_used = 1;
	    CCDSanity (sts->detector, "BIASCORR");
	    if ((status = GetNewRef (phdr, "BIASFILE", sciref)))
		return (status);
	}

	if (sci_sw->darkcorr == PERFORM) {
	    sts->sci_basic_2d = PERFORM;   refimage_used = 1;
	    if ((status = GetNewRef (phdr, "DARKFILE", sciref)))
		return (status);
	}

	if (sci_sw->flatcorr == PERFORM) {
	    sts->sci_basic_2d = PERFORM;   refimage_used = 1;
	    if ((status = GetNewRef (phdr, "PFLTFILE", sciref)))
		return (status);
	    if ((status = GetNewRef (phdr, "DFLTFILE", sciref)))
		return (status);
	    if ((status = GetNewRef (phdr, "LFLTFILE", sciref)))
		return (status);
	}

	if (sci_sw->shadcorr == PERFORM) {
	    sts->sci_basic_2d = PERFORM;   refimage_used = 1;
	    CCDSanity (sts->detector, "SHADCORR");
	    if ((status = GetNewRef (phdr, "SHADFILE", sciref)))
		return (status);
	}

	/* Note that we set sci_basic_2d_a rather than sci_basic_2d. */
	if (sci_sw->blevcorr == PERFORM) {
	    sts->sci_basic_2d_a = PERFORM;
	    CCDSanity (sts->detector, "BLEVCORR");
	} else if (sci_sw->blevcorr == OMIT && refimage_used &&
			sts->detector == CCD_DETECTOR) {
	    /* Dark, flat, etc., assume the overscan has been subtracted. */
	    printf (
"Warning  For science file, should do BLEVCORR to remove overscan \\\n");
	    printf (
"         before doing other steps that use reference images.\n");
	}

	if (sci_sw->photcorr == PERFORM) {
	    sts->sci_basic_2d = PERFORM;
	    if (sts->obstype != IMAGING_TYPE)
		printf (
	"Warning  PHOTCORR = PERFORM, but OBSTYPE is not IMAGING.\n");
	    if ((status = GetNewRef (phdr, "IMPHTTAB", sciref)))
		return (status);
	    if ((status = GetNewRef (phdr, "APERTAB", sciref)))
		return (status);
	}

	/* These are the tables we need for wavecal processing. */
	if (sci_sw->wavecorr == PERFORM) {
	    sts->sci_wavecorr = PERFORM;
	    if ((status = GetNewRef (phdr, "WCPTAB", sciref)))
		return (status);
	    if ((status = GetNewRef (phdr, "LAMPTAB", sciref)))
		return (status);
	    if ((status = GetNewRef (phdr, "APDESTAB", sciref)))
		return (status);
	}

	if (sts->obstype == SPECTROSCOPIC_TYPE) {
	    if (sci_sw->x2dcorr == PERFORM) {
		sts->sci_2d_rect = PERFORM;
		if ((status = GetNewRef (phdr, "SDCTAB", sciref)))
		    return (status);
		if (sci_sw->fluxcorr == PERFORM) {
		    /* We look for PHOTTAB and APERTAB later. */
		    if ((status = GetNewRef (phdr, "PCTAB", sciref)))
			return (status);
		}
	    }
	} else {
	    if (sci_sw->geocorr == PERFORM) {
		sts->sci_geocorr = PERFORM;
		if ((status = GetNewRef (phdr, "IDCTAB", sciref)))
		    return (status);
	    }
	}

	if (sci_sw->x1dcorr == PERFORM) {
	    sts->sci_1d_extract = PERFORM;
	    if ((status = GetNewRef (phdr, "XTRACTAB", sciref)))
		return (status);
	    if ((status = GetNewRef (phdr, "SDCTAB", sciref)))
		return (status);
	    if (sci_sw->sc2dcorr == PERFORM) {
		if ((status = GetNewRef (phdr, "CDSTAB", sciref)))
		    return (status);
		if ((status = GetNewRef (phdr, "ECHSCTAB", sciref)))
		    return (status);
		if ((status = GetNewRef (phdr, "EXSTAB", sciref)))
		    return (status);
		if ((status = GetNewRef (phdr, "HALOTAB", sciref)))
		    return (status);
		if ((status = GetNewRef (phdr, "TELTAB", sciref)))
		    return (status);
		if ((status = GetNewRef (phdr, "RIPTAB", sciref)))
		    return (status);
		if ((status = GetNewRef (phdr, "SRWTAB", sciref)))
		    return (status);
	    }
	}

	/* These reference files are used by both 1-D extraction and
	   2-D spectral rectification.
	*/
	if (sts->sci_2d_rect == PERFORM || sts->sci_1d_extract == PERFORM) {

	    if ((status = GetNewRef (phdr, "SPTRCTAB", sciref)))
		return (status);

	    if (sci_sw->dispcorr == PERFORM || sts->sci_2d_rect == PERFORM) {
		if ((status = GetNewRef (phdr, "APDESTAB", sciref)))
		    return (status);
		if ((status = GetNewRef (phdr, "DISPTAB", sciref)))
		    return (status);
		if ((status = GetNewRef (phdr, "INANGTAB", sciref)))
		    return (status);
	    }

	    if (sci_sw->fluxcorr == PERFORM) {
		if ((status = GetNewRef (phdr, "PHOTTAB", sciref)))
		    return (status);
		if ((status = GetNewRef (phdr, "APERTAB", sciref)))
		    return (status);
	    }
	}

	/* Also used by geometric correction of images. */
	if (sts->sci_2d_rect == PERFORM || sts->sci_geocorr == PERFORM ||
	    sts->sci_1d_extract == PERFORM) {

	    if (sci_sw->sgeocorr == PERFORM) {
		if ((status = GetNewRef (phdr, "SDSTFILE", sciref)))
		    return (status);
	    }
	}

	/* Reasonableness checks. */
	if (sts->obstype == SPECTROSCOPIC_TYPE) {

	    if (sci_sw->geocorr == PERFORM) {
		printf (
	"Warning  GEOCORR = PERFORM, but OBSTYPE is not IMAGING.\n");
	    }

	} else {

	    if (sci_sw->x2dcorr == PERFORM) {
		printf (
	"Warning  X2DCORR = PERFORM, but OBSTYPE is not SPECTROSCOPIC.\n");
	    }
	}

	return (0);
}

static void MAMASanity (int detector, char *calswitch) {

	if (detector == CCD_DETECTOR)
	    printf ("Warning  %s = PERFORM, but detector is CCD.\n",
			calswitch);
}

static void CCDSanity (int detector, char *calswitch) {

	if (detector != CCD_DETECTOR)
	    printf ("Warning  %s = PERFORM, but detector is MAMA.\n",
			calswitch);
}
