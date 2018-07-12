# include <stdio.h>
# include <stdlib.h>	/* calloc */
# include <string.h>

# include "c_iraf.h"
# include "hstio.h"

# include "stis.h"
# include "calstis1.h"
# include "hstcalerr.h"
# include "stisdef.h"

static int UpdateSwitch (char *, int, Hdr *, int *);

/* These routines check whether pedigree indicates that the calibration
   file is dummy, and log history info accordingly.  From one to three
   history records will be logged:  the status of the calibration (i.e.
   "complete" or "skipped"), the value of the pedigree keyword, and the
   value of the descrip keyword.  The calibration switch (e.g. FLATCORR)
   in the header will be updated from PERFORM to either COMPLETE or SKIPPED.

   Note that since this info is written to the primary header, these
   routines should only be called for EXTVER = 1.

   Phil Hodge, 1997 Oct 28:
	Change blevHistory to include driftcorr.

   Phil Hodge, 1997 Nov 13:
	Change photHistory to include APERTAB.

   Phil Hodge, 1998 Jan 16:
	Write "uncertainty" instead of "error" in noiseHistory.

   Phil Hodge, 1998 July 30:
	In UpdateSwitch, add a section for SKIPPED.

   Paul Barrett, 2003 Sep 25:
        Added TDSTAB to photHistory.

   Phil Hodge, 2011 May 9:
	In photHistory, don't check filtcorr or tdscorr.

   Phil Hodge, 2011 Nov 17:
	In photHistory, check tdscorr.

   Phil Hodge, 2012 Oct 15:
	In darkHistory, include history for the tdctab if detector is NUV-MAMA.
*/

int atodHistory (StisInfo1 *sts, Hdr *phdr) {

	int status;

	int logit;			/* true if we should log file name */

	if (OmitStep (sts->atodcorr))		/* nothing to do */
	    return (0);

	if ((status = UpdateSwitch ("ATODCORR", sts->atodcorr, phdr, &logit)))
	    return (status);

	/* Write history records for the A-to-D table. */
	if (logit) {
	    if (TabHistory (&sts->atod, phdr))
		return (status);
	}

	return (0);
}

int biasHistory (StisInfo1 *sts, Hdr *phdr) {

	int status;

	int logit;			/* true if we should log file name */

	if (OmitStep (sts->biascorr))
	    return (0);

	if ((status = UpdateSwitch ("BIASCORR", sts->biascorr, phdr, &logit)))
	    return (status);

	/* Write history records for the bias image. */
	if (logit) {
	    if ((status = ImgHistory (&sts->bias, phdr)))
		return (status);
	}

	return (0);
}

int blevHistory (StisInfo1 *sts, Hdr *phdr, int done, int driftcorr) {

	int status;

	if (OmitStep (sts->blevcorr))
	    return (0);

	if ((status = Put_KeyS (phdr, "BLEVCORR", "COMPLETE", "")))
	    return (status);
	if (done) {
	    addHistoryKw (phdr,
	"BLEVCORR complete; bias level from overscan was subtracted.");
	} else {
	    addHistoryKw (phdr,
	"BLEVCORR complete, but default bias level was subtracted.");
	}
	if (hstio_err())
	    return (HEADER_PROBLEM);

	if (driftcorr) {
	    addHistoryKw (phdr,
	"BLEVCORR includes correction for drift along lines.");
	} else {
	    addHistoryKw (phdr,
	"BLEVCORR does not include correction for drift along lines.");
	}
	if (hstio_err())
	    return (HEADER_PROBLEM);

	return (0);
}

/* Info about CCD parameters table (no specific calibration step). */

int CCDHistory (StisInfo1 *sts, Hdr *phdr) {

	int status;

	addHistoryKw (phdr, "CCD parameters table ...");
	if ((status = TabHistory (&sts->ccdpar, phdr)))
	    return (status);

	return (0);
}

int darkHistory (StisInfo1 *sts, Hdr *phdr) {

	int status;
	int logit;			/* true if we should log file name */

	if (OmitStep (sts->darkcorr))
	    return (0);

	if ((status = UpdateSwitch ("DARKCORR", sts->darkcorr, phdr, &logit)))
	    return (status);

	if (logit) {
	    if ((status = ImgHistory (&sts->dark, phdr)))
		return (status);
	    if (sts->detector == NUV_MAMA_DETECTOR) {
		if ((status = TabHistory (&sts->tdctab, phdr)))
		    return (status);
	    }
	}

	return (0);
}

int dqiHistory (StisInfo1 *sts, Hdr *phdr) {

	int status;
	int logit;			/* true if we should log file name */
	int flag;

	if (OmitStep (sts->dqicorr))
	    return (0);

	if (sts->bpix.exists == EXISTS_YES) {
	    if (sts->bpix.goodPedigree == GOOD_PEDIGREE)
		flag = PERFORM;
	    else if (sts->bpix.goodPedigree == DUMMY_PEDIGREE)
		flag = DUMMY;
	    else
		flag = IGNORED;
	} else {
	    flag = PERFORM;
	}

	if ((status = UpdateSwitch ("DQICORR", flag, phdr, &logit)))
	    return (status);

	if (logit) {
	    addHistoryKw (phdr, "  values checked for saturation");
	    if (sts->bpix.exists == EXISTS_YES) {
		addHistoryKw (phdr, "  DQ array initialized ...");
		if ((status = TabHistory (&sts->bpix, phdr)))
		    return (status);
	    }
	}

	return (0);
}

int flatHistory (StisInfo1 *sts, Hdr *phdr) {

	int status;

	int logit;			/* true if we should log file name */

	if (OmitStep (sts->flatcorr))
	    return (0);

	if ((status = UpdateSwitch ("FLATCORR", sts->flatcorr, phdr, &logit)))
	    return (status);

	if (logit) {

	    if (GotFileName (sts->pflt.name)) {
		if ((status = ImgHistory (&sts->pflt, phdr))) /* pixel-to-pixel */
		    return (status);
	    }

	    if (GotFileName (sts->dflt.name)) {
		if ((status = ImgHistory (&sts->dflt, phdr))) /* delta flat */
		    return (status);
	    }

	    if (GotFileName (sts->lflt.name)) {
		if ((status = ImgHistory (&sts->lflt, phdr))) /* low-order flat */
		    return (status);
	    }
	}

	return (0);
}

/* There are no reference files for conversion to low-res. */

int loresHistory (StisInfo1 *sts, Hdr *phdr, int done) {

	int status;

	if (OmitStep (sts->lorscorr))
	    return (0);

	if (done) {
	    if ((status = Put_KeyS (phdr, "LORSCORR", "COMPLETE", "")))
		return (status);
	    addHistoryKw (phdr,
	"LORSCORR complete; image converted from high-res to low-res.");
	} else {
	    if ((status = Put_KeyS (phdr, "LORSCORR", "SKIPPED", "")))
		return (status);
	    addHistoryKw (phdr,
	"LORSCORR skipped; image is already low-res.\n");
	}
	if (hstio_err())
	    return (HEADER_PROBLEM);

	return (0);
}

/* There are no reference images for assigning initial error values. */

int noiseHistory (Hdr *phdr) {

	addHistoryKw (phdr, "Uncertainty array initialized.");
	if (hstio_err())
	    return (HEADER_PROBLEM);

	return (0);
}

int nonlinHistory (StisInfo1 *sts, Hdr *phdr) {

	int status;
	int logit;			/* true if we should log file name */

	logit = 0;

	if (!OmitStep (sts->glincorr))
	    if ((status = UpdateSwitch ("GLINCORR", sts->glincorr, phdr, &logit)))
		return (status);

	if (!OmitStep (sts->lflgcorr))
	    if ((status = UpdateSwitch ("LFLGCORR", sts->lflgcorr, phdr, &logit)))
		return (status);

	if (logit)
	    if ((status = TabHistory (&sts->mlin, phdr)))
		return (status);

	return (0);
}

int shadHistory (StisInfo1 *sts, Hdr *phdr) {

	int status;
	int logit;			/* true if we should log file name */

	if (OmitStep (sts->shadcorr))
	    return (0);

	if ((status = UpdateSwitch ("SHADCORR", sts->shadcorr, phdr, &logit)))
	    return (status);

	if (logit) {
	    if ((status = ImgHistory (&sts->shad, phdr)))
		return (status);
	}

	return (0);
}

/* Photometry. */

int photHistory (StisInfo1 *sts, Hdr *phdr) {

	int status;
	int logit;			/* true if we should log file name */

	if (OmitStep (sts->photcorr))
	    return (0);

	if ((status = UpdateSwitch ("PHOTCORR", sts->photcorr, phdr, &logit)))
	    return (status);

	if (logit) {
	    if ((status = TabHistory (&sts->phot, phdr)))
		return (status);
	    if (sts->tdscorr == PERFORM) {
		if ((status = TabHistory (&sts->tdstab, phdr)))
		    return (status);
	    } else {
		addHistoryKw (phdr,
		"  Note:  Time-dependent sensitivity was not included "
                              "with PHOTCORR");
	    }
	}

	return (0);
}

/* Update the calibration switch in the primary header to COMPLETE
   or SKIPPED.
*/

static int UpdateSwitch (char *calSwitch, int flag, Hdr *phdr, int *logit) {

/* arguments:
char *calSwitch   i: name of calibration switch
int flag          i: value of calibration switch
Hdr *phdr         o: primary header
int *logit        o: true if we should log reference file names
*/

	int status;

	char *history;

	if ((history = (char *) calloc (STIS_LINE+1, sizeof (char))) == NULL)
	    return (OUT_OF_MEMORY);

	strcpy (history, calSwitch);

	*logit = 0;
	if (flag == PERFORM) {
	    if ((status = Put_KeyS (phdr, calSwitch, "COMPLETE", "")))
		return (status);
	    strcat (history, " complete ...");
	    addHistoryKw (phdr, history);
	    if (hstio_err())
		return (HEADER_PROBLEM);
	    *logit = 1;
	} else if (flag == DUMMY) {
	    if ((status = Put_KeyS (phdr, calSwitch, "SKIPPED", "")))
		return (status);
	    strcat (history, " skipped due to dummy reference file ...");
	    addHistoryKw (phdr,history);
	    if (hstio_err())
		return (HEADER_PROBLEM);
	    *logit = 1;
	} else if (flag == IGNORED) {
	    if ((status = Put_KeyS (phdr, calSwitch, "SKIPPED", "")))
		return (status);
	    strcat (history, " not performed ...");  /* for some other reason */
	    addHistoryKw (phdr,history);
	    if (hstio_err())
		return (HEADER_PROBLEM);
	    *logit = 1;
	} else if (flag == SKIPPED) {
	    if ((status = Put_KeyS (phdr, calSwitch, "SKIPPED", "")))
		return (status);
	    *logit = 0;		/* don't log this */
	}

	free (history);

	return (0);
}
