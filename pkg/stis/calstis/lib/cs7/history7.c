# include <stdio.h>
# include <string.h>

# include "c_iraf.h"
# include "hstio.h"

# include "stis.h"
# include "calstis7.h"
# include "hstcalerr.h"
# include "stisdef.h"

/* This routine appends history records to a header.  The records list
   each calibration step performed (or skipped), and for each step the
   reference files used are listed, along with their pedigree and descrip
   values.

   Phil Hodge, 1997 Nov 13:
	Include pctab in fluxcorr section.

   Phil Hodge, 1998 Jan 26:
	Change message about PCTAB not being included in FLUXCORR to
	not being included in DIFF2PT.

   Phil Hodge, 1998 Feb 4:
	Set DISPCORR to COMPLETE if x2dcorr was done.

   Phil Hodge, 2000 Aug 9:
	Remove references to MAMA offset table.

   Phil Hodge, 2000 Oct 6:
	Only write history for APDESTAB for spectroscopic data.

   Phil Hodge, 2006 Sept 21:
	If the input file is output from the wx2d task, print a message
	saying that the trace was applied earlier by wx2d.
*/

int History7 (StisInfo7 *sts, Hdr *phdr) {

/* arguments:
StisInfo7 *sts   i: calibration switches and info
Hdr *phdr        io: header to receive history records
*/

	int status;

	int logit;		/* true if we log history info */

	logit = 0;
	if (sts->obstype == SPECTROSCOPIC_TYPE) {
	    if (sts->x2dcorr == PERFORM) {
		Put_KeyS (phdr, "X2DCORR", "COMPLETE", "");
		Put_KeyS (phdr, "DISPCORR", "COMPLETE", "");
		addHistoryKw (phdr, "X2DCORR complete ...");
		logit = 1;
	    } else if (sts->x2dcorr == DUMMY) {
		Put_KeyS (phdr, "X2DCORR", "SKIPPED", "");
		addHistoryKw (phdr,
			"X2DCORR skipped due to dummy reference file ...");
		logit = 1;
	    }
	} else {
	    if (sts->x2dcorr == PERFORM) {
		Put_KeyS (phdr, "GEOCORR", "COMPLETE", "");
		addHistoryKw (phdr, "GEOCORR complete ...");
		logit = 1;
	    } else if (sts->x2dcorr == DUMMY) {
		Put_KeyS (phdr, "GEOCORR", "SKIPPED", "");
		addHistoryKw (phdr,
			"GEOCORR skipped due to dummy reference file ...");
		logit = 1;
	    }
	}
	if (logit) {
	    if (hstio_err())
		return (HEADER_PROBLEM);
	    if ((status = TabHistory (&sts->distntab, phdr)))
		return (status);
	    if (sts->obstype == SPECTROSCOPIC_TYPE) {
		if ((status = TabHistory (&sts->apdestab, phdr)))
		    return (status);
		if ((status = TabHistory (&sts->disptab, phdr)))
		    return (status);
		if ((status = TabHistory (&sts->inangtab, phdr)))
		    return (status);
		if ((status = TabHistory (&sts->sptrctab, phdr)))
		    return (status);
		if (sts->wx2dcorr == COMPLETE) {
		    addHistoryKw (phdr,
			"spectral trace was applied earlier, by wx2d");
		}
	    }
	}

	logit = 0;
	if (sts->sgeocorr == PERFORM) {
	    Put_KeyS (phdr, "SGEOCORR", "COMPLETE", "");
	    addHistoryKw (phdr, "SGEOCORR complete ...");
	    logit = 1;
	} else if (sts->sgeocorr == DUMMY) {
	    Put_KeyS (phdr, "SGEOCORR", "SKIPPED", "");
	    addHistoryKw (phdr,
			"SGEOCORR skipped due to dummy reference file ...");
	    logit = 1;
	}
	if (logit) {
	    if (hstio_err())
		return (HEADER_PROBLEM);
	    if ((status = ImgHistory (&sts->sdstfile, phdr)))
		return (status);
	}

	if (sts->heliocorr == PERFORM) {
	    Put_KeyS (phdr, "HELCORR", "COMPLETE", "");
	    addHistoryKw (phdr, "HELCORR complete");
	    if (hstio_err())
		return (HEADER_PROBLEM);
	}

	logit = 0;
	if (sts->fluxcorr == PERFORM) {
	    Put_KeyS (phdr, "FLUXCORR", "COMPLETE", "");
	    addHistoryKw (phdr, "FLUXCORR complete ...");
	    logit = 1;
	} else if (sts->fluxcorr == DUMMY) {
	    Put_KeyS (phdr, "FLUXCORR", "SKIPPED", "");
	    addHistoryKw (phdr,
			"FLUXCORR skipped due to dummy reference file ...");
	    logit = 1;
	}
	if (logit) {
	    if (hstio_err())
		return (HEADER_PROBLEM);
	    if ((status = TabHistory (&sts->phottab, phdr)))
		return (status);
	    if ((status = TabHistory (&sts->apertab, phdr)))
		return (status);

	    /* pctcorr is not independent; it's associated with fluxcorr. */
	    if (sts->pctcorr == PERFORM) {
		if ((status = TabHistory (&sts->pctab, phdr)))
		    return (status);
	    } else {
		addHistoryKw (phdr,
		"  Note:  PCTAB correction was not included in DIFF2PT");
		if (hstio_err())
		    return (HEADER_PROBLEM);
	    }

	    /* tdscorr is not independent; it's associated with fluxcorr. */
	    if (sts->tdscorr == PERFORM) {
		if ((status = TabHistory (&sts->tdstab, phdr)))
		    return (status);
	    } else {
		addHistoryKw (phdr,
		"  Note:  TDSTAB correction was not included in fluxcorr");
		if (hstio_err())
		    return (HEADER_PROBLEM);
	    }
	}

	if (sts->statcorr == PERFORM) {
	    addHistoryKw (phdr, "Statistics computed");
	    if (hstio_err())
		return (HEADER_PROBLEM);
	}

	return (0);
}
