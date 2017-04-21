# include <stdio.h>
# include "c_iraf.h"
# include "hstio.h"

# include "stis.h"
# include "calstis4.h"
# include "err.h"
# include "stisdef.h"

/* This routine updates the WAVECORR switch to COMPLETE (or SKIPPED,
   if there was a dummy reference table) and appends history records
   to the header.  The history records list the reference files used,
   along with their pedigree and descrip values.

   Phil Hodge, 1998 Dec 11:
	Include WCPTAB table.

   Phil Hodge, 2000 Jan 5:
	Write history for tables used for echelle.

   Phil Hodge, 2001 Mar 7:
	Write history for tables used for prism (same as for echelle).
*/

int History4 (StisInfo4 *sts, Hdr *phdr) {

/* arguments:
StisInfo4 *sts   i: calibration switches and info
Hdr *phdr        io: header to receive history records
*/

	int status;

	int logit;		/* true if we log history info */

	logit = 0;
	if (sts->wavecorr == PERFORM) {
	    if ((status = Put_KeyS (phdr, "WAVECORR", "COMPLETE", "")))
		return (status);
	    addHistoryKw (phdr, "WAVECORR complete ...");
	    logit = 1;
	} else if (sts->wavecorr == DUMMY) {
	    if ((status = Put_KeyS (phdr, "WAVECORR", "SKIPPED", "")))
		return (status);
	    addHistoryKw (phdr,
			"WAVECORR skipped due to dummy reference file ...");
	    logit = 1;
	}
	if (logit) {
	    if (hstio_err())
		return (HEADER_PROBLEM);
	    if (sts->wcptab.exists == EXISTS_YES) {
		if ((status = TabHistory (&sts->wcptab, phdr)))
		    return (status);
	    }
	    if ((status = TabHistory (&sts->lamptab, phdr)))
		return (status);
	    if ((status = TabHistory (&sts->apdestab, phdr)))
		return (status);
	    if (sts->disp_type == ECHELLE_DISP ||
		sts->disp_type == PRISM_DISP) {
		if ((status = TabHistory (&sts->disptab, phdr)))
		    return (status);
		if ((status = TabHistory (&sts->inangtab, phdr)))
		    return (status);
		if ((status = TabHistory (&sts->sptrctab, phdr)))
		    return (status);
	    }
	    if (sts->disp_type == PRISM_DISP) {
		if ((status = TabHistory (&sts->sdctab, phdr)))
		    return (status);
	    }
	}

	return (0);
}
