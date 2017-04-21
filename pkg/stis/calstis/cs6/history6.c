# include <stdio.h>
# include <stdlib.h>
# include <string.h>

# include "xtables.h"
# include "hstio.h"

# include "stis.h"
# include "stisdef.h"
# include "calstis6.h"
# include "hstcalerr.h"

static int outHistory (int, char *, Hdr *);

/*
   Append history records to a header. The records list each calibration
   step performed (or skipped), and for each step the reference files used
   are listed, along with their pedigree and descrip values.

   Revision history:
   ----------------
   20 Feb 97  -  Adapted from similar routine in calstis7 (I.Busko)
   10 Apr 97  -  Changes after code review (IB):
                 - added function outHistory, removed logit flag.
   25 Apr 97  -  Removed esptab (IB)
   28 Jan 98  -  Add PCTAB support (IB)
   07 Aug 00  -  Add IDT algorithm support (IB)
   19 Jul 04  -  Add records for TDS and CTE in the FLUXCORR section (PEH)
   08 Apr 05  -  Add records for GAC in the FLUXCORR section (PEH)
   12 Sep 05  -  Add a record for xoffset (PEH)
*/

int History6 (StisInfo6 *sts, Hdr *phdr, int idt) {

/* arguments:
StisInfo6 *sts   i: calibration switches and info
Hdr *phdr        io: header to receive history records
*/

	char *newString;
	int status;

	/* This is used by standard extraction. */

	if (!idt) {

	    if (outHistory (sts->x1d, "X1DCORR", phdr)) {
	        if (hstio_err())
	            return (HEADER_PROBLEM);
	        if ((status = TabHistory (&sts->sptrctab, phdr)))
		    return (status);
	        if ((status = TabHistory (&sts->xtrctab, phdr)))
		    return (status);
	    }

	    if (outHistory (sts->backcorr, "BACKCORR", phdr)) {
	        if (hstio_err())
		    return (HEADER_PROBLEM);
	        if ((status = TabHistory (&sts->xtrctab, phdr)))
		    return (status);
	    }

	    if (outHistory (sts->dispcorr, "DISPCORR", phdr)) {
	        if (hstio_err())
		    return (HEADER_PROBLEM);
                if ((status = TabHistory (&sts->disptab, phdr)))
                    return (status);
                if ((status = TabHistory (&sts->inangtab, phdr)))
                    return (status);
                if ((status = TabHistory (&sts->apdestab, phdr)))
                    return (status);
	    }

	    if (outHistory (sts->sgeocorr, "SGEOCORR", phdr)) {
	        if (hstio_err())
		    return (HEADER_PROBLEM);
	        if ((status = ImgHistory (&sts->sdstfile, phdr)))
		    return (status);
	    }

	    if (sts->heliocorr == PERFORM) {
	        addHistoryKw (phdr,
			    "HELCORR performed (no reference file needed)");
	        if (hstio_err())
		     return (HEADER_PROBLEM);
	    }

	    if (outHistory (sts->fluxcorr, "FLUXCORR", phdr)) {
	        if (hstio_err())
		    return (HEADER_PROBLEM);
	        if ((status = TabHistory (&sts->phottab, phdr)))
		    return (status);
	        if ((status = TabHistory (&sts->pctab, phdr)))
		    return (status);
	        if ((status = TabHistory (&sts->apertab, phdr)))
		    return (status);
		if (sts->gaccorr == PERFORM) {
		    if ((status = TabHistory (&sts->gactab, phdr)))
			return (status);
		}
		if (sts->ctecorr == PERFORM) {
		    if ((status = TabHistory (&sts->ccdtab, phdr)))
			return (status);
		}
		if (sts->tdscorr == PERFORM) {
		    if ((status = TabHistory (&sts->tdstab, phdr)))
			return (status);
		}
	    }
	    if (sts->xoffset != 0.) {
		newString = (char *) calloc (STIS_LINE, sizeof (char));
		sprintf (newString,
	"Offset of %g low-res pixels added in dispersion direction",
			sts->xoffset);
		addHistoryKw (phdr, newString);
		if (hstio_err())
		    return (HEADER_PROBLEM);
	    }
	} else {
	    if (outHistory (sts->idt, "SC2DCORR", phdr)) {
	        if (hstio_err())
		    return (HEADER_PROBLEM);
	        if ((status = TabHistory (&sts->echsctab, phdr)))
		    return (status);
	        if ((status = TabHistory (&sts->exstab, phdr)))
		    return (status);
	        if ((status = TabHistory (&sts->cdstab, phdr)))
		    return (status);
	        if ((status = TabHistory (&sts->riptab, phdr)))
		    return (status);
	        if ((status = TabHistory (&sts->srwtab, phdr)))
		    return (status);
	        if ((status = TabHistory (&sts->halotab, phdr)))
		    return (status);
	        if ((status = TabHistory (&sts->psftab, phdr)))
		    return (status);
	    }
	}

	return (0);
}


static int outHistory (int cal_switch, char *switch_name, Hdr *phdr) {

	char *newString;

	newString = (char *) calloc (STIS_LINE, sizeof (char));

	if (cal_switch == PERFORM) {
	    strcpy (newString, switch_name);
	    strcat (newString, " performed ...");
	    addHistoryKw (phdr, newString);
	    free (newString);
	    return (1);

	} else if (cal_switch == DUMMY) {
	    strcpy (newString, switch_name);
	    strcat (newString, " skipped due to dummy reference file ...");
	    addHistoryKw (phdr, newString);
	    free (newString);
	    return (1);
	}

	free (newString);
	return (0);
}
