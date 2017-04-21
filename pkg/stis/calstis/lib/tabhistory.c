# include <stdio.h>
# include <string.h>
# include "c_iraf.h"
# include "hstio.h"

# include "stis.h"
# include "err.h"
# include "stisdef.h"

/* This routine writes history records for a reference table, including
   the name of the file and the pedigree and descrip values.
   Note:  This is the same as imgHistory except for the name and the
   declaration of ref, and both descrip (from the header) and descrip2
   (from the row) can be written.
*/

int TabHistory (RefTab *ref, Hdr *phdr) {

/* arguments:
RefImage *ref     i: info about reference image
Hdr *phdr         io: header to receive history records
*/

	char history[STIS_LINE];

	strcpy (history, "  reference table ");
	strcat (history, ref->name);
	addHistoryKw (phdr, history);
	if (hstio_err())
	    return (HEADER_PROBLEM);

	if (ref->pedigree[0] != '\0') {
	    strcpy (history, "    ");
	    strcat (history, ref->pedigree);
	    addHistoryKw (phdr, history);
	    if (hstio_err())
		return (HEADER_PROBLEM);
	}

	if (ref->descrip[0] != '\0') {		/* descrip from the header */
	    strcpy (history, "    ");
	    strcat (history, ref->descrip);
	    addHistoryKw (phdr, history);
	    if (hstio_err())
		return (HEADER_PROBLEM);
	}

	if (ref->descrip2[0] != '\0') {		/* descrip from the row */
	    strcpy (history, "    ");
	    strcat (history, ref->descrip2);
	    addHistoryKw (phdr, history);
	    if (hstio_err())
		return (HEADER_PROBLEM);
	}

	return (0);
}
