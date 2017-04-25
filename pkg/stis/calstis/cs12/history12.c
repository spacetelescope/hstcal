# include <stdio.h>
# include <stdlib.h>
# include <string.h>

# include "c_iraf.h"
# include "hstio.h"

# include "stis.h"
# include "calstis12.h"
# include "hstcalerr.h"
# include "stisdef.h"

/* This routine updates the WAVECORR switch to COMPLETE and appends a
   history record to the header.

   Phil Hodge, 1998 Mar 18:
	Remove references to the apdestab, and write wavecal name instead.

   Phil Hodge, 2000 June 16:
	Update the comment.
*/

int History12 (Hdr *phdr, char *wavecal_name) {

/* arguments:
Hdr *phdr            i: header to receive history records
char *wavecal_name   i: as above, but from wavecal row
*/

	int status;

	char *history;

	if ((history = calloc (STIS_LINE+1, sizeof(char))) == NULL)
	    return (OUT_OF_MEMORY);

	if ((status = Put_KeyS (phdr, "WAVECORR", "COMPLETE", "")))
	    return (status);
	addHistoryKw (phdr, "WAVECORR complete ...");
	if (hstio_err())
	    return (HEADER_PROBLEM);

	strcpy (history, "  wavecal = ");
	strcat (history, wavecal_name);
	addHistoryKw (phdr, history);
	if (hstio_err())
	    return (HEADER_PROBLEM);

	free (history);
	return (0);
}
