# include <stdio.h>
# include <stdlib.h>
# include "hstio.h"
# include "acs.h"
# include "hstcalerr.h"

/* Check calibration switch, and set flag to PERFORM, COMPLETE or
   OMIT.  The other options (SKIPPED, DUMMY) are not set here.
   It is not an error for a switch to be absent; the step will just
   be set to OMIT.

   Phil Hodge, 1998 Mar 9:
	Make the test on the switch value case insensitive, and return
	an error if the value is not one of PERFORM, COMPLETE, SKIPPED,
	or OMIT.
   Warren Hack, 1998 May 26:
	Revised for ACS.
*/

int GetSwitch (Hdr *phdr, char *calswitch, int *flag) {

/* arguments:
Hdr *phdr         i: primary header
char *calswitch   i: name of keyword (e.g. FLATCORR)
int *flag         o: value of switch:  PERFORM, OMIT, or COMPLETE
*/

	extern int status;

	FitsKw key;		/* keyword location in header */
	char *word;		/* scratch space for header keyword value */
	int streq_ic (char *, char *);	/* strings equal? (case insensitive) */

	key = findKw (phdr, calswitch);
	if (key == NotFound) {
	    *flag = OMIT;
	    return (status);
	}

	if ((word = (char *) calloc (ACS_FNAME+1, sizeof(char))) == NULL)
	    return (status = OUT_OF_MEMORY);

	getStringKw (key, word, ACS_FNAME);
	if (hstio_err()) {
	    free (word);
	    sprintf (MsgText, "Error getting keyword `%s'.", calswitch);
	    trlerror (MsgText);
		return (status = HEADER_PROBLEM);
	}

	if (streq_ic (word, "perform")) {
	    *flag = PERFORM;
	} else if (streq_ic (word, "complete")) {
	    *flag = COMPLETE;
	} else if (streq_ic (word, "skipped")) {
	    *flag = OMIT;
	} else if (streq_ic (word, "omit")) {
	    *flag = OMIT;
	} else {
	    *flag = OMIT;
	    sprintf (MsgText, "Keyword %s = %s is invalid.", calswitch, word);
	    trlerror (MsgText);
		free (word);
	    return (status = HEADER_PROBLEM);
	}

	free (word);

	return (status);
}
