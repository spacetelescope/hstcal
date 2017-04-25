# include <stdio.h>
# include <stdlib.h>
# include "c_iraf.h"
# include "hstio.h"

# include "stis.h"
# include "hstcalerr.h"
# include "stisdef.h"

/* Check calibration switch, and set flag to PERFORM, COMPLETE or
   OMIT.  The other options (SKIPPED, DUMMY) are not set here.
   It is not an error for a switch to be absent; the step will just
   be set to OMIT.

   Phil Hodge, 1998 Mar 9:
	Make the test on the switch value case insensitive, and return
	an error if the value is not one of PERFORM, COMPLETE, SKIPPED,
	or OMIT.
*/

int GetSwitch (Hdr *phdr, char *calswitch, int *flag) {

/* arguments:
Hdr *phdr         i: primary header
char *calswitch   i: name of keyword (e.g. FLATCORR)
int *flag         o: value of switch:  PERFORM, OMIT, or COMPLETE
*/

	FitsKw key;		/* keyword location in header */
	char *word;		/* scratch space for header keyword value */
	int streq_ic (char *, char *);	/* strings equal? (case insensitive) */

	key = findKw (phdr, calswitch);
	if (key == NotFound) {
	    *flag = OMIT;
	    return (0);
	}

	if ((word = (char *) calloc (STIS_FNAME+1, sizeof(char))) == NULL)
	    return (OUT_OF_MEMORY);

	getStringKw (key, word, STIS_FNAME);
	if (hstio_err()) {
	    free (word);
	    printf ("ERROR    Error getting keyword `%s'.\n", calswitch);
	    return (HEADER_PROBLEM);
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
	    printf ("ERROR    Keyword %s = %s is invalid.\n",
			calswitch, word);
	    free (word);
	    return (HEADER_PROBLEM);
	}

	free (word);

	return (0);
}
