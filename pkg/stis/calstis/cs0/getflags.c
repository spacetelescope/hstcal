# include <stdio.h>
# include "c_iraf.h"
# include "hstio.h"

# include "stis.h"
# include "calstis0.h"
# include "hstcalerr.h"
# include "stisdef.h"

static int GetSw (Hdr *, char *, int *);
static int optionalSwitch (Hdr *, char *, int *);

/* This routine attempts to read the values of all known calibration
switches.  It returns these values as either OMIT or PERFORM (0 or 1)
in the CalSwitch struct.  OMIT and PERFORM are the only possible values
returned by this function.  Missing switches will be set to OMIT.  The
switch will be set to PERFORM only if the calibration switch value in
the header exists and is "PERFORM".

   Phil Hodge, 2000 Oct 5:
	Also get SC2DCORR.

   Phil Hodge, 2004 Mar 1:
	Also get CTECORR, using optionalSwitch.
*/

int GetFlags (CalSwitch *sw, Hdr *phdr) {

/* arguments:
CalSwitch *sw   o: values (0 or 1) of calibration switches
Hdr *phdr       i: primary header
*/

	int status;

	int use_default = 1;	/* use default if keyword is missing */

	if ((status = GetSw (phdr, "ATODCORR", &sw->atodcorr)))
	    return (status);
	if ((status = GetSw (phdr, "BACKCORR", &sw->backcorr)))
	    return (status);
	if ((status = GetSw (phdr, "BIASCORR", &sw->biascorr)))
	    return (status);
	if ((status = GetSw (phdr, "BLEVCORR", &sw->blevcorr)))
	    return (status);
	if ((status = GetSw (phdr, "CRCORR",   &sw->crcorr)))
	    return (status);
	if ((status = GetSw (phdr, "DARKCORR", &sw->darkcorr)))
	    return (status);
	if ((status = GetSw (phdr, "DISPCORR", &sw->dispcorr)))
	    return (status);
	if ((status = GetSw (phdr, "DOPPCORR", &sw->doppcorr)))
	    return (status);
	if ((status = GetSw (phdr, "DQICORR",  &sw->dqicorr)))
	    return (status);
	if ((status = GetSw (phdr, "EXPSCORR", &sw->expscorr)))
	    return (status);
	if ((status = GetSw (phdr, "FLATCORR", &sw->flatcorr)))
	    return (status);
	if ((status = GetSw (phdr, "FLUXCORR", &sw->fluxcorr)))
	    return (status);
	if ((status = optionalSwitch (phdr, "CTECORR", &sw->ctecorr)))
	    return (status);
	if ((status = GetSw (phdr, "GEOCORR",  &sw->geocorr)))
	    return (status);
	if ((status = GetSw (phdr, "GLINCORR", &sw->glincorr)))
	    return (status);
	if ((status = GetSw (phdr, "HELCORR",  &sw->helcorr)))
	    return (status);
	if ((status = GetSw (phdr, "LFLGCORR", &sw->lflgcorr)))
	    return (status);
	if ((status = GetSw (phdr, "LORSCORR", &sw->lorscorr)))
	    return (status);
	if ((status = GetSw (phdr, "PHOTCORR", &sw->photcorr)))
	    return (status);
	if ((status = GetSw (phdr, "RPTCORR",  &sw->rptcorr)))
	    return (status);
	if ((status = GetSw (phdr, "SC2DCORR", &sw->sc2dcorr)))
	    return (status);
	if ((status = GetSw (phdr, "SGEOCORR", &sw->sgeocorr)))
	    return (status);
	if ((status = GetSw (phdr, "SHADCORR", &sw->shadcorr)))
	    return (status);
	if ((status = GetSw (phdr, "WAVECORR", &sw->wavecorr)))
	    return (status);
	if ((status = GetSw (phdr, "X1DCORR",  &sw->x1dcorr)))
	    return (status);
	if ((status = GetSw (phdr, "X2DCORR",  &sw->x2dcorr)))
	    return (status);

	if ((status = Get_KeyI (phdr, "STATFLAG", use_default, 1, &sw->statcorr)))
	    return (status);

	return (0);
}

/* This routine just calls GetSwitch to get the value of the switch,
   but then the value returned is limited to OMIT (0) and PERFORM (1).
*/

static int GetSw (Hdr *phdr, char *calswitch, int *flag) {

/* arguments:
Hdr *phdr        i: primary header
char *calswitch  i: name of keyword (e.g. FLATCORR)
int *flag        o: value (0 or 1) of calibration switch
*/

	int status;

	if ((status = GetSwitch (phdr, calswitch, flag)))
	    return (status);

	if (*flag != PERFORM)
	    *flag = OMIT;

	return (0);
}

/* Check calibration switch, and set flag to PERFORM or OMIT.  This differs
   from GetSwitch in that the flag will be set to PERFORM if the keyword is
   missing, and it will be set to OMIT if the keyword is "complete".

   Phil Hodge, 2004 Mar 1:
	Copy from GetSwitch.
*/

static int optionalSwitch (Hdr *phdr, char *calswitch, int *flag) {

/* arguments:
Hdr *phdr         i: primary header
char *calswitch   i: name of keyword (e.g. CTECORR)
int *flag         o: value of switch:  PERFORM or OMIT
*/

	FitsKw key;		/* keyword location in header */
	char *word;		/* scratch space for header keyword value */
	int streq_ic (char *, char *);	/* strings equal? (case insensitive) */

	key = findKw (phdr, calswitch);
	if (key == NotFound) {
	    *flag = PERFORM;
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
	    *flag = OMIT;
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
