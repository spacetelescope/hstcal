# include <stdio.h>
# include "hstio.h"
# include "acs.h"
# include "acscorr.h"
# include "err.h"
# include "calacs.h"

static int GetSw (Hdr *, char *, int *);

/* This routine attempts to read the values of all known calibration
switches.  It returns these values as either OMIT or PERFORM (0 or 1)
in the CalSwitch struct.  OMIT and PERFORM are the only possible values
returned by this function.  Missing switches will be set to OMIT.  The
switch will be set to PERFORM only if the calibration switch value in
the header exists and is "PERFORM".

5 Mar 99 WJH - Revised to produce error when both RPTCORR and CRCORR
                are set to PERFORM.
11 Sept 2000  WJH - Revised to add POSTFLASH processing keyword.
23 Apr 2001   WJH - Revised to support both 'POSTFLSH' and 'FLSHCORR'
                    for STSDAS distribution and pipeline support.
16 Jan 2002   WJH - Corrected logic error in defaulting to FLSHCORR keyword.
17 Apr 2002   WJH - Removed reference to 'STATFLAG'.  Always perform 'doStat'.
*/

int GetFlags (CalSwitch *sw, Hdr *phdr) {

/* arguments:
CalSwitch *sw   o: values (0 or 1) of calibration switches
Hdr *phdr       i: primary header
*/

	extern int status;
	FitsKw key;		/* keyword location in header */
    char flashkey[ACS_CBUF];

	int GetKeyInt (Hdr *, char *, int, int, int *);

	if (GetSw (phdr, "ATODCORR", &sw->atodcorr))
	    return (status);
	if (GetSw (phdr, "BIASCORR", &sw->biascorr))
	    return (status);
	if (GetSw (phdr, "BLEVCORR", &sw->blevcorr))
	    return (status);
	if (GetSw (phdr, "CRCORR",   &sw->crcorr))
	    return (status);
	if (GetSw (phdr, "DARKCORR", &sw->darkcorr))
	    return (status);
	if (GetSw (phdr, "DQICORR",  &sw->dqicorr))
	    return (status);
	if (GetSw (phdr, "FLATCORR", &sw->flatcorr))
	    return (status);
//if (GetSw (phdr, "PCTECORR", &sw->pctecorr))
//    return (status);
  if (GetSw (phdr, "PCTECORR", &sw->pctecorr) == HEADER_PROBLEM)
      sw->pctecorr = OMIT;


    sprintf(flashkey,"FLSHCORR");
	key = findKw (phdr, flashkey);
	if (key == NotFound) {
        sprintf(MsgText,"FLSHCORR keyword not found...");
        trlwarn(MsgText);

	    key = findKw (phdr, "POSTFLSH");
        if (key != NotFound) {
            sprintf(flashkey,"POSTFLSH");

            /* Now warn the user to change keyword name to FLSHCORR. */
            sprintf(MsgText,"Using old keyword POSTFLSH!");
            trlwarn(MsgText);
            sprintf(MsgText,"Please rename keyword to FLSHCORR in header.");
            trlwarn(MsgText);
        }
	}

	if (GetSw (phdr, flashkey, &sw->flashcorr))
        return(status);

	if (GetSw (phdr, "GLINCORR", &sw->glincorr))
	    return (status);
	if (GetSw (phdr, "LFLGCORR", &sw->lflgcorr))
	    return (status);
	if (GetSw (phdr, "PHOTCORR", &sw->photcorr))
	    return (status);
	if (GetSw (phdr, "RPTCORR",  &sw->rptcorr))
	    return (status);
	if (GetSw (phdr, "SHADCORR", &sw->shadcorr))
	    return (status);
	if (GetSw (phdr, "EXPSCORR", &sw->expscorr))
	    return (status);

    /* Check to insure that only one of these switches are
        set, since they are exclusive options.
    */
    if (sw->rptcorr == PERFORM && sw->crcorr == PERFORM) {
            trlerror("RPTCORR and CRCORR both set to PERFORM!");
            trlerror("    One switch needs to be set to OMIT!");
            status = HEADER_PROBLEM;
            return (status);
    }
	return (status);
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

	extern int status;

	int GetSwitch (Hdr *, char *, int *);

	if (GetSwitch (phdr, calswitch, flag))
	    return (status);

	if (*flag != PERFORM) {
        if (*flag != COMPLETE)
	        *flag = OMIT;
    }

	return (status);
}
