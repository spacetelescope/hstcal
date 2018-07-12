# include <stdio.h>
# include "hstio.h"
# include "wf3.h"
# include "wf3info.h"

static int GetSw (Hdr *, char *, int *);

/* This routine attempts to read the values of all known calibration
switches.  It returns these values as either OMIT or PERFORM (0 or 1)
in the CCDSwitch struct.  OMIT and PERFORM are the only possible values
returned by this function.  Missing switches will be set to OMIT.  The
switch will be set to PERFORM only if the calibration switch value in
the header exists and is "PERFORM".
*/

int GetccdSw (WF3Info *wf3, Hdr *phdr) {


	extern int status;

	if (GetSw (phdr, "ATODCORR", &wf3->atodcorr))
	    return (status);
	if (GetSw (phdr, "BIASCORR", &wf3->biascorr))
	    return (status);
	if (GetSw (phdr, "BLEVCORR", &wf3->blevcorr))
	    return (status);
	if (GetSw (phdr, "DQICORR",  &wf3->dqicorr))
	    return (status);
	if (GetSw (phdr, "FLSHCORR", &wf3->flashcorr))
	    return (status);
    if (GetSw (phdr, "PCTECORR", &wf3->pctecorr))
        return(status);

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

	if (*flag != PERFORM)
	    *flag = OMIT;

	return (status);
}
