# include <stdio.h>
# include  "hstio.h"
# include "wf3.h"
# include "wf3info.h"

static int GetSw (Hdr *, char *, int *);

/* This routine attempts to read the values of all known calibration
switches.  It returns these values as either OMIT or PERFORM (0 or 1)
in the CCDSwitch struct.  OMIT and PERFORM are the only possible values
returned by this function.  Missing switches will be set to OMIT.  The
switch will be set to PERFORM only if the calibration switch value in
the header exists and is "PERFORM".

Revision History:

  H.Bushouse	19-Apr-2001	Initial version copied from equivalent
				WFC3 UVIS routine for WFC3 IR use.

  H.Bushouse	21-Oct-2010	Modified GetSw to not reset switch values to
				OMIT if they are set to something other than
				PERFORM - just return the actual value. This
				is to support re-entrant processing where some
				switches may be COMPLETE. (PR 66081)
*/

int GetirSw (WF3Info *wf3, Hdr *phdr) {


	extern int status;

	if (GetSw (phdr, "ZSIGCORR", &wf3->zsigcorr))
	    return (status);
	if (GetSw (phdr, "ZOFFCORR", &wf3->zoffcorr))
	    return (status);
	if (GetSw (phdr, "DQICORR",  &wf3->dqicorr))
	    return (status);
	if (GetSw (phdr, "DARKCORR", &wf3->darkcorr))
	    return (status);
	if (GetSw (phdr, "BLEVCORR", &wf3->blevcorr))
	    return (status);
	if (GetSw (phdr, "NLINCORR", &wf3->nlincorr))
	    return (status);
	if (GetSw (phdr, "FLATCORR", &wf3->flatcorr))
	    return (status);
	if (GetSw (phdr, "UNITCORR", &wf3->unitcorr))
	    return (status);
	if (GetSw (phdr, "PHOTCORR", &wf3->photcorr))
	    return (status);
	if (GetSw (phdr, "CRCORR",   &wf3->crcorr))
	    return (status);

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

/*	Return the actual value; don't reset
	if (*flag != PERFORM)
	    *flag = OMIT;
*/

	return (status);
}
