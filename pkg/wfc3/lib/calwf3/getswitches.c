# include <stdio.h>

# include "hstio.h"

# include "wf3.h"
# include "wf3corr.h"
# include "hstcalerr.h"
# include "calwf3.h"

static int GetSw (Hdr *, char *, int *);

/* This routine attempts to read the values of all known calibration
switches.  It returns these values as either OMIT or PERFORM (0 or 1)
in the CCD_Switch struct.  OMIT and PERFORM are the only possible values
returned by this function.  Missing switches will be set to OMIT.  The
switch will be set to PERFORM only if the calibration switch value in
the header exists and is "PERFORM". Produces an error when both RPTCORR
and CRCORR are set to PERFORM for CCD images.

H.Bushouse, 2001 May 8:
	Revised to add FLSHCORR processing keyword.
H.Bushouse, 2002 June 20:
	Removed use of STATFLAG switch (always perform doStat by default).
M.Sosey, 2013 Dec 4:
    revised to add FLUXCORR processing keyword
    
M.sosey 2014 August
    revised for CTE correction in UVIS, PCTECORR
*/

int GetCCDSws (CCD_Switch *sw, Hdr *phdr) {

/* arguments:
CCD_Switch *sw  o: values (0 or 1) of calibration switches
Hdr *phdr       i: primary header
*/

	extern int status;

	if (GetSw (phdr, "ATODCORR", &sw->atodcorr))
	    return (status);
	if (GetSw (phdr, "BIASCORR", &sw->biascorr))
	    return (status);
	if (GetSw (phdr, "BLEVCORR", &sw->blevcorr))
	    return (status);
	if (GetSw (phdr, "CRCORR",   &sw->crcorr))
	    return (status);
	if (GetSw (phdr, "RPTCORR",  &sw->rptcorr))
	    return (status);
	if (GetSw (phdr, "DARKCORR", &sw->darkcorr))
	    return (status);
	if (GetSw (phdr, "DQICORR",  &sw->dqicorr))
	    return (status);
	if (GetSw (phdr, "FLATCORR", &sw->flatcorr))
	    return (status);
	if (GetSw (phdr, "FLSHCORR", &sw->flashcorr))
	    return (status);
	if (GetSw (phdr, "PHOTCORR", &sw->photcorr))
	    return (status);
	if (GetSw (phdr, "SHADCORR", &sw->shadcorr))
	    return (status);
	if (GetSw (phdr, "EXPSCORR", &sw->expscorr))
	    return (status);
	if (GetSw (phdr, "FLUXCORR", &sw->fluxcorr))
	    return (status);
	if (GetSw (phdr, "PCTECORR", &sw->pctecorr))
	    return (status);
    
        
	/* Check to ensure that only one of these switches are
	** set, since they are exclusive options. */
	if (sw->rptcorr == PERFORM && sw->crcorr == PERFORM) {
	    trlerror("RPTCORR and CRCORR both set to PERFORM!");
	    trlerror("    One switch needs to be set to OMIT!");
	    status = HEADER_PROBLEM;
	    return (status);  
	}
    
    if (sw->photcorr == OMIT && sw->fluxcorr == PERFORM){
        trlerror("PHOTCORR and FLUXCORR both need to be set to PERFORM!");
        trlerror("     FLUXCORR needs PHOTCORR to set required keywords");
        status = HEADER_PROBLEM;
        return (status);
    }
	return (status);
}

int GetIRSws (IR_Switch *sw, Hdr *phdr) {

/* arguments:
IR_Switch *sw  o: values (0 or 1) of calibration switches
Hdr *phdr       i: primary header
*/

	extern int status;

	if (GetSw (phdr, "ZSIGCORR", &sw->zsigcorr))
	    return (status);
	if (GetSw (phdr, "ZOFFCORR", &sw->zoffcorr))
	    return (status);
	if (GetSw (phdr, "DQICORR",  &sw->dqicorr))
	    return (status);
	if (GetSw (phdr, "BLEVCORR", &sw->blevcorr))
	    return (status);
	sw->noiscorr = PERFORM;
	if (GetSw (phdr, "DARKCORR", &sw->darkcorr))
	    return (status);
	if (GetSw (phdr, "NLINCORR", &sw->nlincorr))
	    return (status);
	if (GetSw (phdr, "FLATCORR", &sw->flatcorr))
	    return (status);
	if (GetSw (phdr, "UNITCORR", &sw->unitcorr))
	    return (status);
	if (GetSw (phdr, "PHOTCALC", &sw->photcorr))
	    return (status);
	if (GetSw (phdr, "CRCORR"  , &sw->crcorr))
	    return (status);
	if (GetSw (phdr, "RPTCORR",  &sw->rptcorr))
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

	if (*flag != PERFORM)
	    *flag = OMIT;

	return (status);
}

