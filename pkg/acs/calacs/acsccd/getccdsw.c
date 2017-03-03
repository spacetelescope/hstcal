# include <stdio.h>
# include "hstio.h"
# include "acs.h"
# include "acsinfo.h"
# include "acserr.h"


static int GetSw (Hdr *, char *, int *);


/* This routine attempts to read the values of all known calibration
 switches.  It returns these values as either OMIT or PERFORM (0 or 1)
 in the CalSwitch struct.  OMIT and PERFORM are the only possible values
 returned by this function.  Missing switches will be set to OMIT.  The
 switch will be set to PERFORM only if the calibration switch value in
 the header exists and is "PERFORM".

 07 May 2001 - Added code to read in FLSHCORR keyword when running ACSCCD.
 08 Mar 2011 MRD - Added code to read PCTECORR keyword.
 12 Dec 2012 PLL - Moved FLSHCORR to ACS2D.
 12 Aug 2013 PLL - Separated PCTECORR from ACSCCD.
 21 Feb 2017 PLL - Added SINKCORR.
 */
int GetccdSw (ACSInfo *acs, Hdr *phdr) {

    extern int status;
    /*FitsKw key;*/        /* keyword location in header */

    if (GetSw (phdr, "DQICORR",  &acs->dqicorr))
        return (status);
    if (GetSw (phdr, "ATODCORR", &acs->atodcorr))
        return (status);
    if (GetSw (phdr, "BIASCORR", &acs->biascorr))
        return (status);
    if (GetSw (phdr, "BLEVCORR", &acs->blevcorr))
        return (status);
    if (GetSw (phdr, "SINKCORR", &acs->sinkcorr))
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
