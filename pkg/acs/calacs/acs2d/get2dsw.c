# include <stdio.h>
# include "hstio.h"
# include "acs.h"
# include "acsinfo.h"
# include "err.h"
# include "acscorr.h"


static int GetSw (Hdr *, char *, int *);


/* This routine attempts to read the values of all known calibration
 switches.  It returns these values as either OMIT or PERFORM (0 or 1)
 in the CalSwitch struct.  OMIT and PERFORM are the only possible values
 returned by this function.  Missing switches will be set to OMIT.  The
 switch will be set to PERFORM only if the calibration switch value in
 the header exists and is "PERFORM".

 9 Oct 2014 PLL - Created based on getccdsw.c
 */
int Get2dSw (CalSwitch *acs2d, Hdr *phdr) {

    extern int status;
    /*FitsKw key;*/      /* keyword location in header */

    if (GetSw (phdr, "DQICORR",  &acs2d->dqicorr))
        return (status);
    if (GetSw (phdr, "GLINCORR",  &acs2d->glincorr))
        return (status);
    if (GetSw (phdr, "LFLGCORR",  &acs2d->lflgcorr))
        return (status);
    if (GetSw (phdr, "DARKCORR",  &acs2d->darkcorr))
        return (status);
    if (GetSw (phdr, "FLSHCORR",  &acs2d->flashcorr))
        return (status);
    if (GetSw (phdr, "FLATCORR",  &acs2d->flatcorr))
        return (status);
    if (GetSw (phdr, "SHADCORR",  &acs2d->shadcorr))
        return (status);
    if (GetSw (phdr, "PHOTCORR",  &acs2d->photcorr))
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
