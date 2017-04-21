# include <stdio.h>
# include "hstio.h"
# include "wf3.h"
# include "wf3info.h"
# include "hstcalerr.h"
# include "wf3corr.h"


static int GetSw (Hdr *, char *, int *);


/* This routine attempts to read the values of all known calibration
 switches.  It returns these values as either OMIT or PERFORM (0 or 1)
 in the CalSwitch struct.  OMIT and PERFORM are the only possible values
 returned by this function.  Missing switches will be set to OMIT.  The
 switch will be set to PERFORM only if the calibration switch value in
 the header exists and is "PERFORM".

 9 Oct 2014 PLL - Created based on getccdsw.c
 */
int Get2dSw (CalSwitch *wf32, Hdr *phdr) {

    extern int status;
    FitsKw key;        /* keyword location in header */

    if (GetSw (phdr, "DQICORR",  &wf32d->dqicorr))
        return (status);
    if (GetSw (phdr, "DARKCORR",  &wf32d->darkcorr))
        return (status);
    if (GetSw (phdr, "FLSHCORR",  &wf32d->flashcorr))
        return (status);
    if (GetSw (phdr, "FLATCORR",  &wf32d->flatcorr))
        return (status);
    if (GetSw (phdr, "SHADCORR",  &wf32d->shadcorr))
        return (status);
    if (GetSw (phdr, "PHOTCORR",  &wf32d->photcorr))
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
