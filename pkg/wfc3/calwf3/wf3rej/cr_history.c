# include   <stdio.h>
# include "hstio.h"
# include   "wf3.h"
# include   "wf3rej.h"

/*  cr_history -- write the history of crrej to the output file

    Description:
    ------------
 
    Date            Author            Description
    ----            ------            -----------
    02-May-1996     J.-C. Hsu         Adapted from the SPP code cr_history.x
    16-May-2008     H. Bushouse       Added detector as function argument and
				      added the logic to set RPTCORR=COMPLETE
				      when processing WFC3 IR images.
    13-Aug-2009     H. Bushouse       Force NEXTEND=3 for IR images.
*/

void cr_history (SingleGroup *sg, clpar *par, int nextend, int detector) {

    int     PutKeyBool (Hdr *, char *, Bool, char *);
    void    UCalVer (Hdr *);
    int     PutKeyFlt (Hdr *, char *, float, char *);	
    int     PutKeyInt (Hdr *, char *, int, char *);
    int     PutKeyStr (Hdr *, char *, char *, char *);

    /* record parameters in the header */
    PutKeyStr (sg->globalhdr, "INITGUES", par->initgues, "");
    PutKeyStr (sg->globalhdr, "SKYSUB", par->sky, "");
    PutKeyStr (sg->globalhdr, "CRSIGMAS", par->sigmas, "");
    PutKeyFlt (sg->globalhdr, "MEANEXP", par->meanexp, "");
    PutKeyFlt (sg->globalhdr, "CRRADIUS", par->radius, "");
    PutKeyFlt (sg->globalhdr, "CRTHRESH", par->thresh, "");
    PutKeyFlt (sg->globalhdr, "SCALENSE", par->scalense, "");
    PutKeyInt (sg->globalhdr, "BADINPDQ", par->badinpdq, "");
    PutKeyBool (sg->globalhdr, "CRMASK", par->mask, "");

    if (detector == IR_DETECTOR) {
        PutKeyStr (sg->globalhdr, "RPTCORR", "COMPLETE", "");
        PutKeyInt (sg->globalhdr, "NEXTEND", 3, "");
    } else {
        PutKeyStr (sg->globalhdr, "CRCORR", "COMPLETE", "");
        PutKeyInt (sg->globalhdr, "NEXTEND", nextend, "");
    }

    UCalVer (sg->globalhdr);
}
