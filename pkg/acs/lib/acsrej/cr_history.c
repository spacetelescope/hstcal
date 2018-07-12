# include   <stdio.h>
# include   "hstio.h"

# include   "acs.h"
# include   "acsrej.h"

/*  cr_history -- write the history of crrej to the output file

    Description:
    ------------
 
    Date            Author            Description
    ----            ------            -----------
    02-May-1996     J.-C. Hsu         Adapted from the SPP code cr_history.x
*/

void cr_history (SingleGroup *sg, clpar *par, int nextend)
{
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

    PutKeyInt (sg->globalhdr, "NEXTEND", nextend, "");
    PutKeyStr (sg->globalhdr, "CRCORR", "COMPLETE", "");

    UCalVer (sg->globalhdr);
}
