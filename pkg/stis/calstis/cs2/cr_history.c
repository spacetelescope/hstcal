# include	<stdio.h>
# include	"c_iraf.h"
# include	"hstio.h"

# include	"stis.h"
# include	"cs2.h"
# include	"stisdef.h"

/*  cr_history -- write the history of crrej to the output file

  Description:
  ------------
  
  Date		Author			Description
  ----		------			-----------
  02-May-1996  J.-C. Hsu		Adapted from the SPP code cr_history.x
  10-Feb-2000  Phil Hodge		Replace put_key routines with putKey[].
*/

void cr_history (SingleGroup *sg, clpar *par)
{
	
	/* record parameters in the header */
        putKeyS (sg->globalhdr, "INITGUES", par->initial, "");
        putKeyS (sg->globalhdr, "SKYSUB", par->sky, "");
        putKeyS (sg->globalhdr, "CRSIGMAS", par->sigmas, "");
        putKeyF (sg->globalhdr, "MEANEXP", par->meanexp, "");
        putKeyF (sg->globalhdr, "CRRADIUS", par->rej, "");
        putKeyF (sg->globalhdr, "CRTHRESH", par->psigma, "");
        putKeyF (sg->globalhdr, "SCALENSE", par->scalenoise, "");
        putKeyI (sg->globalhdr, "BADINPDQ", par->badbits, "");
        putKeyB (sg->globalhdr, "CRMASK", par->mask, "");

        putKeyI (sg->globalhdr, "NEXTEND", EXT_PER_GROUP, "");
        putKeyS (sg->globalhdr, "CRCORR", "COMPLETE", "");

	UCalVer (sg->globalhdr);
}
