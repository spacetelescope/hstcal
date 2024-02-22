# include "hstio.h"
# include "acs.h"
# include "acsinfo.h"

/* This routine resets switches that are inappropriate for the current
   detector.
*/

void Sanity2d (ACSInfo *acs2d) {

	if (acs2d->detector != MAMA_DETECTOR) {

	    acs2d->glincorr = OMIT;
	    acs2d->lflgcorr = OMIT;

	} else {

	    acs2d->shadcorr = OMIT;
        acs2d->flashcorr = OMIT;

	}
}
