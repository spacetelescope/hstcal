# include <stdio.h>

# include "stis.h"
# include "calstis1.h"

/* This routine resets switches that are inappropriate for the current
   detector.
*/

void Sanity1 (StisInfo1 *sts) {

	if (sts->detector == CCD_DETECTOR) {

	    sts->lorscorr = OMIT;
	    sts->glincorr = OMIT;
	    sts->lflgcorr = OMIT;
	    sts->doppcorr = OMIT;

	} else {

	    sts->atodcorr = OMIT;
	    sts->blevcorr = OMIT;
	    sts->biascorr = OMIT;
	    sts->shadcorr = OMIT;
	}
}
