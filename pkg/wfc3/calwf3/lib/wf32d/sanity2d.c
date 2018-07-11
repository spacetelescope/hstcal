# include "hstio.h"

# include "wf3.h"
# include "wf3info.h"

/* This routine resets switches that are inappropriate for the current
   detector.
*/

void Sanity2d (WF3Info *wf32d) {

	if (wf32d->detector == IR_DETECTOR) {

	    wf32d->shadcorr = OMIT;
	}
}
