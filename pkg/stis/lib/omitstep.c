# include "stis.h"

/* This routine returns one if the flag is OMIT or COMPLETE.  This means
   that the calibration step was not even attempted.

   If the flag is not OMIT or COMPLETE (i.e. if it is PERFORM, DUMMY,
   SKIPPED, or IGNORED), we have at least tried (or will try) to perform
   the step.  We could therefore log a message and/or add history records
   to the output header.
*/

int OmitStep (int flag) {

/* argument:
int *flag         i: value of switch:  PERFORM, OMIT, COMPLETE, etc.
*/

	if (flag == OMIT || flag == COMPLETE)
	    return (1);
	else
	    return (0);
}
