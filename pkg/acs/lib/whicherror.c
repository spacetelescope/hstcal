# include <stdio.h>
# include <string.h>
#include "hstcal.h"
# include "hstio.h"
# include <c_iraf.h>
# include "hstcalerr.h"
# include "acs.h"	/* for message output */

/* This routine prints the status value and message for HSTIO errors,
   IRAF error, or it prints a generic message depending on the value of
   status.

   Phil Hodge, 1997 Nov 7:
	Include INTERNAL_ERROR, and print status number if no other error.

   Warren Hack, 1998 May 26:
	Revised to use ACS header files
*/

void WhichError (int status) {

	if (status == ACS_OK)
	    return;

	if (hstio_err())
	    sprintf (MsgText,"         HSTIO error %d:  %s", status, hstio_errmsg());
	else if (c_iraferr())
	    sprintf (MsgText, "         IRAF error %d:  %s", c_iraferr(), c_iraferrmsg());
	trlerror (MsgText);
	
	/* Now, check value of STATUS */
	if (status == OUT_OF_MEMORY) {
	    trlerror("Out of memory.");

	} else if (status == OPEN_FAILED) {
	    trlerror("Open failed.");

	} else if (status == CAL_FILE_MISSING) {
	    trlerror("Calibration file(s) missing.");

	} else if (status == CAL_STEP_NOT_DONE) {
	    trlerror("Previous required calibration processing not done.");

	} else if (status == NOTHING_TO_DO) {
	    trlerror("No output from current processing step.");

	} else if (status == KEYWORD_MISSING) {
	    trlerror("Required keyword missing.");

	} else if (status == ALLOCATION_PROBLEM) {
	    trlerror("Allocation problem.");

	} else if (status == HEADER_PROBLEM) {
	    trlerror("Header problem.");

	} else if (status == SIZE_MISMATCH) {
	    trlerror("Size mismatch.");

	} else if (status == TABLE_ERROR) {
	    trlerror("Table error.");

	} else if (status == COLUMN_NOT_FOUND) {
	    trlerror("Column not found.");

	} else if (status == NO_GOOD_DATA) {
	    trlerror("All data were flagged as bad.");

	} else if (status == REF_TOO_SMALL) {
	    trlerror("Reference image is binned more than science image.");

	} else if (status == INTERNAL_ERROR) {
	    trlerror("Internal error.");

	} else if (status < 0) {
	    sprintf (MsgText,"         Internal error; status = %d.", status);
		trlerror (MsgText);

	} else {
	    sprintf (MsgText,"         status = %d", status);
		trlerror (MsgText);
	}

}
