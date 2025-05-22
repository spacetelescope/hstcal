# include <stdio.h>
# include "c_iraf.h"
# include "hstio.h"

# include "stis.h"
# include "hstcalerr.h"

/* This routine prints the status value and message for HSTIO errors,
   IRAF error, or it prints a generic message depending on the value of
   status.

   Phil Hodge, 1997 Nov 7:
	Include INTERNAL_ERROR, and print status number if no other error.

   Phil Hodge, 1998 Oct 5:
	Don't print status number if no other error.
*/

void WhichError (int status) {

	if (status == 0)
	    return;

	if (hstio_err())
	    trlerror("         HSTIO error %d:  %s", status, hstio_errmsg());
	else if (c_iraferr())
	    trlerror("         IRAF error %d:  %s", c_iraferr(), c_iraferrmsg());

	if (status == OUT_OF_MEMORY) {
	    trlerror("         Out of memory.");

	} else if (status == OPEN_FAILED) {
	    trlerror("         Open failed.");

	} else if (status == CAL_FILE_MISSING) {
	    trlerror("         Calibration file(s) missing.");

	} else if (status == NOTHING_TO_DO) {
	    trlerror("         No output from current processing step.");

	} else if (status == KEYWORD_MISSING) {
	    trlerror("         Required keyword missing.");

	} else if (status == ALLOCATION_PROBLEM) {
	    trlerror("         Allocation problem.");

	} else if (status == HEADER_PROBLEM) {
	    trlerror("         Header problem.");

	} else if (status == SIZE_MISMATCH) {
	    trlerror("         Size mismatch.");

	} else if (status == TABLE_ERROR) {
	    trlerror("         Table error.");

	} else if (status == COLUMN_NOT_FOUND) {
	    trlerror("         Column not found.");

	} else if (status == NO_GOOD_DATA) {
	    trlerror("         All data were flagged as bad.");

	} else if (status == REF_TOO_SMALL) {
	    trlerror("         Reference image is binned more than science image.");

	} else if (status == INTERNAL_ERROR) {
	    trlerror("         Internal error.");

	} else if (status < 0) {
	    trlerror("         Internal error; status = %d.", status);
	}
}
