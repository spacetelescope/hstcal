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
	    trlmessage("         HSTIO error %d:  %s\n", status, hstio_errmsg());
	else if (c_iraferr())
	    printf (
	"         IRAF error %d:  %s\n", c_iraferr(), c_iraferrmsg());

	if (status == OUT_OF_MEMORY) {
	    trlmessage("         Out of memory.\n");

	} else if (status == OPEN_FAILED) {
	    trlmessage("         Open failed.\n");

	} else if (status == CAL_FILE_MISSING) {
	    trlmessage("         Calibration file(s) missing.\n");

	} else if (status == NOTHING_TO_DO) {
	    trlmessage("         No output from current processing step.\n");

	} else if (status == KEYWORD_MISSING) {
	    trlmessage("         Required keyword missing.\n");

	} else if (status == ALLOCATION_PROBLEM) {
	    trlmessage("         Allocation problem.\n");

	} else if (status == HEADER_PROBLEM) {
	    trlmessage("         Header problem.\n");

	} else if (status == SIZE_MISMATCH) {
	    trlmessage("         Size mismatch.\n");

	} else if (status == TABLE_ERROR) {
	    trlmessage("         Table error.\n");

	} else if (status == COLUMN_NOT_FOUND) {
	    trlmessage("         Column not found.\n");

	} else if (status == NO_GOOD_DATA) {
	    trlmessage("         All data were flagged as bad.\n");

	} else if (status == REF_TOO_SMALL) {
	    printf (
	"         Reference image is binned more than science image.\n");

	} else if (status == INTERNAL_ERROR) {
	    trlmessage("         Internal error.\n");

	} else if (status < 0) {
	    trlmessage("         Internal error; status = %d.\n", status);
	}
}
