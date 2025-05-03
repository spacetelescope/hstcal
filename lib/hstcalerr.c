# include <stdio.h>
# include <string.h>
#include "hstcal.h"
# include "hstio.h"
# include <c_iraf.h>
# include "trlbuf.h"
# include "hstcalerr.h"

extern int status;

/**
 *
 * @param filename Path to file containing the call to REPORT_ERROR
 * @param line_number The line number the call to REPORT_ERROR appears in `filename`
 * @param function_name The name of the function that called REPORT_ERROR
 * @return a populated `hstcal_error_state` structure when global `status` is non-zero
 * @return an empty `hstcal_error_state` structure when global `status` is zero
 */
struct hstcal_error_state hstcal_error_state_populate(const char *filename, const int line_number, const char *function_name) {
	if (!status) {
		return (struct hstcal_error_state) {0};
	}

#ifndef NDEBUG
	fprintf(stderr, "DEBUG: %s called in %s() at %s:%d\n", __FUNCTION__, function_name, filename, line_number);
#endif

	struct hstcal_error_state e;
	e.origin.file_name = filename;
	e.origin.line_number = line_number;
	e.origin.function_name = function_name;
	e.hstio.code = hstio_err();
	e.hstio.message = hstio_errmsg();
	e.iraf.code = c_iraferr();
	e.iraf.message = c_iraferrmsg();
	e.runtime.code = status;
	e.runtime.message = NULL;

	/* Map the status code to an appropriate runtime error message */
	switch (e.runtime.code) {
		case OUT_OF_MEMORY:
			e.runtime.message = "Out of memory.";
			break;
		case OPEN_FAILED:
			e.runtime.message = "Open failed.";
			break;
		case CAL_FILE_MISSING:
			e.runtime.message = "Calibration file(s) missing.";
			break;
		case NOTHING_TO_DO:
			e.runtime.message = "No output from current processing step.";
			break;
		case ALLOCATION_PROBLEM:
			e.runtime.message = "Allocation problem.";
			break;
		case KEYWORD_MISSING:
			e.runtime.message = "Required keyword missing.";
			break;
		case HEADER_PROBLEM:
			e.runtime.message = "Header problem.";
			break;
		case SIZE_MISMATCH:
			e.runtime.message = "Size mismatch.";
			break;
		case IO_ERROR:
			e.runtime.message = "I/O error.";
			break;
		case CAL_STEP_NOT_DONE:
			e.runtime.message = "Previous required calibration processing not done.";
			break;
		case PHOTTABLE_ERROR:
			e.runtime.message = "Phottable error.";
			break;
		case TABLE_ERROR:
			e.runtime.message = "Table error.";
			break;
		case COLUMN_NOT_FOUND:
			e.runtime.message = "Column not found.";
			break;
		case ELEMENT_NOT_FOUND:
			e.runtime.message = "Element not found.";
			break;
		case ROW_NOT_FOUND:
			e.runtime.message = "Row not found.";
			break;
		case NO_GOOD_DATA:
			e.runtime.message = "All data were flagged as bad.";
			break;
		case NO_CHIP_FOUND:
			e.runtime.message = "No chip found.";
			break;
		case REF_TOO_SMALL:
			e.runtime.message = "Reference image is binned more than science image.";
			break;
		case INTERNAL_ERROR:
			e.runtime.message = "Internal error.";
			break;
		case INVALID_EXPTIME:
			e.runtime.message = "Invalid exposure time.";
			break;
		case INVALID_FILENAME:
			e.runtime.message = "Invalid filename.";
			break;
		case WRITE_FAILED:
			e.runtime.message = "Write failed.";
			break;
		case INVALID_TEMP_FILE:
			e.runtime.message = "Invalid temporary file.";
			break;
		case FILE_NOT_READABLE:
			e.runtime.message = "File not readable.";
			break;
		case COPY_NOT_POSSIBLE:
			e.runtime.message = "Copy not possible.";
			break;
		case INVALID_VALUE:
			e.runtime.message = "Invalid value.";
			break;
		case GENERIC_ERROR_CODE:
			e.runtime.message = "Generic error.";
			break;
		case UNSUPPORTED_APERTURE:
			e.runtime.message = "Unsupported aperture.";
			break;
		case UNKNOWN_MAGIC_CODE_1991:
			e.runtime.message = "PLACEHOLDER for unknown magic code 1991.";
			break;
		case UNKNOWN_MAGIC_CODE_2011:
			e.runtime.message = "PLACEHOLDER for unknown magic code 2011.";
			break;
		default:
			e.runtime.message = "Unknown error.";
			break;
	}

	return e;
}

/* This routine prints the status value and message for HSTIO errors,
   IRAF error, or it prints a generic message depending on the value of
   status.

   Phil Hodge, 1997 Nov 7:
	Include INTERNAL_ERROR, and print status number if no other error.

   Warren Hack, 1998 May 26:
	Revised to use ACS header files

   Joseph Hunkeler, 2025 May 3:
    Rename from WhichError to hstcal_error_state_show, and hook into hstcal_error_state_populate.
    See: REPORT_ERROR and WhichError macros in hstcalerr.h
*/
void hstcal_error_state_show(const char *filename, const int line_number, const char *function_name) {
	const struct hstcal_error_state e = hstcal_error_state_populate(filename, line_number, function_name);
	if (e.hstio.code) {
		trlerror("HSTIO error %d: %s\n", e.hstio.code, e.hstio.message);
	} else if (e.iraf.code) {
		trlerror("IRAF error %d: %s\n", e.iraf.code, e.iraf.message);
	}
	if (e.runtime.message) {
		trlerror("%s (status = %d)", e.runtime.message, e.runtime.code);
	}

}
