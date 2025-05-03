#ifndef HSTCAL_ERR_INCL
#define HSTCAL_ERR_INCL

struct hstcal_error_state {
    struct hstcal_error_hstio {
        int code;
        char *message;
    } hstio;

    struct hstcal_error_state_iraf {
        int code;
        char *message;
    } iraf;

    struct hstcal_error_state_runtime {
        int code;
        const char *message;
    } runtime;

    struct hstcal_error_state_origin {
        const char *file_name;
        const char *function_name;
        int line_number;
    } origin;
};

struct hstcal_error_state hstcal_error_state_populate(const char *filename, const int line_number, const char *function_name);
void hstcal_error_state_show(const char *filename, const int line_number, const char *function_name);

// Report the error associated with the last status
#define REPORT_ERROR_STATE() hstcal_error_state_populate(__FILE__, __LINE__, __FUNCTION__)

// Handle WhichError deprecation.
// WhichError took the global status variable as an argument.
// This retains the function signature, however the status is ignored
#define WhichError(NOP) REPORT_ERROR_STATE()

// Error codes

#define ACS_OK                 0
#define STIS_OK                0
#define WF3_OK                 0
#define PHOT_OK                0
#define HSTCAL_OK              0

#define ERROR_RETURN           2

#define OUT_OF_MEMORY          111

#define OPEN_FAILED            114
#define CAL_FILE_MISSING       115
#define NOTHING_TO_DO          116
#define KEYWORD_MISSING        117
#define ALLOCATION_PROBLEM     118
#define HEADER_PROBLEM         119
#define SIZE_MISMATCH          120
#define IO_ERROR               121

#define CAL_STEP_NOT_DONE      130

#define PHOTTABLE_ERROR        140
#define TABLE_ERROR            141
#define COLUMN_NOT_FOUND       142
#define ELEMENT_NOT_FOUND      143
#define ROW_NOT_FOUND          144

#define NO_GOOD_DATA           151
#define NO_CHIP_FOUND          152

#define REF_TOO_SMALL          171

#define INTERNAL_ERROR         999

#define INVALID_EXPTIME        1001
#define INVALID_FILENAME       1011
#define WRITE_FAILED           1020
#define INVALID_TEMP_FILE      1021
#define FILE_NOT_READABLE      1023
#define COPY_NOT_POSSIBLE      1025

#define INVALID_VALUE          1111
#define GENERIC_ERROR_CODE     1112 //used in STIS

#define UNSUPPORTED_APERTURE   1030

#define UNKNOWN_MAGIC_CODE_1010   1010  // ./pkg/acs/lib/acsccd/doatod.c:124
                                        // ./pkg/wfc3/lib/wf3ccd/doatod.c:122

#define UNKNOWN_MAGIC_CODE_1991   1991  // ./pkg/acs/lib/binupdate.c:43
                                        // ./pkg/wfc3/lib/binupdate.c:43

#define UNKNOWN_MAGIC_CODE_2011   2011  // ./pkg/stis/lib/mkoutname.c:193

#endif
