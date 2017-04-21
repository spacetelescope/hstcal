#ifndef ERR_INCL
#define ERR_INCL

// Error codes

#define ACS_OK                 0
#define STIS_OK                0
#define WF3_OK                 0
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
#define GENERIC_ERROR_CODE     1111 //used in STIS

#define UNSUPPORTED_APERTURE   1030

#endif
