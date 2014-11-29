# include "c_iraf.h"  /* For Bool type */

/* Constants and Definitions for use with IMPHTTAB library */
/* Definitions based on those defined in acs.h */

#define YES         1
#define NO          0

# define SZ_LINE      255
# define SZ_FITS_REC   82
# define SZ_FNAME     255

/* Standard string for use in Error Messages */
char MsgText[SZ_LINE+1];
void errchk ();                 /* HSTIO error check */
int status;

/* Codes to specify whether a reference file exists or not. */
# define EXISTS_YES       1
# define EXISTS_NO        0
# define EXISTS_UNKNOWN (-1)

/* For a reference file name, this string means that a name was
   intentionally not given.
*/
# define NOT_APPLICABLE   "n/a" /* Upper case in headers, converted to
                                    lower case by CALACS */

/* nearest integer function */
# define NINT(x)  ((x >= 0.) ? (int) (x + 0.5) : (int) (x - 0.5))

/* Codes for goodPedigree in RefImage and RefTab to specify whether
   pedigree is good, dummy, or not specified.
*/
# define GOOD_PEDIGREE      1
# define DUMMY_PEDIGREE     0
# define PEDIGREE_UNKNOWN (-1)


/* Error codes for IMPHTTAB code. [based on acserr.h] */

# define PHOT_OK                 0

/* set in stxtools/errxit.x */
# define ERROR_RETURN           2

# define OUT_OF_MEMORY          111

# define OPEN_FAILED            114
# define CAL_FILE_MISSING       115
# define NOTHING_TO_DO          116
# define KEYWORD_MISSING        117
# define ALLOCATION_PROBLEM     118
# define HEADER_PROBLEM         119
# define SIZE_MISMATCH          120

# define CAL_STEP_NOT_DONE      130

# define PHOTTABLE_ERROR        140
# define TABLE_ERROR            141
# define COLUMN_NOT_FOUND       142
# define ELEMENT_NOT_FOUND      143
# define ROW_NOT_FOUND          144

# define NO_GOOD_DATA           151
# define NO_CHIP_FOUND          152

# define REF_TOO_SMALL          171

# define INTERNAL_ERROR         999

# define INVALID_EXPTIME        1001
# define INVALID_FILENAME       1011
# define WRITE_FAILED           1020
# define INVALID_TEMP_FILE      1021
# define FILE_NOT_READABLE      1023
# define COPY_NOT_POSSIBLE      1025

# define INVALID_VALUE          1111
# define UNSUPPORTED_APERTURE   1030


# define MAXPARS 3          /* Max number of parameterizations supported + 1 */

/* The following definition needs to be kept updated with MAXPARS */
static char colnames[9][12] = {"OBSMODE", "DATACOL", "RESULT", "NELEM1", "NELEM2", "PAR1VALUES", "PAR2VALUES", "RESULT1", "RESULT2"};

static char *photnames[6] = {"PHOTZPT","PHOTFLAM", "PHOTPLAM", "PHOTBW","PHTFLAM1","PHTFLAM2"};


typedef struct {
    char name[SZ_LINE+1];            /* name of table */
    int goodPedigree;               /* DUMMY_PEDIGREE if dummy */
    char pedigree[SZ_FITS_REC];    /* value of pedigree (header or row) */
    char descrip[SZ_FITS_REC];     /* value of descrip from header */
    char descrip2[SZ_FITS_REC];    /* value of descrip from row */
    int exists;                     /* does reference table exist? */

    char obsmode[SZ_FITS_REC];    /* obsmode of science data */
    char photmode[SZ_FITS_REC]; /* obsmode used for comparison with IMPHTTAB */

    /* Extrapolation flag */
    Bool extrap;

    /* parsed out value of any parameterized values */
    double *parvalues;
    char **parnames;
    int npar;

    /* Output values derived from table */
    double photflam;
    double photplam;
    double photbw;
    double photzpt;
    double phtflam1;
    double phtflam2;

} PhotPar;


/* Prototypes for IMPHTTAB functions */
int GetPhotTab (PhotPar *obs, char *photmode);
void InitPhotPar(PhotPar *obs, char *name, char *pedigree);
void FreePhotPar(PhotPar *obs);
