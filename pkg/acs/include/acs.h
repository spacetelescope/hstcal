#ifndef ACS_INCL
#define ACS_INCL

/* acs.h generic header for calacs */
/*      */
# include <stdio.h>             /* To insure that FILE is defined for TrlPtr */
# include "imphttab.h"
# include "hstcal.h"
#include "trlbuf.h"

# define ACS_CBUF           24  /* small buffer for e.g. rootname */
# define ACS_FITS_REC       82
# define SZ_STRKWVAL        68

#define ACS_WFC_N_COLUMNS_PER_CHIP_INCL_OVERSCAN 4144
#define ACS_WFC_N_COLUMNS_PER_CHIP_EXCL_OVERSCAN 4096

/* Macros for dusing GetKey/PutKey functions.... */
# define USE_DEFAULT    1       /* Use default if keyword is missing */
# define NO_DEFAULT     0       /* missing keyword is fatal error */

typedef unsigned char Byte;

#define SIZE_BYTE   8

# define MAX_DQ     65535

/* Number of lines to extract from binned images for unbinning */
# define SECTLINES  2

/* Three extensions per SingleGroup. */
# define EXT_PER_GROUP 3

void errchk ();                 /* HSTIO error check */

/* Integer codes for string-valued header keywords. */

/* Flag values for calibration switches. */
# define DUMMY  (-1)
# define OMIT     0
# define PERFORM  1
# define COMPLETE 2
# define SKIPPED  3
# define IGNORED  4         /* e.g. SHADCORR when the exposure time is zero */

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

/* Possible values for detector. */
# define UNKNOWN_DETECTOR   (-1)
# define WFC_CCD_DETECTOR   1
# define HRC_CCD_DETECTOR   2
# define MAMA_DETECTOR      3

/* Definitions for interpreting CCDAMP, and using multiamp struct */
# define        NAMPS           4
# define        AMPSTR1         "CD"
# define        AMPSTR2         "AB"
# define        AMPSORDER       "ABCD"
/* Define array indices for each amp for clarity of code. 22Mar99, WJH */
# define        AMP_A           0
# define        AMP_B           1
# define        AMP_C           2
# define        AMP_D           3

# define       DEFAULT_OFFSET   3

# define        SM4MJD          54967

/* October 01, 2016: Date of first observation in Cycle 24. The new ACS subarray configurations validated for Cycle 24. */
# define        CYCLE24         57662

/* A reference image. */
typedef struct {
    char name[CHAR_FNAME_LENGTH];            /* name of image */
    char pedigree[ACS_FITS_REC];    /* value of pedigree keyword */
    char descrip[ACS_FITS_REC];     /* value of descrip keyword */
    int exists;                     /* does reference image exist? */
    int goodPedigree;               /* DUMMY_PEDIGREE if dummy */
} RefImage;

/* A reference table. */
typedef struct {
    char name[CHAR_FNAME_LENGTH];            /* name of table */
    char pedigree[ACS_FITS_REC];    /* value of pedigree (header or row) */
    char descrip[ACS_FITS_REC];     /* value of descrip from header */
    char descrip2[ACS_FITS_REC];    /* value of descrip from row */
    int exists;                     /* does reference table exist? */
    int goodPedigree;               /* DUMMY_PEDIGREE if dummy */
} RefTab;

/* This section is for saving the names of reference files. */
# define ACS_KEYWORD        10
typedef struct ref *RefFileInfoPtr;
typedef struct ref {
	char keyword[ACS_KEYWORD];
	char filename[ACS_FITS_REC];
	RefFileInfoPtr next;
} RefFileInfo;

typedef struct {
    float val[NAMPS];
    int colx;       /* If colx == 0, single amp readout */
    int coly;       /* if coly == 0, NOT 4-amp readout */
    int chip;       /* Chip being processed */
    int detector;   /* Which detector was used */
} multiamp;

/* The following function definitions handle the messages created
	by CALACS during operations.  These will have counterparts which
	send output both the STDOUT and a TRL file.
*/
void asnwarn (char *message);
void asnerror (char *message);

#endif
