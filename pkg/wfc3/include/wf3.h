#ifndef INCL_WF3_H
#define INCL_WF3_H

/* wf3.h generic header for calwf3 */

# include <stdio.h>             /* To insure that FILE is defined for TrlPtr */
# include "msg.h"
# include "imphttab.h"
# include "hstcal.h"
#include "trlbuf.h"

/* Macros for dusing GetKey/PutKey functions.... */
# define USE_DEFAULT    1       /* Use default if keyword is missing */
# define NO_DEFAULT     0       /* missing keyword is fatal error */

typedef unsigned char Byte;

#define SIZE_BYTE   8

# define MAX_DQ     65535

# define ATOD_SATURATE 65534

/* Number of lines to extract from binned images for unbinning */
# define SECTLINES  2

/* Three extensions per SingleGroup. */
# define EXT_PER_GROUP 3

/* Integer codes for string-valued header keywords. */

/* Flag values for calibration switches. */
enum SwitchVals_ {BLANK=(-2), DUMMY, OMIT, PERFORM, COMPLETE, SKIPPED,
		  IGNORED, SKIP, OMITTED, PERFORMED};
typedef enum SwitchVals_ SwitchVals;

/* Codes to specify whether a reference file exists or not. */
# define EXISTS_YES       1
# define EXISTS_NO        0
# define EXISTS_UNKNOWN (-1)

/* For a reference file name, this string means that a name was
   intentionally not given.
*/
# define NOT_APPLICABLE   "n/a" /* Upper case in headers, converted to
                                    lower case by CALWF3 */

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
# define CCD_DETECTOR   1
# define IR_DETECTOR    2

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

# define	DEFAULT_OFFSET	3

/* used for the CTE correction in UVIS where the amps are stacked
   in the order they are read out. Defined here so that the rest of
   the pipeline has access to them, for use mostly in the sink pixel
   mask creation which happens in wf3ccd. This is the full extent size
	 of a raw file which has been rotated in amp order. */
# define RAZ_COLS 8412
# define RAZ_ROWS 2070


/* A reference image. */
typedef struct {
    char name[CHAR_FNAME_LENGTH+1];            /* name of image */
    char type[SZ_FITS_REC+1];        /* value of filetype */
    char pedigree[SZ_FITS_REC+1];    /* value of pedigree keyword */
    char descrip[SZ_FITS_REC+1];     /* value of descrip keyword */
    int exists;                     /* does reference image exist? */
    int goodPedigree;               /* DUMMY_PEDIGREE if dummy */
} RefImage;


/* A reference table. */
typedef struct {
    char name[CHAR_FNAME_LENGTH+1];            /* name of table */
    char type[SZ_FITS_REC+1];        /* value of filetype */
    char pedigree[SZ_FITS_REC+1];    /* value of pedigree (header or row) */
    char descrip[SZ_FITS_REC+1];     /* value of descrip from header */
    char descrip2[SZ_FITS_REC+1];    /* value of descrip from row */
    int exists;                     /* does reference table exist? */
    int goodPedigree;               /* DUMMY_PEDIGREE if dummy */
} RefTab;

/* This section is for saving the names of reference files. */
typedef struct ref *RefFileInfoPtr;

typedef struct ref {
	char keyword[SZ_KEYWORD+1];
	char filename[SZ_FITS_REC+1];
	RefFileInfoPtr next;
} RefFileInfo;

typedef struct {
    float val[NAMPS];
    int colx;       /* If colx == 0, single amp readout */
    int coly;       /* if coly == 0, NOT 4-amp readout */
    int chip;       /* Chip being processed */
    int detector;   /* Which detector was used */
} multiamp;

#endif /* INCL_WF3_H */
