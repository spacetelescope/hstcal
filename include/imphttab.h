#ifndef IMPHTTAB_INCL
#define IMPHTTAB_INCL

# include "c_iraf.h"  /* For Bool type */
# include "xtables.h" /* for SZ_COLNAME*/
# include "hstcalerr.h"

/* Constants and Definitions for use with IMPHTTAB library */
/* Definitions based on those defined in acs.h */

#define YES         1
#define NO          0

# define SZ_LINE      255
# define SZ_FITS_REC   82
# define SZ_FNAME     255

/* Standard string for use in Error Messages */
/*char MsgText[SZ_LINE+1];*/
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

# define MAXPARS 3          /* Max number of parameterizations supported + 1 */

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


typedef struct {
    IRAFPointer tp;            /* pointer to table descriptor */
    IRAFPointer cp_obsmode;        /* column descriptors */
    IRAFPointer cp_datacol;
    IRAFPointer cp_result[MAXPARS+1];
    IRAFPointer cp_nelem[MAXPARS];
    IRAFPointer cp_parnames[MAXPARS];
    IRAFPointer cp_parvalues[MAXPARS];
    IRAFPointer cp_pedigree;
    IRAFPointer cp_descrip;
    char photvar[SZ_COLNAME]; /* photometric parameter in this table */
    int nrows;            /* number of rows in table */
    int ncols;          /* number of columns in table */
    int parnum;         /* number of parameterized variables in tables */
} PhtCols;

typedef struct {
    char obsmode[SZ_FITS_REC];  /* obsmode string read from table row */
    char datacol[SZ_FITS_REC];
    char **parnames; /* record the par. variable names for comparison with obsmode string */
    int parnum;      /* number of parameterized variables */
    double *results;  /* results[telem] or results[nelem1*nelem2*...] */
    int telem;     /* total number of parameterized variable values */
    int *nelem;    /* multiple paramerized variables will each N values */
    double **parvals; /* need to support multiple parameterized variables */
} PhtRow;


typedef struct {
    int    ndim;   /* number of dimensions for each position */
    double *index; /* array index along each axis,
              used to determine which axis is being interpolated */
    double *pos;   /* value of each axis at position given by index */
    double value;  /* value at position */
} BoundingPoint;



/* Prototypes for IMPHTTAB functions */
int GetPhotTab (PhotPar *obs, char *photmode);
void InitPhotPar(PhotPar *obs, char *name, char *pedigree);
void FreePhotPar(PhotPar *obs);

#endif
