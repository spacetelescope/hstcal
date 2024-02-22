#ifndef INCL_WF3REJ_H
#define INCL_WF3REJ_H

/* The header file wf3.h must be called prior to calling this
    file, in order to define the macros used here.
*/

/* define the user input parameter constants */
# define    TOTAL       0
# define    SCALENSE    1
# define    INITGUES    2
# define    SKYSUB      3
# define    CRSIGMAS    4
# define    CRRADIUS    5
# define    CRTHRESH    6
# define    BADINPDQ    7
# define    CRMASK      8
# define    MAX_PAR     8

/*  define the parameter structure */
typedef struct {
    char    tbname[CHAR_FNAME_LENGTH+1];      /* Name of CCDTAB to be read */
    float   radius;
    float   thresh;
    int     nexpnames;
    float   meanexp;
    float   fillval;
    char    sigmas[SZ_FITS_VAL+1];
    char    initgues[SZ_FITS_REC+1];
    char    sky[SZ_FITS_REC+1];
    float   scalense;
    char    expname[SZ_FITS_REC+1];
    char    skyname[SZ_FITS_REC+1];
    int     mask;
    short   crval;
    short   badinpdq;
    int     shadcorr;
    int     printtime;
    int     verbose;
} clpar;

#endif /* INCL_WF3REJ_H */
