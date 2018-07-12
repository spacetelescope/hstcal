#ifndef INCL_ACSREJ_H
#define INCL_ACSREJ_H

/* The header file acs.h must be called prior to calling this
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
    char    tbname[CHAR_FNAME_LENGTH];      /* Name of CCDTAB to be read */
    float   radius;
    float   thresh;
    int     nexpnames;
    float   meanexp;
    float   fillval;
    char    sigmas[SZ_STRKWVAL];
    char    initgues[ACS_FITS_REC];
    char    sky[ACS_FITS_REC];
    float   scalense;
    char    expname[ACS_FITS_REC];
    char    skyname[ACS_FITS_REC];
    int     mask;
    short   crval;
    short   badinpdq;
    int     shadcorr;
    int     printtime;
    int     verbose;
    int     readnoise_only;
} clpar;

#endif /* INCL_ACSREJ_H */
