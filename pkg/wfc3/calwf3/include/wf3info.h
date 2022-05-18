#ifndef INCL_WF3INFO_H
#define INCL_WF3INFO_H

/* wf3info.h

    Warren Hack, 1998 June 10:
    Initial ACS version.

    Howard Bushouse, 2000 August 22:
    Initial WFC3 version.

    H. Bushouse, 2001 May 8:
    Updated to include post-flash correction entries;
    Updated to change filter1, filter2 to just filter.

    H. Bushouse, 2001 Nov 15:
    Updated to support new PHOTCORR methods: added 'graph' and 'comp' to
    replace 'phottab' and 'apertab' and removed 'filtcorr'; simplified
    PhotInfo to only contain a single set of arrays for wave and thru.
    Also, replaced 'ccdbias' with an array 'ccdbias[4]' to be applied
    for each AMP. All changes in accordance with changes in CALACS.

    H. Bushouse, 2002 Jan 28:
    Increased dimension of trimx arrays from 2 to 4 to accomodate extra
    serial virtual overscan regions in WFC3 UVIS CCD's. Added definitions
    for biassectc and biassectd to accomodate extra serial virtual
    overscan regions in WFC3 UVIS CCD's.

    H. Bushouse, 2002 Mar 1:
    Increased dimensions of vx and vy arrays from 2 to 4 to accomodate use
    of separate parallel overscan regions for each WFC3 CCD amp.

    H. Bushouse, 2002 May 8:
    Added NlinData structure definition for use with IR data (copied from
    original calnica.h).

    H. Bushouse, 2002 June 17:
    Removed "statcorr" variable, because statistics are now always computed
    by default. Added expstart and expend. Changed obsmode to obstype.
    (tracking CALACS changes).

    H. Bushouse, 2003 Oct 16:
    Changed datatype of ccdgain from int to float to accomodate WFC3
    non-integer gain settings. Also added expscorr switch and blev array
    (CALACS changes).

    H. Bushouse, 2007 Feb 14:
    Added "subtype" to WF3Info structure for use with IR subarrays.

    H. Bushouse, 2008 Aug 29:
    Added "ncoeff" and "nerr" to NlinData structure.

    H. Bushouse, 2009 Feb 25:
    Added "crrej" to WF3Info structure for the CRREJTAB ref table.

    H. Bushouse, 2009 Oct 21:
    Added "mean_gain" to WF3Info structure for use in FLATCORR.

    H. Bushouse, 2011 Sep 7:
    Replaced graphtab and comptab RefTab entries with new phot table
    in WF3Info struct.

    M. Sosey, 2013 Dec 4:
    Added new FLUXCORR switch for UVIS data to scale chip photometry
    
    M. Sosey 2014 August:
    Added new keywords and step for UVIS CTE correction
    
    M. Sosey  ARG, this is where the wf3info structure is defined for the
    calibration steps. It has all the reference file names from the science
    image. There is ANOTHER wf3info struct, yes, called exactly the same
    thing which is defined in calwf3.h and exists one level from from the
    proccd calling code. It contains only the calibration switches. It would
    be much clearer to call them different things.

    M. Sosey 2015 May:
    Updates for UVIS 2.0
    
    M. sosey 2015 June:
    Updates to fluxcorr for subarray images necessitate adding
    a FLAM variable to the info array so that I can access the 
    value to scale subarray images which have only 1 set of science
    images. When the chip being processed is chip2, then the flam
    for chip1 will be saved in the structure for use in fluxcorr.

    M. De La Pena 2020 March:
    Read the PCTERNOI value from the raw image header for possible 
    use in the CTE reduction.

    M. De La Pena 2022 February:
    Added the satmap variable which rendered the "saturate" variable
    obsolete.  Removed "saturate" as a cleanup operation.
 
*/

# include "msg.h"

# define NAMPS  4    /* Maximum number of amps for a single readout */
# define MAX_MAREADS	26	/* Max number of MultiAccum reads */

/* Dark image interpolation type */
enum DarkTypes_ {MATCH, INTERP, EXTRAP};
typedef enum DarkTypes_ DarkTypes;

/* Data units */
enum DataUnits_ {COUNTS, COUNTRATE};
typedef enum DataUnits_ DataUnits;

/* Structure describing SINGLE CHIP exposure and its reference files */
typedef struct {
    /* input and output image names */
    char input[CHAR_LINE_LENGTH+1];	   /* input image to be calibrated */
    char output[CHAR_LINE_LENGTH+1];        /* output calibrated image */

    char rootname[CHAR_LINE_LENGTH+1];      /* root name for set of obs */

    /* command-line flags */
    int printtime;                  /* print time after each step? */
    int verbose;                    /* print additional info? */

    /* keywords and file names for reference files */
    RefFileInfo *refnames;

    /* info about input image */
    char det[SZ_CBUF+1];	   /* name of detector */
    char aperture[SZ_CBUF+1];      /* aperture name */
    char filter[SZ_CBUF+1];        /* name of filter used */
    char obstype[SZ_CBUF+1];       /* e.g. Imaging, internal */
    int detector;                   /* integer code for detector */
    int chip;
    double chip1_flam;              /*the flam for chip1 from imphttab*/
    double chip2_flam;              /*the flam for chip2 from the imphttab*/
    int ncombine;                   /* number previously summed together */
    int nimsets;                    /* number of "groups" in file */
    int members;                /* # of members associated with this exposure */
    char mtype[SZ_FITS_VAL+1];      /* Role of exposure in association */

    /* Exposure time keywords */
    double expstart, expend;	    /* exposure start and end times (MJD) */

    short sdqflags;                 /* serious data quality values */

    /* coordinate info */
    int subarray;                   /* is current image a subarray? */
    int bin[2];                     /* size of pixel in detector coordinates */
    double offsetx, offsety;        /* LTV1,LTV2 for readout */

    /* CCD-specific info */
    char ccdamp[NAMPS+1];           /* CCD amplifier read out (A,B,C,D) */
    float ccdgain;                  /* commanded gain of CCD */
    int binaxis[2];                 /* BINAXIS1, BINAXIS2 from header */
    int ccdoffset[NAMPS];           /* commanded offset for each amp of CCD */
    float ccdbias[NAMPS];           /* calibrated offset for each amp of CCD */
    float atodgain[NAMPS];     /* actual gain for each amp used to read CCD */
    float readnoise[NAMPS];    /* readout noise for each amp used to read chip*/
    double blev[NAMPS];	/* bias level value fit for each amp from overscan */
    float mean_gain;    /* mean actual gain of all amps */
    int ampx;           /* first column affected by amps on 2/4amp readout*/
    int ampy;           /* first row affected by amps on 2/4amp readout*/
    int trimx[4];       /* Width of overscan to trim off ends of each line */
    int trimy[2];       /* Amount of overscan to trim off ends of each col */
    int vx[4];
    int vy[4];          /* Coordinates of virtual overscan region to use */
    int biassecta[2];   /* Columns to use for leading overscan region amp 1 */
    int biassectb[2];   /* Columns to use for leading overscan region amp 2 */
    int biassectc[2];   /* Columns to use for trailing overscan region amp 1 */
    int biassectd[2];   /* Columns to use for trailing overscan region amp 2 */
    float flashdur;	/* duration of post-flash (in seconds) */
    char flashstatus[SZ_CBUF+1]; /* status of post-flash exposure */
    float pcternoi;     /* Read noise amplitude */

    /* IR-specific info */
    int group;		/* current group being processed */
    int ngroups;	/* total number of groups */
    int nsamp;		/* number of samples */
    char sampseq[SZ_FITS_VAL+1];	/* sample sequence name */
    char subtype[SZ_FITS_VAL+1];	/* IR subarray type */
    double sampzero;	/* sample zero exptime */
    DarkTypes DarkType; /* dark interpolation type */
    int darkframe1;     /* 1st dark ref frame used */
    int darkframe2;     /* 2nd dark ref frame used */
    int ndarks;		/* number of dark ref images */
    double dtimes[MAX_MAREADS]; /* dark image exposure times */
    double exptime[MAX_MAREADS]; /* exposure time for each readout */
    DataUnits bunit[MAX_MAREADS]; /* BUNIT value for each readout */
    int tdftrans[MAX_MAREADS];    /* number of TDF transitions */
    float crthresh;               /* CR rejection threshold */
    float zsthresh;               /* ZSIG signal threshold */
    int samp_rej;                 /* Number of samples to reject */

    /* calibration flags (switches) */
    int dqicorr;        /* data quality initialization */
    int atodcorr;       /* analog to digital correction */
    int blevcorr;       /* subtract bias from overscan */
    int biascorr;       /* subtract bias image */
    int flashcorr;	/* subtract post-flash image */
    int noiscorr;       /* initialize error array? (yes) */
    int darkcorr;       /* subtract dark image */
    int nlincorr;	/* IR detector non-linearity correction */
    int flatcorr;       /* apply flat field */
    int pfltcorr;       /* apply pixel-to-pixel flat */
    int dfltcorr;       /* apply delta flat */
    int lfltcorr;       /* apply low-order flat */
    int shadcorr;       /* correct short exposures for shutter time */
    int photcorr;       /* compute photometry header keyword values */
    int expscorr;	/* calibrate blv_tmp products? */
    int unitcorr;	/* convert units to count rate */
    int zsigcorr;	/* IR zero-read signal correction */
    int zoffcorr;	/* IR zero-read subtraction */
    int crcorr;		/* IR CR rejection */
    int fluxcorr;   /*Uvis chip flux correction*/
    int pctecorr;  /*UVIS chip CTE correction */

    /* calibration images and tables */
    RefImage bias;      /* bias image */
    RefImage sink;     /*the sink pixel file*/
    RefImage biac;     /*biasc image for cte correction*/ 
    RefImage flash;  	/* post-flash image */
    RefTab bpix;        /* bad pixel table */
    RefTab ccdpar;      /* CCD parameters table */
    RefTab oscn;        /* Overscan parameters table */
    RefTab atod;        /* analog to digital correction table */
    RefTab crrej;       /* CR rejection parameters table */
    RefImage dark;      /* dark image */
    RefImage pflt;      /* pixel-to-pixel flat field */
    RefImage dflt;      /* delta flat */
    RefImage lflt;      /* low-order flat */
    RefImage shad;      /* shutter shading correction image */
    RefTab phot;    	/* photometry table */
    RefImage nlin;	    /* non-linearity coefficients image */
    RefImage zoff;      /* super zero  */
    RefTab pctetab;     /*uvis CTE parameter table*/
    RefImage satmap;    /* full-well saturation image */
    
} WF3Info;

/* This contains the throughput curve (from _pht table) and the aperture
    (filter) throughput curve (from _apt table), which we use for computing
    PHOTFLAM, etc.
*/

typedef struct {
    /* These three are for filter throughputs. */
    int nelem;          /* size of arrays */
    float *f_thru;      /* array of filter throughputs */
    float *f_wl;        /* array of wavelengths from _apt */
    float *f_err;	/* array of throughput errors */
} PhotInfo;

/* Expanded form of SingleGroupLine to represent an image section */
typedef struct {
    char *filename;     /* Name of input file	*/
    Hdr *globalhdr;     /* header information */
    int group_num;      /* EXTVER number		*/
    int line_num;       /* line number of FIRST line loaded	*/
    int nlines;         /* Number of lines in section */
    int npix;           /* Number of pixels per line */
    Bool phdr_loaded;   /* Global Header already loaded? */
    SciHdrLine *sci;    /* Array of science data lines */
    ErrHdrLine *err;    /* Array of error data lines */
    DQHdrLine *dq;      /* Array of DQ data lines */
} WF3sect;

/* IR Non-linearity reference data structure */
typedef struct {
	int	      ncoeff;
	int	      nerr;
	Hdr	     *globalhdr;
	FloatHdrData *coeff;
	FloatHdrData *error;
	ShortHdrData *dqual;
	FloatHdrData *nodes;
	FloatHdrData *zsci;
	FloatHdrData *zerr;
} NlinData;


#endif /* INCL_WF3INFO_H */
