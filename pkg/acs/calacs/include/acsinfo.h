/* acsinfo.h

    Warren Hack, 1998 June 10:
    Initial ACS version.

    2001-10-29 WJH: Replaced 'ccdbias' with an array 'ccdbias[4]' to
        be applied for each AMP.  Added 'graph' and 'comp' to
        replace 'phottab' and 'apertab' and removed 'filtcorr'. Simplified
        PhotInfo to only contain a single set of arrays for wave and thru.
    2001-12-04 WJH: Added expstart and expend.
    2017-02-21 PLL: Added SINKCORR varibles.
*/

# define NAMPS  4    /* Maximum number of amps for a single readout */

/* Structure describing SINGLE CHIP exposure and its reference files */
typedef struct {
    /* input and output image names */
    char input[ACS_LINE+1];       /* input image to be calibrated */
    char output[ACS_LINE+1];        /* output calibrated image */

    char rootname[ACS_LINE+1];      /* root name for set of obs */

    /* command-line flags */
    int printtime;                  /* print time after each step? */
    int verbose;                    /* print additional info? */
    int onecpu;                     /* turn off OpenMP usage if True */

    /* keywords and file names for reference files */
    RefFileInfo *refnames;

    /* info about input image */
    char det[ACS_CBUF+1];	        /* name of detector */
    char aperture[ACS_CBUF+1];      /* aperture name */
    char filter1[ACS_CBUF+1];       /* name of filter used */
    char filter2[ACS_CBUF+1];       /* name of filter used */
    char obstype[ACS_CBUF+1];   /* e.g. IMAGING, CORONOGRAPHIC, INTERNAL */
    char jwrotype[ACS_CBUF+1];      /* for WFC: DS_int or CLAMP */
    int detector;                   /* integer code for detector */
    int chip;
    int ncombine;                   /* number previously summed together */
    int nimsets;                    /* number of "groups" in file */
    int members;                /* # of members associated with this exposure */
    char mtype[SZ_STRKWVAL+1];      /* Role of exposure in association */

    /* Exposure time keywords */
    double exptime;                 /* exposure time */
    double expstart,expend;         /* exposure start and end times (MJD)*/

    short sdqflags;                 /* serious data quality values */

    /* coordinate info */
    int subarray;                   /* is current image a subarray? */
    int bin[2];                     /* size of pixel in detector coordinates */
    double offsetx, offsety;        /* LTV1,LTV2 for readout */

    /* MAMA-specific info (regarding linearity) */
    double global_limit;            /* count rate for 10% global nonlinearity */
    double tau;                     /* time constant for global nonlinearity */
    double local_limit;              /* count rate for 10% local nonlinearity */
    float expand;                   /* radius for flagging local nonlinearity */

    /* MAMA-specific info */
    double globrate;                /* global count rate */

    /* CCD-specific info */
    char ccdamp[NAMPS+1];           /* CCD amplifier read out (A,B,C,D) */
    float ccdgain;                    /* commanded gain of CCD */
    int binaxis[2];                 /* BINAXIS1, BINAXIS2 from header */
    int ccdoffset[NAMPS];           /* commanded offset for each amp of CCD */
    float ccdbias[NAMPS];           /* calibrated offset for each amp of CCD */
    float atodgain[NAMPS];      /* actual gain for each amp used to read CCD */
    float readnoise[NAMPS];     /* readout noise for each amp used to read chip*/
    double blev[NAMPS];  /* bias level value fit for each amp from overscan*/
    int ampx;           /* first column affected by amps on 2/4amp readout*/
    int ampy;           /* first row affected by amps on 2/4amp readout*/
    float saturate;     /* CCD saturation level */
    int trimx[2];       /* Width of overscan to trim off ends of each line */
    int trimy[2];       /* Amount of overscan to trim off ends of each col */
    int vx[2];
    int vy[2];          /* Coordinates of virtual overscan region to use */
    int biassecta[2];   /* Columns to use for leading overscan region */
    int biassectb[2];   /* Columns to use for trailing overscan region */
    float flashdur; 	/* duration of post-flash (in seconds) */
    char flashstatus[ACS_CBUF+1];		/* status of post-flash exposure */

    /* calibration flags (switches) for ACSCCD */
    int dqicorr;        /* data quality initialization */
    int atodcorr;       /* analog to digital correction */
    int blevcorr;       /* subtract bias from overscan */
    int biascorr;       /* subtract bias image */
    int flashcorr;      /* subtract post-flash image */
    int noisecorr;      /* initialize error array? (yes) */
    int sinkcorr;       /* flag sink pixels */

    /* calibration flags (switches) for ACSCTE */
    int pctecorr;       /* perform pixel CTE correction */

    /* calibration flags (switches) for ACS2D */
    int glincorr;       /* global nonlinearity correction */
    int lflgcorr;       /* flag local nonlinearity */
    int darkcorr;       /* subtract dark image */
    int flatcorr;       /* apply flat field */
    int pfltcorr;       /* apply pixel-to-pixel flat */
    int dfltcorr;       /* apply delta flat */
    int lfltcorr;       /* apply low-order flat */
    int cfltcorr;       /* apply spot flat */
    int shadcorr;       /* correct short exposures for shutter time */
    int photcorr;       /* compute photometry header keyword values */
    int expscorr;       /* calibrate blv_tmp products?  */

    /* calibration images and tables for ACSCCD */
    RefImage bias;      /* bias image */
    RefTab bpix;        /* bad pixel table */
    RefTab ccdpar;      /* CCD parameters table */
    RefTab oscn;        /* Overscan parameters table */
    RefTab atod;        /* analog to digital correction table */
    RefTab spot;        /* Spotflat offset table */
    RefImage sink;      /* sink pixel image */

    /* calibration images and tables for ACSCTE */
    RefTab pcte;        /* Pixel CTE parameters table */

    /* calibration images and tables for ACS2D */
    RefImage dark;      /* dark image */
    RefImage darkcte;   /* cte corrected dark */
    RefImage flash;     /* post-flash image */
    RefImage flashcte;  /* cte corrected post-flash */
    RefImage pflt;      /* pixel-to-pixel flat field */
    RefImage dflt;      /* delta flat */
    RefImage lflt;      /* low-order flat */
    RefImage cflt;      /* coronographic (spot) flat */
    RefImage shad;      /* shutter shading correction image */
    RefTab mlin;        /* MAMA nonlinearity info table */
    RefTab phot;        /* photometry table (processed from pysynphot) */
} ACSInfo;

/* This contains the throughput curve (from _pht table) and the aperture
    (filter) throughput curve (from _apt table), which we use for computing
    PHOTFLAM, etc.
*/

typedef struct {
    /* These three are for filter throughputs. */
    int nelem;        /* size of arrays */
    float *f_thru;     /* array of filter throughputs */
    float *f_wl;       /* array of wavelengths*/
    float *f_err;       /* array of throughput errors */
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
} ACSsect;
