#ifndef INCL_CALSTIS6_H
#define INCL_CALSTIS6_H

# include <float.h>

/* calstis6.h 1-D spectral extraction */

/* requires "../stis.h" and <xtables.h> */

/* default values for command-line parameters */

# define NO_RANGE	-1
# define NO_POSITION	-1.0
# define NO_SIZE	-1.0
# define NO_ORDER	-1
# define NO_TILT	-100.

# define ALL_BITS	0x7F	/* default mask: all bits set */
# define PI		3.1415926535897932384626433

# define MAX_PRISM_WAVELENGTH  (6000.)

/* return codes. */

# define FATAL			0
# define NO_FATAL		1
# define INVALID		1001
# define WARN			1
# define NO_WARN		0

/* control for reference file info printout */

# define XTRAC_INFO	1
# define DISP_INFO	2
# define FLUX_INFO	3
# define SGEO_INFO	4
# define CCD_INFO       5

/* constants used in photometry calibration */

# define  HST_AREA  45238.93416		/* cm**2, including obscured areas */
# define  H_PLANCK  (6.6260755e-27)	/* Planck's constant, erg * sec */
# define  C_LIGHT   (2.9979245800e10)	/* cm / sec */
# define  CM_PER_ANGSTROM  (1.e-8)	/* cm / angstrom */

/* heliocentric correction factor */

# define  HCFACTOR 	299792.458

/* pixel fill values and generic rejection DQ flag */

# define PIX_FILL	0.0	/* fill values. Revised OPRs demand that */
# define ERR_FILL	0.0	/* the pixel fill value must be zero     */
# define REJECT		1

/* algorithm selection strings. These must match command line parameters. */

# define	UNWEIGHTED	"UNWEIGHTED"
# define	OPTIMAL		"OPTIMAL"
# define	SCATTER_CORR	"SCATTER_CORR"
# define	IDT_CORR	" + IDT"

/* degree of background fitting polynomial in scatter_corr algorithm. */

# define	BACKP	7

/* Safety margin for profile builder in optimal extraction. This constant
   is used by both the profile builder and the optimal extractor, thus it
   creates a cross dependency between the two tasks. Profile reference
   files must always be created and used under the same version of the
   code.
*/

# define	PROF_MARGIN	4

/* column names in _x1d table. Only the common columns used by both
   standard and idt extractions are defined here.
*/

# define	SPORDER		"SPORDER"
# define	NELEM		"NELEM"
# define	WAVELENGTH	"WAVELENGTH"
# define	GROSS		"GROSS"
# define	BACKGROUND	"BACKGROUND"
# define	NET		"NET"
# define	EXTRLOCY	"EXTRLOCY"

/* Image description and reference files for calstis6. */

typedef struct {

	/* input and output image names */
	char input[STIS_LINE];		/* input image to be calibrated */
	char output[STIS_LINE];		/* output calibrated image */
	char rootname[STIS_CBUF];	/* root name for set of obs */
	char outw[STIS_LINE];		/* output weights image */

	/* command-line parameters */
	double cl_a2center;		/* spectrum center */
	double xoffset;			/* ref star position for slitless data*/
	int maxsearch;			/* cross correlation range */
	double extrsize;		/* size of spectrum box */
	double bksize[2];		/* size of background boxes */
	double bkoffset[2];		/* offset of background boxes */
	double bktilt;			/* angle of background boxes */
	double blazeshift;		/* blaze shift (in pixels) */
	int bkord;			/* order of backr. polyn. fit */
	int sporder;			/* spectral order */
	char xtracalg[STIS_CBUF];	/* extraction algorithm */
	char profilefile[STIS_FNAME];	/* file with profiles for optm. extr. */
	char fluxfile[STIS_FNAME];	/* file with fluxes for optm. extr. */
	int printtime;			/* print time after each step? */
	int verbose;			/* print additional info? */
	int extrloc;			/* output extraction location info? */
	int pipeline;			/* running from pipeline ? */
	int optimal;			/* optimal extraction ? */
	int scatter;			/* scattered light correction ? */
	int variance_image;		/* variance instead of weight image ? */
	int fflux;			/* flux instead of net in ref spec  ? */
	int lfilter;			/* Lee filter window size */
	int imset;			/* selected IMSET (0 -> process all) */
	int bks_mode;			/* bkgr. smoothing mode */
	int bks_order;			/* bkgr. smoothing polynomial order */

	/* info about input image */
	char obsmode[STIS_FITS_REC];
	char aperture[STIS_CBUF];
	char opt_elem[STIS_CBUF];
	char det[STIS_CBUF];	/* name of detector */
	int detector;		/* NUV-MAMA, FUV-MAMA, or CCD */
	int wavecal;		/* true if ASN_MTYP is "WAVECAL" */
	int nimages;		/* number of IMSETs in file */
	int echelle;		/* echelle mode? */
	short sdqflags;		/* serious data quality bits */
	short cc_sdqflags;	/* DQ flag used in crosscor */
	short sdqflags_orig;	/* as read originaly from image header */

	/* CCD-specific info */
	char ccdamp[STIS_CBUF]; /* CCD amplifier read out (A,B,C,D) */
	double readnse;		/* CCD readout noise */
	double rn2;		/* CCD readout variance */
        int ccdoffset;          /* CCD on-board pixel offset */
        int binaxis[2];         /* CCD on-board pixel bin size */
	double atodgain;	/* calibrated CCD gain */
	int ccdgain;		/* commanded CCD gain */
	int crsplit;		/* CCD CRSPLIT value */
	double gain;		/* used in optimal extraction */

	/* values specific to each IMSET */
	double exptime;		/* exposure time (sec) */
	double expstart;	/* exposure start time (MJD) */
	double expend;		/* exposure stop time (MJD) */
	double doppzero;	/* Doppler shift zero phase time (MJD) */
	double doppmag;		/* Doppler shift magnitude (high-res pixels) */
	double orbitper;	/* assumed HST orbital period (seconds) */
	double hfactor;		/* heliocentric correction factor */
        double meandark;        /* mean dark count rate */
        int ncombine;           /* number of combined images */

	/* Detector temperature (or housing temperature for side-2 CCD)
	   is used for the temperature dependence of the sensitivity.
	*/
	double detector_temp;	/* degrees Celsius */

	/* coordinate info */
	double ltm[2];		/* matrix part of transformation to det. */
	double ltv[2];		/* offset to detector coordinates */
	double cd[2];		/* derivatives in CD matrix */
	double crpix[2];	/* WCS reference pixel (physical) */
	double crval[2];	/* WCS value at reference pixel */
	int cenwave;		/* central wavelength */
	int dispaxis;		/* dispersion axis, 1 ==> X, 2 ==> Y */
	double ra_targ;		/* right ascension of target, degrees */
	double dec_targ;	/* declination of target, degrees */
	double pos_targ[2];	/* POS TARG, arcsec */

	/* various offsets */
	double dither_offset[2];/* deliberate MAMA offset */
	double msm_offset[2];	/* MSM slop offset */
	double mama_offset[2];	/* deliberate MAMA "dither" offset */
	double ap_offset[2];	/* aperture offset in pixels */
	double total_offset[2];	/* sum of above offsets */

	/* cross correlation - related quantities */
	double cc_a2center;	/* a2center after cross correlation */
	double crscroff;	/* offset from cross correlation */
	double gcrscroff;	/* global offset from cross correlation */
	double nm_a2center;	/* nominal a2center */
	double cc_peakf;	/* peak flux in cc function */
	double cc_peak2f;	/* 2nd below peak flux */
	double cc_peak3f;	/* 3rd below peak flux */
	double cc_avef;		/* average (peak + 2 neighbors) flux */
	double cc_minf;		/* minimum flux */
	double cc_thresh;	/* threshold for accepting crosscor */
	int cc_global;		/* use global fit everywhere ? */
	double *cc_off;		/* array with individual offsets */
	int *cc_spord;		/* array with individual sporders */
	int *cc_rej;		/* array with flags on rejected sporders */
	int cc_size;		/* size of above arrays */
	int avoid1a;		/* Lya avoidance window (in physical pixels)*/
	int avoid2a;
	int avoid1b;		/* 1300 line avoidance window */
	int avoid2b;

	/* background-related quantities */
	double backval;		/* value from command line */
	double backerr;		/* error from command line */
	double bck[BACKP+3];	/* values */
	double vbck[BACKP+3];	/* variances */
	double ebck;		/* total error (rms around fit) */
	short  dqbck;		/* DQ flags */
	double sin_bktilt;	/* tilt */
	double cos_bktilt;

	/* profile-related quantities */
	int do_profile;		/* build profile ? */
	double **profile;	/* profile array */
	short **profile_dq;	/* array with DQs associated with profile */
	short *profile_rej;	/* array with profile rejection flags */
	double *profile_offset;	/* array with profile offsets */
	double *profile_centroid; /* array with centroid corrections */
	double *profile_rejf;	/* array with flux rejection factors */
	short profile_x;	/* profile array dimension */
	short profile_y;
	char rejranges[STIS_LINE];/* rejection ranges */
	double sclip;		/* sigma-clip used in opt. extr. cleaning */
	double psclip;		/* sigma-clip used in profile building */

	/* subsampled profile */
	int subprof_size;	/* array size */
	double subscale;	/* drop -> drizzled pixel scale factor */
	double **subprof;	/* drizzled profile array */

	/* compressed profile-related quantities */
	int profile_pstep;	/* pixel step used to compress profiles */
	double profile_wstep;	/* wavelength step used to compress profiles */
	double profile_minsn;	/* minimum acceptable S/N */
	double *profile_sn;	/* signal-to-noise in profile section */
	double *profile_minw;	/* array with minimum wavelengths */
	double *profile_maxw;	/* array with maximum wavelengths */
	short *profile_minp;	/* array with minimum pixels */
	short *profile_maxp;	/* array with maximum pixels */
	int profile_msize;	/* X size of above arrays */

	/* calibration switches */
	int wavecorr;		/* has wavecal processing been completed? */
	int x1d;		/* 1-D spectral extraction? (fatal flag) */
	int backcorr;		/* extract and subtract background? */
	int dispcorr;		/* generate wavelength array? */
	int sgeocorr;		/* correct for small-scale distortions? */
	int heliocorr;		/* convert to heliocentric wavelengths? */
	int fluxcorr;		/* convert to absolute flux units? */
        int ctecorr;            /* correct for Charge Transfer Inefficiency */
	/* these are not set from calibration switches; their values are
	   set based on the presence of the associated reference file.
	 */
	int pctcorr;		/* include PCT correction with fluxcorr? */
	int gaccorr;		/* include grating-aperture correction? */
	int tdscorr;		/* include TDS correction with fluxcorr? */
        int wecfcorr;           /* use the stis echelle model? */

	/* control flag for skipping dummy reference entries */
	int x1d_o;

	/* contol flags for printing reference file info only once */
	int dispinfo;		/* these are used by routine Message6 to */
	int fluxinfo;		/* print reference file info only when the */
	int sgeoinfo;		/* file is referenced for the first time */
        int ccdinfo;

	/* calibration images and tables */
        RefTab distntab;        /* SDC or IDC, coordinate & dist info table */
	RefImage sdstfile;	/* SSD, MAMA small-scale distortion file */
	RefTab apdestab;	/* APD, aperture description table */
	RefTab apertab;		/* APT, relative aperture throughput table */
	RefTab phottab;		/* PHT, photometric throughput table */
	RefTab tdstab;		/* TDS, time-dependent sensitivity table */
	RefTab disptab;		/* DSP, dispersion coefficients table */
	RefTab inangtab;	/* IAC, incidence-angle correction table */
	RefTab sptrctab;	/* 1DT, 1-D spectrum trace table */
	RefTab xtrctab;		/* 1DX, 1-D extraction info table */
	RefTab pctab;		/* PCT, photometric correction table */
	RefTab gactab;		/* GAC, grating-aperture correction table */
	RefTab pftab;		/* optimal extraction profile table */
	RefTab pxtab;		/* optimal extraction flux table */
        RefTab ccdtab;          /* CCD, CCD reference table */
	RefTab cfgtab;          /* STIS echelle configuration table */

	/* calibration swicth and files associated with the IDT algorithm */
	int idt;
        RefTab echsctab;
        RefTab exstab;
        RefTab cdstab;
        RefTab riptab;
        RefTab srwtab;
        RefTab halotab;
        RefTab psftab;

	/* parameters used by IDT algorithm. */
	double ap_xsize;	/* aperture size */
	double ap_ysize;

	double trace_rotation;  /* trace rotation angle */

} StisInfo6;



/* CoordInfo describes the coordinates for the output image. Caltsis6
   needs just the plate scale in the cross-dispersion (slit) direction
   to correct for POSTARG2. The value is read from the _sdc table.
   Revised 07Jan00: the scattered light correction algorithm is based
   on the calstis7 output images. Thus we need to store, for each spectral
   order, the number of pixels in the A2 direction spanned by each
   rectified image.
*/

typedef struct coo *CoordPtr;

typedef struct coo {
        int sporder;            /* order number */
        double a2center;        /* Y detector location */
        double cdelt2;          /* increment per pixel */
	int npix;		/* size of rctified image in A2 direction */
        CoordPtr next;          /* pointer to next struct in list */
} CoordInfo;                    /* for info from SDC */



/* XtractInfo stores extraction parameters, read from the _1dx table. */
# define MAX_SLIT_COEFF 8	/* length of array of slit tilt coeffs */
# define MAX_BACK_COEFF 8	/* length of array of background tilt coeffs */

typedef struct xtract *XtractPtr;

typedef struct xtract {
	short sporder;		/* order number */
	float extrsize;		/* size of spectrum extraction box */
	float bksize[2];	/* sizes of background extraction boxes */
	float bkoffset[2];	/* offsets of back. boxes from spect. box */
	double bktcoeff[MAX_BACK_COEFF]; /* coeffs. for background tilt */
	int   ncoeffbk;		/* size of background tilt coeff. array */
	double sltcoeff[MAX_SLIT_COEFF]; /* coeffs. for spectrum slit tilt*/
	int   ncoeffsl;		/* size of spectrum tilt coeff. array */
	short backord;		/* order of polynomial to fit background */
	char xtracalg[STIS_CBUF]; /* extraction algorithm */
	short maxsearch;	/* maximum search range for cross correlation */
	XtractPtr next;		/* pointer to next struct in list */
} XtractInfo;



/* Spectrum trace info, the displacement in the cross-dispersion
   direction at each reference pixel, read from the _1dt table.
*/
# define MAX_SP_TRACE  1024	/* max. length of array for spectrum trace */
# define NO_TRACE      -2	/* no spectrum trace for current order */
# define ERROR_TRACE   -2	/* error in trace table */

typedef struct sptrc *TracePtr;

typedef struct sptrc {
	double a2center;	/* Y location on detector */
	double a1center;	/* X location on detector */
	int nelem;		/* actual size of array */
	double a2displ[MAX_SP_TRACE];	/* spectrum trace */
	TracePtr next;		/* pointer to next struct in list */
} SpTrace;



/* This stores profile arrays used by the optimal extraction algorithm. */

typedef struct prof *ProfilePtr;

typedef struct prof {
	int npts;		/* array size as read from file */
	int nptsoff;		/* offsets array size as read from file */
	double minw;		/* minimum wavelength for this range */
	double maxw;		/* maximum wavelength for this range */
	int minp;		/* minimum pixel for this range */
	int maxp;		/* maximum pixel for this range */
	int xpix;		/* auxiliary storage */
	float sn;		/* signal-to-noise */
	double *prof;		/* profile array */
	double *profoff;	/* offsets array */
	ProfilePtr next;	/* pointer to next struct in list */
} ProfileArray;



/* This contains a description of the slit, including throughput.  Values
   are gotten from the aperture description table (_apd) and aperture
   throughput table (_apt).
*/

typedef struct {
	/* aperture throughput _apt */
	int allocated;		/* true if memory has been allocated */
	int nelem;		/* size of wavelength and throughput arrays */
	double *wl;		/* array of wavelengths */
	double *thr;		/* array of fraction of light passed by slit */
	/* geometric description _apd */
	float ap_offset[2];	/* offset from nominal location in X and Y */
	double width[2];	/* width of slit in X and Y directions */
	/* grating-aperture correction information */
	int gac_allocated;	/* true if memory has been allocated */
	int gac_nelem;		/* size of wavelength and throughput arrays */
	double *gac_wl;		/* array of wavelengths for GAC factors */
	double *gac_thr;	/* array of GAC factors */
} ApInfo;



/* Dispersion relation at A2CENTER, a Y location on the detector.  These
   values are read from the _dsp table.
*/
# define MAX_DISP_COEFF  10	/* length of array of dispersion coeffs */

typedef struct dsp *DispPtr;

typedef struct dsp {
	double a2center;		/* Y location on detector */
	int ncoeff;			/* number of coefficients */
	double coeff[MAX_DISP_COEFF];	/* array of coefficients */
	char ref_aper[STIS_CBUF];	/* ref. aperture for disp. relation */
	DispPtr next;			/* pointer to next struct in list */
} DispRelation;



/* This contains the incidence-angle correction coefficients, which are
   read from the _iac table.
*/

typedef struct {
	int allocated;		/* true if memory has been allocated */
	int ncoeff1;		/* sizes of coefficient arrays */
	int ncoeff2;
	double *coeff1;		/* arrays of coefficients */
	double *coeff2;
} InangInfo;



/* This structure contains the Charge Transfer Inefficiency data for the
   CTI correction algorithm.
 */

typedef struct {
	double ctinorm;         /* gross counts scaling for zero background */
	double ctigpower;       /* gross counts power relation */
	double ctibgfac;        /* exponential roll-off scaling */
	double ctibgpower;      /* exponential roll-off power relation */
	double ctitimfac;       /* time dependence scaling (per year) */
	double ctirefmjd;       /* time dependence reference MJD */
	double spurcharge;      /* spurious charge per pixel (depends on
				   CCD gain & binsize) */
	double halofac;		/* relevant for G750L or G750M */
	double halominfrac;	/* relevant for G750L or G750M */
	int nelem;              /* size of haloterm */
	double *haloterm;       /* array of halo factors (for G750L & G750M) */
	int allocated;          /* true if haloterm has been allocated */
} CTICorrInfo;



/* This stores the intensity/wavelength arrays used by the optimal
   extraction algorithm.
*/

typedef struct {
	int allocated;		/* true if memory has been allocated */
	int nelem;		/* array size as read from file */
	double *wave;		/* wavelength array */
	double *intens;		/* intensity array (in counts) */
} IntensArray;



/* This stores the table descriptors for the main cs6 output table,
   including the main table pointer and pointers to each column
   descriptor.
*/

typedef struct {
	IRAFPointer tp;			/* pointer to table descriptor */
	int nrows;			/* number of rows in table */
	int array_size;			/* number of elements in array */
	IRAFPointer sporder;
	IRAFPointer npts;
	IRAFPointer wave;
	IRAFPointer gross;
	IRAFPointer back;
	IRAFPointer net;
	IRAFPointer flux;
	IRAFPointer error;
	IRAFPointer net_error;	/* this was added 21Oct11 */
	IRAFPointer dq;
        IRAFPointer a2center;   /* these were added in 10Apr98 */
	IRAFPointer extrsize;
	IRAFPointer bk1size;
	IRAFPointer bk2size;
	IRAFPointer bk1offset;	
	IRAFPointer bk2offset;	
	IRAFPointer maxsearch;
	IRAFPointer extrlocy;
	IRAFPointer cc_offset;   /* this was added 23Jun98 */
} TblDesc;



/* This stores the table descriptors for the profile generator output
   table, including the main table pointer and pointers to each column
   descriptor.
*/

typedef struct {
	IRAFPointer tp;			/* pointer to table descriptor */
	int nrows;			/* number of rows in table */
	int array_size;			/* number of elements in array */
	int array_size_off;		/* number of elements in offset array */
	IRAFPointer sporder;
	IRAFPointer npts;
	IRAFPointer nptsoff;
	IRAFPointer subscale;
	IRAFPointer minwave;
	IRAFPointer maxwave;
	IRAFPointer minpix;		/* reference */
	IRAFPointer maxpix;		/* reference */
	IRAFPointer s_n;
	IRAFPointer profoff;
	IRAFPointer profcent;
	IRAFPointer prof;
} ProfTblDesc;



/* This stores the output row contents that go into each row of the
   output table, that is, the actual product of calstis6.

   The same structure is used by the IDT algorithm.
*/

typedef struct {
	short	sporder;
	short	npts;
	double	scale;		/* used by idt algorithm */
	double	*wave;
	float	*gross;
	float	*back;
	float	*net;
	float	*flux;
	float	*error;
	float   *net_error;	/* this was added 21Oct11 */
	short	*dq;
        float	a2center;	/* these were added 10Apr98 */
	float	extrsize;
	float	bk1size;
	float	bk2size;
	float	bk1offset;	
	float	bk2offset;	
	short	maxsearch;
	float	*extrlocy;
	float	cc_offset;	/* this was added 23Jun98 */
} RowContents;

#endif /* INCL_CALSTIS6_H */
