#ifndef INCL_STIS_H
#define INCL_STIS_H

# include <float.h>

/* stis.h generic header for calstis */

# define STIS_CBUF      25	/* small buffer for e.g. rootname */
# define STIS_FNAME   1025
# define STIS_LINE    1025
# define STIS_FITS_REC  81

/* This signals that "no value" was passed as input parameter in the 
   command line to the blazeshift variable in calstis6. This definition
   was moved from calstis6 to here since it is now shared by calstis6, 
   calstis7, calstis0, and a library function (IB, 4/16/02).
*/

# define NO_VALUE	DBL_MAX

/* Background smoothing algorithm modes. This definition was moved from
   calstis6.h to here since now calstis0 must set the background smoothing
   switch (IB, 10/29/02).
*/

# define	BKS_OFF		0
# define	BKS_MEDIAN	1
# define	BKS_AVERAGE	2

/* Three extensions per SingleGroup. */
# define EXT_PER_GROUP 3

/* Integer codes for string-valued header keywords. */

/* Flag values for calibration switches. */
# define DUMMY  (-1)
# define OMIT     0
# define PERFORM  1
# define COMPLETE 2
# define SKIPPED  3
# define IGNORED  4	/* e.g. SHADCORR when the exposure time is zero */

/* Codes to specify whether a reference file exists or not. */
# define EXISTS_YES       1
# define EXISTS_NO        0
# define EXISTS_UNKNOWN (-1)

/* For a reference file name, this string means that a name was
   intentionally not given.
*/
# define NOT_APPLICABLE   "n/a"

/* nearest integer function */
# define NINT(x)  ((x >= 0.) ? (int) (x + 0.5) : (int) (x - 0.5))

/* Codes for goodPedigree in RefImage and RefTab to specify whether
   pedigree is good, dummy, or not specified.
*/
# define GOOD_PEDIGREE      1
# define DUMMY_PEDIGREE     0
# define PEDIGREE_UNKNOWN (-1)

/* Possible values for detector. */
# define UNKNOWN_DETECTOR (-1)
# define NUV_MAMA_DETECTOR  1
# define FUV_MAMA_DETECTOR  2
# define CCD_DETECTOR       3

/* These four are possible values of disp_type. */
# define GRATING_DISP  1	/* first-order grating */
# define PRISM_DISP    2
# define ECHELLE_DISP  3
# define RECTIFIED     4	/* already 2-D rectified */

/* prism dispersion diverges near this wavelength */
# define MAX_PRISM_WAVELENGTH  (6000.)

/* These two are for flagging wavecal shifts as undefined.
   UNDEFINED_SHIFT is the value that CalStis4 puts in the header if
   the shift was not determined.
   Any shift as large as UNREASONABLE_SHIFT must be bad.  Note that
   this is smaller in absolute value than UNDEFINED_SHIFT, so it can be
   used to check for good shifts using:
     if (fabs (shift) < UNREASONABLE_SHIFT) ...
*/
# define UNDEFINED_SHIFT    (-9999.)
# define UNREASONABLE_SHIFT  (1000.)

/* These two are the possible values for the algorithm for
   2-D interpolation of the error array (ERR extension).
   WGT_VARIANCE (weight the variance, this is the default) is:
	sqrt (sum ( w[i] * err[i]^2 )) / sum (w[i])
   while WGT_ERROR (weight the error) is:
	sqrt (sum ( (w[i] * err[i])^2 )) / sum (w[i])
*/
# define WGT_VARIANCE   0
# define WGT_ERROR      1

/* A reference image. */
typedef struct {
	char name[STIS_LINE];		/* name of image */
	char pedigree[STIS_FITS_REC];	/* value of pedigree keyword */
	char descrip[STIS_FITS_REC];	/* value of descrip keyword */
	int exists;			/* does reference image exist? */
	int goodPedigree;		/* DUMMY_PEDIGREE if dummy */
} RefImage;

/* A reference table. */
typedef struct {
	char name[STIS_LINE];		/* name of table */
	char pedigree[STIS_FITS_REC];	/* value of pedigree (header or row) */
	char descrip[STIS_FITS_REC];	/* value of descrip from header */
	char descrip2[STIS_FITS_REC];	/* value of descrip from row */
	int exists;			/* does reference table exist? */
	int goodPedigree;		/* DUMMY_PEDIGREE if dummy */
} RefTab;

/* This section is for saving the names of reference files. */
# define STIS_KEYWORD    10
typedef struct ref *RefFileInfoPtr;
typedef struct ref {
	char keyword[STIS_KEYWORD];
	char filename[STIS_FITS_REC];
	RefFileInfoPtr next;
} RefFileInfo;

/* calibration switches for calstis1. Note that last entry is not a switch,
   but a parameter, required to support the darkscale requirement.
 */

typedef struct {
	int dqicorr;
	int atodcorr;
	int blevcorr;
	int doppcorr;
	int lorscorr;
	int glincorr;
	int lflgcorr;
	int biascorr;
	int darkcorr;
	int flatcorr;
	int shadcorr;
	int photcorr;
	int statcorr;
	char darkscale_string[STIS_LINE];
} cs1_switch;

/*  define the parameter structure for calstis2 */
typedef struct {
	char	tbname[132];
	char	obstype[32];
	float	rej;
	float	psigma;
	int	nexpnames;
	float	meanexp;
	float	fillval;
	char	sigmas[64];
	char	initial[10];
	char	sky[10];
	float	scalenoise;
	char	expname[10];
	char	skyname[10];
	int	mask;
	short	crval;
	short	badbits;
	int	printtime;
	int	verbose;
} clpar;

/* Prototypes for calstisN functions. */

int CalStis0 (char *rawfile, char *wavfile, char *outroot,
	int printtime, int save_tmp, int verbose);

int CalStis1 (char *input, char *output, char *outblev,
	cs1_switch *cs1_sw, RefFileInfo *refnames,
	int printtime, int verbose);

int CalStis2 (char *input, char *fout, clpar *par, int newpar[]);

int CalStis4 (char *input, char *dbgfile,
	RefFileInfo *refnames, int printtime, int verbose, double slit_angle);

int CalStis6 (char *input, char *output,
	int backcorr, int dispcorr, int fluxcorr, int helcorr, int sgeocorr,
        int cticorr, int sc2dcorr,
	double cl_a2center, int maxsearch, double extrsize,
	double bk1size, double bk2size, double bk1offset, double bk2offset,
	double bktilt, int bkord,
	int sporder, char *xtracalg, int printtime, int verbose,
	int extrloc,
	int ccglobal, double ccthresh,
	int do_profile, int pstep, double wstep, double minsn,
	char *rejranges,
	char *profilefile, char *fluxfile, char *outw,
	double backval, double backerr, int variance, int fflux,
	double psclip, double sclip,
	int lfilter, char *idtfile,
	double subscale, double blazeshift, int bks_mode, int bks_oder,
	double stpos, int pipeline);

int CalStis7 (char *input, char *output,
	int sgeocorr, int helcorr, int fluxcorr, int statcorr,
	RefFileInfo *refnames, int printtime, int verbose,
	int center_target, double blazeshift, int err_algorithm);

int CalStis8 (char *input, char *output, int printtime, int verbose);

int CalStis11 (char *inwav, char *insci, char *output,
	int printtime, int verbose);

int CalStis12 (char *inwav, char *insci,
	int w_option, int printtime, int verbose);

/* Prototypes from functions in the lib subdirectory. */

int BinUpdate (double *block, double *offset,
	double *ltm, double *ltv, double *cd, double *crpix);

int checkImsetOK (char *input, int extver, int *imset_ok);

int DefSwitch (char *keyword);

double evalDisp (double coeff[], int ncoeff, double m, double wl);
double prismDisp (double coeff[], double ix_r);
double evalDerivDisp (double coeff[], int ncoeff, double m, double wl);
int evalInvDisp (double coeff[], int ncoeff, double m, double pixel,
		double wl_estimate, double tolerance, double *wl);

int FileExists (char *fname);

int GotFileName (char *filename);

void InitRefTab (RefTab *);

double interp1d (double x, double wl[], double f[], int n, int *starti);
double extrap1d (double x, double wl[], double f[], int n, int *starti);

double MedianDouble (double *v, int n, int inplace);
float MedianFloat (float *v, int n, int inplace);

int MkName (char *input, char *isuffix, char *osuffix,
		char *output, int maxch);

int MkOutName (char *input, char **isuffix, char **osuffix, int nsuffix,
		char *output, int maxch);

int OmitStep (int flag);

double orbitalDopp (double t1, double t2,
	double doppzero, double orbitper, double doppmag);

void PrBegin (int csnumber);
void PrEnd (int csnumber);
void PrFullVersion (void);
void PrVersion (void);
void PrFileName (char *label, char *filename);
void PrHdrInfo (char *obsmode, char *aperture, char *opt_elem, char *detector);
void PrSwitch (char *keyword, int value);
void PrGrpBegin (char *label, int n);
void PrGrpEnd (char *label, int n);
void PrRefInfo (char *keyword, char *filename,
	char *pedigree, char *descrip, char *descrip2);

void InitRefFile (RefFileInfo *ref);
int NewRefFile (RefFileInfo *ref, char *keyword, char *filename);
void FindRefFile (RefFileInfo *ref, char *keyword,
	char *filename, int *foundit);
void FreeRefFile (RefFileInfo *ref);

double RadialVel (double, double, double);

double rotatetrace(double, double, double, double *, int);

int SameInt (int rowvalue, int value);
int SameString (char *rowvalue, char *value);

int splint_nr (double *xa, double *ya, int n,
		double *x, double *y, int nelem);

void pseudoap (char *propaper, char *aperture, int verbose);
int strcmptail (char *s1, char *s2);

int streq_ic (char *s1, char *s2);

void TimeStamp (char *message, char *rootname);

void WhichError (int status);

/* An array of complex values is stored as an array of pairs of floats,
   each pair consisting of the real part followed by the imaginary part.
   The array can be 1-D or 2-D, and the lengths of the first and second
   axes are given by nx and ny respectively.
*/
typedef struct {
	int allocated;		/* true if memory has been allocated */
	float *data;		/* 1-D or 2-D array of data */
	int nx, ny;		/* size of array of data */
	float *workx, *worky;	/* scratch space used by FFT routines */
} CmplxArray;

/* real and imaginary parts of a complex 1-D array */
# define RPIX1D(z,i) ((z)->data[2*i])
# define IPIX1D(z,i) ((z)->data[2*i+1])

/* real and imaginary parts of a complex 2-D array */
# define RPIX2D(z,i,j) ((z)->data[2 * (i + (z)->nx * j)])
# define IPIX2D(z,i,j) ((z)->data[2 * (i + (z)->nx * j) + 1])

/* complex array memory management */
void InitCmplxArray (CmplxArray *);
int AllocCmplxArray (CmplxArray *, int, int);
void FreeCmplxArray (CmplxArray *);

/* forward and inverse Fourier transforms of 2-D complex arrays */
int fft2d (CmplxArray *);
int ifft2d (CmplxArray *);

/* utility functions for complex arrays */
void FFTShift (CmplxArray *);
void CpyCmplx (CmplxArray *, CmplxArray *);

/* These are the Fourier transform functions written by Paul Swarztrauber. */
# if defined(NO_UNDERSCORE)
# define FTDREC_U ftdrec
# define CFFTI cffti
# define CFFTF cfftf
# define CFFTB cfftb
# else
# define CFFTI cffti_
# define CFFTF cfftf_
# define CFFTB cfftb_
# endif

void CFFTI (int *, float *);
void CFFTF (int *, float *, float *);
void CFFTB (int *, float *, float *);

#endif /* INCL_STIS_H */
