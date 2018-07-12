# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <math.h>

# include "xtables.h"
# include "ximio.h"
# include "hstio.h"

# include "stis.h"
# include "stisdef.h"
# include "hstcalerr.h"
# include "calstis6.h"
# include "idtalg.h"

static int AddGhost (Hdr *, float **, int, int);
static void BuildTempNames (char *, char *, char *, char *);
static int CheckIDT (Hdr *, StisInfo6 *, int *);
static void FillArray (float **, float **, float **, int, int,
                       double, double, double, double, double *, int);
static int GetMatchingOrder (int, ScatterFunctions *);
static int GetScatterOrder (int, ScatterFunctions *);
static double Interpolate (double, int, double *, double *, int, int *);
static double Linear (double, double, double, double, double);
static double MedianRowsPerAngstrom (RowContents **, int);
static int RedoX1DFile (char *, RowContents **, int);
static double Select (unsigned long, unsigned long, double *);

/* Not used */
/*
static double BLog (double);
static int Debug (char *, float **, int, int);
static int DDebug (char *, double *, int);
static int CDebug (char *, CmplxArray *);
static void CReport (CmplxArray *);
static void Report (float **, int, int);
*/

/*
   calstis6idt -- IDT 1-D extraction with scattered light removal.

   Revision history:
   ----------------
   16 Feb 00  -  Implemented (I.Busko)
   25 Jul 00  -  Got it to work ! (IB)
   19 Apr 01  -  Fix verbosity (IB)
    8 Mar 02  -  Again, fix verbosity (IB)
   16 Apr 02  -  Blaze shift from command line (IB)
   24 Jul 02  -  Background smoothing parameters (IB)
   21 Aug 02  -  Fix temporary file names in subdirectories (IB)
    5 Dec 03  -  Fixed bug in IDT status messages (PB)
   28 Jul 04  -  Replace UpdHdrSwitch was a call to Put_KeyS (PEH)
   01 Jun 05  -  Initialize slit.gac_allocated to 0 (PEH)
   15 Aug 06  -  In the last call to CalStis6Std, pass verbose+1
                 instead of verbose (PEH)
    6 Oct 06  -  Add xoffset to the calling sequence, pass to CalStis6Std.
    9 Oct 08  -  Also copy sporder to gx1d; modify RedoX1DFile to match
                 original gross with final output row using spectral order
                 instead of row number.
    3 Nov 08  -  Call checkImsetOK to get imset_ok, and skip imset if F (PEH).
*/

int CalStis6IDT (char *input, char *output, Hdr *phdr, int backcorr,
                 int dispcorr, int fluxcorr, int helcorr, int sgeocorr,
                 int ctecorr,
                 double cl_a2center, int maxsearch, double extrsize,
                 double bk1size, double bk2size, double bk1offset,
                 double bk2offset, double bktilt, int bkord, int sporder,
                 char *xtracalg, int printtime, int verbose, int extrloc,
                 int ccglobal, double ccthresh, int do_profile, int pstep,
                 double wstep, double minsn, char *rejranges,
                 char *profilefile, char *fluxfile, char *outw, double backval,
                 double backerr, int variance, int fflux, double psclip,
                 double sclip, int lfilter, int imset, int pipeline,
                 char *idtfile, double blazeshift,
                 int bks_mode, int bks_order, double xoffset) {

/* arguments
char *input;		i: input file name
char *output;		i: output file name
Hdr *phdr       	i: primary header
int  backcorr;          i: calibration switch
int  dispcorr;          i: calibration switch
int  fluxcorr;          i: calibration switch
int  helcorr;           i: calibration switch
int  sgeocorr;          i: calibration switch
double cl_a2center;	i: nominal A2 spectrum center
int maxsearch;		i: cross correlation range
double extrsize;	i: size of spectrum box
double bk1size;		i: size of background box 1
double bk2size;		i: size of background box 2
double bk1offset;	i: offset of background box 1
double bk2offset;	i: offset of background box 2
double bktilt;		i: angle of background boxes
int bkord;		i: order of backr. polyn. fit
int sporder;		i: spectral order
char *xtracalg;		i: extraction algorithm
int printtime;		i: print time after each step?
int verbose;		i: print additional info ?
int extrloc;		i: print extraction location info ?
int ccglobal;		i: use global crosscor everywhere ?
double ccthresh;	i: crosscor theshold
int do_profile;		i: use calstis6 as profile generator ?
int pstep;		i: pixel step used to compress profiles
double wstep;		i: wavelength step used to compress profiles
double minsn;		i: minimum acceptable S/N
char *rejranges;	i: rejection ranges
char *profilefile;      i: profiles for optimal extractions
char *fluxfile;         i: fluxes for optimal extractions
char *outw;		i: output weights file name
double backval;		i: background value from command line
double backerr;		i: background error from command line
int variance;		i:variance image instead of weight image ?
int *fflux;		i: use FLUX instead of NET in opt. extraction
double psclip;          i: sigma-clip in profile builder
double sclip;           i: sigma-clip in extraction
int lfilter;		i: Lee filter window size
int imset;		i: selected IMSET (0 -> process entire file)
int pipeline;		i: is calstis6 being run from the pipeline ?
char *idtfile		i: file with IDT final deconvolved image
double blazeshift;	i: blaze shift (in pixels) from command line
int bks_mode;		i: backgr. smoothing mode
int bks_order;		i: backgr. smoothing polynomial order
*/
	StisInfo6 sts_idt;		/* stores calibration file names */
	IODescPtr tdesc;		/* x1d table primary header descrip. */
	Hdr tphdr;			/* x1d table primary header */
	SingleGroup in;			/* input image */
	SingleGroup win;		/* working input image */
	SingleGroup wout;		/* working output image */
	TblDesc tabptr;			/* _x1d temporary table */
	RowContents **x1d;		/* array with table data. */
	RowContents **wx1d;		/* working array with table data. */
	RowContents **wx1d2;		/* working array with table data. */
	RowContents ***gx1d;		/* array with GROSS values. */
	int *gx1d_nrows;		/* size of each gx1d array */
	char temp_x1d[STIS_FNAME];	/* _x1d temporary data file */
	char temp_ima[STIS_FNAME];	/* temporary image file */
	char temp_ima1[STIS_FNAME];	/* temporary image file */
	char temp_out[STIS_FNAME];	/* temporary output table */
	char out_ima[STIS_FNAME];	/* output image file */
	float **bpos;			/* positions of extr. background */
	char opt_elem[STIS_CBUF];	/* grating */
	double ltm[2],ltv[2];		/* LTM values */
	int use_def = 1;                /* use default if missing keyword */
	int no_default = 0;             /* missing keyword is fatal error */
	int imset_ok;			/* value of header keyword IMSET_OK */
	int first_extver;		/* true when processing first imset */
	double exptime;			/* exposure time */
	ApInfo slit;			/* dummy slit */
	double ap_xsize, ap_ysize;	/* aperture size */
	int cenwave;			/* CENWAVE keyword value */
	double linepos[3];
	double mpsfpos[3];
	double *mline;
	double dydw;			/* median rows per Angstrom */
	float **wbig;			/* extended wavelength array */
	float **blaze;			/* extended ripple array */
	float **fbig;			/* extended flux (net) array */
	double inc_scale;		/* algorithm parameters */
	double yyamp;
	ScatterFunctions scf;		/* array with scattering functions */
	Spliced merge;			/* spliced spectrum */
	int iter;			/* current iteration */
	float **im_mod;			/* scattering model */
	float **im_mod1;
	float **im_mod2;
	float **im_mod3;
	float **o_mod;			/* model image without scattering */
	float **o_mod1;
	float **o_mod2;
	float **o_mod3;
	double *scale_lsf;		/* scattering function for curr. order*/
	double *scale_lsf1;		/* above with clipped peak */
	double *scale_lsf2;		/* and wings only */
	int lsf_off;			/* pixel offset in profile */
	double *echelle_scat_offset;
	int status, missing;
	int nextend, nimages, extver, extver0, i, j, i1, i2, k, k1, k2;
	int iorder, ipix, npix, mrip;
	double disp, dw, y;
	double *hold1, *hold2, *hold3;
	float *eblaze, *eonimage, *worder, *forder, *yorder_onimage;
	double *worder_onimage;
	double f, fr2; /* clumsy names ! kept to match IDT code */
	double fr1;
	int y1, ixx, iy1, iy2, ilsf;
	int istart, istop, image1pos, image2pos, ilsf1, ilsf2;
	int nx, ny;
	CmplxArray zhold;		/* buffers for convolution */
        IODescPtr tim;           	/* temporary primary header */
        int prtimestamp;
	double dummy;

	float **Alloc2DArrayF (int, int);
	RowContents **AllocX1DTable (int);
	int EchScatRead (Hdr *, double, double, ScatterFunctions *, int);
	int FFTConvolve2 (CmplxArray *, CmplxArray *);
	void FreeScatter (ScatterFunctions *);
	void Free2DArrayF (float **, int);
	void FreeThroughput6 (ApInfo *);
	void FreeX1DTable (RowContents **, int);
	int GetApDes6 (StisInfo6 *, ApInfo *);
	int GetFlags6 (StisInfo6 *, Hdr *);
	int GetKeyInfo6 (StisInfo6 *, Hdr *);
	int GetX1DTable (TblDesc *, RowContents **);
	int History6 (StisInfo6 *, Hdr *, int);
	int OpenX1DTable (char *, int, TblDesc *);
	void Float2Cmplx (float **, int, int, CmplxArray *);
	int RebinData (SingleGroup *, SingleGroup *, RowContents **,
                       RowContents **, int, int);
	int Splice (RowContents **, int, double, Spliced *);
	void StisInit6 (StisInfo6 *);
	void Cmplx2Float (CmplxArray *, float **, int, int);
	int CalStis6Std (char *, char *, int, int, int, int, int, int,
                         double, int,
                         double, double, double, double, double, double,
                         int, int, char *, int, int, int, int, double, int,
                         int, double, double, char *, char *, char *, char *,
                         double, double, int, int, double, double, int, int,
                         int, double *, double *, double, double, int, int,
                         double);

	/* Print ID greeting. */
	prtimestamp = 0;

	if (verbose == 1 || verbose == 2 || verbose == 3) {
	    PrBegin (6);
	    if (verbose == 1) {
	        TimeStamp ("CALSTIS-6 started", "");
	        prtimestamp = 1;
	    }
	    printf ("\n2-D scattering correction algorithm activated.\n\n");
	    fflush (stdout);
	}

	/* Check whether output file already exist. */

	if ((status = FileExists (output)))
	    return (status);

	/* Initialize. */

	StisInit6 (&sts_idt);
	sts_idt.pipeline  = pipeline;
	sts_idt.backcorr  = backcorr;
	sts_idt.dispcorr  = dispcorr;
	sts_idt.fluxcorr  = fluxcorr;
	sts_idt.heliocorr = helcorr;
	sts_idt.sgeocorr  = sgeocorr;
	InitCmplxArray (&zhold);

	/* Find out how many extensions there are in this file. */

	if ((status = Get_KeyI (phdr, "NEXTEND", use_def, EXT_PER_GROUP,
                                &nextend)))
	    return (status);

	/* Convert number of extensions to number of SingleGroups. */

	nimages = nextend / EXT_PER_GROUP;
	if (nextend != nimages * EXT_PER_GROUP) {
	    printf (
            "ERROR    NEXTEND must be a multiple of %d\n", EXT_PER_GROUP);
	    return (HEADER_PROBLEM);
	}

	/*  Alloc structure to hold the original GROSS data. */

	gx1d = (RowContents ***) malloc (nimages * sizeof (RowContents **));
	if (gx1d == NULL)
	    return (OUT_OF_MEMORY);
	gx1d_nrows = (int *) malloc (nimages * sizeof (int));
	if (gx1d_nrows == NULL)
	    return (OUT_OF_MEMORY);

	/* Get CENWAVE */

	if ((status = Get_KeyI (phdr, "CENWAVE", no_default, 0, &cenwave)))
	    return (status);

	/* Check status of reference files specific of the IDT algorithm.
           This also fills up the sts structure with file name strings.
           These are used later to generate HISTORY records.
        */
	missing = 0;
	sts_idt.idt = PERFORM;
	if ((status = CheckIDT (phdr, &sts_idt, &missing))) {
	    if (status == KEYWORD_MISSING) {
	        printf (
"ERROR    Input file header may not have keywords specific to SC2DCORR.\n");
	        printf (
"         Process it with task 'sc2dref'to install the required keywords.\n");
	        fflush (stdout);
	    }
	    return (status);
	}
	if (missing) {
	    return (CAL_FILE_MISSING);
	}

	/* Get aperture size. This code re-does operations already
           done inside CalStis6Std, but must be executed here so the kernel
           FT computations in EchScatRead can be performed before the
           first call to CalStis6Std, which is inside the loop over IMSETs.
           Since the Get functions also check for consistency and do other
           useful stuff, it is OK to use then here so the execution
           aborts before the time-consuming FT computation even starts.
        */
	if ((status = Get_KeyS (phdr, "OPT_ELEM", no_default, "",
                                opt_elem, STIS_CBUF)))
	    return (status);
	strcpy (scf.opt_elem, opt_elem);
	if ((status = GetKeyInfo6 (&sts_idt, phdr)))
	    return (status);
	if ((status = GetFlags6 (&sts_idt, phdr)))
	    return (status);
	slit.allocated  = 0;
	slit.gac_allocated  = 0;
	if ((status = GetApDes6 (&sts_idt, &slit)))
	    return (status);
	ap_xsize = sts_idt.ap_xsize;
	ap_ysize = sts_idt.ap_ysize;

	/* Read scattering functions and build kernels. */

	if (verbose) {
	    printf ("Read scattering functions and build kernels.\n");
	    fflush (stdout);
	}
	if ((status = EchScatRead (phdr, ap_xsize, ap_ysize, &scf, verbose)))
	    return (status);

	/* Loop thru IMSETS in input file. */

	first_extver = 1;		/* true */
	for (extver = 1;  extver <= nimages;  extver++) {
	    extver0 = extver - 1;

	    if (verbose && (nimages > 1)) {
	        printf ("Begin processing IMSET %d\n", extver);
	        fflush (stdout);
	    }

	    x1d           = NULL;
	    gx1d[extver0] = NULL;
	    tabptr.tp     = NULL;

	    status = checkImsetOK (input, extver, &imset_ok);
	    if (status < 0)
		continue;
	    else if (status > 0)
                return status;
	    if (!imset_ok) {
		printf ("Warning  imset %d skipped (IMSET_OK = F)\n",
			extver);
		continue;
	    }

	    /* Build names for temporary files. */

	    BuildTempNames (output, temp_x1d, temp_ima, temp_ima1);

	    /* Create _x1d file for the current IMSET. Make sure the
               EXTRLOCY column gets written and heliocentric correction
               is turned off. Fluxcor is also turned off to force extraction
               of orders not matched by the photometric reference tables.
            */
	    if (verbose) {
	        printf ("Begin extraction of 1-D data.\n");
	        fflush (stdout);
	    }
	    if ((status = CalStis6Std (input, temp_x1d, backcorr, dispcorr,
                         OMIT, OMIT, sgeocorr, ctecorr, cl_a2center, maxsearch,
                         extrsize, bk1size, bk2size, bk1offset, bk2offset,
                         bktilt, bkord, sporder, xtracalg, printtime,
                         verbose, 1, ccglobal, ccthresh, do_profile, pstep,
                         wstep, minsn, rejranges, profilefile, fluxfile, outw,
                         backval, backerr, variance, fflux, psclip, sclip,
                         lfilter, extver, pipeline, &ap_xsize, &ap_ysize, 1.0,
                         blazeshift, bks_mode, bks_order, xoffset)))
	        return (status);

	    if (verbose) {
	        printf ("End extraction of 1-D data.\n");
	        fflush (stdout);
	    }

	    /* Get image data from current IMSET. */

	    if (verbose) {
	        printf ("Read input image.\n");
	        fflush (stdout);
	    }
	    initSingleGroup (&in);
	    getSingleGroup (input, extver, &in);
	    if (hstio_err())
		return (OPEN_FAILED);

	    /* Open table file created by CalStis6Std */

	    if (verbose) {
	        printf ("Read 1-D data.\n");
	        fflush (stdout);
	    }
	    if ((status = OpenX1DTable (temp_x1d, 0, &tabptr)))
	        return (status);

	    /* Alloc memory for storing table contents. */

	    if ((x1d = AllocX1DTable (tabptr.nrows)) == NULL)
	        return (OUT_OF_MEMORY);

	    /* Get data from table and store in RowContents array. */

	    if ((status = GetX1DTable (&tabptr, x1d))) {
	        FreeX1DTable (x1d, tabptr.nrows);
	        c_tbtclo (tabptr.tp);
	        remove (temp_x1d);
	        return (status);
	    }

	    /* Copy original GROSS array into safe storage. */

	    if ((gx1d[extver0] = AllocX1DTable (tabptr.nrows)) == NULL)
	        return (OUT_OF_MEMORY);
	    gx1d_nrows[extver0] = tabptr.nrows;
	    for (j = 0; j < tabptr.nrows; j++) {
	        gx1d[extver0][j]->sporder = x1d[j]->sporder;
	        gx1d[extver0][j]->npts = x1d[j]->npts;
	        gx1d[extver0][j]->gross = (float *) calloc (
                                                   gx1d[extver0][j]->npts,
                                                   sizeof(float));
	        if (gx1d[extver0][j]->gross == NULL) {
	            printf ("Not enough memory to allocate data arrays.\n");
	            return (OUT_OF_MEMORY);
	        }
	        for (i = 0; i < gx1d[extver0][j]->npts; i++)
	            gx1d[extver0][j]->gross[i] = x1d[j]->gross[i];

	    }

	    /* Get image sampling and exposure time. */

            if ((status = GetLT0 (&in.sci.hdr, ltm, ltv)))
                return (status);
	    if ((status = Get_KeyD (&in.sci.hdr, "EXPTIME", no_default, 0.,
                                    &exptime)))
                return (status);

	    /* Rebin image and x1d data to handle hires data. */

	    if (verbose) {
	        printf ("Rebin/copy image to work area.\n");
	        fflush (stdout);
	    }
	    initSingleGroup (&win);
	    if ((wx1d = AllocX1DTable (tabptr.nrows)) == NULL)
	        return (OUT_OF_MEMORY);
	    if (ltm[0] > 1.0) {
	        allocSingleGroup (&win, in.sci.data.nx / 2,
                                        in.sci.data.ny / 2, True);
	        if (RebinData (&in, &win, x1d, wx1d, 2, tabptr.nrows))
	            return (ERROR_RETURN);
	    } else {
	        allocSingleGroup (&win, in.sci.data.nx,
                                        in.sci.data.ny, True);
	        if (RebinData (&in, &win, x1d, wx1d, 1, tabptr.nrows))
	            return (ERROR_RETURN);
	    }
	    copyHdr (&(win.sci.hdr), &(in.sci.hdr));

            /* Rebinned image header must be updated so calstis6
               doesn't get confused...
            */
	    if (ltm[0] > 1.0) {
	        if ((status = Put_KeyD (&(win.sci.hdr), "LTV1", 0.0, "")))
	            return (status);
	        if ((status = Put_KeyD (&(win.sci.hdr), "LTV2", 0.0, "")))
	            return (status);
	        if ((status = Put_KeyD (&(win.sci.hdr), "LTM1_1", 1.0, "")))
	            return (status);
	        if ((status = Put_KeyD (&(win.sci.hdr), "LTM2_2", 1.0, "")))
	            return (status);
	    }

	    /* Original extraction data is not needed anymore, since the
               GROSS arrays are already secure.
            */
	    FreeX1DTable (x1d, tabptr.nrows);

	    if (verbose) {
	        printf ("Compute algorithm parameters.\n");
	        fflush (stdout);
	    }

	    /* BPOS (FROM THE IDT CODE) IS NEVER USED ! */

 	    /* Compute background positions (midpoints between extrlocy).
               To make life simple, this code assumes that all rows in
               the x1d table store arrays of the same size.
            */
	    if ((bpos = Alloc2DArrayF (wx1d[0]->npts,
                                       tabptr.nrows + 1)) == NULL)
	        return (OUT_OF_MEMORY);
	    for (j = 1; j < tabptr.nrows-1; j++) {
	        for (i = 0; i < wx1d[0]->npts; i++)
	            bpos[j][i] = (wx1d[j-1]->extrlocy[i] +
                                  wx1d[j+1]->extrlocy[i]) / 2.0;
	    }
	    /* Lines at both extremes are computed by linear extrapolation. */
	    for (i = 0; i < wx1d[0]->npts; i++) {
	        bpos[0][i] = 2.0 * wx1d[0]->extrlocy[i] -
                                   wx1d[1]->extrlocy[i];
	        bpos[tabptr.nrows-1][i] = 2.0 *
                                          wx1d[tabptr.nrows-2]->extrlocy[i] -
                                          wx1d[tabptr.nrows-3]->extrlocy[i];
	    }

	    /* Define internal algorithm parameters based on
               OPT_ELEM from header.
            */

	    if (streq_ic (opt_elem, "E140M")) {
	        for (j = 0; j < tabptr.nrows; j++) {
	            wx1d[j]->scale = 14.5 - 0.275 * wx1d[j]->sporder +
                                     0.001435 *
                                     wx1d[j]->sporder * wx1d[j]->sporder;
                    wx1d[j]->scale = (wx1d[j]->scale < 2.5) ?
                                     2.5 : wx1d[j]->scale;
	        }
	        inc_scale = 1.4;
	        yyamp = 1.0;
	    } else if (streq_ic (opt_elem, "E140H")) {
	        for (j = 0; j < tabptr.nrows; wx1d[j++]->scale = 1.5);
	        inc_scale = 1.3;
	        yyamp = 1.25;
	    } else if (streq_ic (opt_elem, "E230M")) {
	        for (j = 0; j < tabptr.nrows; wx1d[j++]->scale = 1.0);
	        inc_scale = 1.0;
	        yyamp = 0.0;
	    } else if (streq_ic (opt_elem, "E230H")) {
	        for (j = 0; j < tabptr.nrows; wx1d[j++]->scale = 1.5);
	        inc_scale = 1.0;
	        yyamp = 0.0;
	    } else {
	        printf ("Non supported grating.\n");
	        return (ERROR_RETURN);
	    }

	    /* Get rid of temporary x1d table. */

	    c_tbtclo (tabptr.tp);
	    remove (temp_x1d);

	    /* Find integer row and fractional order where scattering
               kernels are centered.
            */

	    hold1 = (double *) malloc (tabptr.nrows * sizeof (double));
	    hold2 = (double *) malloc (tabptr.nrows * sizeof (double));
	    hold3 = (double *) malloc (tabptr.nrows * sizeof (double));
	    if (hold1 == NULL || hold2 == NULL || hold3 == NULL)
	        return (OUT_OF_MEMORY);
	    for (j = 0; j < tabptr.nrows; j++) {
	        hold1[j] = (double)(wx1d[j]->wave[(wx1d[j]->npts)/2]);
	        hold2[j] = (double)(wx1d[j]->extrlocy[(wx1d[j]->npts)/2]);
	        hold3[j] = (double)(wx1d[j]->sporder);
	    }
	    k1 = 0;
	    k2 = 0;
	    for (i = 0; i < 3; i++) {
	        if (scf.kernw[i] > 0.0) {
	            linepos[i] = Interpolate (scf.kernw[i], k1, hold1, hold2,
                                              tabptr.nrows, &k1);
	            linepos[i] = NINT(linepos[i]);
	            if (linepos[i] < 0.0)
	                linepos[i] = 0.0;
	            if (linepos[i] >= (double)in.sci.data.ny)
	                linepos[i] =  (double)(in.sci.data.ny - 1);
	            mpsfpos[i] = Interpolate (linepos[i], k2, hold2, hold3,
                                              tabptr.nrows, &k2);
	        } else {
	            linepos[i] = 0.0;
	            mpsfpos[i] = 0.0;
	        }
	    }
	    free (hold1);

	    /* Find effective spectral order number for each image line. */

	    mline = (double *) malloc (win.sci.data.ny * sizeof (double));
	    k1 = 0;
	    for (i = 0; i < win.sci.data.ny; i++)
	        mline[i] = Interpolate ((double)i, k1, hold2, hold3,
                                        tabptr.nrows, &k1);
	    free (hold3);
	    free (hold2);

	    /* Compute extended wavelength array. */

	    if (verbose) {
	        printf ("Compute extended arrays.\n");
	        fflush (stdout);
	    }
	    if ((wbig = Alloc2DArrayF (NSBIG, tabptr.nrows)) == NULL)
	        return (OUT_OF_MEMORY);

	    for (j = 0; j < tabptr.nrows; j++) {
	        /* Limits of existing wavelength range */
	        i1 = (NSBIG - wx1d[j]->npts) / 2;
	        i2 = i1 + wx1d[j]->npts - 1;
	        /* Average dispersion in existing wavelength range. */
	        disp = wx1d[j]->wave[wx1d[j]->npts-1] - wx1d[j]->wave[0];
	        disp /= (double)wx1d[j]->npts;
	        /* Extrapolate lower end. */
	        for (i = 0; i < i1; i++)
	            wbig[j][i] = (i - i1) * disp + wx1d[j]->wave[0];
	        /* Copy existing range into array. */
	        for (i = i1, k = 0; i <= i2; i++, k++)
	            wbig[j][i] = wx1d[j]->wave[k];
	        /* Extrapolate high end. */
	        for (i = i2 + 1; i < NSBIG; i++)
	            wbig[j][i] = (i - i2) * disp +
                                 wx1d[j]->wave[wx1d[j]->npts-1];
	    }

	    /* Compute extended ripple array. Note arbitrary extrapolation.
               Tests were performed with several extraoplation schemes.
            */
	    if ((blaze = Alloc2DArrayF (NSBIG, tabptr.nrows)) == NULL)
	        return (OUT_OF_MEMORY);
	    for (j = 0; j < tabptr.nrows; j++) {
	        if ((mrip = GetMatchingOrder (wx1d[j]->sporder, &scf)) < 0) {
	            printf ("No matching spectral order in ripple table.\n");
	            return (ERROR_RETURN);
	        }
	        if ((hold1 = (double *) malloc (scf.rpfunc[mrip].nelem *
                                               sizeof (double))) == NULL)
	            return (OUT_OF_MEMORY);
	        for (i = 0; i < scf.rpfunc[mrip].nelem; i++)
	            hold1[i] = scf.rpfunc[mrip].wavelengths[i] *
                               ((double)scf.rpfunc[mrip].sporder /
                                (double)wx1d[j]->sporder);
	        /* Limits of existing wavelength range */
	        i1 = (NSBIG - wx1d[j]->npts) / 2;
	        i2 = i1 + wx1d[j]->npts - 1;
	        /* Interpolate existing range into array. */
	        k1 = 0;
	        for (i = i1; i <= i2; i++)
	            blaze[j][i] = Interpolate ((double)wbig[j][i], k1, hold1,
                                               scf.rpfunc[mrip].values,
                                               scf.rpfunc[mrip].nelem, &k1);
	        /* Extrapolate lower end. */

	        disp = (scf.rpfunc[mrip].values[1] -
                        scf.rpfunc[mrip].values[0]) /
                       (hold1[1] - hold1[0]);
	        for (i = 0; i < i1; i++)
	            blaze[j][i] = (wbig[j][i] - wbig[j][i1]) * disp +
                                   blaze[j][i1];

	        /* Extrapolate higher end. */

	        disp = (scf.rpfunc[mrip].values[scf.rpfunc[mrip].nelem-1]
                  -     scf.rpfunc[mrip].values[scf.rpfunc[mrip].nelem-2]) /
                    (hold1[scf.rpfunc[mrip].nelem-1] -
                     hold1[scf.rpfunc[mrip].nelem-2]);
	        for (i = i2 + 1; i < NSBIG; i++)
	            blaze[j][i] = (wbig[j][i] - wbig[j][i2]) * disp +
                                   blaze[j][i2];
	        free (hold1);
	    }

	    /* Compute median of number of rows per Angstrom (cross dispersion).
            */
	    if ((dydw = MedianRowsPerAngstrom (wx1d, tabptr.nrows)) == 0.0)
	        return (OUT_OF_MEMORY);

	    /* Create model spectrum. */

	    if (verbose) {
	        printf ("Create model 1-D spectrum.\n");
	        fflush (stdout);
	    }
	    /* Multiply by scale and divide by on-image blaze. */
	    for (j = 0; j < tabptr.nrows; j++) {
	        for (i = 0; i < wx1d[j]->npts;
                     wx1d[j]->net[i++] *= wx1d[j]->scale);
	        for (i = 0, i1 = (NSBIG - wx1d[j]->npts) / 2;
                     i < wx1d[j]->npts; i++, i1++) {
	            if (blaze[j][i1] != 0.0)
	                wx1d[j]->net[i] /= blaze[j][i1];
	        }
	    }

	    /* Splice orders together. */
	    if ((status = Splice (wx1d, tabptr.nrows, exptime, &merge)))
	        return (status);

	    /* Apply threshold to counts or c/sec ? */

	    for (i = 0; i < merge.npts; i++)
	        merge.fmerge[i] = (merge.fmerge[i] > (-20. * exptime)) ?
                                   merge.fmerge[i] : (-20. * exptime);

	    /* Now multiply back by blaze and divide by scale so the
               net data extracted from the input image gets restored to
               its original state.
             */
	    for (j = 0; j < tabptr.nrows; j++) {
	        for (i = 0, i1 = (NSBIG - wx1d[j]->npts) / 2;
                     i < wx1d[j]->npts; i++, i1++) {
	            if (blaze[j][i1] != 0.0)
	                wx1d[j]->net[i] *= blaze[j][i1];
	        }
	        for (i = 0; i < wx1d[j]->npts;
                     wx1d[j]->net[i++] /= wx1d[j]->scale);
	    }

	    /* Main loop. */

	    if (verbose) {
	        printf ("Begin main loop.\n");
	        fflush (stdout);
	    }

	    /* These are used after the main loop ends. */

	    o_mod  = Alloc2DArrayF (win.sci.data.nx, win.sci.data.ny);
	    im_mod = Alloc2DArrayF (win.sci.data.nx, win.sci.data.ny);
	    if (o_mod == NULL || im_mod == NULL)
	        return (OUT_OF_MEMORY);

	    for (iter = 1; iter <= NITER; iter++) {

	        if (verbose) {
	            printf ("\nBegin iteration  %d\n", iter);
	            fflush (stdout);
	        }

	        /* Allocate and zero arrays for scattering model. */

	        im_mod1 = Alloc2DArrayF (win.sci.data.nx, win.sci.data.ny);
	        im_mod2 = Alloc2DArrayF (win.sci.data.nx, win.sci.data.ny);
	        if (im_mod1 == NULL || im_mod2 == NULL)
	            return (OUT_OF_MEMORY);
	        for (j = 0; j < win.sci.data.ny; j++)
	            for (i = 0; i < win.sci.data.nx; o_mod[j][i++] = 0.0);

	        /* Create model spectrum on extended wavelength scale by
                   mapping spliced spectrum onto wbig array. */

	        if (verbose) {
	            printf ("Map 1-D spectrum onto 2-D arrays.\n");
	            fflush (stdout);
	        }

	        if ((fbig = Alloc2DArrayF (NSBIG, tabptr.nrows)) == NULL)
	            return (OUT_OF_MEMORY);
	        for (j = 0; j < tabptr.nrows; j++) {
	            k1 = 0;
	            for (i = 0; i < NSBIG; i++)
	                fbig[j][i] = Interpolate ((double)wbig[j][i], k1,
                                                   merge.wmerge,
                                                   merge.fmerge,
                                                   merge.npts, &k1) *
                                     blaze[j][i];
	        }

	        /* Handle end sections. These must be extrapolated by
                   copying the last valid value, which is not the way
                   the Interpolate function works...
                 */

	        i1 = (NSBIG - wx1d[0]->npts) / 2;
	        i2 = i1 + wx1d[0]->npts - 1;
	        for (i = i2; i < NSBIG; fbig[0][i++] = fbig[0][i2-10]);
	        i1 = (NSBIG - wx1d[tabptr.nrows-1]->npts) / 2;
	        for (i = 0; i < i1; fbig[tabptr.nrows-1][i++] =
                                    fbig[tabptr.nrows-1][i1+10]);

	        /* Build image containing scattered light prediction. */

	        if (verbose) {
	            printf ("Build image with scattered light prediction.\n");
	            fflush (stdout);
	        }

	        /* Loop through orders. */

	        for (iorder = 0; iorder < tabptr.nrows; iorder ++) {

	            if ((k = GetScatterOrder (wx1d[iorder]->sporder, &scf)) ==
                        -1) {
	                printf (
                        "Error: No matching order (# %d) in ECHSCTAB table.\n",
                        wx1d[iorder]->sporder);
	                return (ERROR_RETURN);
	            }

	            scale_lsf = scf.scfunc[k].values;
	            lsf_off = scf.scfunc[k].nelem / 2;

	            echelle_scat_offset = (double *) malloc (
                                                     scf.scfunc[k].nelem *
                                                     sizeof (double));
	            hold1 = (double *) malloc (5 * sizeof (double));
	            hold2 = (double *) malloc (5 * sizeof (double));
	            if (echelle_scat_offset == NULL || hold1 == NULL ||
                        hold2 == NULL)
	                return (OUT_OF_MEMORY);
	            hold1[0] = 0.0;
	            hold2[0] = 0.0;
	            hold1[1] = 0.25 * scf.scfunc[k].nelem;
	            hold2[1] = yyamp;
	            hold1[2] = 0.5  * scf.scfunc[k].nelem;
	            hold2[2] = 0.0;
	            hold1[3] = 0.75 * scf.scfunc[k].nelem;
	            hold2[3] = -yyamp;
	            hold1[4] = scf.scfunc[k].nelem;
	            hold2[4] = 0.0;
	            k1 = 0;
	            for (i = 0; i < scf.scfunc[k].nelem; i++)
	                echelle_scat_offset[i] = Interpolate ((double)i, k1,
                                                 hold1, hold2, 5, &k1);
	            free (hold2);
	            free (hold1);

	            /* Define two additional echelle scattering functions,
                       one which gets scattered in y-direction and one
                       that doesn't. Why ?
                    */

	            scale_lsf1 = (double *) malloc (scf.scfunc[k].nelem *
                                                    sizeof (double));
	            scale_lsf2 = (double *) malloc (scf.scfunc[k].nelem *
                                                    sizeof (double));
	            if (scale_lsf1 == NULL || scale_lsf2 == NULL)
	                return (OUT_OF_MEMORY);
	            for (i = 0; i < scf.scfunc[k].nelem; i++) {
	                scale_lsf1[i] = scale_lsf[i];
	                if (i >= (lsf_off - 5) && i <= (lsf_off + 5))
	                    scale_lsf1[i] = scale_lsf[lsf_off - 5];
	                scale_lsf2[i] = scale_lsf[i] - scale_lsf1[i];
	            }

	            /* Define pointers into data for current order. */

	            i1 = (NSBIG - wx1d[iorder]->npts) / 2;
	            i2 = i1 + wx1d[iorder]->npts - 1;

	            eblaze = blaze[iorder];
	            eonimage = eblaze + i1;   /* part actually on image */
	            worder = wbig[iorder];
	            forder = fbig[iorder];
	            yorder_onimage = wx1d[iorder]->extrlocy;
	            worder_onimage = wx1d[iorder]->wave;

	            /* Define indexes for extended arrays where scattering
                       by lsf can affect actual image region i1 to i2.
                    */
	            istart = i1 - lsf_off;
	            istop  = i2 + lsf_off;

	            /* Loop through pixels istart to istop in ext. arrays. */

	            for (ipix = istart; ipix <= istop; ipix++) {

	                /* Define indices used for addressing arrays
                           while processing current pixel:

                           image1pos and image2pos bound pixels in observed
                           image that can receive flux from current pixel by
                           LSF scattering.

                           ilsf1 and ilsf2 bound points in lsf arrays to use
                           in calculation.
                        */
	                image1pos = ipix - lsf_off - i1;
	                image2pos = ipix + lsf_off - i1;

	                ilsf1 = 0;
	                ilsf2 = scf.scfunc[k].nelem - 1;

	                if (image1pos < 0) {
	                    ilsf1 -= image1pos;
	                    image1pos = 0;
	                 }
	                 if (image2pos > (win.sci.data.nx - 1)) {
	                     ilsf2 = (win.sci.data.nx - 1) + ilsf2 - image2pos;
	                     image2pos = win.sci.data.nx - 1;
	                 }
	                 if (image2pos < image1pos) {
	                     printf ("Invalid indices.\n");
	                     return (ERROR_RETURN);
	                 }

	                 /* Loop thru pixels to receive scattered light. */

	                for (j = 0; j <= (image2pos - image1pos); j++) {

	                    /* Translate fbig in ipix to flux elsewhere
                               in echelle blaze function.
                            */
	                    if (eblaze[ipix] != 0.0)
	                        f = forder[ipix] * eonimage[image1pos+j] /
                                    eblaze[ipix];

	                    /* Determine where to place scattered light in
                               the model images.
                            */
	                    dw = worder_onimage[image1pos+j] - worder[ipix];
	                    y  = yorder_onimage[image1pos+j] - dydw * dw +
                                 echelle_scat_offset[ilsf1+j];

	                    /* Add scattered light into model images. Add into
                               each pixel along extrlocy positions the observed
                               flux in pixel ipix, modified for blaze function
                               changes, line spread function, and fractional
                               splitting between two bracketting rows.
                            */
	                    ixx = image1pos + j;
	                    y1 = (int)y;
	                    iy1 = y1;
	                    iy2 = y1 + 1;
	                    fr2 = y - y1;
	                    fr1 = 1.0 - fr2;
	                    ilsf = ilsf1 + j;
	                    if (ixx >= 0 && ixx < win.sci.data.nx &&
                                iy1 >= 0 && iy1 < win.sci.data.ny &&
                                iy2 >= 0 && iy2 < win.sci.data.ny &&
                                ilsf < ilsf2) {
	                        im_mod1[iy1][ixx] += f * scale_lsf1[ilsf] * fr1;
	                        im_mod1[iy2][ixx] += f * scale_lsf1[ilsf] * fr2;
	                        im_mod2[iy1][ixx] += f * scale_lsf2[ilsf] * fr1;
	                        im_mod2[iy2][ixx] += f * scale_lsf2[ilsf] * fr2;

	                        /* If final iteration of main loop, then
                                   build object image.
                                */
	                        if (iter == NITER) {
	                            o_mod[iy1][ixx] += f * scale_lsf2[ilsf] *
                                                           fr1;
	                            o_mod[iy2][ixx] += f * scale_lsf2[ilsf] *
                                                           fr2;
	                        }
	                    }
	                }
	            }

	            free (scale_lsf1);
	            free (scale_lsf2);
	            free (echelle_scat_offset);
	        }

	        /* Convolve scattering component im_mod1 with
                   cross-disperser profile.
                */
	        if (verbose) {
	            printf ("Convolve with x-disperser profile.\n");
	            fflush (stdout);
	        }
	        npix = win.sci.data.ny + scf.nspsf;
	        i1 = scf.nspsf / 2;
	        hold1 = (double *) calloc (npix, sizeof (double));
	        hold2 = (double *) calloc (npix, sizeof (double));
                if (hold1 == NULL || hold2 == NULL)
	            return (OUT_OF_MEMORY);
	        for (i = 0; i < win.sci.data.nx; i++) {

	            /* Get current colum into zero-padded buffer. */

	            for (j = 0, k = i1; j < win.sci.data.ny; j++, k++)
	                hold1[k] = im_mod1[j][i];

	            /* Convolve. Kernel file is norm. to unit area. */

	            for (j = 0; j < npix; j++) {
		        k1 = j - scf.nspsf / 2;
	                for (k = 0; k < scf.nspsf; k++, k1++) {
	                    k1 = (k1 < 0) ? 0 : k1;
	                    k1 = (k1 >= npix) ? npix - 1 : k1;
	                    hold2[j] += hold1[k1] * scf.spsf[k];
	                }
	            }

	            /* Copy column back into image array. */

	            for (j = 0, k = i1; j < win.sci.data.ny; j++, k++)
	                im_mod1[j][i] = hold2[k];

	            /* Clear buffers for next column processing. */

	            for (k = 0; k < npix; k++) {
	                hold1[k] = 0.0;
	                hold2[k] = 0.0;
	            }
	        }
	        free (hold2);
	        free (hold1);

	        /* Sum model components to construct composite model. */

	        if (verbose) {
	            printf ("Build composite model.\n");
	            fflush (stdout);
	        }
	        for (j = 0; j < win.sci.data.ny; j++) {
	            for (i = 0; i < win.sci.data.nx; i++)
	                im_mod[j][i] = im_mod1[j][i] + im_mod2[j][i];
	        }

	        /* Convolve model image with telescope PSF and
                   detector halo profile.
                */

	        if (verbose) {
	            printf ("Convolve model with PSF and detector halo.\n");
	            fflush (stdout);
	        }

	        nx = win.sci.data.nx;
	        ny = win.sci.data.ny;
	        if (AllocCmplxArray (&zhold, 2 * nx, 2 * ny))
	            return (OUT_OF_MEMORY);

	        Float2Cmplx (im_mod, nx, ny, &zhold);

	        if ((status = FFTConvolve2 (&zhold, &(scf.ft1))))
	            return (status);

	        Cmplx2Float (&zhold, im_mod1, nx, ny);

	        if (scf.nwave > 1) {
	            if (verbose) {
	                printf ("Convolve at 2nd wavelength.\n");
	                fflush (stdout);
	            }
	            Float2Cmplx (im_mod, nx, ny, &zhold);
	            if ((status = FFTConvolve2 (&zhold, &(scf.ft2))))
	                return (status);
	            Cmplx2Float (&zhold, im_mod2, nx, ny);
	        }
	        if (scf.nwave > 2) {
	            if (verbose) {
	                printf ("Convolve at 3rd wavelength.\n");
	                fflush (stdout);
	            }
	            Float2Cmplx (im_mod, nx, ny, &zhold);
	            if ((status = FFTConvolve2 (&zhold, &(scf.ft3))))
	                return (status);
	            im_mod3 = Alloc2DArrayF (win.sci.data.nx, win.sci.data.ny);
	            if (im_mod3 == NULL)
	                return (OUT_OF_MEMORY);
	            Cmplx2Float (&zhold, im_mod3, nx, ny);
	        }

	        /* If final iteration, convolve object model o_mod
                   with PSF and halo.
                */

	        if (iter == NITER) {

	            if (verbose) {
	                printf ("Convolve final object image.\n");
	                fflush (stdout);
	            }

	            o_mod1 = Alloc2DArrayF (win.sci.data.nx, win.sci.data.ny);
	            if (o_mod1 == NULL)
	                return (OUT_OF_MEMORY);
	            Float2Cmplx (o_mod, nx, ny, &zhold);
	            if ((status = FFTConvolve2 (&zhold, &(scf.fto1))))
	                return (status);
	            Cmplx2Float (&zhold, o_mod1, nx, ny);

	            if (scf.nwave > 1) {
	                o_mod2 = Alloc2DArrayF (win.sci.data.nx,
                                                win.sci.data.ny);
	                if (o_mod2 == NULL)
	                    return (OUT_OF_MEMORY);
	                if (verbose) {
	                    printf ("Convolve at 2nd wavelength.\n");
	                    fflush (stdout);
	                }
	                Float2Cmplx (o_mod, nx, ny, &zhold);
	                if ((status = FFTConvolve2 (&zhold, &(scf.fto2))))
	                    return (status);
	                Cmplx2Float (&zhold, o_mod2, nx, ny);
	            }
	            if (scf.nwave > 2) {
	                o_mod3 = Alloc2DArrayF (win.sci.data.nx,
                                                win.sci.data.ny);
	                if (o_mod3 == NULL)
	                    return (OUT_OF_MEMORY);
	                if (verbose) {
	                    printf ("Convolve at 3rd wavelength.\n");
	                    fflush (stdout);
	                }
	                Float2Cmplx (o_mod, nx, ny, &zhold);
	                if ((status = FFTConvolve2 (&zhold, &(scf.fto3))))
	                    return (status);
	                Cmplx2Float (&zhold, o_mod3, nx, ny);
	            }
	        }

	        FreeCmplxArray (&zhold);

	        /* Coadd convolved model images with row weights based
                   on order number. This scheme is designed to handle
                   variations in PSF shape with row (wavelength). Up to
                   3 model images are constructed using constant PSF kernels.
                */

	        if (verbose) {
	            printf ("Co-add individual model images.\n");
	            fflush (stdout);
	        }
	        for (j = 0; j < win.sci.data.ny; j++) {
	            for (i = 0; i < win.sci.data.nx; i++)
	                im_mod[j][i] = im_mod1[j][i];
	        }

	        if (scf.nwave > 1) {
	            FillArray (im_mod, im_mod1, im_mod2, win.sci.data.nx,
                               win.sci.data.ny,
                               linepos[0], linepos[1], mpsfpos[0], mpsfpos[1],
                               mline, in.sci.data.ny);
	            for (j = (int)linepos[1]; j < win.sci.data.ny; j++) {
	                for (i = 0; i < win.sci.data.nx; i++)
	                    im_mod[j][i] = im_mod2[j][i];
	            }
	        }
	        if (scf.nwave > 2) {
	            FillArray (im_mod, im_mod2, im_mod3, win.sci.data.nx,
                               win.sci.data.ny,
                               linepos[1], linepos[2], mpsfpos[1], mpsfpos[2],
                               mline, in.sci.data.ny);
	            for (j = (int)linepos[2]; j < win.sci.data.ny; j++) {
	                for (i = 0; i < win.sci.data.nx; i++)
	                    im_mod[j][i] = im_mod3[j][i];
	            }
	        }

	        /* If final iteration, coadd object models with weights. */

	        if (iter == NITER) {

	            if (verbose) {
	                printf ("Co-add individual object images.\n");
	                fflush (stdout);
	            }
	            for (j = 0; j < win.sci.data.ny; j++) {
	                for (i = 0; i < win.sci.data.nx; i++)
	                    o_mod[j][i] = o_mod1[j][i];
	            }
	            if (scf.nwave > 1) {
	                FillArray (o_mod, o_mod1, o_mod2, win.sci.data.nx,
                                   win.sci.data.ny,
                                   linepos[0], linepos[1], mpsfpos[0],
                                   mpsfpos[1], mline, in.sci.data.ny);
	                for (j = (int)linepos[1]; j < win.sci.data.ny; j++) {
	                    for (i = 0; i < win.sci.data.nx; i++)
	                        o_mod[j][i] = o_mod2[j][i];
	                }
	            }

	            if (scf.nwave > 2) {
	                FillArray (o_mod, o_mod2, o_mod3, win.sci.data.nx,
                                   win.sci.data.ny,
                                   linepos[1], linepos[2], mpsfpos[1],
                                   mpsfpos[2], mline, in.sci.data.ny);
	                for (j = (int)linepos[2]; j < win.sci.data.ny; j++) {
	                    for (i = 0; i < win.sci.data.nx; i++)
	                        o_mod[j][i] = o_mod3[j][i];
	                }
	            }
	        }

	        /* Add ghosts due to reflections from window above NUV MAMA. */

	        if ((status = AddGhost (phdr, im_mod, win.sci.data.nx,
                                        win.sci.data.ny)))
	            return (status);

	        /* If final iteration, clean up memory and leave main loop. */

	        if (iter == NITER) {
	            free (merge.fmerge);
	            free (merge.wmerge);
	            Free2DArrayF (o_mod1, win.sci.data.ny);
	            Free2DArrayF (im_mod1, win.sci.data.ny);
	            Free2DArrayF (im_mod2, win.sci.data.ny);
	            if (scf.nwave > 1)
	                Free2DArrayF (o_mod2,  win.sci.data.ny);
	            if (scf.nwave > 2) {
	                Free2DArrayF (o_mod3,  win.sci.data.ny);
	                Free2DArrayF (im_mod3, win.sci.data.ny);
	            }
	            Free2DArrayF (fbig, tabptr.nrows);
	            break;
	        }

	        /* Write temporary image to disk so it can be accessed by
                   the standard 1-D extraction code.
                */
	        if (verbose) {
	            printf ("Write current model image to disk.\n");
	            fflush (stdout);
	        }
	        remove (temp_ima);
	        initSingleGroup (&wout);
	        allocSingleGroup (&wout, win.sci.data.nx, win.sci.data.ny, True);
	        for (j = 0; j < wout.sci.data.ny; j++) {
	            for (i = 0; i < wout.sci.data.nx; i++)
	                Pix (wout.sci.data, i, j) = im_mod[j][i];
	        }
	        copyHdr (&(wout.sci.hdr), &(win.sci.hdr));

                tim = openOutputImage (temp_ima, "", 0, phdr, 0, 0, FITSBYTE);
                if (hstio_err())
                    return (OPEN_FAILED);
                closeImage (tim);
                if (putSingleGroup (temp_ima, 1, &wout, 0))
                     return (OPEN_FAILED);
	        freeSingleGroup (&wout);

	        /* Create _x1d file from temporary image. Make sure
                   heliocentric and photo corrections are turned off.
                   Note that extver is set identically to 1 since this
                   is a temporary file with one IMSET.
                */

	        if (verbose) {
	            printf ("Begin extraction of 1-D data.\n");
	            fflush (stdout);
	        }
	        if ((status = CalStis6Std (temp_ima, temp_x1d, backcorr,
                             dispcorr, OMIT, OMIT, sgeocorr, ctecorr,
                             cl_a2center, maxsearch, extrsize, bk1size,
                             bk2size, bk1offset, bk2offset, bktilt, bkord,
                             sporder, xtracalg, printtime, verbose, 1,
                             ccglobal, ccthresh, do_profile, pstep, wstep,
                             minsn, rejranges, profilefile, fluxfile, outw,
                             backval, backerr, variance, fflux, psclip, sclip,
                             lfilter, 1, pipeline, &dummy, &dummy, 1.0,
                             blazeshift, bks_mode, bks_order, xoffset)))
	            return (status);
	        remove (temp_ima);

	        if (verbose) {
	            printf ("End extraction of 1-D data.\n");
	            fflush (stdout);
	        }

	        /* Open table file just created by CalStis6Std */

	        if (verbose) {
	            printf ("Read 1-D data.\n");
	            fflush (stdout);
	        }
	        if ((status = OpenX1DTable (temp_x1d, 0, &tabptr)))
	            return (status);

	        /* Alloc memory for storing table contents. */

	        if ((wx1d2 = AllocX1DTable (tabptr.nrows)) == NULL)
	            return (OUT_OF_MEMORY);

	        /* Get data from table and store in RowContents array. */

	        if ((status = GetX1DTable (&tabptr, wx1d2))) {
	            FreeX1DTable (wx1d2, tabptr.nrows);
	            c_tbtclo (tabptr.tp);
	            remove (temp_x1d);
	            return (status);
	        }

	        /* Copy original extrlocy positions into new data
                   structure. This is to avoid effects from residual
                   shifts caused by calstis6 recentering abilities.
                */

	        for (j = 0; j < tabptr.nrows; j++) {
	            for (i = 0; i < wx1d2[j]->npts; i++)
                        wx1d2[j]->extrlocy[i] = wx1d[j]->extrlocy[i];
	        }

	        /* Difference between net arrays (now - original data)
                   is used to update the model spectrum. Remember that
                   spectrum estimate is expressed in counts, not c/s !
                */
	        if (verbose) {
	            printf ("Update with correction term.\n");
	            fflush (stdout);
	        }
	        for (j = 0; j < tabptr.nrows; j++) {
	            i1 = (NSBIG - wx1d2[j]->npts) / 2;
	            i2 = i1 + wx1d2[j]->npts - 1;
	            for (i = i1, k = 0; i <= i2; i++, k++)
	                wx1d2[j]->net[k] = fbig[j][i] / exptime +
                            (wx1d[j]->net[k] - wx1d2[j]->net[k]) * inc_scale;
	        }

	        /* Normalize by the blaze so the result can be spliced
                   together.
                */
	        for (j = 0; j < tabptr.nrows; j++) {
	            for (i = 0, i1 = (NSBIG - wx1d2[j]->npts) / 2;
                         i < wx1d2[j]->npts; i++, i1++) {
	                if (blaze[j][i1] != 0.0)
                             wx1d2[j]->net[i] /= blaze[j][i1];
	            }
	        }

	        /* x1d table is not necessary anymore. */

	        c_tbtclo (tabptr.tp);
	        remove (temp_x1d);

	        /* Update model. */

	        if (verbose) {
	            printf ("Build new 1-D spectrum model.\n");
	            fflush (stdout);
	        }
	        /* Must be fred now so they can be reused with arrays
                   of potentially different size.
                */
	        free (merge.fmerge);
	        free (merge.wmerge);

	        /* Splice orders together and transform to counts. */

	        if ((status = Splice (wx1d2, tabptr.nrows, exptime, &merge)))
	            return (status);

	        /* Free memory in current iteration. */

	        FreeX1DTable (wx1d2, tabptr.nrows);
	        Free2DArrayF (im_mod1, win.sci.data.ny);
	        Free2DArrayF (im_mod2, win.sci.data.ny);
	        if (scf.nwave == 3)
	            Free2DArrayF (im_mod3, win.sci.data.ny);
	        Free2DArrayF (fbig, tabptr.nrows);
	    }

	    if (verbose) {
	        printf ("End main loop.\n");
	        fflush (stdout);
	    }

	    /* Free some more memoery. */

	    Free2DArrayF (blaze, tabptr.nrows);
	    Free2DArrayF (wbig, tabptr.nrows);
	    free (mline);
	    Free2DArrayF (bpos, tabptr.nrows+1);

	    /* Construct scattered light image. */

	    for (j = 0; j < win.sci.data.ny; j++) {
	        for (i = 0; i < win.sci.data.nx; i++)
	            Pix (win.sci.data, i, j) = im_mod[j][i] - o_mod[j][i];
	    }

	    Free2DArrayF (o_mod,  win.sci.data.ny);
	    Free2DArrayF (im_mod, win.sci.data.ny);

	    /* Rebin image back to hi-res if necessary. */

	    if (verbose) {
	        printf ("Rebin/copy work image to original sampling.\n");
	        fflush (stdout);
	    }
	    initSingleGroup (&wout);
	    if ((wx1d2 = AllocX1DTable (tabptr.nrows)) == NULL)
	        return (OUT_OF_MEMORY);
	    allocSingleGroup (&wout, in.sci.data.nx, in.sci.data.ny, True);
	    if (ltm[0] > 1.0) {
	        if (RebinData (&win, &wout, wx1d, wx1d2, -2, tabptr.nrows))
	            return (ERROR_RETURN);
	    } else {
	        if (RebinData (&win, &wout, wx1d, wx1d2, 1, tabptr.nrows))
	            return (ERROR_RETURN);
	    }
	    /* Copy ERR and DQ arrays to output, as RebinData doesn't do that for
	       binned data */
	    if (ltm[0] > 1.0) {
	        for (j = 0; j < in.sci.data.ny; j++) {
		    for (i = 0; i < in.sci.data.nx; i++) {
		        Pix(wout.err.data, i, j) = Pix(in.err.data, i, j);
		        DQPix(wout.dq.data, i, j) = DQPix(in.dq.data, i, j);
		    }
		}
	    }
	    copyHdr (&(wout.sci.hdr), &(in.sci.hdr));
	    FreeX1DTable (wx1d2, tabptr.nrows);

	    /* Subtract scattered light image from original raw image. */

	    if (verbose) {
	        printf ("Subtract scattered light from input image.\n");
	        fflush (stdout);
	    }
	    for (j = 0; j < in.sci.data.ny; j++) {
	        for (i = 0; i < in.sci.data.nx; i++)
	            Pix (wout.sci.data, i, j) = Pix (in.sci.data, i, j) -
                                                Pix (wout.sci.data, i, j);
	    }

	    /* Write image data to disk so it can be accessed by
               the standard 1-D extraction code.
            */
	    if (verbose) {
	        printf ("Write corrected input image to disk.\n");
	        fflush (stdout);
	    }
	    if (first_extver) {
	        remove (temp_ima1);
                tim = openOutputImage (temp_ima1, "", 0, phdr, 0, 0, FITSBYTE);
                if (hstio_err())
                    return (OPEN_FAILED);
                closeImage (tim);
		first_extver = 0;		/* false */
	    }
            if (putSingleGroup (temp_ima1, extver, &wout, 0))
                return (OPEN_FAILED);
	    freeSingleGroup (&wout);

	    if (verbose && (nimages > 1)) {
	        printf ("End processing IMSET %d\n\n", extver);
	        fflush (stdout);
	    }

	    FreeX1DTable (wx1d, tabptr.nrows);
	    freeSingleGroup (&win);
	    freeSingleGroup (&in);
	}

	/* Free convolution kernel memory. */

	FreeCmplxArray (&(scf.ft1));
	FreeCmplxArray (&(scf.ft2));
	FreeCmplxArray (&(scf.ft3));
	FreeCmplxArray (&(scf.fto1));
	FreeCmplxArray (&(scf.fto2));
	FreeCmplxArray (&(scf.fto3));
	FreeScatter (&scf);
	FreeThroughput6 (&slit);

	/* Create final _x1d file. CalStis6Std must be run over
           all IMSETs in the file ! Don't forget to set helcorr and
           fluxcorr properly.
        */
	if (verbose) {
	    printf ("Begin final extraction of 1-D data.\n");
	    fflush (stdout);
	}

	/* Reset the verbosity level. */

	if ((status = CalStis6Std (temp_ima1, output, backcorr, dispcorr,
                     fluxcorr, helcorr, sgeocorr, ctecorr, cl_a2center,
                     maxsearch, extrsize, bk1size, bk2size, bk1offset,
                     bk2offset, bktilt, bkord, sporder, xtracalg, printtime,
                     verbose == 3? 0: verbose+1, 1, ccglobal, ccthresh,
                     do_profile, pstep, wstep, minsn, rejranges, profilefile,
                     fluxfile, outw, backval, backerr, variance, fflux,
                     psclip, sclip, lfilter, 0, pipeline, &dummy, &dummy, 1.0,
                     blazeshift, bks_mode, bks_order, xoffset)))
	    return (status);

	if (idtfile[0] != '\0') {
	    sprintf (out_ima, "%s.fits", idtfile);
	    rename (temp_ima1, out_ima);
	    if (verbose) {
	        printf ("File %s contains deconvolved image.\n", out_ima);
	        fflush (stdout);
	    }
	} else
	    remove (temp_ima1);

	/* Output table must have some of its columns re-computed.
           Notice that the very original 1-D extraction must be
           supplied, not one of the possibly rebinned later versions.
         */
	if (verbose) {
	    printf ("Update columns in final output table.\n");
	    fflush (stdout);
	}
	for (extver = 1; extver <= nimages; extver++) {
	    extver0 = extver - 1;

	    /* This will still be NULL if the imset was skipped due to
		exptime being zero.
	    */
	    if (gx1d[extver0] == NULL)
		continue;

	    sprintf (temp_out, "%s[SCI,%d]", output, extver);

	    if ((status = RedoX1DFile (temp_out,
                                       gx1d[extver0], gx1d_nrows[extver0]))) {
	        return (status);
	    }

	    FreeX1DTable (gx1d[extver0], gx1d_nrows[extver0]);
	}
	free (gx1d);
	free (gx1d_nrows);

	/* Now the primary header of the final output table must be
           updated (SC2DCORR switch and HISTORY keywords).
        */
	if (verbose) {
	    printf ("Update header in final output table.\n");
	    fflush (stdout);
	}
	initHdr (&tphdr);
	tdesc = openUpdateImage (output, "", 0, &tphdr);
	if (hstio_err())
	    return (OPEN_FAILED);
	Put_KeyS (&tphdr, "sc2dcorr", "COMPLETE", "");
	if ((status = History6 (&sts_idt, &tphdr, 1)))
	    return (status);
	putHeader (tdesc);
	if (hstio_err())
	    return (OPEN_FAILED);
	closeImage (tdesc);
	freeHdr (&tphdr);

        if (verbose == 1 || verbose == 2 || verbose == 3) {
	    PrEnd (6);
            if (prtimestamp)
	        TimeStamp ("CALSTIS-6 completed", sts_idt.rootname);
	}

	return (STIS_OK);
}



/*
   Gets scattering function index given the spectral order.
*/

static int GetScatterOrder (int sporder, ScatterFunctions *scf) {

/* arguments:
int           sporder;	i: sporder number to match
ScatterFunctions *scf;  i: scattering functions

return: the index of the matching entry, or -1 if error.
*/
	int i, index;

	index = -1;
	for (i = 0; i < scf->nsc; i++) {
	    if (sporder == scf->scfunc[i].sporder) {
	        index = i;
	    }
	}
	return (index);
}



/*
   Gets ripple function index for given spectral order. This function
   returns the index of the closest matching spectral order entry in the
   ripple function array.
*/

static int GetMatchingOrder (int sporder, ScatterFunctions *scf) {

/* arguments:
int           sporder;	i: sporder number to match
ScatterFunctions *scf;  i: scattering functions

return: the index of the matching entry, or -1 if error.
*/
	int i, index, dist, mindist;

	mindist = 100000;
	index = -1;
	for (i = 0; i < scf->nrp; i++) {
	    dist = abs (sporder - scf->rpfunc[i].sporder);
	    if (dist < mindist) {
	        mindist = dist;
	        index = i;
	    }
	}
	return (index);
}



/*
   Linear interpolation.

   This function is optimized to work within loops that provide
   independent variable values that differ by small amounts from
   iteration to iteration. It starts searching in the data arrays
   from a supplied index value, and outputs the index of the data
   point immediately before the sought data value. This index can
   then be used in the next call to start the search at a point
   already close to the sought one.

   The ordering (ascending or descending) of the independent variable
   array is not relevant.

   For data values outside the range the extrapolation based on the
   two first/last data points is returned. If the beginning index is
   out of bounds, the first/last value of the dependent variable is
   returned.
*/

static double Interpolate (double x1, int begin, double *x, double *y,
                           int n, int *index) {

/* arguments:
double  x1;	i: the independent variable value where to interpolate
int  begin;	i: begin search from this index
double  *x;	i: array of independent values
double  *y;	i: array of dependent values
int      n;	i: array size
int *index;	o: index of the bin immediately before the sought one

return: the interpolated/extrapolated value.
*/
	int i, i1, i2;

	/* Check if beginning index is out of range. */

	if (begin < 0)  {
	    *index = 0;
	    return (y[0]);
	}
	if (begin >= n)  {
	    *index = n-2;
	    return (y[n-1]);
	}

	/* Scan array, beginning from supplied index and looking incrementaly
           in both directions. This takes care of ascending and descending
           arrays, as well as independent variable values provided in any
           sequence.
        */

	for (i = 0; i < n; i++) {

	    i1 = begin + i;
	    i2 = begin + i + 1;
	    if (i2 < n) {
	        if ((x1 >= x[i1] && x1 <= x[i2]) ||
                    (x1 <= x[i1] && x1 >= x[i2])) {
	            *index = i1;
	            return (Linear (x[i1], y[i1] , x[i2], y[i2], x1));
	        }
	    }

	    i1 = begin - i;
	    i2 = begin - i + 1;
	    if (i1 >= 0) {
	        if ((x1 >= x[i1] && x1 <= x[i2]) ||
                    (x1 <= x[i1] && x1 >= x[i2])) {
	            *index = i1;
	            return (Linear (x[i1], y[i1] , x[i2], y[i2], x1));
	        }
	    }
	}

	/* Take care of out-of-range condition. */

	if (x[1] > x[0]) {
	    if (x1 < x[0]) {
	        *index = 0;
	        return (Linear (x[0], y[0] , x[1], y[1], x1));
	    } else {
	        *index = n-2;
	        return (Linear (x[n-1], y[n-1] , x[n-2], y[n-2], x1));
	    }
	} else {
	    if (x1 > x[0]) {
	        *index = 0;
	        return (Linear (x[0], y[0] , x[1], y[1], x1));
	    } else {
	        *index = n-2;
	        return (Linear (x[n-1], y[n-1] , x[n-2], y[n-2], x1));
	    }
	}
}



/*
    Linear interpolation/extrapolation between two points.
*/

static double Linear (double x1, double y1 , double x2, double y2, double x) {

	double a, b;

	if (x2 != x1) {
	    a = (y2 - y1) / (x2 - x1);
	    b = y2 - a * x2;
	    return (a * x + b);
	} else
	    return (0.0);
}



/*
   Compute median of number of rows per Angstrom (cross dispersion).
*/

static double MedianRowsPerAngstrom (RowContents **x1d, int n) {

/* arguments:
RowContents **x1d;	i: extracted 1-D data array
int n;			i: number of RowContents elements in array

return: the median.
*/
	int j;
	double *temp, median = 0.0;

	temp = (double *) malloc ((n - 1) * sizeof (double));
	if (temp == NULL)
	    return (0.0);

	for (j = 1; j < n; j++)
	    temp[j-1] = (x1d[j]->extrlocy[(x1d[j]->npts)/2] -
                         x1d[j-1]->extrlocy[(x1d[j-1]->npts)/2]) /
                        (x1d[j]->wave[(x1d[j]->npts)/2] -
                         x1d[j-1]->wave[(x1d[j-1]->npts)/2]);

	median = Select ((long)(n)/2, (long)(n-1), temp-1);
	free (temp);

	return (median);
}


#define D_SWAP(a,b) { double temp=(a);(a)=(b);(b)=temp; }

/* Algorithm in the public domain (http://blog.beamng.com/a-faster-selection-algorithm/)
 * with modifications to maintain the HSTCAL interface and added some error checking.
 * Note: The code needs more than 2 elements to work.
 *
 * This routine maintains the same interface as the original routine.
 * Arrays are 1-indexed, thus the calling sequence should be
 * something as:
 *  med = select ((long)(n+1)/2, (long)n, array-1);
 *
 * k: the K-th smallest value in array, range = [1:length]
 * length: length of the array
 * array: the array to search
 */
static double Select (unsigned long k, unsigned long length, double *array) {
    unsigned long l=0, m=length-1, i=l, j=m;
    double x;

    /* The original routine did no checking internally, but if its calling
     * routine returns 0.0, it is an "OUT_OF_MEMORY" error.  This is kludgy,
     * but it is the only method to return an error in the current context.
     */
    if ((length < 2) || (k <= 0) || (k > length)) {
        //sprintf (MsgText, "Requested value %d is out of range (1 - %d}", k, length);
        //trlerror (MsgText)
        return (0.0);
    }

    /* Adjust k to be for a zero-indexed array */
    array += 1;
    k--;

    while (l < m) {
        if(array[k] < array[i]) D_SWAP(array[i], array[k]);
        if(array[j] < array[i]) D_SWAP(array[i], array[j]);
        if(array[j] < array[k]) D_SWAP(array[k], array[j]);

        x=array[k];
        while (j > k && i < k) {
            do i++; while (array[i] < x);
            do j--; while (array[j] > x);

            D_SWAP(array[i], array[j]);
        }
        i++; j--;

        if (j < k) {
            while (array[i] < x) i++;
            l = i; j = m;
        }
        if (k < i) {
            while (x < array[j]) j--;
            m = j; i = l;
        }
    }
    return array[k];
}
# undef  D_SWAP


/*
   This function is used to fill up parts of the model array when
   combinng data from different wavelengths (PSFs).
*/

static void FillArray (float **out, float **in1, float **in2, int nx, int ny,
                       double linepos1, double linepos2, double mpsfpos1,
                       double mpsfpos2, double *mline, int msize) {
/* arguments:
float **out;		o: output array
float **in1;		i: input array # 1
float **in2;		i: input array # 2
int nx;			i: # of columns in input/output arrays
int ny;			i: # of rows in input/output arrays
double linepos1;	i: image lines where to update the output array
double linepos2;
double mpsfpos1;	i: order # corrersponding to lines above
double mpsfpos2;
double* mline;		i: fractional order # associated with each image row
int msize;		i: size of mline array
*/
	int i, j;
	float frac1, frac2;

	/* Scan array section. */

	for (j = (int)linepos1; j < (int)linepos2; j++) {

	    if (j >=0 && j < ny) {

	        /* Compute weights for current row. */

	        frac2 = (mline[j] - mpsfpos1) / (mpsfpos2 - mpsfpos1);
	        frac1 = 1.0 - frac2;

	        /* Fill pixels in array row. */

	        for (i = 0; i < nx; i++)
	            out[j][i] = frac1 * in1[j][i] + frac2 * in2[j][i];
	    }
	}
}


/*
   Adds ghost to NUV MAMA image.
*/

static int AddGhost (Hdr *phdr, float **im, int nx, int ny) {

	char opt_elem[STIS_CBUF];	/* grating */
	double kx[2][2], ky[2][2];	/* ghost image warp matrices */
	double xx, yy;
	float **ghost, **cghost;
	int i, j, ii, jj, status;

	float **Alloc2DArrayF (int, int);
	void Free2DArrayF (float **, int);

	if ((status = Get_KeyS (phdr, "OPT_ELEM", 0, "", opt_elem, STIS_CBUF)))
	    return (status);
	if (!(streq_ic (opt_elem, "E230M") || streq_ic (opt_elem, "E230H")))
	    return (STIS_OK);

	if (streq_ic (opt_elem, "E230M")) {
	    kx[0][0] = -33.9909;
	    kx[0][1] = -0.000844246;
	    kx[1][0] =  0.998378;
	    kx[1][1] = 3.51736e-06;
	    ky[0][0] = 5.46450;
	    ky[0][1] = 1.00869;
	    ky[1][0] = 9.48275e-05;
	    ky[1][1] = -1.10486e-06;
	} else if (streq_ic (opt_elem, "E230H")) {
	    kx[0][0] = -16.4834;
	    kx[0][1] = -0.000280014;
	    kx[1][0] =  0.998752;
	    kx[1][1] = 1.51807e-06;
	    ky[0][0] = 5.59971;
	    ky[0][1] = 1.00791;
	    ky[1][0] = 1.58386e-05;
	    ky[1][1] = -7.14570e-07;
	} else
	    return (STIS_OK);

	if ((ghost = Alloc2DArrayF (nx, ny)) == NULL)
	    return (OUT_OF_MEMORY);

	for (j = 0; j < ny; j++) {
	    for (i = 0; i < nx; i++) {
	        xx = kx[0][0] + j * kx[0][1] + i * kx[1][0] + j * i * kx[1][1];
	        yy = ky[0][0] + j * ky[0][1] + i * ky[1][0] + j * i * ky[1][1];
	        ii = (int)NINT(xx);
	        jj = (int)NINT(xx);
	        ii = (ii < 0) ? 0 : ii;
	        jj = (jj < 0) ? 0 : jj;
	        ii = (ii >= nx) ? (nx - 1) : ii;
	        jj = (jj >= ny) ? (ny - 1) : jj;
	        ghost[j][i] = 0.003 * im[jj][ii];
	    }
	}

	if ((cghost = Alloc2DArrayF (nx, ny)) == NULL)
	    return (OUT_OF_MEMORY);

	for (j = 1; j < ny-1; j++) {
	    for (i = 1; i < nx-1; i++) {
	        cghost[j][i] = ghost[j-1][i-1] * 0.2 +
                               ghost[j-1][i]   * 0.5 +
                               ghost[j-1][i+1] * 0.2 +
                               ghost[j][i-1]   * 0.5 +
                               ghost[j][i]           +
                               ghost[j][i+1]   * 0.5 +
                               ghost[j+1][i-1] * 0.2 +
                               ghost[j+1][i]   * 0.5 +
                               ghost[j+1][i+1] * 0.2;
	        cghost[j][i] /= 3.8;
	    }
	}

	Free2DArrayF (ghost, ny);

	for (j = 0; j < ny; j++) {
	    for (i = 0; i < nx; i++)
	        im[j][i] += cghost[j][i];
	}

	Free2DArrayF (cghost, ny);
	return (STIS_OK);
}


/* Not used */
/*
static double BLog (double arg) {

	double hold;

	hold = (arg > 1.E-5) ? arg : 1.E-5;
	return (log (arg));
}
*/


/*
   This function re-computes the GROSS and BACKGROUND columns in the
   final output x1d table. The new values are defined as:

   gross      = non-idt gross (original extraction from calstis6std)
   background = non-idt gross - net
*/

static int RedoX1DFile (char *output, RowContents **x1d, int x1d_nrows) {

	TblDesc tabptr;
	float *new_gross, *new_back, *new_net;
	short npts, sporder;
	int i, irow, jrow, nnet, status;

	int OpenX1DTable (char *, int, TblDesc *);

	tabptr.tp = NULL;
	if ((status = OpenX1DTable (output, 1, &tabptr)))
	    return (status);

	/* Get the required columns. */

	c_tbcfnd1 (tabptr.tp, SPORDER,    &(tabptr.sporder));
	c_tbcfnd1 (tabptr.tp, GROSS,      &(tabptr.gross));
	c_tbcfnd1 (tabptr.tp, BACKGROUND, &(tabptr.back));
	c_tbcfnd1 (tabptr.tp, NET,        &(tabptr.net));
	c_tbcfnd1 (tabptr.tp, NELEM,      &tabptr.npts);
	if ((tabptr.gross   == 0) || (tabptr.back == 0) ||
            (tabptr.npts    == 0) || (tabptr.net  == 0) ||
            (tabptr.sporder == 0)) {
	    c_tbtclo (tabptr.tp);
	    printf ("Column not found in output table.\n");
	    return (TABLE_ERROR);
	}

	/* Process each row. */

	for (irow = 0; irow < tabptr.nrows; irow++) {

	    c_tbegts (tabptr.tp, tabptr.npts, irow+1, &npts);
	    new_gross = (float *) calloc (npts, sizeof(float));
	    new_back  = (float *) calloc (npts, sizeof(float));
	    new_net   = (float *) calloc (npts, sizeof(float));
	    if ((new_gross == NULL) || (new_back == NULL) ||
                (new_net == NULL)) {
	        c_tbtclo (tabptr.tp);
                printf ("Not enough memory to allocate data arrays.\n");
	        return (OUT_OF_MEMORY);
	    }

	    c_tbegts (tabptr.tp, tabptr.sporder, irow+1, &sporder);
	    nnet   = c_tbagtr (tabptr.tp, tabptr.net, irow+1,
                               new_net, 1, npts);
	    if (nnet != npts) {
	        c_tbtclo (tabptr.tp);
	        printf ("Unexpected number of elements read from table.\n");
	        return (TABLE_ERROR);
	    }

	    /* Find the row in x1d with the same spectral order as the
	       current row of the output table.
	    */
	    jrow = -1;
	    for (i = 0;  i < x1d_nrows;  i++) {
		if (x1d[i]->sporder == sporder) {
		    jrow = i;
		    break;
		}
	    }
	    if (jrow == -1) {
		printf ("ERROR    sporder %d not found in gx1d array.\n",
			sporder);
		return (INTERNAL_ERROR);
	    }

	    for (i = 0; i < npts; i++) {
	        new_gross[i] = x1d[jrow]->gross[i];
	        new_back[i]  = x1d[jrow]->gross[i] - new_net[i];
	    }

	    c_tbaptr (tabptr.tp, tabptr.gross, irow+1, new_gross, 1, npts);
	    c_tbaptr (tabptr.tp, tabptr.back,  irow+1, new_back,  1, npts);

	    free (new_gross);
	    free (new_back);
	    free (new_net);
	}

	c_tbtclo (tabptr.tp);
	return (STIS_OK);
}



/* Checks SC2DCORR reference files. */

static int CheckIDT (Hdr *phdr, StisInfo6 *sts, int *missing) {

/* arguments:
Hdr *phdr       i: primary header
StisInfo6 *sts  o: file names
int *missing    io: incremented if the file is missing
*/

	int status;
 	int GetCheckRef (Hdr *, char *, RefTab *, int *, int *, int);

	if ((status = GetCheckRef (phdr, "ECHSCTAB", &sts->echsctab,
                                   &sts->idt, missing, FATAL)))
	    return (status);
	if ((status = GetCheckRef (phdr, "EXSTAB", &sts->exstab,
                                   &sts->idt, missing, FATAL)))
	    return (status);
	if ((status = GetCheckRef (phdr, "CDSTAB", &sts->cdstab,
                                   &sts->idt, missing, FATAL)))
	    return (status);
	if ((status = GetCheckRef (phdr, "RIPTAB", &sts->riptab,
                                   &sts->idt, missing, FATAL)))
	    return (status);
	if ((status = GetCheckRef (phdr, "SRWTAB", &sts->srwtab,
                                   &sts->idt, missing, FATAL)))
	    return (status);
	if ((status = GetCheckRef (phdr, "HALOTAB", &sts->halotab,
                                   &sts->idt, missing, FATAL)))
	    return (status);
	if ((status = GetCheckRef (phdr, "TELTAB", &sts->psftab,
                                   &sts->idt, missing, FATAL)))
	    return (status);

	if (sts->idt != PERFORM) {
	    printf (
            "Warning  SC2DCORR skipped due to dummy reference file.\n");
	    return (NOTHING_TO_DO);
	}
	return (0);
}


static void BuildTempNames (char *basename, char *out1, char *out2, char *out3) {

	char file_path[STIS_FNAME];
	char file_name[STIS_FNAME];
	char separator[2];
	char *slash;
	char *aname;

	aname = strdup (basename);
	slash = strrchr (aname, '/');

	if (slash != NULL) {
	    *slash = '\0';
	    strcpy (file_path, aname);
	    strcpy (file_name, slash+1);
	    strcpy (separator, "/");
	} else {
	    strcpy (file_path, "");
	    strcpy (file_name, basename);
	    strcpy (separator, "");
	}

	strcpy (out1, file_path);
	strcat (out1, separator);
	strcat (out1, "TEMPX1D_");
	strcat (out1, file_name);

	strcpy (out2, file_path);
	strcat (out2, separator);
	strcat (out2, "TEMPIMA_");
	strcat (out2, file_name);

	strcpy (out3, file_path);
	strcat (out3, separator);
	strcat (out3, "TEMPIMA1_");
	strcat (out3, file_name);

	free(aname);
}




/***************************************************************************/




/* These functions dump arrays (float and complex) as IMSET images and
   write at stdout. They were used for debugging during development and
   are kept in here just in case. (Not used)
*/
/*
static int Debug (char *name, float **array, int nx, int ny) {

	SingleGroup out;
	int i, j;

	initSingleGroup (&out);
	if (allocSingleGroup (&out, nx, ny))
	    return (OUT_OF_MEMORY);
	for (j = 0; j < ny; j++) {
	    for (i = 0; i < nx; i++)
	        Pix (out.sci.data, i, j) = array[j][i];
	}

	printf ("Writing %s image.\n", name);

	if (putSingleGroup (name, 1, &out, 0))
	    return (OPEN_FAILED);
	freeSingleGroup (&out);

	printf ("Done writing.\n");

	return (STIS_OK);
}

static int DDebug (char *name, double *array, int nx) {

	SingleGroup out;
	int i;

	initSingleGroup (&out);
	if (allocSingleGroup (&out, nx, 1))
	    return (OUT_OF_MEMORY);
	for (i = 0; i < nx; i++)
	    Pix (out.sci.data, i, 0) = array[i];

	printf ("Writing %s image.\n", name);

	if (putSingleGroup (name, 1, &out, 0))
	    return (OPEN_FAILED);
	freeSingleGroup (&out);

	printf ("Done writing.\n");

	return (STIS_OK);
}

static int CDebug (char *name, CmplxArray *z) {

	SingleGroup out;
	int i, j;

	initSingleGroup (&out);
	if (allocSingleGroup (&out, z->nx, z->ny) == -1)
	    return (OUT_OF_MEMORY);
	for (j = 0; j < z->ny; j++) {
	    for (i = 0; i < z->nx; i++) {
	        Pix (out.sci.data, i, j) = RPIX2D (z, i, j);
	        Pix (out.err.data, i, j) = IPIX2D (z, i, j);
	    }
	}

	printf ("Writing %s image.\n", name);

	if (putSingleGroup (name, 1, &out, 0))
	    return (OPEN_FAILED);
	freeSingleGroup (&out);

	printf ("Done writing.\n");

	return (STIS_OK);
}

static void CReport (CmplxArray *z) {

	int i, j;
	double sumr, sumi;

	sumr = 0.0;
	sumi = 0.0;
	for (j = 0; j < z->ny; j++) {
	    for (i = 0; i < z->nx; i++) {
	        sumr += RPIX2D (z, i, j);
	        sumi += IPIX2D (z, i, j);
	    }
	}

	sumr /= (double)z->nx * (double)z->ny;
	sumi /= (double)z->nx * (double)z->ny;

	printf ("real = %g  imag = %g  n = %d\n", sumr, sumi, z->nx * z->ny);
	fflush (stdout);
}

static void Report (float **a, int nx, int ny) {

	int i, j;
	double sum;

	sum = 0.0;
	for (j = 0; j < ny; j++) {
	    for (i = 0; i < nx; i++)
	        sum += a[j][i];
	}

	sum /= (double)nx * (double)ny;

	printf ("aver = %g  n = %d\n", sum, nx * ny);
	fflush (stdout);
}
*/
