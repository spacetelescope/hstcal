# include <stdio.h>
# include <string.h>

# include "xtables.h"
# include "hstio.h"

# include "stis.h"
# include "stisdef.h"
# include "calstis6.h"
# include "hstcalerr.h"

/*
   calstis6std -- Standard 1-D spectral extraction and profile generator.

   The "standard" single pass algorithm is executed.


   Revision history:
   ----------------
   20 Feb 97  -  Implemented (I.Busko)
   11 Apr 97  -  Changes after code review (IB):
                 - added TimeStamp6 function.
   28 Apr 97  -  Add UFilename and UCalVer functions (IB)
   02 May 97  -  Do1Dx routine physically writes the primary header (IB)
   09 May 97  -  Conform to new _trl standard, remove TimeStamp (IB)
   13 May 97  -  Command line calibration switches (IB)
   14 May 97  -  Flag signals CalStis6 is being run from the pipeline (IB)
   27 Jan 98  -  PCTAB support (IB)
   10 Apr 98  -  Replace debug swicth by extrloc, remove debug file (IB)
   22 Jun 98  -  Global croscor fit mode, crosscor threshold (IB)
   24 Jun 98  -  Profile generator (IB)
   29 Jun 98  -  Rejection ranges in profile generator (IB)
   14 Sep 98  -  Support for optimal extractions (IB)
   23 Sep 98  -  Free header here, not in Do1Dx (IB)
   23 Sep 98  -  Weights image (IB)
   18 Nov 98  -  Background value and its error from command line (IB)
   04 May 99  -  FLUX or NET reference spectrum in opt. extraction (IB)
   09 Dec 99  -  Sigma clip in profile builder and extraction (IB)
   05 Jan 00  -  Lee filter size (IB)
   17 Feb 00  -  Old CalStis6 reanmed to CalStis6Std (IB)
   17 Feb 00  -  IMSET selection (IB)
   23 Feb 00  -  Return aperture size (IB)
   14 Apr 00  -  Verbosity control (IB)
   04 Dec 00  -  Subsampling factor for profile builder (IB)
   08 Mar 02  -  Add back TimeStamp removed on 09 May 97 (IB)
   16 Apr 02  -  Blaze shift from command line (IB)
   24 Jul 02  -  Background smoothing (IB)
   16 Dec 02  -  Reference star A1 position for slitless data (IB)
   18 Jun 03  -  Added CTE correction switch (PB)
    4 Feb 04  -  Initialize tdscorr (PEH)
    8 Apr 05  -  Initialize gaccorr (PEH)
   21 Apr 05  -  Rename variable stpos to xoffset (PEH)
*/

int CalStis6Std (char *input, char *output, int backcorr, int dispcorr,
                 int fluxcorr, int helcorr, int sgeocorr, int ctecorr,
                 double cl_a2center, int maxsearch, double extrsize,
                 double bk1size, double bk2size, double bk1offset,
                 double bk2offset, double bktilt, int bkord, int sporder,
                 char *xtracalg, int printtime, int verbose, int extrloc,
                 int ccglobal, double ccthresh, int do_profile, int pstep,
                 double wstep, double minsn, char *rejranges,
                 char *profilefile, char *fluxfile, char *outw, double backval,
                 double backerr, int variance, int fflux, double psclip,
                 double sclip, int lfilter, int imset, int pipeline,
                 double *ap_xsize, double *ap_ysize, double subscale,
                 double blazeshift, int bks_mode, int bks_order,
                 double xoffset) {

/* arguments
char *input;		i: input file name
char *output;		i: output file name
int  backcorr;          i: calibration switch
int  dispcorr;          i: calibration switch
int  fluxcorr;          i: calibration switch
int  helcorr;           i: calibration switch
int  sgeocorr;          i: calibration switch
int  ctecorr;           i: calibration switch
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
double *ap_xsize;       o: aperture size
double *ap_ysize;       o: aperture size
double subscale;	i: subsampling factor in profile builder
double blazeshift;	i: blaze shift (in pixels) from command line
int bks_mode;		i: backgr. smoothing mode
int bks_order;		i: backgr. smoothing polynomial order
double xoffset;		i: for slitless data, an offset in dispersion direction
*/
	int status;

	StisInfo6 sts;	/* calibration switches, reference files, etc. */

	IODescPtr im;		/* descriptor for primary header */
	Hdr phdr;		/* primary header from input image */

	int CheckOptimal (StisInfo6 *);
	int Do1Dx (StisInfo6 *, Hdr *);
	int GetFlags6 (StisInfo6 *, Hdr *);
	int GetKeyInfo6 (StisInfo6 *, Hdr *);
	int GetRefCommLine (StisInfo6 *);
	int History6 (StisInfo6 *, Hdr *, int);
	void StisInit6 (StisInfo6 *);

	sts.pipeline = pipeline;

	/* Initialize structure containing calstis6 information. */
	StisInit6 (&sts);

	/* These are set from the argument list. */
	sts.backcorr  = backcorr;
	sts.dispcorr  = dispcorr;
	sts.fluxcorr  = fluxcorr;
	sts.heliocorr = helcorr;
	sts.sgeocorr  = sgeocorr;

	/* Copy command-line arguments into sts. */
	strcpy (sts.input, input);
	strcpy (sts.output, output);
	strcpy (sts.outw, outw);
	sts.cl_a2center = cl_a2center;
	sts.xoffset     = xoffset;
	sts.maxsearch   = maxsearch;
	sts.extrsize    = extrsize;
	sts.bksize[0]   = bk1size;
	sts.bksize[1]   = bk2size;
	sts.bkoffset[0] = bk1offset;
	sts.bkoffset[1] = bk2offset;
	sts.bktilt      = bktilt;
	sts.bkord       = bkord;
	sts.sporder     = sporder;
	sts.lfilter	= lfilter;
	sts.blazeshift  = blazeshift;
	strcpy (sts.xtracalg, xtracalg);
	strcpy (sts.profilefile, profilefile);
	strcpy (sts.fluxfile, fluxfile);
	sts.printtime      = printtime;
	sts.verbose        = verbose;
	sts.variance_image = variance;
	sts.fflux          = fflux;
	sts.extrloc        = extrloc;
	sts.cc_global      = ccglobal;
	sts.cc_thresh      = ccthresh;
	strcpy (sts.rejranges, rejranges);
	if (pipeline) {
	    sts.backval = NO_VALUE;
	    sts.backerr = NO_VALUE;
	} else {
	    sts.backval = backval;
	    sts.backerr = backerr;
	}

	sts.do_profile    = do_profile;
	sts.profile_pstep = pstep;
	sts.profile_wstep = wstep;
	sts.profile_minsn = minsn;
	sts.psclip        = psclip;
	sts.sclip         = sclip;
	sts.imset         = imset;
	sts.ap_xsize      = 0.0;
	sts.ap_ysize      = 0.0;
	sts.subscale      = subscale;
	sts.bks_mode      = bks_mode;
	sts.bks_order     = bks_order;
        sts.ctecorr       = ctecorr;

	/* Print greeting. */
	if (!sts.do_profile && (verbose == 1 || verbose == 2)) {
	    PrBegin (6);
	    if (verbose == 2)
	        TimeStamp ("CALSTIS-6 started", "");
	    printf ("\n");
	}

	/* PCT, GAC & TDS corrections are performed only if fluxcorr is
	   selected.
	*/
	if (sts.fluxcorr == PERFORM) {
	    sts.pctcorr = PERFORM;
	    sts.gaccorr = PERFORM;
	    sts.tdscorr = PERFORM;
	} else {
	    sts.pctcorr = OMIT;
	    sts.gaccorr = OMIT;
	    sts.tdscorr = OMIT;
	}

	/* Open input image in order to read its primary header. */
	im = openInputImage (sts.input, "", 0);
	if (hstio_err())
	    return (OPEN_FAILED);
	initHdr (&phdr);
	getHeader (im, &phdr);
	if (hstio_err())
	    return (OPEN_FAILED);
	closeImage (im);

	/* Get keyword values from primary header. */
	if ((status = GetKeyInfo6 (&sts, &phdr)))
	    return (status);

	/* Check whether output file already exist. */
	if ((status = FileExists (sts.output)))
	    return (status);

	/* Get calibration flags and file names from input image header.
           This function will not read switches from the header if cs6
           is being run in standalone mode, but will still check for
           inconsistencies in between the several possible switch
           combinations. Nothing related to optimal extraction is
           handled yet.
        */
	if ((status = GetFlags6 (&sts, &phdr)))
	    return (status);

	/* See if reference file names were supplied in the command line.
           This routine must be called _after_ GetFlags6 so that names
           provided in the command line supersede names taken from header.
           Currently this handles only the file names associated with
           optimal extraction.
        */
	if ((status = GetRefCommLine (&sts)))
	    return (status);

	/* Check consistency of all information associated with optimal
           extraction available at this point. This may or may not
           generate an abort. Notice that more information will be
           available only at the point where extraction information is
           retrieved from XTRACTAB for each particular spectral order.
        */
	if ((status = CheckOptimal (&sts)))
	    return (status);

	/* Create output's primary header, update FILENAME and CAL_VER
           keywords, and and add history to it. The actual output is
           delayed until the first output extension becomes ready. In
           this way no incomplete file is created if the program aborts
           prematurely.
        */
	UFilename (sts.output, &phdr);
	UCalVer (&phdr);
	if ((status = History6 (&sts, &phdr, 0)))
	    return (status);

	/* Print information about the input file. */
	if (verbose == 1 || verbose == 2) {
	    PrFileName ("Input",    sts.input);
	    PrFileName ("Output",   sts.output);
	    PrFileName ("Rootname", sts.rootname);
	    PrHdrInfo (sts.obsmode, sts.aperture, sts.opt_elem, sts.det);
	    printf ("\n");
	}

	/* Do 1-D spectral extraction. */
	if ((status = Do1Dx (&sts, &phdr)))
	    return (status);

	/* These are used by the IDT algorithm. */
	*ap_xsize = sts.ap_xsize;
	*ap_ysize = sts.ap_ysize;

	freeHdr (&phdr);
	if (!sts.do_profile && (verbose == 1 || verbose == 2)) {
	    PrEnd (6);
	    if (verbose == 2)
	        TimeStamp ("CALSTIS-6 completed", sts.rootname);
	}

	return (0);
}
