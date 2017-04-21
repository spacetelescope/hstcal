# include <stdio.h>
# include <string.h>

# include "xtables.h"
# include "hstio.h"

# include "stis.h"
# include "stisdef.h"
# include "calstis6.h"
# include "err.h"

/* 
   calstis6 -- This is the entry point for both pipeline and standalone 
               versions. This function basically gets the algorithm 
               selected by either the command line of header keyword, 
               and calls either the  "standard" single pass algorithm 
               or the IDT scattered light correction iterative algorithm.



   Revision history:
   ----------------
   17 Feb 00  -  Implemented (I.Busko)
   18 Jun 03  -  Added CTE correction switch (PB)
    5 Dec 04  -  Corrected bug in IDT status messages (PB)
    6 Oct 06  -  Rename stpos to xoffset, add to the calling sequence of
                 CalStis6IDT.
*/

int CalStis6 (char *input, char *output, int backcorr, int dispcorr,
              int fluxcorr, int helcorr, int sgeocorr, int ctecorr,
              int sc2dcorr, double cl_a2center, int maxsearch, 
              double extrsize, double bk1size, double bk2size, 
              double bk1offset, double bk2offset, double bktilt, int bkord, 
              int sporder, char *xtracalg, int printtime, int verbose, 
              int extrloc, int ccglobal, double ccthresh, int do_profile, 
              int pstep, double wstep, double minsn, char *rejranges, 
              char *profilefile, char *fluxfile, char *outw, double backval,
              double backerr, int variance, int fflux, double psclip,
              double sclip, int lfilter, char *idtfile, double subscale,
              double blazeshift, int bks_mode, int bks_order, double xoffset,
              int pipeline) {

/* arguments
char *input;		i: input file name
char *output;		i: output file name
int  backcorr;          i: calibration switch
int  dispcorr;          i: calibration switch
int  fluxcorr;          i: calibration switch
int  helcorr;           i: calibration switch
int  sgeocorr;          i: calibration switch
int  ctecorr;           i: calibration switch
int  sc2dcorr;		i: calibration switch
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
char *idtfile		i: file with IDT final deconvolved image
double subscal;		i: subsampling factor in profile builder
double blazeshift;	i: blaze shift (in pixels) from command line
int bks_mode;		i: backgr. smoothing mode
int bks_order;		i: backgr. smoothing polynomial order
double xoffset;		i: an offset in dispersion direction
int pipeline;		i: is calstis6 being run from the pipeline ?
*/
	int status;

	IODescPtr im;		/* descriptor for primary header */
	Hdr phdr;		/* primary header from input image */
	int x2dcorr;
	double dummy;

	int CalStis6Std (char *, char *, int, int, int, int, int, int, 
                         double, int, 
                         double, double, double, double, double, double, 
                         int, int, char *, int, int, int, int, double, int, 
                         int, double, double, char *, char *, char *, char *,
                         double, double, int, int, double, double, int, int, 
                         int, double *, double *, double, double, int, int,
                         double);
	int CalStis6IDT (char *, char *, Hdr *,
                         int, int, int, int, int, int, double, int, 
                         double, double, double, double, double, double, 
                         int, int, char *, int, int, int, int, double, int, 
                         int, double, double, char *, char *, char *, char *,
                         double, double, int, int, double, double, int, int,
                         int, char *, double, int, int, double);

	/* Open input image in order to read its primary header. 
           This step will be repeated later, but is necessary in 
           here in order to get the keyword that controls the
           execution of the IDT algorithm.
        */
	im = openInputImage (input, "", 0);
	if (hstio_err())
	    return (OPEN_FAILED);
	initHdr (&phdr);
	getHeader (im, &phdr);
	if (hstio_err())
	    return (OPEN_FAILED);
	closeImage (im);

	/* If it is an _x2d image, bail out. */
	if (GetSwitch (&phdr, "X2DCORR", &x2dcorr) == 0) {
	    if (x2dcorr == COMPLETE) {
	        freeHdr (&phdr);
	        printf ("ERROR    Input file was already rectified (x2d).\n");
	        return (ERROR_RETURN);
	    }
	}

	if (sc2dcorr) {

	    /* IDT extraction. Zero IMSET is meaningless in this 
               context (for now ?). If the H&S scattering correction
               algorithm was selected, revert to unweighted.
            */

	    if (streq_ic (xtracalg, SCATTER_CORR))
	        strcpy (xtracalg, UNWEIGHTED);

	    /* This setting is undone at the last call inside CalStis6IDT */

	    status = CalStis6IDT (input, output, &phdr, backcorr, dispcorr,
                     fluxcorr, helcorr, sgeocorr, ctecorr, cl_a2center,
                     maxsearch, extrsize, bk1size, bk2size, bk1offset,
                     bk2offset, bktilt, bkord, sporder, xtracalg, printtime,
                     verbose == 0? 3: verbose,
                     extrloc, ccglobal, ccthresh, do_profile, pstep, wstep,
                     minsn, rejranges, profilefile, fluxfile, outw, backval,
                     backerr, variance, fflux, psclip, sclip, lfilter, 0,
                     pipeline, idtfile, blazeshift, BKS_OFF,
		     bks_order, xoffset);

	} else {

	    /* Standard extraction. Note zeroed IMSET to force extraction
               on entire file. 
            */

	    status = CalStis6Std (input, output, backcorr, dispcorr,
                     fluxcorr, helcorr, sgeocorr, ctecorr, cl_a2center, 
                     maxsearch, extrsize, bk1size, bk2size, bk1offset,
                     bk2offset, bktilt, bkord, sporder, xtracalg, printtime,
                     verbose+1,
                     extrloc, ccglobal, ccthresh, do_profile, pstep, wstep,
                     minsn, rejranges, profilefile, fluxfile, outw, backval,
                     backerr, variance, fflux, psclip, sclip, lfilter, 0,
                     pipeline, &dummy, &dummy, subscale, blazeshift, bks_mode,
                     bks_order, xoffset);
	}

	freeHdr (&phdr);
	return (status);
}
