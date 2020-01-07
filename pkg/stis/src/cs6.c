# include <stdio.h>
# include <stdlib.h>
# include <string.h>

#include "hstcal_memory.h"
# include "c_iraf.h"		/* for c_irafinit, IRAFPointer */
# include "ximio.h"		/* for the imt routines */

# include "stis.h"
# include "hstcalerr.h"
# include "calstis6.h"

/*
   Main entry point for the standalone version of calstis6.

   This module calls the Commline function to decode the command line
   parameters, and then expands the input and output file name lists
   and calls the CalStis6 function in a loop to process each entry in
   turrn.


   Revision history:
   ----------------
   20 Feb 97  -  Implemented (I.Busko)
   13 May 97  -  Calibration switches are input from command line (IB)
   14 May 97  -  Flag tells CalStis6 it is being run from the pipeline (IB)
   22 Jan 98  -  Check output file name extension (moved from Commline, IB)
   22 Jan 98  -  Expand file template/lists (IB)
   10 Feb 98  -  Add MkOutName call (IB)
   10 Apr 98  -  Replace debug switch by extrloc, remove debug file (IB)
   22 Jun 98  -  Global croscor fit mode, crosscorr theshold (IB)
   26 Jun 98  -  Profile generator (IB)
   29 Jun 98  -  Rejection ranges in profile generator (IB)
   14 Sep 98  -  File names to support optimal extractions (IB)
   23 Sep 98  -  Output weights image name (IB)
   18 Nov 98  -  Background value and its error read from command line (IB)
   11 Dec 98  -  Weigth or variance image in optimal extraction (IB)
   04 May 99  -  FLUX or NET reference spectrum in opt. extraction (IB)
   10 Dec 99  -  Sigma clip in profile builder and extraction (IB)
   05 Jan 00  -  Lee filter size (IB)
   17 Feb 00  -  IDT's scattered light algorithm (IB)
   05 Apr 00  -  Output deconvolved image (IB)
   05 Oct 00  -  sc2d switch (IB)
   04 Dec 00  -  Subsampling factor for profile builder (IB)
   16 Apr 02  -  Blaze shift from command line (IB)
   16 Dec 02  -  reference star A1 position for slitless data (IB)
   18 Jun 03  -  Add CTE correction switch (PB)
    6 Oct 06  -  Rename stpos to xoffset.
   16 Dec 11  -  Include <stdlib.h> for the declaration of exit().
*/

int main (int argc, char **argv) {

	int status = 0;			/* zero is OK */
	char input[STIS_LINE];	/* list names */
	char output[STIS_LINE];
	char outw[STIS_LINE];
	char finput[STIS_FNAME];/* file names */
	char foutput[STIS_FNAME];
	char foutw[STIS_FNAME];
	IRAFPointer inlist;	/* input file list */
	IRAFPointer outlist;	/* output file list */
	IRAFPointer outwlist;	/* output weights file list */
	int innf, outnf, outwnf;/* number of files in lists */
	int  backcorr;		/* calibration switch */
	int  dispcorr;		/* calibration switch */
	int  fluxcorr;		/* calibration switch */
	int  helcorr;		/* calibration switch */
	int  sgeocorr;		/* calibration switch */
        int  ctecorr;           /* calibration switch */
	double cl_a2center;	/* nominal A2 spectrum center */
	int maxsearch;		/* cross correlation range */
	double extrsize;	/* size of spectrum box */
	double bk1size;		/* size of background box 1 */
	double bk2size;		/* size of background box 2 */
	double bk1offset;	/* offset of background box 1 */
	double bk2offset;	/* offset of background box 2 */
	double bktilt;		/* angle of background boxes */
	double blazeshift;	/* blaze shift (in pixels) from comm. line */
	int bkord;		/* order of backr. polyn. fit */
	int sporder;		/* spectral order */
	char xtracalg[STIS_CBUF];/* extraction algorithm */
	int printtime;		/* print time after each step ? */
	int verbose;		/* print additional info ? */
	int extrloc;		/* output extraction location info ? */
	int ccglobal;		/* use global crosscor everywhere ? */
	double ccthresh;	/* crosscorr threshold */
	int do_profile;		/* use calstis6 as profile generator ? */
	int pstep;		/* pixel step used to compress profiles */
	double wstep;		/* wavelength step used to compress profiles*/
	double minsn;		/* minimum acceptable S/N */
	double psclip;		/* sigma clip in profile builder */
	double sclip;		/* sigma clip in extraction */
	char rejranges[STIS_LINE];/* rejection ranges */
	char profilefile[STIS_FNAME]; /* profiles for optimal extractions */
	char fluxfile[STIS_FNAME]; /* fluxes for optimal extractions */
	double backval;	        /* background value from command line */
	double backerr;		/* background error from command line */
	int variance;		/* variance image instead of weight image ? */
	int fflux;		/* FLUX instead of NET in opt. extraction ? */
	int lfilter;		/* Lee filter window size */
	int idt;		/* IDT's scattered light algorithm */
	char idtfile[STIS_FNAME]; /* file with IDT final deconvolved image */
	double subscale;	/* subsampling factor for profile builder */
	double xoffset;		/* an offset in dispersion direction */
	int bks_mode;		/* back. smooth mode */
	int bks_order;		/* back. smooth polynomila order */
	/*int bks_size;*/	/* back. smooth box size */
	int ifile;

	char *isuffix[]  = {"_flt"};  /* default suffixes */
	char *osuffix[]  = {"_x1d"};
	char *owsuffix[] = {"_owt"};
	int nsuffix = 1;

	int CalStis6 (char *, char *, int, int, int, int, int, int, int,
                      double, int, double, double, double, double, double,
                      double, int, int, char *, int, int, int, int, double,
                      int, int, double, double, char *, char *, char *,
                      char *, double, double, int, int, double, double, int,
                      char *, double, double, int, int, double, int);
	int CommLine (int, char**, char *, char *, int *, int *, int *, int *,
                      int *, int *, double *, int *, double *, double *,
                      double *, double *, double *, double *, int *, int *,
                      char *, int *, int *, int *, int *, double *, int *,
                      int *, double *, double *, char *, char *, char *,
                      char *, double *, double *, int *, int *, double *,
                      double *, int *, int *, char *, double *, double *,
                      int *, int *, double *);

	c_irafinit (argc, argv);

	/* Get command line parameters. */
	if (CommLine (argc, argv, input, output, &backcorr, &dispcorr,
                      &fluxcorr, &helcorr, &sgeocorr, &ctecorr, &cl_a2center,
                      &maxsearch, &extrsize, &bk1size, &bk2size, &bk1offset,
                      &bk2offset, &bktilt, &bkord, &sporder, xtracalg,
                      &printtime, &verbose, &extrloc, &ccglobal, &ccthresh,
                      &do_profile, &pstep, &wstep, &minsn, rejranges,
                      profilefile, fluxfile, outw, &backval, &backerr,
                      &variance, &fflux, &psclip, &sclip, &lfilter, &idt,
                      idtfile, &subscale, &blazeshift, &bks_mode, &bks_order,
                      &xoffset))
	    exit (ERROR_RETURN);

	/* Expand input and output file lists. */
	PtrRegister ptrReg;
	initPtrRegister(&ptrReg);
	inlist   = c_imtopen (input);
	addPtr(&ptrReg, inlist, &c_imtclose);
	outlist  = c_imtopen (output);
    addPtr(&ptrReg, outlist, &c_imtclose);
	outwlist = c_imtopen (outw);
    addPtr(&ptrReg, outwlist, &c_imtclose);
	innf     = c_imtlen (inlist);
	outnf    = c_imtlen (outlist);
	outwnf   = c_imtlen (outwlist);
	if (outnf != 0 && innf != outnf) {
	    printf ("ERROR: Input and output lists have different sizes.\n");
	    freeOnExit(&ptrReg);
	    exit (ERROR_RETURN);
	}
	if (outwnf != 0 && innf != outwnf) {
	    printf ("ERROR: Input and weights lists have different sizes.\n");
	    freeOnExit(&ptrReg);
	    exit (ERROR_RETURN);
	}
	for (ifile = 0; ifile < innf; ifile++) {
	    c_imtgetim (inlist,  finput,  STIS_FNAME);
	    if (outnf != 0)
	        c_imtgetim (outlist, foutput, STIS_FNAME);
	    else
	        strcpy (foutput, "");
	    if (outwnf != 0)
	        c_imtgetim (outwlist, foutw, STIS_FNAME);
	    else
	        strcpy (foutw, "");

	    /* Check or build output file names. */

	    if (MkOutName (finput, isuffix, osuffix, nsuffix,
                foutput, STIS_FNAME))
	    {
	        freeOnExit(&ptrReg);
	        exit (ERROR_RETURN);
	    }
	    if (MkOutName (finput, isuffix, owsuffix, nsuffix,
                foutw, STIS_FNAME))
	    {
	        freeOnExit(&ptrReg);
	        exit (ERROR_RETURN);
	    }

	    /* Extract 1-D STIS data. */

	    if ((status = CalStis6 (finput, foutput,
			  backcorr, dispcorr, fluxcorr, helcorr, sgeocorr,
                          ctecorr, idt, cl_a2center, maxsearch, extrsize,
                          bk1size, bk2size, bk1offset, bk2offset, bktilt,
                          bkord, sporder, xtracalg, printtime, verbose,
                          extrloc, ccglobal, ccthresh, do_profile, pstep,
                          wstep, minsn, rejranges, profilefile, fluxfile,
                          foutw, backval, backerr, variance, fflux, psclip,
                          sclip, lfilter, idtfile, subscale, blazeshift,
                          bks_mode, bks_order, xoffset, 0))) {
	        printf ("Error processing %s.\n", finput);
	        WhichError (status);
	    }
	    if (status)
	    {
	        freeOnExit(&ptrReg);
	        exit (ERROR_RETURN);
	    }
	}
	freeOnExit(&ptrReg);
	exit (0);
}
