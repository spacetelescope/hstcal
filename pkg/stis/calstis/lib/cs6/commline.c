# include <stdio.h>
# include <stdlib.h>		/* atoi, atof */
# include <string.h>

# include "xtables.h"

# include "stis.h"
# include "calstis6.h"
# include "hstcalerr.h"

static int syntax_error (char *);
static int getArgI (char **, int, int *, int *);
static int getArgD (char **, int, int *, double *);

/*
   Processes command-line parameters.

   Revision history:
   ----------------
   20 Feb 97  -  Implemented (I.Busko)
   09 Apr 97  -  Changes after code review (IB):
                 - added functions getArgI, getArgD
                 - included return error code in syntax_error function
   01 May 97  -  Check output file name extension (IB)
   08 May 97  -  Conform to new _trl standard (IB)
   13 May 97  -  Add calibration switches to command line (IB)
   02 Jun 97  -  Add -x1d switch to command line (IB)
   22 Jan 98  -  Moved checkOutName to main cs6 module (IB)
   10 Feb 98  -  Optional output name (IB)
   10 Apr 98  -  Replace debug swicth by -e, remove debug file (IB)
   22 Jun 98  -  Global croscor fit mode, croscor threshold (IB)
   26 Jun 98  -  Profile generator (IB)
   29 Jun 98  -  Rejection ranges in profile generator (IB)
   15 Sep 98  -  File names to support optimal extractions (IB)
   23 Sep 98  -  Output weights file name (IB)
   18 Nov 98  -  Background value and its error (IB)
   11 Dec 98  -  Weigth or variance image in optimal extraction (IB)
   04 May 99  -  FLUX or NET reference spectrum in opt. extraction (IB)
   10 Dec 99  -  Sigma clip in profile builder and extraction (IB)
   05 Jan 00  -  Lee filter size (IB)
   17 Feb 00  -  IDT's scattered light algorithm (IB)
   05 Apr 00  -  Output deconvolved image (IB)
   04 Dec 00  -  Subsampling factor for profile builder (IB)
   16 Apr 02  -  Blaze shift (IB)
   24 Jul 02  -  Background smoothing (IB)
   16 Dec 02  -  Reference star A1 position for slitless data (IB)
   17 Jun 03  -  Add CTE correction switch (PB)
   21 Apr 05  -  Rename variable stpos to xoffset (PEH)
   06 Feb 06  -  Change the default bks_mode (background smoothing mode)
		from BKS_OFF to BKS_MEDIAN, consistent with calstis0 (PEH)
   27 May 11  -  Add new command line option -bn for bks_mode=off (PEH)
   06 Jul 11  -  Add new command line option --version (PEH)
   10 Feb 12  -  If -r is the first word on the command line (following cs6.e),
		print the full version string and exit 0 (PEH)
*/

int CommLine (int argc, char **argv, char *input, char *output,
              int *backcorr, int *dispcorr, int *fluxcorr, int *helcorr,
              int *sgeocorr, int *ctecorr, double *cl_a2center,
              int *maxsearch, double *extrsize, double *bk1size,
              double *bk2size, double *bk1offset, double *bk2offset,
              double *bktilt, int *bkord, int *sporder, char *xtracalg,
              int *printtime, int *verbose, int *extrloc, int *ccglobal,
              double *ccthresh, int *do_profile, int *pstep, double *wstep,
              double *minsn, char *rejranges, char *profilefile,
              char *fluxfile, char *outw, double *backval, double *backerr,
              int *variance, int *fflux, double *psclip, double *sclip,
              int *lfilter, int *idt, char *idtfile, double *subscale,
              double *blazeshift, int *bks_mode, int *bks_order,
              double *xoffset) {

/* arguments
int argc;               i: input command-line parameters
char **argv;
char *input;		o: input file name
char *output;		o: output file name
int  *backcorr;         o: calibration switch
int  *dispcorr;         o: calibration switch
int  *fluxcorr;         o: calibration switch
int  *helcorr;          o: calibration switch
int  *sgeocorr;         o: calibration switch
int  *ctecorr;          o: calibration switch
double *cl_a2center;	o: nominal A2 spectrum center
int *maxsearch;		o: cross correlation range
double *extrsize;	o: size of spectrum box
double *bk1size;	o: size of background box 1
double *bk2size;	o: size of background box 2
double *bk1offset;	o: offset of background box 1
double *bk2offset;	o: offset of background box 2
double *bktilt;		o: angle of background boxes
int *bkord;		o: order of backr. polyn. fit
int *sporder;		o: spectral order
char *xtracalg;		o: extraction algorithm
int *printtime;		o: print time after each step?
int *verbose;		o: print additional info ?
int *extrloc;		o: output extraction info ?
int *ccglobal;		o: use global crosscor everywhere ?
double *ccthresh;	o: croscorr threshold
int *do_profile;	o: use calstis6 as profile generator ?
int *pstep;		o: pixel step used to compress profiles
double *wstep;		o: wavelength step used to compress profiles
double *minsn;		o: minimum acceptable S/N
char *rejranges         o: string with rejection ranges
char *profilefile	o: file with profiles for optimal extraction
char *fluxfile		o: _x1d-type file with flux for optimal extractions
char *outw;		o: output weights file name
double *backval;	o: background value
double *backerr;	o: background error
int *variance		o: variance image instead of weights image ?
int *fflux;		o: use FLUX instead of NET in opt. extraction
double *psclip;		o: sigma-clip in profile builder
double *sclip;		o: sigma-clip in optimal extraction
int *lfilter;		o: scattered light correction Lee filter size
int *idt;		o: IDT's scattered light algorithm
char *idtfile		o: file with IDT final deconvolved image
double *subscale;	o: subsampling factor for profile builder
double *blazeshift;	o: blaze shift value in pixels
int *bks_mode;		o: backgr. smoothing mode
int *bks_order;		o: backgr. smoothing polynomial order
double *xoffset;	o: offset in dispersion direction for slitless data
*/
	int ctoken;	/* current command-line token being processed */
	int switches;   /* are there calibration switches ? */
	int back;	/* these hold what is read from the command     */
	int disp;	/* line, before setting the output. This is to  */
	int flux;	/* handle the logic of double defaults.         */
	int hel;
	int sgeo;
        int cte;

	/* First, set all output from this routine to default values. */
	input[0]       = '\0';
	output[0]      = '\0';
	*cl_a2center   = NO_POSITION;
	*maxsearch     = NO_RANGE;
	*extrsize      = NO_SIZE;
	*bk1size       = NO_SIZE;
	*bk2size       = NO_SIZE;
	*bk1offset     = NO_SIZE;
	*bk2offset     = NO_SIZE;
	*bktilt        = NO_TILT;
	*bkord         = NO_ORDER;
	*sporder       = NO_ORDER;
	xtracalg[0]    = '\0';
	rejranges[0]   = '\0';
	*printtime     = 0;
	*verbose       = 0;
	*extrloc       = 1;
	*ccglobal      = 0;
	*ccthresh      = 0.0;
	profilefile[0] = '\0';
	fluxfile[0]    = '\0';
	outw[0]        = '\0';
	*backval       = NO_VALUE;
	*backerr       = 0.0;
        *variance      = 0;
	*fflux         = 1;
	*do_profile    = 0;
	*lfilter       = 17;
	*blazeshift    = NO_VALUE;
	*idt           = 0;
	idtfile[0]     = '\0';
	*bks_mode      = BKS_MEDIAN;
	*bks_order     = 3;
	*xoffset       = 0.;
	*psclip        = 5.0;
	*sclip         = 6.0;
	*minsn         = 0.0;
	*pstep         = 0.0;   /* Zero in here disables wavel./pixel scan */
	*wstep         = 0.0;   /* and force profile summing over default  */
                                /* ranges.                                 */

	/* Calibration switches are set to their defaults, but if
           at least one appears explictly in the command line, the
           default becomes OMIT, so just the explictly set ones are
           activated.
        */
	*backcorr  = DefSwitch ("backcorr");
	*dispcorr  = DefSwitch ("dispcorr");
	*fluxcorr  = DefSwitch ("fluxcorr");
	*helcorr   = DefSwitch ("helcorr");
	*sgeocorr  = DefSwitch ("sgeocorr");
	*ctecorr   = DefSwitch ("ctecorr");

	back = 0;
	disp = 0;
	flux = 0;
	hel  = 0;
	sgeo = 0;
        cte  = 0;
	switches = 0;

	if (argc > 1 && strcmp (argv[1], "-r") == 0) {
	    PrFullVersion();
	    exit (0);
	}
	for (ctoken = 1;  ctoken < argc;  ctoken++) {
	    if (strcmp (argv[ctoken], "--version") == 0) {
		PrVersion();
		exit (0);
	    }
	}

	if (argc == 1) {
printf ("cs6 input [output] [-t] [-v] [-e] [-back] [-disp] [-flux] [-hel]\n");
printf ("    [-cte] [-idt] [-g] [-c a2center] [-r maxsearch] [-p ccthresold]\n");
printf ("    [-x extrsize] [-if sc2dfile] [-bs blazeshift]\n");
printf ("    [-b1 bk1size] [-b2 bk2size] [-o1 bk1offset] [-o2 bk2offset]\n");
printf ("    [-k bktilt] [-n bkord] [-s sporder] [-a xtracalg]\n");
printf ("    [-bk backval] [-be backerr] [-sp sclip] [-l lfiltersize]\n");
printf ("    [-pf profilefile] [-px fluxfile] [-pa] [-wf weightsfile]\n");
printf ("    [-va] [-bb|-bm|-bn] [-bo bsorder] [-st xoffset]\n");
printf ("    Box sizes and offsets are in REFERENCE pixel units, \n");
printf ("    extraction postion a2center is in IMAGE pixel units. \n");
	    return (0);
	}
	/* Additional command-line arguments not included above are:
	   -pr  --> do_profile = 1
	   -bp  pstep
	   -bw  wstep
	   -re  rejranges
	   -sc  sclip
	   -sn  minsn
	   -ss  subscale
	*/

	/* Get name of input file/list/template. This is mandatory. */
	strcpy (input,  argv[1]);

	/* If second argument starts with a '-', this means that
           there is no output file name and it must be built
           elsewhere. An empty string signals this condition.
        */
	if (argc > 2) {
	    if (argv[2][0] != '-') {
	        strcpy (output, argv[2]);  /* output name exists */
	        if (argc > 3)
	            ctoken = 3;            /* additional parameters */
	        else
	            ctoken = argc;         /* no additional parameters */
	    } else {
	        strcpy (output, "");       /* output name does not exist */
	        ctoken = 2;
	    }
	} else
	    ctoken = argc;                 /* no remaining parameters */

	/* Scan eventual remaining parameters. Notice that two-character
           switches must be scanned before one-character ones, to avoid
           confusion between switches that begin with the same first
           character, as in -p and -pr.
        */
	while (ctoken < argc) {

	    /* Parameters must begin with a '-' ! */
	    if (argv[ctoken][0] != '-')
	        return (syntax_error (argv[ctoken]));

	    else {

	        /* These do not require additional arguments. */

	        if (argv[ctoken][1] == 't') {
	            *printtime = 1;
	            ctoken++;

	        } else if (argv[ctoken][1] == 'v' &&
	                   argv[ctoken][2] == 'a') {
	            *variance = 1;
	            ctoken++;

	        } else if (argv[ctoken][1] == 'v') {
	            *verbose = 1;
	            ctoken++;

	        } else if (argv[ctoken][1] == 'g') {
	            *ccglobal = 1;
	            ctoken++;

	        } else if (argv[ctoken][1] == 'e') {
	            *extrloc = 0;
	            ctoken++;

	        } else if (argv[ctoken][1] == 'x' &&
	                   argv[ctoken][2] == '1' &&
	                   argv[ctoken][3] == 'd') {
	            switches = 1;
	            ctoken++;

	        } else if (argv[ctoken][1] == 'b' &&
	                   argv[ctoken][2] == 'a' &&
	                   argv[ctoken][3] == 'c' &&
	                   argv[ctoken][4] == 'k') {
	            back = 1;
	            switches = 1;
	            ctoken++;

	        } else if (argv[ctoken][1] == 'd' &&
	                   argv[ctoken][2] == 'i' &&
	                   argv[ctoken][3] == 's' &&
	                   argv[ctoken][4] == 'p') {
	            disp = 1;
	            switches = 1;
	            ctoken++;

	        } else if (argv[ctoken][1] == 'f' &&
	                   argv[ctoken][2] == 'l' &&
	                   argv[ctoken][3] == 'u' &&
	                   argv[ctoken][4] == 'x') {
	            flux = 1;
	            switches = 1;
	            ctoken++;

	        } else if (argv[ctoken][1] == 'h' &&
	                   argv[ctoken][2] == 'e' &&
	                   argv[ctoken][3] == 'l') {
	            hel = 1;
	            switches = 1;
	            ctoken++;

	        } else if (argv[ctoken][1] == 'i' &&
	                   argv[ctoken][2] == 'd' &&
	                   argv[ctoken][3] == 't') {

	            *idt = 1;
	            ctoken++;

	        } else if (argv[ctoken][1] == 's' &&
	                   argv[ctoken][2] == 'g' &&
	                   argv[ctoken][3] == 'e' &&
	                   argv[ctoken][4] == 'o') {
	            sgeo = 1;
	            switches = 1;
	            ctoken++;

	        } else if (argv[ctoken][1] == 'c' &&
	                   argv[ctoken][2] == 't' &&
	                   argv[ctoken][3] == 'e') {
	            cte = 1;
                    switches = 1;
	            ctoken++;

	        } else if (argv[ctoken][1] == 'b' &&
	                   argv[ctoken][2] == 'b') {
	            *bks_mode = BKS_AVERAGE;
	            ctoken++;

	        } else if (argv[ctoken][1] == 'b' &&
	                   argv[ctoken][2] == 'm') {
	            *bks_mode = BKS_MEDIAN;
	            ctoken++;

	        } else if (argv[ctoken][1] == 'b' &&
	                   argv[ctoken][2] == 'n') {
	            *bks_mode = BKS_OFF;
	            ctoken++;

	        } else if (argv[ctoken][1] == 'p' &&
	                   argv[ctoken][2] == 'r') {
	            *do_profile = 1;
	            ctoken++;

	        } else if (argv[ctoken][1] == 'p' &&
	                   argv[ctoken][2] == 'a') {
	            *fflux = 0;
	            ctoken++;

	       /* These require one additional argument, which is
                   handled either by the getArg functions or is just
                   a string to be copied.
               */

	        } else if (argv[ctoken][1] == 'a') {
	            if (ctoken <= (argc-2)) {
	                strcpy (xtracalg, argv[++ctoken]);
	                ctoken++;
	            } else
	                return (syntax_error (argv[ctoken]));

	        } else if (argv[ctoken][1] == 'b' &&
                           argv[ctoken][2] == 'o') {
	            if (getArgI (argv, argc, &ctoken, bks_order))
	                return (1);

	        } else if (argv[ctoken][1] == 'b' &&
                           argv[ctoken][2] == 'p') {
	            if (getArgI (argv, argc, &ctoken, pstep))
	                return (1);

	        } else if (argv[ctoken][1] == 'b' &&
                           argv[ctoken][2] == 'w') {
	            if (getArgD (argv, argc, &ctoken, wstep))
	                return (1);

	        } else if (argv[ctoken][1] == 'b' &&
                           argv[ctoken][2] == '1') {
	            if (getArgD (argv, argc, &ctoken, bk1size))
	                return (1);

	        } else if (argv[ctoken][1] == 'b' &&
                           argv[ctoken][2] == '2') {
	            if (getArgD (argv, argc, &ctoken, bk2size))
	                return (1);

	        } else if (argv[ctoken][1] == 'b' &&
                           argv[ctoken][2] == 'k') {
	            if (getArgD (argv, argc, &ctoken, backval))
	                return (1);

	        } else if (argv[ctoken][1] == 'b' &&
                           argv[ctoken][2] == 'e') {
	            if (getArgD (argv, argc, &ctoken, backerr))
	                return (1);

	        } else if (argv[ctoken][1] == 'b' &&
                           argv[ctoken][2] == 's') {
	            if (getArgD (argv, argc, &ctoken, blazeshift))
	                return (1);

	        } else if (argv[ctoken][1] == 'c') {
	            if (getArgD (argv, argc, &ctoken, cl_a2center))
	                return (1);
	            (*cl_a2center) -= 1.0;  /* 0-indexing */

	        } else if (argv[ctoken][1] == 'i' &&
                           argv[ctoken][2] == 'f') {
	            if (ctoken <= (argc-2)) {
	                strcpy (idtfile, argv[++ctoken]);
	                ctoken++;
	            } else
	                return (syntax_error (argv[ctoken]));

	        } else if (argv[ctoken][1] == 'k') {
	            if (getArgD (argv, argc, &ctoken, bktilt))
	                return (1);

	        } else if (argv[ctoken][1] == 'l') {
	            if (getArgI (argv, argc, &ctoken, lfilter))
	                return (1);

	        } else if (argv[ctoken][1] == 'n') {
	            if (getArgI (argv, argc, &ctoken, bkord))
	                return (1);

	        } else if (argv[ctoken][1] == 'o' &&
                           argv[ctoken][2] == '1') {
	            if (getArgD (argv, argc, &ctoken, bk1offset))
	                return (1);

	        } else if (argv[ctoken][1] == 'o' &&
                           argv[ctoken][2] == '2') {
	            if (getArgD (argv, argc, &ctoken, bk2offset))
	                return (1);

	        } else if (argv[ctoken][1] == 'p' &&
                           argv[ctoken][2] == 'f') {
	            if (ctoken <= (argc-2)) {
	                strcpy (profilefile, argv[++ctoken]);
	                ctoken++;
	            } else
	                return (syntax_error (argv[ctoken]));

	        } else if (argv[ctoken][1] == 'p' &&
                           argv[ctoken][2] == 'x') {
	            if (ctoken <= (argc-2)) {
	                strcpy (fluxfile, argv[++ctoken]);
	                ctoken++;
	            } else
	                return (syntax_error (argv[ctoken]));

	        } else if (argv[ctoken][1] == 'p') {
	            if (getArgD (argv, argc, &ctoken, ccthresh))
	                return (1);

	        } else if (argv[ctoken][1] == 'r' &&
                           argv[ctoken][2] == 'e') {
	            if (ctoken <= (argc-2)) {
	                strcpy (rejranges, argv[++ctoken]);
	                ctoken++;
	            } else
	                return (syntax_error (argv[ctoken]));

	        } else if (argv[ctoken][1] == 'r') {
	            if (getArgI (argv, argc, &ctoken, maxsearch))
	                return (1);

	        } else if (argv[ctoken][1] == 's' &&
                           argv[ctoken][2] == 'c') {
	            if (getArgD (argv, argc, &ctoken, sclip))
	                return (1);

	        } else if (argv[ctoken][1] == 's' &&
                           argv[ctoken][2] == 'n') {
	            if (getArgD (argv, argc, &ctoken, minsn))
	                return (1);

	        } else if (argv[ctoken][1] == 's' &&
                           argv[ctoken][2] == 'p') {
	            if (getArgD (argv, argc, &ctoken, psclip))
	                return (1);

	        } else if (argv[ctoken][1] == 's' &&
                           argv[ctoken][2] == 's') {
	            if (getArgD (argv, argc, &ctoken, subscale))
	                return (1);

	        } else if (argv[ctoken][1] == 's' &&
                           argv[ctoken][2] == 't') {
	            if (getArgD (argv, argc, &ctoken, xoffset))
	                return (1);

	        } else if (argv[ctoken][1] == 's') {
	            if (getArgI (argv, argc, &ctoken, sporder))
	                return (1);

	        } else if (argv[ctoken][1] == 'w' &&
                           argv[ctoken][2] == 'f') {
	            if (ctoken <= (argc-2)) {
	                strcpy (outw, argv[++ctoken]);
	                ctoken++;
	            } else
	                return (syntax_error (argv[ctoken]));

	        } else if (argv[ctoken][1] == 'x') {
	            if (getArgD (argv, argc, &ctoken, extrsize))
	                return (1);

	        /* No match. */
	        } else
	            return (syntax_error (argv[ctoken]));
	    }
	}

	/* Here we change the calibration switch defaults to OMIT and
           set just the chosen ones. x1d is always set to PERFORM.
        */
	if (switches) {
	    *backcorr = OMIT;
	    *dispcorr = OMIT;
	    *fluxcorr = OMIT;
	    *helcorr  = OMIT;
	    *sgeocorr = OMIT;
            *ctecorr  = OMIT;
	    if (back) *backcorr = PERFORM;
	    if (disp) *dispcorr = PERFORM;
	    if (flux) *fluxcorr = PERFORM;
	    if (hel)  *helcorr  = PERFORM;
	    if (sgeo) *sgeocorr = PERFORM;
	    if (cte)  *ctecorr  = PERFORM;
	}

	return (0);
}


/*  Get a int parameter from the argv array. */

static int getArgI (char **argv, int argc, int *ctoken, int *value) {
	if (*ctoken <= (argc-2)) {
	    *value = atoi (argv[++(*ctoken)]);
	    (*ctoken)++;
	    return (0);
	} else
	    return (syntax_error (argv[*ctoken]));
}



/*  Get a double parameter from the argv array. */

static int getArgD (char **argv, int argc, int *ctoken, double *value) {
	if (*ctoken <= (argc-2)) {
	    *value = atof (argv[++(*ctoken)]);
	    (*ctoken)++;
	    return (0);
	} else
	    return (syntax_error (argv[*ctoken]));
}


static int syntax_error (char *msg) {

	printf ("ERROR  Syntax error: %s\n", msg);
	return (1);
}
