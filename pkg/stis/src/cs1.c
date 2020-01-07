/* calstis1 -- basic 2-D image reduction */

# include <stdio.h>
# include <stdlib.h>		/* calloc */
# include <time.h>
# include <string.h>

#include "hstcal_memory.h"
# include "c_iraf.h"		/* for c_irafinit */
# include "ximio.h"

# include "stis.h"
# include "calstis1.h"
# include "hstcalerr.h"
# include "hstcalversion.h"

static int CompareNumbers (int, int, char *);
static void printSyntax(void)
{
    printf ("syntax:  cs1.e [--help] [-t] [-v] [--version] [--gitinfo] input output [outblev]\n");
    printf ("  command-line switches:\n");
    printf ("       -dqi  -atod -blev\n");
    printf ("       -dopp -lors -glin -lflg\n");
    printf ("       -bias -dark -flat -shad -phot -stat\n");
}
static void printHelp(void)
{
    printSyntax();
}

/* This is the main module for calstis1.  It gets the input and output
   file names, calibration switches, and flags, and then calls CalStis1.

   Phil Hodge, 1997 Dec 11:
	Define RefFileInfo refnames, and pass it (an empty list) to CalStis1.

   Phil Hodge, 1998 Jan 27:
	Accept lists of input and output file names, and loop over
	individual file names.

   Phil Hodge, 1998 Feb 5:
	Call MkOutName instead of MkName.

   Phil Hodge, 1998 June 4:
	Add ("_wav", "_fwv") to the list of suffix pairs for MkOutName.

   Phil Hodge, 1998 Oct 5:
	Change "j = c_imtgetim" to "junk = c_imtgetim".
	Change status value 1 to ERROR_RETURN.

   Ivo Busko, 2002 Mar 19:
	Add support for darkscale command line parameter.

   Phil Hodge, 2011 July 6:
	Include command-line option '--version'.

   Phil Hodge. 2012 Jan 23:
	Allocate memory for inlist, outlist, blevlist based on the size
	of the string on the command line, rather than using a fixed size.

   Phil Hodge, 2012 Feb 10:
	Include command-line option '-r'.
*/

int main (int argc, char **argv) {

	int status;			/* zero is OK */

	char *inlist;		/* list of input file names */
	char *outlist;		/* list of output file names */
	char *blevlist;		/* list of output blev file names */
	int switch_on = 0;	/* was any switch specified? */
	int printtime = 0;	/* print time after each step? */
	int verbose = 0;	/* print additional info? */
	int too_many = 0;	/* too many command-line arguments? */
	int i, j;		/* loop indexes */
	int junk;

	IRAFPointer i_imt, o_imt, b_imt;	/* imt list pointers */
	char *input;		/* name of input science file */
	char *output;		/* optional name of output file */
	char *outblev;		/* optional file for blev values */
	int n_in, n_out, n_blev;	/* number of files in each list */
	int n;

	/* Input and output suffixes. */
	char *isuffix[] = {"_raw", "_blv_tmp", "_crj_tmp", "_wav"};
	char *osuffix[] = {"_flt", "_flt",     "_crj",     "_fwv"};
	int nsuffix = 4;

	/* A structure to pass the calibration switches to CalStis1 */
	cs1_switch cs1_sw;

	/* reference file keywords and names */
	RefFileInfo refnames;

	c_irafinit (argc, argv);

    PtrRegister ptrReg;
    initPtrRegister(&ptrReg);
	/* Allocate space for file names. */
	inlist = calloc (1, sizeof (char));	/* allocated later */
	addPtr(&ptrReg, inlist, &free);
	outlist = calloc (1, sizeof (char));
	addPtr(&ptrReg, outlist, &free);
	blevlist = calloc (1, sizeof (char));
	addPtr(&ptrReg, blevlist, &free);
	input = calloc (STIS_LINE+1, sizeof (char));
	addPtr(&ptrReg, input, &free);
	output = calloc (STIS_LINE+1, sizeof (char));
	addPtr(&ptrReg, output, &free);
	outblev = calloc (STIS_LINE+1, sizeof (char));
	addPtr(&ptrReg, outblev, &free);
	if (!inlist || !outlist || !blevlist || !input || !output || !outblev) {
	    printf ("ERROR:  Can't even begin; out of memory.\n");
	    freeOnExit(&ptrReg);
	    exit (ERROR_RETURN);
	}

	/* Initialize the lists of reference file keywords and names. */
	InitRefFile (&refnames);
	addPtr(&ptrReg, &refnames, &FreeRefFile);

	/* Initial values. */
	cs1_sw.dqicorr = OMIT;
	cs1_sw.atodcorr = OMIT;
	cs1_sw.blevcorr = OMIT;
	cs1_sw.doppcorr = OMIT;
	cs1_sw.lorscorr = OMIT;
	cs1_sw.glincorr = OMIT;
	cs1_sw.lflgcorr = OMIT;
	cs1_sw.biascorr = OMIT;
	cs1_sw.darkcorr = OMIT;
	cs1_sw.flatcorr = OMIT;
	cs1_sw.shadcorr = OMIT;
	cs1_sw.photcorr = OMIT;
	cs1_sw.statcorr = OMIT;

	strcpy (cs1_sw.darkscale_string, "");

	for (i = 1;  i < argc;  i++) {

	    if (strcmp (argv[i], "--version") == 0) {
		PrVersion();
		freeOnExit(&ptrReg);
		exit (0);
	    }
        if (!(strcmp(argv[i],"--gitinfo")))
        {
            printGitInfo();
            freeOnExit(&ptrReg);
            exit(0);
        }
        if (!(strcmp(argv[i],"--help")))
        {
            printHelp();
            freeOnExit(&ptrReg);
            exit(0);
        }
	    if (strcmp (argv[i], "-r") == 0) {
		PrFullVersion();
		freeOnExit(&ptrReg);
		exit (0);
	    }
	    if (strcmp (argv[i], "-dqi") == 0) {	/* turn on */
		cs1_sw.dqicorr = PERFORM;
		switch_on = 1;
	    } else if (strcmp (argv[i], "-atod") == 0) {
		cs1_sw.atodcorr = PERFORM;
		switch_on = 1;
	    } else if (strcmp (argv[i], "-blev") == 0) {
		cs1_sw.blevcorr = PERFORM;
		switch_on = 1;

	    } else if (strcmp (argv[i], "-dopp") == 0) {
		cs1_sw.doppcorr = PERFORM;
		switch_on = 1;
	    } else if (strcmp (argv[i], "-lors") == 0) {
		cs1_sw.lorscorr = PERFORM;
		switch_on = 1;
	    } else if (strcmp (argv[i], "-glin") == 0) {
		cs1_sw.glincorr = PERFORM;
		switch_on = 1;
	    } else if (strcmp (argv[i], "-lflg") == 0) {
		cs1_sw.lflgcorr = PERFORM;
		switch_on = 1;
	    } else if (strcmp (argv[i], "-bias") == 0) {
		cs1_sw.biascorr = PERFORM;
		switch_on = 1;
	    } else if (strcmp (argv[i], "-dark") == 0) {
		cs1_sw.darkcorr = PERFORM;
		switch_on = 1;
	    } else if (strcmp (argv[i], "-flat") == 0) {
		cs1_sw.flatcorr = PERFORM;
		switch_on = 1;
	    } else if (strcmp (argv[i], "-shad") == 0) {
		cs1_sw.shadcorr = PERFORM;
		switch_on = 1;
	    } else if (strcmp (argv[i], "-phot") == 0) {
		cs1_sw.photcorr = PERFORM;
		switch_on = 1;
	    } else if (strcmp (argv[i], "-stat") == 0) {
		cs1_sw.statcorr = PERFORM;
		switch_on = 1;
	    } else if (strcmp (argv[i], "-dscl") == 0) {
		strcpy (cs1_sw.darkscale_string, argv[++i]);
		switch_on = 1;
	    } else if (argv[i][0] == '-') {
		for (j = 1;  argv[i][j] != '\0';  j++) {
		    if (argv[i][j] == 't') {
			printtime = 1;
		    } else if (argv[i][j] == 'v') {
			verbose = 1;
		    } else {
			printf ("ERROR:  Unrecognized option %s\n", argv[i]);
			printSyntax();
			freeOnExit(&ptrReg);
			exit (1);
		    }
		}
	    } else if (inlist[0] == '\0') {
	    freePtr(&ptrReg, inlist);
		if ((inlist = calloc (strlen(argv[i])+1, sizeof(char)))
			== NULL) {
		    printf ("ERROR:  Out of memory.\n");
		    freeOnExit(&ptrReg);
		    exit (ERROR_RETURN);
		}
		strcpy (inlist, argv[i]);
	    } else if (outlist[0] == '\0') {
	    freePtr(&ptrReg, outlist);
		if ((outlist = calloc (strlen(argv[i])+1, sizeof(char)))
			== NULL) {
		    printf ("ERROR:  Out of memory.\n");
		    freeOnExit(&ptrReg);
		    exit (ERROR_RETURN);
		}
		strcpy (outlist, argv[i]);
	    } else if (blevlist[0] == '\0') {
	    freePtr(&ptrReg, blevlist);
		if ((blevlist = calloc (strlen(argv[i])+1, sizeof(char)))
			== NULL) {
		    printf ("ERROR:  Out of memory.\n");
		    freeOnExit(&ptrReg);
		    exit (ERROR_RETURN);
		}
		strcpy (blevlist, argv[i]);
	    } else {
		too_many = 1;
	    }
	}
	if (inlist[0] == '\0' || too_many) {
	    printSyntax();
	    freeOnExit(&ptrReg);
	    exit (ERROR_RETURN);
	}

	/* Was no calibration switch specified on command line? */
	if (!switch_on) {	/* default values (mostly PERFORM) */
	    cs1_sw.dqicorr  = DefSwitch ("dqicorr");
	    cs1_sw.atodcorr = DefSwitch ("atodcorr");
	    cs1_sw.blevcorr = DefSwitch ("blevcorr");
	    cs1_sw.doppcorr = DefSwitch ("doppcorr");
	    cs1_sw.lorscorr = DefSwitch ("lorscorr");
	    cs1_sw.glincorr = DefSwitch ("glincorr");
	    cs1_sw.lflgcorr = DefSwitch ("lflgcorr");
	    cs1_sw.biascorr = DefSwitch ("biascorr");
	    cs1_sw.darkcorr = DefSwitch ("darkcorr");
	    cs1_sw.flatcorr = DefSwitch ("flatcorr");
	    cs1_sw.shadcorr = DefSwitch ("shadcorr");
	    cs1_sw.photcorr = DefSwitch ("photcorr");
	    cs1_sw.statcorr = DefSwitch ("statcorr");
	}

	/* Expand the templates. */
	i_imt = c_imtopen (inlist);
	addPtr(&ptrReg, i_imt, &c_imtclose);
	o_imt = c_imtopen (outlist);
    addPtr(&ptrReg, o_imt, &c_imtclose);
	b_imt = c_imtopen (blevlist);
    addPtr(&ptrReg, b_imt, &c_imtclose);
	n_in = c_imtlen (i_imt);
	n_out = c_imtlen (o_imt);
	n_blev = c_imtlen (b_imt);

	/* The number of input and output files must be the same. */
	status = 0;
	if (CompareNumbers (n_in, n_out, "output"))
	    status = ERROR_RETURN;
	if (CompareNumbers (n_in, n_blev, "outblev"))
	    status = ERROR_RETURN;
	if (status) {
	    freeOnExit(&ptrReg);
	    exit (ERROR_RETURN);
	}

	/* Loop over the list of input files. */
	for (n = 0;  n < n_in;  n++) {

	    junk = c_imtgetim (i_imt, input, STIS_LINE);
	    if (n_out > 0)
		junk = c_imtgetim (o_imt, output, STIS_LINE);
	    else
		output[0] = '\0';
	    if (n_blev > 0)
		junk = c_imtgetim (b_imt, outblev, STIS_LINE);
	    else
		outblev[0] = '\0';

	    status = 0;
	    if ((status = MkOutName (input, isuffix, osuffix, nsuffix,
                                     output, STIS_LINE))) {
		WhichError (status);
		printf ("Skipping %s\n", input);
		continue;
	    }

	    /* Calibrate the current input file. */
	    if ((status = CalStis1 (input, output, outblev,
                                    &cs1_sw, &refnames, printtime, verbose))) {
		printf ("Error processing %s.\n", input);
		WhichError (status);
	    }
	}

    freeOnExit(&ptrReg);

	if (status)
	    exit (ERROR_RETURN);
	else
	    exit (0);
}

/* This function checks that the number of input and output files are
   the same, unless no output file was specified.  If the numbers match,
   zero will be returned; if they aren't the same, one will be returned.
   That is, this function returns true for the error condition.
*/

static int CompareNumbers (int n_in, int n_out, char *str_out) {

	if (n_out > 0 && n_out != n_in) {
	    printf ("You specified %d input file", n_in);
	    if (n_in > 1)
		printf ("s");
	    printf (" but %d %s file", n_out, str_out);
	    if (n_out > 1)
		printf ("s");
	    printf (".\n");
	    return (1);
	}

	return (0);
}
