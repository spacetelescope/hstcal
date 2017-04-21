/* calstis7 -- 2-D spectral rectification */

# include <stdio.h>
# include <stdlib.h>		/* calloc, atof */
# include <string.h>


# include "c_iraf.h"		/* for c_irafinit */
# include "ximio.h"

# include "stis.h"
# include "calstis7.h"
# include "hstcalerr.h"

static int CompareNumbers (int, int, char *);
static void FreeNames (char *, char *, char *, char *);

/* This is the main module for calstis7.  It gets the input and output
   file names, calibration switches, and flags, and then calls CalStis7.

   Phil Hodge, 1997 Dec 11:
	Define RefFileInfo refnames, and pass it (an empty list) to CalStis7.

   Phil Hodge, 1997 Dec 24:
	Initialize and free refnames.

   Phil Hodge, 1998 Jan 27:
	Accept lists of input and output file names, and loop over
	individual file names.

   Phil Hodge, 1998 Feb 5:
	Call MkOutName instead of MkName.

   Phil Hodge, 2000 Mar 11:
	Add _fwv --> _w2d and _cwv --> _w2d to isuffix and osuffix.

   Phil Hodge, 2000 Apr 5:
	Get -c (center_target) from command line; add center_target
	to the calling sequence of CalStis7.

   Ivo Busko, 2002 Apr 24:
   	Add -b command line argument.

   Phil Hodge, 2006 Feb 7:
	Add -wgt_err command line argument and err_algorithm CalStis7 argument.

   Phil Hodge, 2011 July 6:
	Include command-line option '--version'.

   Phil Hodge, 2012 Feb 10:
	Include command-line option '-r'.
*/

int main (int argc, char **argv) {

	int status;		/* zero is OK */
	char *inlist;		/* list of input file names */
	char *outlist;		/* list of output file names */
	int switch_on = 0;	/* was any switch specified? */
	int sgeocorr = OMIT;	/* calibration switches */
	int helcorr = OMIT;
	int fluxcorr = OMIT;
	int statcorr = OMIT;
	int err_algorithm = WGT_VARIANCE;
	int printtime = 0;	/* print time after each step? */
	int verbose = 0;	/* print additional info? */
	int center_target = 0;	/* center target in output image? */
	int too_many = 0;	/* too many command-line arguments? */
	int i, j;		/* loop indexes */
	double blazeshift = NO_VALUE;

	IRAFPointer i_imt, o_imt;	/* imt list pointers */
	char *input;		/* name of input science file */
	char *output;		/* optional name of output file */
	int n_in, n_out;	/* number of files in each list */
	int n;

	/* Input and output suffixes. */
	char *isuffix[] =
		{"_flt", "_crj", "_fwv", "_cwv", "_fwv_tmp", "_cwv_tmp"};
	char *osuffix[] =
		{"_x2d", "_sx2", "_w2d", "_w2d", "_w2d_tmp", "_w2d_tmp"};
	int nsuffix = 6;

	/* reference file keywords and names */
	RefFileInfo refnames;

	c_irafinit (argc, argv);

	inlist = calloc (STIS_LINE+1, sizeof (char));
	outlist = calloc (STIS_LINE+1, sizeof (char));
	input = calloc (STIS_LINE+1, sizeof (char));
	output = calloc (STIS_LINE+1, sizeof (char));
	if (inlist == NULL || outlist == NULL ||
		input == NULL || output == NULL) {
	    printf ("ERROR:  Can't even begin:  out of memory.\n");
	    exit (ERROR_RETURN);
	}

	/* Get command-line arguments. */
	for (i = 1;  i < argc;  i++) {
	    if (strcmp (argv[i], "-x2d") == 0 ||
		strcmp (argv[i], "-geo") == 0) {
		switch_on = 1;
	    } else if (strcmp (argv[i], "-sgeo") == 0) {	/* turn on */
		sgeocorr = PERFORM;
		switch_on = 1;
	    } else if (strcmp (argv[i], "-hel") == 0) {
		helcorr = PERFORM;
		switch_on = 1;
	    } else if (strcmp (argv[i], "-flux") == 0) {
		fluxcorr = PERFORM;
		switch_on = 1;
	    } else if (strcmp (argv[i], "-stat") == 0) {
		statcorr = PERFORM;
		switch_on = 1;
	    } else if (strcmp (argv[i], "-wgt_err") == 0) {
		err_algorithm = WGT_ERROR;
	    } else if (argv[i][0] == '-') {
		if (strcmp (argv[i], "--version") == 0) {
		    PrVersion();
		    exit (0);
		}
		if (strcmp (argv[i], "-r") == 0) {
		    PrFullVersion();
		    exit (0);
		}
		for (j = 1;  argv[i][j] != '\0';  j++) {
		    if (argv[i][j] == 't') {
			printtime = 1;
		    } else if (argv[i][j] == 'v') {
			verbose = 1;
		    } else if (argv[i][j] == 'c') {
			center_target = 1;	/* yes, center target */
		    } else if (argv[i][j] == 'b') {
	                blazeshift = (double) atof (argv[++i]);
	                if (i == argc-1)
	                    break;
		    } else {
			printf ("ERROR:  Unrecognized option %s\n", argv[i]);
			exit (1);
		    }
		}
	    } else if (inlist[0] == '\0') {
		strcpy (inlist, argv[i]);
	    } else if (outlist[0] == '\0') {
		strcpy (outlist, argv[i]);
	    } else {
		too_many = 1;
	    }
	}
	if (inlist[0] == '\0' || too_many) {
	    printf (
"syntax:  cs7.e [-t] [-v] [-c] [-wgt_err] [-b blazeshift] input output\n");
	    printf ("  command-line switches:  -x2d -sgeo -hel -flux -stat\n");
	    FreeNames (inlist, outlist, input, output);
	    exit (ERROR_RETURN);
	}

	/* Was no calibration switch specified on command line? */
	if (!switch_on) {	/* default values (mostly PERFORM) */
	    sgeocorr = DefSwitch ("sgeocorr");
	    helcorr  = DefSwitch ("helcorr");
	    fluxcorr = DefSwitch ("fluxcorr");
	    statcorr = DefSwitch ("statcorr");
	}

	/* Initialize the list of reference file keywords and names. */
	InitRefFile (&refnames);

	/* Expand the templates. */
	i_imt = c_imtopen (inlist);
	o_imt = c_imtopen (outlist);
	n_in = c_imtlen (i_imt);
	n_out = c_imtlen (o_imt);

	/* The number of input and output files must be the same. */
	if (CompareNumbers (n_in, n_out, "output")) {
	    FreeNames (inlist, outlist, input, output);
	    exit (ERROR_RETURN);
	}

	/* Loop over the list of input files. */
	for (n = 0;  n < n_in;  n++) {

	    j = c_imtgetim (i_imt, input, STIS_LINE);
	    if (n_out > 0)
		j = c_imtgetim (o_imt, output, STIS_LINE);
	    else
		output[0] = '\0';

	    status = 0;
	    if ((status = MkOutName (input, isuffix, osuffix, nsuffix,
                                     output, STIS_LINE))) {
		WhichError (status);
		printf ("Skipping %s\n", input);
		continue;
	    }

	    /* Calibrate the current input file. */
	    if ((status = CalStis7 (input, output,
			sgeocorr, helcorr, fluxcorr, statcorr,
			&refnames, printtime, verbose, center_target,
                                    blazeshift, err_algorithm))) {
		printf ("Error processing %s.\n", input);
		WhichError (status);
	    }
	}

	/* Close lists of file names, and free name buffers. */
	c_imtclose (i_imt);
	c_imtclose (o_imt);
	FreeRefFile (&refnames);
	FreeNames (inlist, outlist, input, output);

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

static void FreeNames (char *inlist, char *outlist,
		char *input, char *output) {

	free (output);
	free (input);
	free (outlist);
	free (inlist);
}
