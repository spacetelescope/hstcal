/* calstis4 -- determine offsets from wavecal */

# include <stdio.h>
# include <stdlib.h>		/* calloc */
# include <string.h>

# include "c_iraf.h"		/* for c_irafinit */
# include "ximio.h"

# include "stis.h"
# include "calstis4.h"
# include "hstcalerr.h"

static int CompareNumbers (int, int);
static void FreeNames (char *, char *, char *, char *);

/* This is the main module for calstis4.  It gets the input file name
   and command-line flags, and then calls CalStis4.

   Phil Hodge, 1997 Dec 11:
	Define RefFileInfo refnames, and pass it (an empty list) to CalStis4.

   Phil Hodge, 1997 Dec 24:
	Initialize and free refnames.

   Phil Hodge, 1998 Jan 27:
	Accept list of input file names, and loop over them.

   Phil Hodge, 1998 Dec 11:
	Get dbgfile, and add it to the calling sequence of CalStis4.

   Phil Hodge, 2004 July 23:
	Get slit angle, and add it to the calling sequence of CalStis4.

   Phil Hodge, 2011 July 6:
	Include command-line option '--version'.

   Phil Hodge, 2012 Feb 10:
	Include command-line option '-r'.
*/

int main (int argc, char **argv) {

	int status;		/* zero is OK */

	char *inlist;		/* list of input file names */
	char *dbglist;		/* list of files for debug output */
	int printtime = 0;	/* print time after each step? */
	int verbose = 0;	/* print additional info? */
	int too_many = 0;	/* too many command-line arguments? */
	int dbg_next = 0;	/* debug file name must be next argument? */
	int angle_next = 0;	/* slit angle must be next argument? */
	double slit_angle = 0.;	/* (degrees) angle of long slit for echelle */
	int i, j;		/* loop indexes */

	IRAFPointer i_imt;	/* imt list pointer for input files */
	IRAFPointer d_imt;	/* imt list pointer for output debug files */
	char *input;		/* name of input science file */
	char *dbgfile;		/* name of file for debug output */
	int n_in, n_dbg;	/* number of names in lists */
	int n;
	int junk;

	/* reference file keywords and names */
	RefFileInfo refnames;

	c_irafinit (argc, argv);

	inlist = calloc (STIS_LINE+1, sizeof (char));
	input = calloc (STIS_LINE+1, sizeof (char));
	dbglist = calloc (STIS_LINE+1, sizeof (char));
	dbgfile = calloc (STIS_LINE+1, sizeof (char));
	if (inlist == NULL || input == NULL ||
	    dbglist == NULL || dbgfile == NULL) {
	    printf ("ERROR:  Can't even begin; out of memory.\n");
	    exit (ERROR_RETURN);
	}

	/* Get command-line arguments. */
	for (i = 1;  i < argc;  i++) {
	    if (dbg_next) {
		if (argv[i][0] == '-') {
		    printf (
	"`%s' encountered when debug file name was expected\n", argv[i]);
		    FreeNames (inlist, input, dbglist, dbgfile);
		    exit (ERROR_RETURN);
		}
		strcpy (dbglist, argv[i]);
		dbg_next = 0;
	    } else if (angle_next) {
		slit_angle = atof (argv[i]);
		angle_next = 0;
	    } else if (argv[i][0] == '-') {
		if (strcmp (argv[i], "--version") == 0) {
		    PrVersion();
		    exit (0);
		}
		if (strcmp (argv[i], "-r") == 0) {
		    PrFullVersion();
		    exit (0);
		}
		if (strcmp (argv[i]+1, "angle") == 0) {
		    /* next argument must be slit angle */
		    angle_next = 1;
		} else {
		    for (j = 1;  argv[i][j] != '\0';  j++) {
			if (argv[i][j] == 't') {
			    printtime = 1;
			} else if (argv[i][j] == 'v') {
			    verbose = 1;
			} else if (argv[i][j] == 'd') {
			    /* next argument must be debug file name */
			    dbg_next = 1;
			} else {
			    printf ("ERROR:  Unrecognized option %s\n",
				argv[i]);
			    exit (1);
			}
		    }
		}
	    } else if (inlist[0] == '\0') {
		strcpy (inlist, argv[i]);
	    } else {
		too_many = 1;
	    }
	}
	if (inlist[0] == '\0' || too_many || dbg_next) {
	    printf (
	"syntax:  cs4.e [-t] [-v] input [-angle slit_angle] [-d debugfile]\n");
	    FreeNames (inlist, input, dbglist, dbgfile);
	    exit (ERROR_RETURN);
	}

	/* Initialize the list of reference file keywords and names. */
	InitRefFile (&refnames);

	/* Expand the templates. */
	i_imt = c_imtopen (inlist);
	d_imt = c_imtopen (dbglist);
	n_in = c_imtlen (i_imt);
	n_dbg = c_imtlen (d_imt);

	/* The number of input and debug files must be the same. */
	if (CompareNumbers (n_in, n_dbg)) {
	    FreeNames (inlist, input, dbglist, dbgfile);
	    exit (ERROR_RETURN);
	}

	/* Loop over the list of input files. */
	for (n = 0;  n < n_in;  n++) {

	    junk = c_imtgetim (i_imt, input, STIS_LINE);
	    if (n_dbg > 0)
		junk = c_imtgetim (d_imt, dbgfile, STIS_LINE);
	    else
		dbgfile[0] = '\0';

	    /* Calibrate the current input file. */
	    status = 0;
	    if ((status = CalStis4 (input, dbgfile, &refnames,
                                    printtime, verbose, slit_angle))) {
		printf ("Error processing %s.\n", input);
		WhichError (status);
	    }
	}

	/* Close lists of file names, and free name buffers. */
	c_imtclose (d_imt);
	c_imtclose (i_imt);
	FreeRefFile (&refnames);
	FreeNames (inlist, input, dbglist, dbgfile);

	if (status)
	    exit (ERROR_RETURN);
	else
	    exit (0);
}

/* This function checks that the number of input and debug files are
   the same, unless no debug file was specified.  If the numbers match,
   zero will be returned; if they aren't the same, one will be returned.
   That is, this function returns true for the error condition.
*/

static int CompareNumbers (int n_in, int n_dbg) {

	if (n_dbg > 0 && n_dbg != n_in) {
	    printf ("You specified %d input file", n_in);
	    if (n_in > 1)
		printf ("s");
	    printf (" but %d debug file", n_dbg);
	    if (n_dbg > 1)
		printf ("s");
	    printf (".\n");
	    return (1);
	}

	return (0);
}

static void FreeNames (char *inlist, char *input,
		char *dbglist, char *dbgfile) {

	free (dbgfile);
	free (dbglist);
	free (input);
	free (inlist);
}
