/* calstis12 -- update SHIFTAi keywords in science file extensions */

# include <stdio.h>
# include <stdlib.h>		/* calloc */
# include <string.h>
# include <ctype.h>		/* isupper, tolower */

static int WavOption (char *, int *);

# include "c_iraf.h"		/* for c_irafinit */
# include "ximio.h"

# include "stis.h"
# include "calstis12.h"
# include "cs12.h"		/* for interpolation options */
# include "err.h"

static int CompareNumbers (int, int);
static void FreeNames (char *, char *, char *, char *, char *);

/* This is the main module for calstis12.  It gets the input file names
   and command-line flags, and then calls CalStis12.

   Phil Hodge, 1997 Dec 24:
	Define RefFileInfo refnames, and pass it (an empty list) to CalStis12.

   Phil Hodge, 1998 Jan 27:
	Accept lists of input and output file names, and loop over
	individual file names.

   Phil Hodge, 1998 Mar 18:
	Remove RefFileInfo refnames from calling sequence of CalStis12.

   Phil Hodge, 1998 Oct 5:
	Change status value 122 to ERROR_RETURN.

   Phil Hodge, 2000 Jan 19:
	Add linear interpolation option, and make that the default.

   Phil Hodge, 2000 June 16:
	Update a comment.

   Phil Hodge, 2011 July 6:
	Include command-line option '--version'.

   Phil Hodge, 2012 Feb 10:
	Include command-line option '-r'.
*/

int main (int argc, char **argv) {

	int status;		/* zero is OK */
	char *wavlist;		/* list of input wavecal file names */
	char *scilist;		/* list of input science file names */
	char *which_wavecal;	/* option for selecting wavecal */
	int printtime = 0;	/* print time after each step? */
	int verbose = 0;	/* print additional info? */
	int too_many = 0;	/* too many command-line arguments? */
	int w_option;		/* which option for selecting wavecal */
	int i, j;		/* loop indexes */

	IRAFPointer w_imt, s_imt;	/* imt list pointers */
	char *inwav, *insci;	/* wavecal and science file names */
	int n_wav, n_sci;	/* number of files in each list */
	int n;

	c_irafinit (argc, argv);

	/* Allocate space for file names. */
	wavlist = calloc (STIS_LINE+1, sizeof (char));
	scilist = calloc (STIS_LINE+1, sizeof (char));
	which_wavecal = calloc (STIS_LINE+1, sizeof (char));
	inwav = calloc (STIS_LINE+1, sizeof (char));
	insci = calloc (STIS_LINE+1, sizeof (char));
	if (wavlist == NULL || scilist == NULL || which_wavecal == NULL ||
		inwav == NULL || insci == NULL) {
	    printf ("ERROR:  Can't even begin; out of memory.\n");
	    exit (ERROR_RETURN);
	}

	/* Get names of input files. */
	for (i = 1;  i < argc;  i++) {
	    if (strcmp (argv[i], "--version") == 0) {
		PrVersion();
		exit (0);
	    }
	    if (strcmp (argv[i], "-r") == 0) {
		PrFullVersion();
		exit (0);
	    }
	    if (argv[i][0] == '-') {
		for (j = 1;  argv[i][j] != '\0';  j++) {
		    if (argv[i][j] == 't') {
			printtime = 1;
		    } else if (argv[i][j] == 'v') {
			verbose = 1;
		    } else {
			printf ("ERROR:  Unrecognized option %s\n", argv[i]);
			exit (1);
		    }
		}
	    } else if (wavlist[0] == '\0') {
		strcpy (wavlist, argv[i]);
	    } else if (scilist[0] == '\0') {
		strcpy (scilist, argv[i]);
	    } else if (which_wavecal[0] == '\0') {
		strcpy (which_wavecal, argv[i]);
	    } else {
		too_many = 1;
	    }
	}
	if (scilist[0] == '\0' || too_many) {
	    printf ("ERROR:  syntax:  cs12.e [-t] [-v] wavecal science\n");
	    exit (ERROR_RETURN);
	}
	if (which_wavecal[0] != '\0') {
	    if ((status = WavOption (which_wavecal, &w_option)))
		exit (status);
	} else {
	    w_option = STIS_LINEAR;
	}

	/* Expand the templates. */
	w_imt = c_imtopen (wavlist);
	s_imt = c_imtopen (scilist);
	n_wav = c_imtlen (w_imt);
	n_sci = c_imtlen (s_imt);

	/* The number of wavecal and science files must be the same. */
	if (CompareNumbers (n_wav, n_sci)) {
	    FreeNames (wavlist, scilist, which_wavecal, inwav, insci);
	    exit (ERROR_RETURN);
	}

	/* Loop over the list of input files. */
	for (n = 0;  n < n_wav;  n++) {

	    j = c_imtgetim (w_imt, inwav, STIS_LINE);
	    j = c_imtgetim (s_imt, insci, STIS_LINE);

	    status = 0;
	    if ((status = CalStis12 (inwav, insci,
                                     w_option, printtime, verbose))) {
		printf ("Error processing %s.\n", insci);
		WhichError (status);
	    }
	}

	/* Close lists of file names, and free name buffers. */
	c_imtclose (w_imt);
	c_imtclose (s_imt);
	FreeNames (wavlist, scilist, which_wavecal, inwav, insci);

	if (status)
	    exit (ERROR_RETURN);
	else
	    exit (0);
}

static int WavOption (char *which_wavecal, int *w_option) {

	char *which_lc;		/* which_wavecal, lower case */
	int len;		/* length of which_wavecal string */
	int i, c;

	len = strlen (which_wavecal);
	if (len < 1)
	    return (0);

	if ((which_lc = calloc (len+1, sizeof (char))) == NULL)
	    return (OUT_OF_MEMORY);

	i = 0;
	while ((c = which_wavecal[i]) != '\0') {
	    if (isupper (c))
		which_lc[i] = tolower (c);
	    else
		which_lc[i] = c;
	    i++;
	}
	which_lc[i] = '\0';

	if (strncmp (which_lc, "nearest", len) == 0) {
	    *w_option = STIS_NEAREST;
	} else if (strncmp (which_lc, "linear", len) == 0) {
	    *w_option = STIS_LINEAR;
	} else {
	    printf ("ERROR:  Don't understand wavecal option = `%s'.\n",
		which_wavecal);
	    return (ERROR_RETURN);
	}

	free (which_lc);

	return (0);
}

/* This function checks that the number of file names in the two lists
   are the same.  Note that this version requires an exact match; it does
   not allow n_out to be zero.  If the numbers match, zero will be
   returned; if they aren't the same, one will be returned.  That is,
   this function returns true for the error condition.
*/

static int CompareNumbers (int n_in, int n_out) {

	if (n_out != n_in) {
	    printf ("You specified %d wavecal file", n_in);
	    if (n_in > 1)
		printf ("s");
	    printf (" but %d science file", n_out);
	    if (n_out > 1)
		printf ("s");
	    printf (".\n");
	    return (1);
	}

	return (0);
}

static void FreeNames (char *wavlist, char *scilist, char *which_wavecal,
		char *inwav, char *insci) {

	free (insci);
	free (inwav);
	free (which_wavecal);
	free (scilist);
	free (wavlist);
}
