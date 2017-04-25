/* calstis11 -- subtract science from wavecal */

# include <stdio.h>
# include <stdlib.h>		/* calloc */
# include <string.h>

# include "c_iraf.h"		/* for c_irafinit */
# include "ximio.h"

# include "stis.h"
# include "calstis11.h"
# include "hstcalerr.h"

static int CompareNumbers (int, char *, int, char *);
static void FreeNames (char *, char *, char *, char *, char *, char *);

/* This is the main module for calstis11.  It gets the file names,
   calibration switches, and flags, and then calls CalStis11.

   Phil Hodge, 1998 Jan 27:
	Accept lists of input and output file names, and loop over
	individual file names.

   Phil Hodge, 1998 Feb 5:
	Call MkOutName instead of MkName.

   Phil Hodge, 1998 Oct 5:
	Change status value 1 to ERROR_RETURN.

   Phil Hodge, 2011 July 6:
	Include command-line option '--version'.

   Phil Hodge, 2012 Feb 10:
	Include command-line option '-r'.
*/

int main (int argc, char **argv) {

	int status;		/* zero is OK */
	char *wavlist;		/* list of input wavecal file names */
	char *scilist;		/* list of input science file names */
	char *outlist;		/* list of output file names */
	int printtime = 0;	/* print time after each step? */
	int verbose = 0;	/* print additional info? */
	int too_many = 0;	/* too many command-line arguments? */
	int i, j;		/* loop indexes */

	IRAFPointer w_imt, s_imt, o_imt;	/* imt list pointers */
	char *inwav, *insci;	/* input wavecal and science file names */
	char *output;		/* output file name */
	int n_wav, n_sci, n_out;	/* number of files in each list */
	int n;

	/* Input and output suffixes. */
	char *isuffix[] = {"_wav", "_fwv_tmp"};
	char *osuffix[] = {"_cwv", "_cwv_tmp"};
	int nsuffix = 2;

	c_irafinit (argc, argv);

	/* Allocate space for file names. */
	wavlist = calloc (STIS_LINE+1, sizeof (char));
	scilist = calloc (STIS_LINE+1, sizeof (char));
	outlist = calloc (STIS_LINE+1, sizeof (char));
	inwav = calloc (STIS_LINE+1, sizeof (char));
	insci = calloc (STIS_LINE+1, sizeof (char));
	output = calloc (STIS_LINE+1, sizeof (char));
	if (wavlist == NULL || scilist == NULL || outlist == NULL ||
		inwav == NULL || insci == NULL || output == NULL) {
	    printf ("ERROR:  Can't even begin:  out of memory.\n");
	    exit (ERROR_RETURN);
	}

	/* Get names of input and output files. */
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
	    } else if (outlist[0] == '\0') {
		strcpy (outlist, argv[i]);
	    } else {
		too_many = 1;
	    }
	}
	if (scilist[0] == '\0' || too_many) {
	    printf (
	"ERROR:  syntax:  cs11.e [-t] [-v] wavecal science output\n");
	    FreeNames (wavlist, scilist, outlist, inwav, insci, output);
	    exit (ERROR_RETURN);
	}

	/* Expand the templates. */
	w_imt = c_imtopen (wavlist);
	s_imt = c_imtopen (scilist);
	o_imt = c_imtopen (outlist);
	n_wav = c_imtlen (w_imt);
	n_sci = c_imtlen (s_imt);
	n_out = c_imtlen (o_imt);

	/* The number of wavecal and science files must be the same.
	   The output may have been omitted, but if any output was
	   specified, the number must match the number of wavecals.
	*/
	status = 0;
	if (CompareNumbers (n_wav, "wavecal", n_sci, "science"))
	    status = ERROR_RETURN;
	if (n_out > 0 && CompareNumbers (n_wav, "wavecal", n_out, "output"))
	    status = ERROR_RETURN;
	if (status) {
	    FreeNames (wavlist, scilist, outlist, inwav, insci, output);
	    exit (ERROR_RETURN);
	}

	/* Loop over the list of input files. */
	for (n = 0;  n < n_wav;  n++) {

	    j = c_imtgetim (w_imt, inwav, STIS_LINE);
	    j = c_imtgetim (s_imt, insci, STIS_LINE);
	    if (n_out > 0)
		j = c_imtgetim (o_imt, output, STIS_LINE);
	    else
		output[0] = '\0';

	    status = 0;
	    if ((status = MkOutName (inwav, isuffix, osuffix, nsuffix,
                                     output, STIS_LINE))) {
		WhichError (status);
		printf ("Skipping %s\n", inwav);
		continue;
	    }

	    if ((status = CalStis11 (inwav, insci, output, printtime, verbose))) {
		printf ("Error processing %s.\n", inwav);
		WhichError (status);
	    }
	}

	/* Close lists of file names, and free name buffers. */
	c_imtclose (w_imt);
	c_imtclose (s_imt);
	c_imtclose (o_imt);
	FreeNames (wavlist, scilist, outlist, inwav, insci, output);

	if (status)
	    exit (ERROR_RETURN);
	else
	    exit (0);
}

/* This function checks that the number of file names in the two lists
   are the same.  Note that this version requires an exact match; it does
   not allow n_out to be zero.  If the numbers match, zero will be
   returned; if they aren't the same, one will be returned.  That is,
   this function returns true for the error condition.
*/

static int CompareNumbers (int n_in, char *str_in,
		int n_out, char *str_out) {

	if (n_out != n_in) {
	    printf ("You specified %d %s file", n_in, str_in);
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

static void FreeNames (char *wavlist, char *scilist, char *outlist,
		char *inwav, char *insci, char *output) {

	free (output);
	free (insci);
	free (inwav);
	free (outlist);
	free (scilist);
	free (wavlist);
}
