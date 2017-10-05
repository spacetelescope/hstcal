/* calstis0 -- integrated calstis processing */

# include <stdio.h>
# include <stdlib.h>		/* calloc */
# include <string.h>

#include "hstcal_memory.h"
# include "c_iraf.h"		/* for c_irafinit */
# include "ximio.h"

# include "stis.h"
# include "calstis0.h"
# include "hstcalerr.h"
# include "hstcalversion.h"

static void printSyntax(void)
{
    printf("syntax:  cs0.e [--help] [-t] [-s] [-v] [--version] [--gitinfo] input [outroot] [-w wavecal]\n");
}
static void printHelp(void)
{
    printSyntax();
}

/* This is the main module for calstis0.  It gets the input and output
   file names, calibration switches, and flags, and then calls CalStis0.

   Phil Hodge, 1998 Jan 21:
	Accept lists of input and outroot file names, and loop over
	individual file names.

   Phil Hodge, 2011 July 6:
	Include command-line option '--version'.

   Phil Hodge, 2011 July 26:
	Add a check that n_raw > 0 after calling c_imtopen.

   Phil Hodge, 2012 Jan 23:
	Allocate memory for rawlist, wavlist, outlist based on the size
	of the string on the command line, rather than using a fixed size.

   Phil Hodge, 2012 Feb 10:
	Include command-line option '-r'.
*/

int main (int argc, char **argv) {

	int status;			/* zero is OK */

	char *rawlist;		/* list of input file names */
	char *wavlist;		/* list of input wavecal file names */
	char *outlist;		/* list of output root names */
	int printtime = 0;	/* print time after each step? */
	int save_tmp = 0;	/* save temporary files? */
	int verbose = 0;	/* print info in individual calstisN? */
	int too_many = 0;	/* too many command-line arguments? */
	int wavecal_next = 0;	/* next argument is wavecal name? */
	int i, j;		/* loop indexes */

	IRAFPointer r_imt, w_imt, o_imt;	/* imt list pointers */
	char *rawfile;		/* name of input science file */
	char *wavfile;		/* optional name of input wavecal */
	char *outroot;		/* optional output root name */
	int many_wavfiles;	/* do we have more than one wavfile? */
	int many_outroots;	/* do we have more than one outroot? */
	int n_raw, n_wav, n_out;	/* number of files in each list */
	int n;

	c_irafinit (argc, argv);

    PtrRegister ptrReg;
    initPtrRegister(&ptrReg);
	/* Allocate space for file names. */
	rawlist = calloc (1, sizeof (char));	/* allocated later */
	addPtr(&ptrReg, rawlist, &free);
	wavlist = calloc (1, sizeof (char));
    addPtr(&ptrReg, wavlist, &free);
	outlist = calloc (1, sizeof (char));
    addPtr(&ptrReg, outlist, &free);
	rawfile = calloc (STIS_LINE+1, sizeof (char));
    addPtr(&ptrReg, rawfile, &free);
	wavfile = calloc (STIS_LINE+1, sizeof (char));
    addPtr(&ptrReg, wavfile, &free);
	outroot = calloc (STIS_LINE+1, sizeof (char));
    addPtr(&ptrReg, outroot, &free);

	if (!rawlist || !wavlist || !outlist || !rawfile || !wavfile || !outroot) {
	    printf ("Can't even begin; out of memory.\n");
	    freeOnExit(&ptrReg);
	    exit (ERROR_RETURN);
	}

	for (i = 1;  i < argc;  i++) {
	    if (wavecal_next) {
		if (argv[i][0] == '-') {
		    printf (
	"`%s' encountered when wavecal file name was expected\n", argv[i]);
		    freeOnExit(&ptrReg);
		    exit (ERROR_RETURN);
		}
		freePtr(&ptrReg, wavlist);
		if ((wavlist = calloc (strlen(argv[i])+1, sizeof(char)))
			== NULL) {
		    printf ("ERROR:  Out of memory.\n");
		    freeOnExit(&ptrReg);
		    exit (ERROR_RETURN);
		}
		addPtr(&ptrReg, wavlist, &free);
		strcpy (wavlist, argv[i]);
		wavecal_next = 0;
	    } else if (argv[i][0] == '-') {
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
		for (j = 1;  argv[i][j] != '\0';  j++) {
		    if (argv[i][j] == 't') {
			printtime = 1;
		    } else if (argv[i][j] == 's') {
			save_tmp = 1;
		    } else if (argv[i][j] == 'v') {
			verbose = 1;
		    } else if (argv[i][j] == 'w') {
			/* next argument must be wavecal name */
			wavecal_next = 1;
		    } else {
			printf ("Unrecognized option %s\n", argv[i]);
			printSyntax();
			freeOnExit(&ptrReg);
			exit (ERROR_RETURN);
		    }
		}
	    } else if (rawlist[0] == '\0') {
		freePtr(&ptrReg, rawlist);
		if ((rawlist = calloc (strlen(argv[i])+1, sizeof(char)))
			== NULL) {
		    printf ("ERROR:  Out of memory.\n");
		    freeOnExit(&ptrReg);
		    exit (ERROR_RETURN);
		}
		addPtr(&ptrReg, rawlist, &free);
		strcpy (rawlist, argv[i]);
	    } else if (outlist[0] == '\0') {
		freePtr(&ptrReg, outlist);
		if ((outlist = calloc (strlen(argv[i])+1, sizeof(char)))
			== NULL) {
		    printf ("ERROR:  Out of memory.\n");
		    exit (ERROR_RETURN);
		}
		addPtr(&ptrReg, outlist, &free);
		strcpy (outlist, argv[i]);
	    } else {
		too_many = 1;
	    }
	}
	if (rawlist[0] == '\0' || too_many || wavecal_next) {
	    printSyntax();
	    freeOnExit(&ptrReg);
	    exit (ERROR_RETURN);
	}

	/* Expand the templates. */
	r_imt = c_imtopen (rawlist);
	addPtr(&ptrReg, r_imt, &c_imtclose);
	w_imt = c_imtopen (wavlist);
	addPtr(&ptrReg, w_imt, &c_imtclose);
	o_imt = c_imtopen (outlist);
	addPtr(&ptrReg, o_imt, &c_imtclose);
	n_raw = c_imtlen (r_imt);
	n_wav = c_imtlen (w_imt);
	n_out = c_imtlen (o_imt);

	/* If there's only one wavfile name, get it now. */
	if (n_wav > 1) {
	    many_wavfiles = 1;
	    if (n_wav != n_raw) {
		printf ("You specified %d input file", n_raw);
		if (n_raw > 1)
		    printf ("s");
		printf (" but %d wavecals.\n", n_wav);
		freeOnExit(&ptrReg);
		exit (ERROR_RETURN);
	    }
	} else {
	    many_wavfiles = 0;
	    if (n_wav == 1)
		j = c_imtgetim (w_imt, wavfile, STIS_LINE);
	    else
		wavfile[0] = '\0';
	}

	/* If there's only one outroot name, get it now. */
	if (n_out > 1) {
	    many_outroots = 1;
	    if (n_out != n_raw) {
		printf ("You specified %d input file", n_raw);
		if (n_raw > 1)
		    printf ("s");
		printf (" but %d outroots.\n", n_out);
		freeOnExit(&ptrReg);
		exit (ERROR_RETURN);
	    }
	} else {
	    many_outroots = 0;
	    if (n_out == 1) {
		j = c_imtgetim (o_imt, outroot, STIS_LINE);
		/* If we have only one outroot, and there is more than one
			raw files, outroot must be a directory.
		*/
		if (n_raw > 1) {
		    if (outroot[strlen(outroot)-1] != '/') {
			printf (
	"You specified multiple input files and only one outroot;\n");
			printf (
	"to do this outroot must be a directory (i.e. end in '/').\n");
			freeOnExit(&ptrReg);
			exit (ERROR_RETURN);
		    }
		}
	    } else {
		outroot[0] = '\0';
	    }
	}
	if (n_raw < 1) {
	    printf ("File not found.\n");
	    freeOnExit(&ptrReg);
	    exit (ERROR_RETURN);
	}

	/* Loop over the list of input files. */
	for (n = 0;  n < n_raw;  n++) {

	    j = c_imtgetim (r_imt, rawfile, STIS_LINE);
	    if (many_wavfiles)
		j = c_imtgetim (w_imt, wavfile, STIS_LINE);
	    if (many_outroots)
		j = c_imtgetim (o_imt, outroot, STIS_LINE);

	    /* Calibrate the current input file. */
	    status = 0;
	    if ((status = CalStis0 (rawfile, wavfile, outroot,
                                    printtime, save_tmp, verbose))) {
		WhichError (status);
		if ((status == NOTHING_TO_DO))
		    status = 0;			/* not an error for cs0 */
		else
		    printf ("calstis0 failed.\n");
	    }
	}

    freeOnExit(&ptrReg);

	if (status)
	    exit (ERROR_RETURN);
	else
	    exit (0);
}
