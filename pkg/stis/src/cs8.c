/* calstis8 -- Sum repeatobs data */

# include <stdio.h>
# include <stdlib.h>		/* calloc */
# include <string.h>

#include "hstcal_memory.h"
# include "c_iraf.h"		/* for c_irafinit */

# include "../stis.h"
# include "calstis8.h"
# include "hstcalerr.h"
# include "hstcalversion.h"

static void printSyntax(void)
{
    printf("syntax:  cs8.e [--help] [-t] [-v] [--version] [--gitinfo] input output\n");
}
static void printHelp(void)
{
    printSyntax();
}

/* This is the main module for calstis8.  It gets the input and output
   file names, calibration switches, and flags, and then calls CalStis8.

   Phil Hodge, 1998 Apr 10:
	Change _s2d to _sx2 when calling MkName (this was just a typo).

   Phil Hodge, 2011 July 6:
	Include command-line option '--version'.

   Phil Hodge, 2012 Feb 10:
	Include command-line option '-r'.
*/

int main (int argc, char **argv) {

	int status;		/* zero is OK */
	char *input, *output;	/* file names */
	int printtime = 0;	/* print time after each step? */
	int verbose = 0;	/* print additional info? */
	int too_many = 0;	/* too many command-line arguments? */
	int i, j;		/* loop indexes */

	c_irafinit (argc, argv);

    PtrRegister ptrReg;
    initPtrRegister(&ptrReg);
	input = calloc (STIS_LINE+1, sizeof (char));
    addPtr(&ptrReg, input, &free);
	output = calloc (STIS_LINE+1, sizeof (char));
    addPtr(&ptrReg, output, &free);
	if (input == NULL || output == NULL) {
	    printf ("ERROR:  Can't even begin:  out of memory.\n");
	    freeOnExit(&ptrReg);
	    exit (ERROR_RETURN);
	}

	/* Get names of input and output files. */
	for (i = 1;  i < argc;  i++) {
	    if (argv[i][0] == '-') {
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
		    } else if (argv[i][j] == 'v') {
			verbose = 1;
		    } else {
			printf ("ERROR:  Unrecognized option %s\n", argv[i]);
			printSyntax();
			freeOnExit(&ptrReg);
			exit (1);
		    }
		}
	    } else if (input[0] == '\0') {
		strcpy (input, argv[i]);
	    } else if (output[0] == '\0') {
		strcpy (output, argv[i]);
	    } else {
		too_many = 1;
	    }
	}
	if (input[0] == '\0' || too_many) {
	    printSyntax();
	    freeOnExit(&ptrReg);
	    exit (ERROR_RETURN);
	}
	if (output[0] == '\0') {
	    if ((status = MkName (input, "_x2d", "_sx2", output, STIS_LINE)))
	    {
	        freeOnExit(&ptrReg);
	        exit (status);
	    }
	}

	/* Sum imsets. */
	if ((status = CalStis8 (input, output, printtime, verbose))) {
	    printf ("Error processing %s.\n", input);
	    WhichError (status);
	}

	freeOnExit(&ptrReg);

	if (status)
	    exit (ERROR_RETURN);
	else
	    exit (0);
}
