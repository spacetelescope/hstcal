/* AcsSum -- Sum repeatobs data */

# include <stdio.h>
# include <stdlib.h>		/* calloc */
# include <string.h>

int status = 0;			/* zero is OK */

# include <c_iraf.h>		/* for c_irafinit */

# include "acs.h"
# include "acssum.h"
# include "hstcalerr.h"

/* 
    This function will only return either 0 (ACS_OK) if everything
        processed normally or ERROR_RETURN (2) if there was some error.
        WJH 21 April 1999
    Removed 'save_tmp' as a parameter to ACSSUM.  This parameter will 
        only be controlled by CALACS, not by individual tasks.
        WJH 19 April 2001
*/
        
int main (int argc, char **argv) {

	char *input, *output;	/* file names */
	char *mtype;        	/* MEMTYPE of output */
	int printtime = NO;	/* print time after each step? */
	int verbose = NO;	/* print additional info? */
	int quiet = NO;		/* suppress STDOUT messages? */
	int too_many = NO;	/* too many command-line arguments? */
	int i, j;		/* loop indexes */

	int AcsSum (char *, char *, char *, int, int);
	int MkName (char *, char *, char *, char *, char *, int);
	void WhichError (int);

	c_irafinit (argc, argv);

	input = calloc (ACS_LINE+1, sizeof (char));
	output = calloc (ACS_LINE+1, sizeof (char));
	mtype = calloc (SZ_STRKWVAL+1, sizeof (char));
	if (input == NULL || output == NULL) {
	    printf ("Can't even begin:  out of memory.\n");
		exit (ERROR_RETURN);
	}
		
    mtype[0] = '\0';
    
	/* Get names of input and output files. */
	for (i = 1;  i < argc;  i++) {
	    if (argv[i][0] == '-') {
		for (j = 1;  argv[i][j] != '\0';  j++) {
		    if (argv[i][j] == 't') {
			printtime = YES;
		    } else if (argv[i][j] == 'v') {
			verbose = YES;
		    } else if (argv[i][j] == 'q') {
			quiet = YES;
		    } else {
				printf (MsgText, "Unrecognized option %s\n", argv[i]);
				free (input);
				free (output);
				exit (1);
		    }
		}
	    } else if (input[0] == '\0') {
		strcpy (input, argv[i]);
	    } else if (output[0] == '\0') {
		strcpy (output, argv[i]);
	    } else {
		too_many = YES;
	    }
	}
	if (input[0] == '\0' || too_many) {
	    printf ("syntax:  acssum [-t] [-v] [-q] input output\n");
		free (input);
		free (output);
		exit (ERROR_RETURN);
	}

	/* Initialize the structure for managing trailer file comments */
	InitTrlBuf ();
	
	/* Copy command-line value for QUIET to structure */
	SetTrlQuietMode(quiet);

	if (output[0] == '\0') {
	    if (MkName (input, "_asn", "_sfl", "", output, ACS_LINE))
		CloseTrlBuf ();
		free (input);
		free (output);
	    WhichError (status);
		exit (ERROR_RETURN);
	}

	/* Sum imsets. */
	if (AcsSum (input, output, mtype, printtime, verbose)) {
	    sprintf (MsgText, "Error processing %s.", input);
		trlerror (MsgText);
	    WhichError (status);
	}
	free (input);
	free (output);
    free (mtype);

	CloseTrlBuf ();

	if (status)
	    exit (ERROR_RETURN);
	else
	    exit (ACS_OK);
}
