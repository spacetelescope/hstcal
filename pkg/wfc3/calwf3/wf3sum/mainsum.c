/* Wf3Sum -- Sum repeatobs data */

# include <stdio.h>
# include <stdlib.h>		/* calloc */
# include <string.h>

extern int status;

#include "hstcal_memory.h"
#include "hstcal.h"
# include "c_iraf.h"		/* for c_irafinit */

# include "wf3.h"
# include "wf3sum.h"
# include "hstcalerr.h"
# include "wf3version.h"
# include "hstcalversion.h"
# include "trlbuf.h"

/* Standard string buffer for use in messages */
char MsgText[MSG_BUFF_LENGTH]; // Global char auto initialized to '\0'
struct TrlBuf trlbuf = { 0 };

static void printSyntax(void)
{
    printf ("syntax:  wf3sum [--help] [-t] [-v] [-q] [-r] [--version] [--gitinfo] input output\n");
}
static void printHelp(void)
{
    printSyntax();
}

/* 
    This function will only return either 0 (WF3_OK) if everything
    processed normally or ERROR_RETURN (2) if there was some error.

    Howard Bushouse, 2001 Nov 16:
	Updates to track changes in CALACS - Removed 'save_tmp' as a
	parameter to WF3SUM. This parameter will only be controlled by
	CALWF3, not by individual tasks.
*/
        
int main (int argc, char **argv) {

	char *input, *output;	/* file names */
	char *mtype;        	/* MEMTYPE of output */
	int printtime = NO;	/* print time after each step? */
	int verbose = NO;	/* print additional info? */
	int quiet = NO;		/* suppress STDOUT messages? */
	int too_many = NO;	/* too many command-line arguments? */
	int i, j;		/* loop indexes */

	int Wf3Sum (char *, char *, char *, int, int);
	int MkName (char *, char *, char *, char *, char *, int);
	void WhichError (int);

/*===========================================================================*/
    status = 0;

	/* Initialize IRAF interface */
	c_irafinit (argc, argv);

    PtrRegister ptrReg;
    initPtrRegister(&ptrReg);
	/* Allocate local memory for file names */
	input  = calloc (CHAR_LINE_LENGTH+1, sizeof (char));
    addPtr(&ptrReg, input, &free);
	output = calloc (CHAR_LINE_LENGTH+1, sizeof (char));
    addPtr(&ptrReg, output, &free);
	mtype  = calloc (SZ_FITS_VAL+1, sizeof (char));
    addPtr(&ptrReg, mtype, &free);
	if (!input || !output || !mtype) {
	    printf ("Can't even begin:  out of memory.\n");
	    freeOnExit(&ptrReg);
	    exit (ERROR_RETURN);
	}

	mtype[0] = '\0';

	/* Get names of input and output files and all arguments. */
	for (i = 1;  i < argc;  i++) {
	    if (argv[i][0] == '-') {
            if (!(strcmp(argv[i],"--version")))
            {
                printf("%s\n",WF3_CAL_VER);
                freeOnExit(&ptrReg);
                exit(0);
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
		for (j = 1;  argv[i][j] != '\0';  j++) {
		    if (argv[i][j] == 't') {
			printtime = YES;
		    } else if (argv[i][j] == 'v') {
			verbose = YES;
		    } else if (argv[i][j] == 'q') {
			quiet = YES;
            } else if (argv[i][j] == 'r'){
                printf ("Current version: %s\n", WF3_CAL_VER);
                freeOnExit(&ptrReg);
                exit(0);
		    } else {
			printf (MsgText, "Unrecognized option %s\n", argv[i]);
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
		too_many = YES;
	    }
	}
	if (input[0] == '\0' || too_many) {
	    printSyntax();
	    freeOnExit(&ptrReg);
	    exit (ERROR_RETURN);
	}

	/* Initialize the structure for managing trailer file comments */
	InitTrlBuf ();
    addPtr(&ptrReg, &trlbuf, &CloseTrlBuf);
    trlGitInfo();

	/* Copy command-line value for QUIET to structure */
	SetTrlQuietMode(quiet);

	if (output[0] == '\0') {
	    if (MkName (input, "_asn", "_sfl", "", output, CHAR_LINE_LENGTH)) {
		WhichError (status);
		freeOnExit(&ptrReg);
		exit (ERROR_RETURN);
	    }
	}

	/* Sum imsets. */
	if (Wf3Sum (input, output, mtype, printtime, verbose)) {
	    sprintf (MsgText, "Error processing %s.", input);
	    trlerror (MsgText);
	    WhichError (status);
	}

	freeOnExit(&ptrReg);

	if (status)
	    exit (ERROR_RETURN);
	else
	    exit (WF3_OK);
}
