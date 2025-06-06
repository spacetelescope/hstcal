# include <stdio.h>
# include <stdlib.h>
# include <string.h>

#include "hstcal_memory.h"
#include "hstcal.h"
# include "wf3.h"
# include "hstcalerr.h"
# include "wf3version.h"
# include "hstcalversion.h"
# include "trlbuf.h"

/* H. Bushouse	07-Sep-2011	Implemented new "--version" command line argument. */
/* M. Sosey     added a -r to also print the version (-v is used, so warren chose r for revision */

/* CALWF3 driver: retrieves input table name and calls pipeline
 */

extern int	status;

/* Standard string buffer for use in messages */
char MsgText[MSG_BUFF_LENGTH]; // Global char auto initialized to '\0'

static void printSyntax(void)
{
    printf("syntax:  calwf3.e [--help] [-t] [-s] [-v] [-q] [-r] [-1] [--version] [--gitinfo] input \n");
}
static void printHelp(void)
{
    printSyntax();
}

int main (int argc, char **argv) {

	/* Local variables */
	char input[CHAR_FNAME_LENGTH+1];
	int printtime = NO;	/* print time after each step? */
	int save_tmp = DUMMY;	/* save temporary files? */
	int verbose = NO;	/* print info during processing? */
	int debug = NO;		/* print debug statements during processing? */
	int quiet = NO;		/* suppress STDOUT messages? */
	int too_many = NO;	/* too many command-line arguments? */
	int onecpu = NO;  /* suppress openmp usage by using only 1 thread?*/
    int i, j;		/* loop indexes */

	/* Function definitions */
	void c_irafinit (int, char **);
	int  CalWf3Run  (char *, int, int, int, int, int);
	void WhichError (int);

	/* Initialize status to OK and MsgText to null */
	status     = WF3_OK;
	MsgText[0] = '\0';
	input[0]   = '\0';

	/* Initialize IRAF environment */
	c_irafinit(argc, argv);
	PtrRegister ptrReg;
	initPtrRegister(&ptrReg);

	/* Command line arguments: 
	 **	0. Check for --version option
	 **	1. input file name
	 **	2. print time?
	 **	3. save intermediate files?
	 **	4. verbose?
	 */
	for (i = 1;  i < argc;  i++) {
		if (!(strcmp(argv[i],"--version"))) {
			printf("%s\n",WF3_CAL_VER_NUM);
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
		if (argv[i][0] == '-') {
			for (j = 1;  argv[i][j] != '\0';  j++) {
				if (argv[i][j] == 't') {
					printtime = YES;
				} else if (argv[i][j] == 's') {
					save_tmp = YES;
				} else if (argv[i][j] == 'r'){
					printf ("Current version: %s\n", WF3_CAL_VER);
					freeOnExit(&ptrReg);
					exit(0);
				} else if (argv[i][j] == 'v') {
					verbose = YES;
				} else if (argv[i][j] == 'd') {
					debug = NO;
				} else if (argv[i][j] == 'q') {
					quiet = YES;
				} else if (argv[i][j] == '1'){
                    onecpu = YES;
                } else {
					printf ("Unrecognized option %s\n", argv[i]);
					printSyntax();
					freeOnExit(&ptrReg);
					exit (ERROR_RETURN);
				}
			}
		} else if (input[0] == '\0') {
			strcpy (input, argv[i]);
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
	SetTrlQuietMode (quiet);

	/* Call the CALWF3 main program */
	if (CalWf3Run (input, printtime, save_tmp, verbose, debug, onecpu)) {

		if (status == NOTHING_TO_DO) {
			/* If there is just nothing to do, 
			 ** as for ACQ images, just quit. */
			status = 0;
			trlmessage("CALWF3 did NOT process %s", input);
			freeOnExit(&ptrReg);
			exit(0); 
		} else {
			/* Error during processing */
			trlerror("CALWF3 processing NOT completed for %s", input);
			/* Provide interpretation of error for user */
			WhichError (status);
			freeOnExit(&ptrReg);
			exit (ERROR_RETURN);
		}
	}

	/* Successful completion */
	trlmessage("CALWF3 completion for %s", input);

	freeOnExit(&ptrReg);

	/* Exit the program */
	exit (0);
}

