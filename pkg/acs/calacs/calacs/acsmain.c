# include <stdio.h>
# include <stdlib.h>
# include <string.h>

# include "acs.h"
# include "acserr.h"
# include "acsversion.h"

/* CALACS driver: retrieves input table name and calls pipeline
*/

int	status;		/* value of zero indicates OK */

int main(int argc, char **argv) {

	/* Local variables */
	char input[ACS_FNAME+1];
	int printtime = NO;	/* print time after each step? */
	int save_tmp = DUMMY;	/* save temporary files? */
	int verbose = NO;	/* print info during processing? */
	int debug = NO;		/* print debug statements during processing? */
	int quiet = NO;		/* suppress STDOUT messages? */
	int onecpu = NO;		/* suppress OpenMP usage? */
	int too_many = NO;	/* too many command-line arguments? */
	int i, j;		/* loop indexes */

	/* Function definitions */
	void c_irafinit (int, char **);
	int CalAcsRun (char *, int, int, int, int, int);
    void WhichError (int);

	/* Initialize status to OK and MsgText to null */
	status = ACS_OK;
	MsgText[0] = '\0';
	input[0] = '\0';

	/* Initialize IRAF environment */
	c_irafinit(argc, argv);

	/* Command line arguments:
	**       0. Check for --version option
  **
  **       1. input file name
	**		   2. print time?
	**		   3. save intermediate files?
	**		   4. verbose?
	*/
	for (i = 1;  i < argc;  i++) {
		if (!(strcmp(argv[i],"--version"))) {
      printf("%s\n",ACS_CAL_VER_NUM);
      exit(0);
    }

    if (argv[i][0] == '-') {
		  for (j = 1;  argv[i][j] != '\0';  j++) {
		    if (argv[i][j] == 't') {
          printtime = YES;
		    } else if (argv[i][j] == 's') {
          save_tmp = YES;
		    } else if (argv[i][j] == 'r') {
	  printf ("Current version: %s\n", ACS_CAL_VER);
	  exit(0);
		    } else if (argv[i][j] == 'v') {
          verbose = YES;
		    } else if (argv[i][j] == 'd') {
          debug = YES;
		    } else if (argv[i][j] == 'q') {
          quiet = YES;
		    } else if (argv[i][j] == '1') {
          onecpu = YES;
		    } else {
          printf ("Unrecognized option %s\n", argv[i]);
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
        printf ("CALACS Version %s\n",ACS_CAL_VER_NUM);
	    printf ("syntax:  calacs.e [-t] [-s] [-v] [-q] [-r] [-1] input \n");
	    exit (ERROR_RETURN);
	}

	/* Initialize the structure for managing trailer file comments */
	InitTrlBuf ();

	/* Copy command-line value for QUIET to structure */
	SetTrlQuietMode(quiet);

	/* Call the CALACS main program */
	if (CalAcsRun (input, printtime, save_tmp, verbose, debug, onecpu)) {

        if (status == NOTHING_TO_DO){
            /* If there is just nothing to do,
                as for ACQ images, just quit...     WJH 27Apr1999 */
            status = 0;
	        sprintf (MsgText, "CALACS did NOT process %s", input);
	        trlmessage (MsgText);
            exit(0);
        } else {
	        /* Error during processing */
	        sprintf (MsgText, "CALACS processing NOT completed for %s", input);
		    trlerror (MsgText);
            /* Added 19 Mar 1999 - provides interpretation of error for user */
            WhichError (status);
		    CloseTrlBuf ();
	        exit (ERROR_RETURN);
        }
	}

	/* Successful completion */
	sprintf (MsgText, "CALACS completion for %s", input);
	trlmessage (MsgText);

	CloseTrlBuf ();

	/* Exit the program */
	exit(0);
}
