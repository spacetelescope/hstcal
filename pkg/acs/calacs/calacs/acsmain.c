# include <stdio.h>
# include <stdlib.h>
# include <string.h>

# include "acs.h"
# include "hstcalerr.h"
# include "acsversion.h"
# include "hstcalversion.h"

# ifdef _OPENMP
#  include <omp.h>
# endif

/* CALACS driver: retrieves input table name and calls pipeline
*/

int	status;		/* value of zero indicates OK */

static void printSyntax()
{
    printf ("syntax:  calacs.e [-t] [-s] [-v] [-q] [-r] [--version] [--gitinfo] [-1|--nthreads <N>] [--ctegen <1|2>] [--pctetab <path>] input \n");
}

int main(int argc, char **argv) {

	/* Local variables */
	char input[ACS_FNAME+1];
	int printtime = NO;	/* print time after each step? */
	int save_tmp = DUMMY;	/* save temporary files? */
	int verbose = NO;	/* print info during processing? */
	int debug = NO;		/* print debug statements during processing? */
	int quiet = NO;		/* suppress STDOUT messages? */
	int onecpu = NO;		/* suppress OpenMP usage? */
    unsigned cteAlgorithmGen = 0; //Use gen1cte algorithm rather than gen2 (default)
    char pcteTabNameFromCmd[ACS_LINE];
    *pcteTabNameFromCmd = '\0';
	int too_many = NO;	/* too many command-line arguments? */
	int i, j;		/* loop indexes */
    unsigned nThreads = 0;

	/* Function definitions */
	void c_irafinit (int, char **);
	int CalAcsRun (char *, int, int, int, int, const unsigned nThreads, const int gen1cte, const char * pcteTabNameFromCmd);
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

    for (i = 1;  i < argc;  i++)
    {
        if (!(strcmp(argv[i],"--version")))
        {
            printf("%s\n",ACS_CAL_VER);
            exit(0);
        }
        if (!(strcmp(argv[i],"--gitinfo")))
        {
            printGitInfo();
            exit(0);
        }
        if (strncmp(argv[i], "--ctegen", 8) == 0)
        {
            if (i + 1 > argc - 1)
            {
                printf("ERROR: --ctegen - CTE algorithm generation not specified.\n");
                exit(1);
            }
            ++i;
            cteAlgorithmGen = (unsigned)atoi(argv[i]);
            if (cteAlgorithmGen != 1 && cteAlgorithmGen != 2)
            {
                printf("ERROR: --ctegen - value out of range. Please specify either generation 1 or 2.\n");
                exit(1);
            }
            continue;
        }
        else if (strncmp(argv[i], "--pctetab", 9) == 0)
        {
            if (i + 1 > argc - 1)
            {
                printf("ERROR: --pctetab - no file specified\n");
                exit(1);
            }
            strcpy(pcteTabNameFromCmd, argv[i+1]);
            printf("WARNING: using pcteTab file '%s'\n", pcteTabNameFromCmd);
            ++i;
            continue;
        }
        else if (strncmp(argv[i], "--nthreads", 10) == 0)
        {
            if (i + 1 > argc - 1)
            {
                printf("ERROR: --nthreads - number of threads not specified\n");
                exit(1);
            }
            ++i;
            nThreads = (unsigned)atoi(argv[i]);
            if (nThreads < 1)
                nThreads = 1;
#ifndef _OPENMP
            printf("WARNING: '--nthreads <N>' used but OPENMP not found!\n");
            nThreads = 1;
#endif
            continue;
        }
        if (argv[i][0] == '-')
        {
            if (argv[i][1] == '-')
            {
                printf ("Unrecognized option %s\n", argv[i]);
                printSyntax();
                exit (ERROR_RETURN);
            }
            for (j = 1;  argv[i][j] != '\0';  j++)
            {
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
                    printSyntax();
                    exit (ERROR_RETURN);
                }
            }
        }
        else if (input[0] == '\0')
            strcpy (input, argv[i]);
        else
            too_many = YES;
    }

    if (input[0] == '\0' || too_many) {
        printf ("CALACS Version %s\n",ACS_CAL_VER_NUM);
        printSyntax();
        exit (ERROR_RETURN);
    }

	/* Initialize the structure for managing trailer file comments */
	InitTrlBuf ();
    trlGitInfo();

	/* Copy command-line value for QUIET to structure */
	SetTrlQuietMode(quiet);

    if (cteAlgorithmGen)
    {
        sprintf(MsgText, "(pctecorr) Using generation %d CTE algorithm", cteAlgorithmGen);
        trlmessage(MsgText);
    }

    if (*pcteTabNameFromCmd != '\0')
    {
        sprintf(MsgText, "(pctecorr) Using cmd line specified PCTETAB file: '%s'", pcteTabNameFromCmd);
        trlmessage(MsgText);
    }

#ifdef _OPENMP
    unsigned ompMaxThreads = omp_get_num_procs();
#endif
    if (onecpu)
    {
        if (nThreads)
            trlwarn("Option '-1' takes precedence when used in conjunction with '--nthreads <N>'");
        nThreads = 1;
    }
    else if (!nThreads)//unset
    {
#ifdef _OPENMP
        nThreads = ompMaxThreads;
#else
        nThreads = 1;
#endif
    }

#ifdef _OPENMP
    omp_set_dynamic(0);
    if (nThreads > ompMaxThreads)
    {
        sprintf(MsgText, "System env limiting nThreads from %d to %d", nThreads, ompMaxThreads);
        nThreads = ompMaxThreads;
    }
    else
        sprintf(MsgText,"Setting max threads to %d out of %d available", nThreads, ompMaxThreads);

    omp_set_num_threads(nThreads);
    trlmessage(MsgText);
#endif

	/* Call the CALACS main program */
	if (CalAcsRun (input, printtime, save_tmp, verbose, debug, nThreads, cteAlgorithmGen, pcteTabNameFromCmd)) {

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
