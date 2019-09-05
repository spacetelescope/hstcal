/* ACSCTE -- CTE loss correction */

# include <stdio.h>
# include <stdlib.h>		/* calloc */
# include <time.h>
# include <string.h>
#include <stdbool.h>
#include <limits.h>

extern int status;			/* zero is OK */

# include <c_iraf.h>		/* for c_irafinit */
#include "hstcal_memory.h"
#include "hstcal.h"
# include "ximio.h"
# include "hstio.h"

# include "acs.h"
# include "acsinfo.h"
# include "hstcalerr.h"
# include "acscorr.h"		/* calibration switch names for acsccd */
# include "hstcalversion.h"
# include "acsversion.h"
#include "trlbuf.h"

# ifdef _OPENMP
#  include <omp.h>
# endif

static char *program;
struct TrlBuf trlbuf = { 0 };

/* Standard string buffer for use in messages */
char MsgText[MSG_BUFF_LENGTH]; // Global char auto initialized to '\0'

/* This is the main module for ACSCTE.  It gets the input and output
 file names, calibration switches, and flags, and then calls ACScte.

 Pey Lian Lim, 2013 Aug 9:
 Separated PCTECORR step from ACSCCD. It is now a separate sub-module.

 */

int doMainCTE (int argc, char **argv);

static void printSyntax()
{
    printf ("syntax:  %s [--help] [-t] [-v] [-q] [--version] [--gitinfo] [-1|--nthreads <N>] [--ctegen <1|2>] [--pctetab <path>] [--forwardModelOnly] input [output]\n", program);
}
static void printHelp(void)
{
    printSyntax();
}

int doMainCTE (int argc, char **argv) {

    char *inlist;		/* input blv_tmp file name */
    char *outlist;		/* output file name */
    /*int switch_on = 0;*/	/* was any switch specified? */
    int printtime = NO;	/* print time after each step? */
    int verbose = NO;	/* print additional info? */
    int onecpu = NO; /* Use OpenMP (multi vs single CPU mode), if available? */
    int quiet = NO;	/* print additional info? */
    unsigned cteAlgorithmGen = 0; //Use gen1cte algorithm rather than gen2 (default)
    unsigned nThreads = 0;
    bool forwardModelOnly = false;
    char pcteTabNameFromCmd[CHAR_LINE_LENGTH];
    *pcteTabNameFromCmd = '\0';
    int too_many = 0;	/* too many command-line arguments? */
    int i, j;		/* loop indexes */
    int k;

    IRAFPointer i_imt, o_imt;	/* imt list pointers */
    char *input;		/* name of input science file */
    char *output;		/* name of output file */
    int n_in, n_out;	/* number of files in each list */
    int n;

    /* Input and output suffixes. */
    char isuffix[CHAR_FNAME_LENGTH] = "_blv_tmp";
    char osuffix[CHAR_FNAME_LENGTH] = "_blc_tmp"; // "_ctefmod" for forwardModelOnly==true - gets assigned later.

    /* A structure to pass the calibration switches to ACSCTE */
    CalSwitch ccd_sw;

    /* reference file keywords and names */
    RefFileInfo refnames;

    void InitRefFile (RefFileInfo *);
    void FreeRefFile (RefFileInfo *);
    void initSwitch (CalSwitch *);

    int ACScte (char *, char *, CalSwitch *, RefFileInfo *, int, int, int,
            const unsigned cteAlgorithmGen, const char * pcteTabNameFromCmd,
            const bool forwardModelOnly);
    int DefSwitch (char *);
    int MkName (char *, char *, char *, char *, char *, int);
    void WhichError (int);
    int CompareNumbers (int, int, char *);

    /* For image header access */
    Hdr phdr;
    int pctecorr;
    int LoadHdr (char *, Hdr *);
    int GetSwitch (Hdr *, char *, int *);

    /* For program basename */
    char program_buf[PATH_MAX];

    status = 0;

    c_irafinit (argc, argv);

    PtrRegister ptrReg;
    initPtrRegister(&ptrReg);
    /* Allocate space for file names. */
    inlist = calloc (CHAR_LINE_LENGTH+1, sizeof (char));
    addPtr(&ptrReg, inlist, &free);
    outlist = calloc (CHAR_LINE_LENGTH+1, sizeof (char));
    addPtr(&ptrReg, outlist, &free);
    input = calloc (CHAR_LINE_LENGTH+1, sizeof (char));
    addPtr(&ptrReg, input, &free);
    output = calloc (CHAR_LINE_LENGTH+1, sizeof (char));
    addPtr(&ptrReg, output, &free);

    if (!inlist || !outlist || !input || !output) {
        printf ("Can't even begin; out of memory.\n");
        freeOnExit(&ptrReg);
        exit (ERROR_RETURN);
    }
    inlist[0] = '\0';
    outlist[0] = '\0';
    input[0] = '\0';
    output[0] = '\0';

    /* Initialize the lists of reference file keywords and names. */
    InitRefFile (&refnames);
    addPtr(&ptrReg, &refnames, &FreeRefFile);

    /* Initial values. */
    initSwitch (&ccd_sw);

    /* Obtain program basename */
    if ((program = strrchr(argv[0], '/')) != NULL) {

        strcpy(program_buf, program + 1);
        program = program_buf;

    } else {

        program = argv[0];
    }

    if (!strcmp(program, "acscteforwardmodel.e")) {
        forwardModelOnly = true;
    }

    for (i = 1;  i < argc;  i++) {

        if (argv[i][0] == '-') {
            if (!(strcmp(argv[i],"--version")))
            {
                printf("%s\n",ACS_CAL_VER);
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
            if (strncmp(argv[i], "--ctegen", 8) == 0)
            {
                if (i + 1 > argc - 1)
                {
                    printf("ERROR: --ctegen - CTE algorithm generation not specified.\n");
                    freeOnExit(&ptrReg);
                    exit(1);
                }
                ++i;
                cteAlgorithmGen = (unsigned)atoi(argv[i]);
                if (cteAlgorithmGen != 1 && cteAlgorithmGen != 2)
                {
                    printf("ERROR: --ctegen - value out of range. Please specify either generation 1 or 2.\n");
                    freeOnExit(&ptrReg);
                    exit(1);
                }
                continue;
            }
            else if (strncmp(argv[i], "--nthreads", 10) == 0)
            {
                if (i + 1 > argc - 1)
                {
                    printf("ERROR: --nthreads - number of threads not specified\n");
                    freeOnExit(&ptrReg);
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
            else if (strncmp(argv[i], "--pctetab", 9) == 0)
            {
                if (i + 1 > argc - 1)
                {
                    printf("ERROR: --pctetab - no file specified\n");
                    freeOnExit(&ptrReg);
                    exit(1);
                }
                ++i;
                strcpy(pcteTabNameFromCmd, argv[i]);
                continue;
            }
            else if (strncmp(argv[i], "--forwardModelOnly", 18) == 0)
            {
                printf("WARNING: running CTE forward model only\n");
                forwardModelOnly = true;
                continue;
            }
            else
            {
                if (argv[i][1] == '-')
                {
                    printf ("Unrecognized option %s\n", argv[i]);
                    printSyntax();
                    freeOnExit(&ptrReg);
                    exit (ERROR_RETURN);
                }
                for (j = 1;  argv[i][j] != '\0';  j++) {
                    if (argv[i][j] == 't') {
                        printtime = YES;
                    } else if (argv[i][j] == 'v') {
                        verbose = YES;
                    } else if (argv[i][j] == 'q') {
                        quiet = YES;
                    } else if (argv[i][j] == '1') {
                        onecpu = YES;
                    } else {
                        printf ("Unrecognized option %s\n", argv[i]);
                        printSyntax();
                        break;
                    }
                }
            }
        } else if (inlist[0] == '\0') {
            strcpy (inlist, argv[i]);
        } else if (outlist[0] == '\0') {
            strcpy (outlist, argv[i]);
        } else {
            too_many = 1;
        }
    }
    if (inlist[0] == '\0' || too_many) {
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

    if (forwardModelOnly && cteAlgorithmGen == 1)
    {
        trlerror("--forwardModelOnly NOT compatible with 1st generation CTE algorithm");
        freeOnExit(&ptrReg);
        exit (INVALID_VALUE);
    }

    if (cteAlgorithmGen)
    {
        sprintf(MsgText, "(pctecorr) Using generation %d CTE algorithm", cteAlgorithmGen);
        trlmessage(MsgText);
    }

    if (*pcteTabNameFromCmd != '\0')
    {
        sprintf (MsgText, "(pctecorr) Using cmd line specified PCTETAB file: '%s'", pcteTabNameFromCmd);
        trlmessage(MsgText);
    }

#ifdef _OPENMP
    unsigned ompMaxThreads = omp_get_num_procs();
#endif
    if (onecpu)
    {
        if (nThreads )
            trlwarn("WARNING: option '-1' takes precedence when used in conjunction with '--nthreads <N>'");
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

    /* Was no calibration switch specified on command line?
       default values (mostly PERFORM)
    if (!switch_on) {*/
    ccd_sw.pctecorr = DefSwitch ("pctecorr");
    /*}*/

    /* Expand the templates. */
    i_imt = c_imtopen (inlist);
    addPtr(&ptrReg, i_imt, &c_imtclose);
    o_imt = c_imtopen (outlist);
    addPtr(&ptrReg, o_imt, &c_imtclose);
    n_in = c_imtlen (i_imt);
    n_out = c_imtlen (o_imt);

    /* The number of input and output files must be the same. */
    if (CompareNumbers (n_in, n_out, "output"))
        status = 1;
    if (status) {
        freeOnExit(&ptrReg);
        exit (ERROR_RETURN);
    }

    /* Loop over the list of input files. */
    for (n = 0;  n < n_in;  n++) {

        k = c_imtgetim (i_imt, input, CHAR_LINE_LENGTH);

        if (n_out > 0)
            k = c_imtgetim (o_imt, output, CHAR_LINE_LENGTH);
        else
            output[0] = '\0';

        /* Open input image in order to read its primary header. */
        if (LoadHdr (input, &phdr)) {
            WhichError (status);
            sprintf (MsgText, "Skipping %s", input);
            trlmessage (MsgText);
            continue;
        }

        /* Determine osuffix. */
        if (GetSwitch (&phdr, "PCTECORR", &pctecorr)) {
            WhichError (status);
            sprintf (MsgText, "Skipping %s", input);
            trlmessage (MsgText);
            continue;
        }
        if (pctecorr == PERFORM)
        {
            if (forwardModelOnly)
                strcpy(osuffix, "_ctefmod");
            else
                strcpy(osuffix, "_blc_tmp");
        }
        else {
            WhichError (status);
            sprintf (MsgText, "Skipping %s because PCTECORR is not set to PERFORM", input);
            trlmessage (MsgText);
            continue;
        }

        if (MkName (input, isuffix, osuffix, "", output, CHAR_LINE_LENGTH)) {
            WhichError (status);
            sprintf (MsgText, "Skipping %s", input);
            trlmessage (MsgText);
            continue;
        }

        /* Calibrate the current input file. */
        if ((status = ACScte (input, output, &ccd_sw, &refnames, printtime, verbose,
                    nThreads, cteAlgorithmGen, pcteTabNameFromCmd, forwardModelOnly))) {
            sprintf (MsgText, "Error processing %s.", input);
            trlerror (MsgText);
            WhichError (status);
        }
    }


    freeOnExit(&ptrReg);

    if (status)
        exit (ERROR_RETURN);
    else
        exit (0);
}
