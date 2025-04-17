/* ACSCTE -- CTE loss correction */

#include "config.h"
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

#include <argp.h>

#define OPT_GITINFO 1000
#define OPT_VERSION 1001

#define OPT_CTEGEN 2000
#define OPT_SINGLE_THREAD 2001
#define OPT_NTHREADS 2002
#define OPT_FORWARD_MODEL_ONLY 2003
#define OPT_PCETAB 2004

#define OPT_GROUP_INFO 0
#define OPT_GROUP_TOGGLE 0
#define OPT_GROUP_INPUT 0

static char doc[] = "ACS CTE";
static char args_doc[] = "INPUT OUTPUT";
static struct argp_option options[] = {
    {
        .name = "version",
        .key = OPT_VERSION,
        .arg = 0,
        .flags = 0,
        .doc = "Print version",
        .group = OPT_GROUP_INFO
    },
    {
        .name = "gitinfo",
        .key = OPT_GITINFO,
        .arg = 0,
        .flags = 0,
        .doc = "Print git version information",
        .group = OPT_GROUP_INFO
    },
    {
        .name = "verbose",
        .key = 'v',
        .arg = 0,
        .flags = 0,
        .doc = "Turn on verbose mode",
        .group = OPT_GROUP_TOGGLE
    },
    {
        .name = "quiet",
        .key = 'q',
        .arg = 0,
        .flags = 0,
        .doc = "Turn on quiet mode",
        .group = OPT_GROUP_TOGGLE
    },
    {
        .name = "print-timestamps",
        .key = 't',
        .arg = 0,
        .flags = 0,
        .doc = "Print timestamps",
        .group = OPT_GROUP_TOGGLE
    },
    {
        .name = "ctegen",
        .key = OPT_CTEGEN,
        .arg = "NUM",
        .flags = 0,
        .doc = "Use ctegen algorithm (default: 1)",
        .group = OPT_GROUP_INPUT
    },
    {
        .name = "single",
        .key = '1',
        .arg = 0,
        .flags = 0,
        .doc = "Use single-threaded mode (no OpenMP)",
        .group = OPT_GROUP_TOGGLE
    },
    {
        .name = "nthreads",
        .key = OPT_NTHREADS,
        .arg = "NUM",
        .flags = 0,
        .doc = "Use multi-threaded mode (OpenMP)",
        .group = OPT_GROUP_INPUT
    },
    {
        .name = "forwardModelOnly",
        .key = OPT_FORWARD_MODEL_ONLY,
        .arg = 0,
        .flags = 0,
        .doc = "Use CTE forward model",
        .group = OPT_GROUP_TOGGLE
    },
    {
        .name = "pcetab",
        .key = OPT_PCETAB,
        .arg = "FILE",
        .flags = 0,
        .doc = "Path to PCETAB file",
        .group = OPT_GROUP_INPUT
    },
    {.name = 0, .key = 0, .arg = 0, .flags = 0, .doc = 0, .group = 0} // END OF ARGUMENTS
};

struct arguments {
    char *args[2]; /* for input and output files */
    int *printtime;	/* print time after each step? */
    int *verbose;	/* print additional info? */
    int *onecpu; /* Use OpenMP (multi vs single CPU mode), if available? */
    int *quiet;	/* print additional info? */
    unsigned *cteAlgorithmGen; // Use gen1cte algorithm rather than gen2 (default)
    unsigned *nThreads; // number of threads to use
    bool *forwardModelOnly;
    char *pcteTabNameFromCmd;
};

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

static error_t parse_opt(int key, char *arg, struct argp_state *state) {
    struct arguments *arguments = state->input;
    switch (key) {
        case 'v':
            *arguments->verbose = YES;
            break;
        case 'q':
            *arguments->quiet = YES;
            break;
        case 't':
            *arguments->printtime = YES;
            break;
        case OPT_VERSION:
            printf("%s\n", ACS_CAL_VER);
            exit(0);
        case OPT_GITINFO:
            trlGitInfo();
            exit(0);
        case OPT_CTEGEN:
            *arguments->cteAlgorithmGen = strtol(arg, NULL, 10);
            if (*arguments->cteAlgorithmGen < 1 || *arguments->cteAlgorithmGen > 2) {
                status = INVALID_VALUE;
                argp_failure(state, status, errno, "--ctegen - value out of range. Please specify either generation 1 or 2.");
            }
            break;
        case OPT_FORWARD_MODEL_ONLY:
            fprintf(stderr, "WARNING: running CTE forward model only\n");
            *arguments->forwardModelOnly = YES;
            break;
        case OPT_SINGLE_THREAD:
            *arguments->onecpu = YES;
            break;
        case OPT_NTHREADS:
            *arguments->nThreads = strtol(arg, NULL, 10);
            if (*arguments->nThreads < 1) {
                *arguments->nThreads = 1;
            }
#ifndef _OPENMP
            fprintf(stderr, "WARNING: '--nthreads <N>' used but OPENMP not found!\n");
            arguments->nThreads = 1;
#endif
            break;
        case OPT_PCETAB:
            strncpy(arguments->pcteTabNameFromCmd, arg, sizeof(arguments->pcteTabNameFromCmd) - 1);
            break;
        case ARGP_KEY_ARG:
            if (state->arg_num >= 2) {
                argp_failure(state, ERROR_RETURN, 0, "Too many arguments\n");
            }
            arguments->args[state->arg_num] = arg;
            break;
        case ARGP_KEY_END:
            if (state->arg_num < 1) {
                argp_failure(state, ERROR_RETURN, 0, "Not enough arguments\n");
            }
            break;
        default:
            return ARGP_ERR_UNKNOWN;
    }

    return 0;
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


    struct argp argp = {options, parse_opt, args_doc, doc, 0, 0, 0};
    struct arguments arguments = {
        .verbose = &verbose,
        .quiet = &quiet,
        .printtime = &printtime,
        .cteAlgorithmGen = &cteAlgorithmGen,
        .onecpu = &onecpu,
        .nThreads = &nThreads,
        .pcteTabNameFromCmd = pcteTabNameFromCmd,
    };


    argp_parse(&argp, argc, argv, 0, 0, &arguments);
    strcpy(inlist, arguments.args[0]);
    if (arguments.args[1]) {
        strcpy(outlist, arguments.args[1]);
    }

    /* Obtain program basename */
    if ((program = strrchr(argv[0], '/')) != NULL) {
        strcpy(program_buf, strlen(program) ? program + 1 : program);
        program = program_buf;

    } else {
        program = argv[0];
    }

    if (inlist[0] == '\0') {
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
        if (nThreads)
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
