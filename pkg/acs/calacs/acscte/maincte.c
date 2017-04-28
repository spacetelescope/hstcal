/* ACSCTE -- CTE loss correction */

# include <stdio.h>
# include <stdlib.h>		/* calloc */
# include <time.h>
# include <string.h>

int status = 0;			/* zero is OK */

# include <c_iraf.h>		/* for c_irafinit */
# include "ximio.h"
# include "hstio.h"

# include "acs.h"
# include "acsinfo.h"
# include "hstcalerr.h"
# include "acscorr.h"		/* calibration switch names for acsccd */

# ifdef _OPENMP
#  include <omp.h>
# endif

static void FreeNames (char *, char *, char *, char *);

/* This is the main module for ACSCTE.  It gets the input and output
 file names, calibration switches, and flags, and then calls ACScte.

 Pey Lian Lim, 2013 Aug 9:
 Separated PCTECORR step from ACSCCD. It is now a separate sub-module.

 */

int main (int argc, char **argv) {

    char *inlist;		/* input blv_tmp file name */
    char *outlist;		/* output blc_tmp file name */
    /*int switch_on = 0;*/	/* was any switch specified? */
    int printtime = NO;	/* print time after each step? */
    int verbose = NO;	/* print additional info? */
    int onecpu = NO; /* Use OpenMP (multi vs single CPU mode), if available? */
    int quiet = NO;	/* print additional info? */
    unsigned cteAlgorithmGen = 0; //Use gen1cte algorithm rather than gen2 (default)
    unsigned nThreads = 0;
    char pcteTabNameFromCmd[ACS_LINE];
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
    char isuffix[] = "_blv_tmp";
    char osuffix[] = "_blc_tmp";

    /* A structure to pass the calibration switches to ACSCTE */
    CalSwitch ccd_sw;

    /* reference file keywords and names */
    RefFileInfo refnames;

    void InitRefFile (RefFileInfo *);
    void FreeRefFile (RefFileInfo *);
    void initSwitch (CalSwitch *);

    int ACScte (char *, char *, CalSwitch *, RefFileInfo *, int, int, int, const unsigned cteAlgorithmGen, const char * pcteTabNameFromCmd);
    int DefSwitch (char *);
    int MkName (char *, char *, char *, char *, char *, int);
    void WhichError (int);
    int CompareNumbers (int, int, char *);

    /* For image header access */
    Hdr phdr;
    int pctecorr;
    int LoadHdr (char *, Hdr *);
    int GetSwitch (Hdr *, char *, int *);

    c_irafinit (argc, argv);

    /* Allocate space for file names. */
    inlist = calloc (ACS_LINE+1, sizeof (char));
    outlist = calloc (ACS_LINE+1, sizeof (char));
    input = calloc (ACS_LINE+1, sizeof (char));
    output = calloc (ACS_LINE+1, sizeof (char));

    if (inlist == NULL || outlist == NULL ||
        input == NULL || output == NULL) {
        printf ("Can't even begin; out of memory.\n");
        exit (ERROR_RETURN);
    }
    inlist[0] = '\0';
    outlist[0] = '\0';
    input[0] = '\0';
    output[0] = '\0';

    /* Initialize the lists of reference file keywords and names. */
    InitRefFile (&refnames);

    /* Initial values. */
    initSwitch (&ccd_sw);

    for (i = 1;  i < argc;  i++) {

        if (argv[i][0] == '-') {
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
            else if (strncmp(argv[i], "--pctetab", 9) == 0)
            {
                if (i + 1 > argc - 1)
                {
                    printf("ERROR: --pctetab - no file specified\n");
                    exit(1);
                }
                ++i;
                strcpy(pcteTabNameFromCmd, argv[i]);
                continue;
            }
            else
            {
                if (argv[i][1] == '-')
                {
                    printf ("Unrecognized option %s\n", argv[i]);
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
        printf ("syntax:  acscte [-t] [-v] [-q] [-1|--nthreads <N>] [--ctegen <1|2>] [--pctetab <path>] input output\n");
        FreeNames (inlist, outlist, input, output);
        exit (ERROR_RETURN);
    }
    /* Initialize the structure for managing trailer file comments */
    InitTrlBuf ();

    /* Copy command-line value for QUIET to structure */
    SetTrlQuietMode(quiet);

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
    o_imt = c_imtopen (outlist);
    n_in = c_imtlen (i_imt);
    n_out = c_imtlen (o_imt);

    /* The number of input and output files must be the same. */
    if (CompareNumbers (n_in, n_out, "output"))
        status = 1;
    if (status) {
        FreeNames (inlist, outlist, input, output);
        CloseTrlBuf();
        exit (ERROR_RETURN);
    }

    /* Loop over the list of input files. */
    for (n = 0;  n < n_in;  n++) {

        k = c_imtgetim (i_imt, input, ACS_LINE);

        if (n_out > 0)
            k = c_imtgetim (o_imt, output, ACS_LINE);
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
            strcpy(osuffix, "_blc_tmp");
        else {
            WhichError (status);
            sprintf (MsgText, "Skipping %s because PCTECORR is not set to PERFORM", input);
            trlmessage (MsgText);
            continue;
        }

        if (MkName (input, isuffix, osuffix, "", output, ACS_LINE)) {
            WhichError (status);
            sprintf (MsgText, "Skipping %s", input);
            trlmessage (MsgText);
            continue;
        }

        /* Calibrate the current input file. */
        if ((status = ACScte (input, output, &ccd_sw, &refnames, printtime, verbose,
                    nThreads, cteAlgorithmGen, pcteTabNameFromCmd))) {
            sprintf (MsgText, "Error processing %s.", input);
            trlerror (MsgText);
            WhichError (status);
        }
    }

    /* Close lists of file names, and free name buffers. */
    c_imtclose (i_imt);
    c_imtclose (o_imt);
    CloseTrlBuf();
    FreeRefFile (&refnames);
    FreeNames (inlist, outlist, input, output);

    if (status)
        exit (ERROR_RETURN);
    else
        exit (0);
}


static void FreeNames (char *inlist, char *outlist, char *input, char *output) {
    free (output);
    free (input);
    free (outlist);
    free (inlist);
}
