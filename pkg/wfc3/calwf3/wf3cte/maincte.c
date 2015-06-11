/* WFC3CTE -- CTE loss correction */

# include <stdio.h>
# include <stdlib.h>		/* calloc */
# include <time.h>
# include <string.h>

int status = 0;			/* zero is OK */

# include <c_iraf.h>		/* for c_irafinit */
# include "ximio.h"
# include "hstio.h"

# include "wf3.h"
# include "wf3info.h"
# include "wf3err.h"
# include "wf3corr.h"		/* calibration switch names for WFC3ccd */

static void FreeNames (char *, char *, char *, char *);

/* 

This is the main module for WF3CTE.  It gets the input and output
file names, calibration switches, and flags, and then calls WFC3cte.
This is necessary for the task to act standalone.

*/

int main (int argc, char **argv) {

    char *inlist;		/* input blv_tmp file name */
    char *outlist;		/* output blc_tmp file name */
    /*int switch_on = 0;*/	/* was any switch specified? */
    int printtime = NO;	/* print time after each step? */
    int verbose = NO;	/* print additional info? */
    int quiet = NO;	/* print additional info? */
    int onecpu = NO; /* Use OpenMP with onely one thread, if available? */
    int too_many = 0;	/* too many command-line arguments? */
    int i, j;		/* loop indexes */
    int k;

    IRAFPointer i_imt, o_imt;	/* imt list pointers */
    char *input;		/* name of input science file */
    char *output;		/* name of output file */
    int n_in, n_out;	/* number of files in each list */
    int n;

    /* Input and output suffixes. */
    char isuffix[] = "_raw";
    char osuffix[] = "_rac";

    /* A structure to pass the calibration switches to WFC3CTE */
    CalSwitch ccd_sw;

    /* reference file keywords and names */
    RefFileInfo refnames;

    void InitRefFile (RefFileInfo *);
    void FreeRefFile (RefFileInfo *);
    void initSwitch (CalSwitch *);

    int WF3cte (char *, char *, CalSwitch *, RefFileInfo *, int, int, int);
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
    inlist = calloc (SZ_FNAME+1, sizeof (char));
    outlist = calloc (SZ_FNAME+1, sizeof (char));
    input = calloc (SZ_FNAME+1, sizeof (char));
    output = calloc (SZ_FNAME+1, sizeof (char));

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
                    printf (MsgText, "Unrecognized option %s\n", argv[i]);
                    FreeNames (inlist, outlist, input, output);
                    exit (1);
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
        printf ("syntax:  WF3cte [-t] [-v] [-q] [-1] input output\n");
        FreeNames (inlist, outlist, input, output);
        exit (ERROR_RETURN);
    }
    /* Initialize the structure for managing trailer file comments */
    InitTrlBuf ();

    /* Copy command-line value for QUIET to structure */
    SetTrlQuietMode(quiet);

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

        k = c_imtgetim (i_imt, input, SZ_FNAME);

        if (n_out > 0)
            k = c_imtgetim (o_imt, output, SZ_FNAME);
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
            strcpy(osuffix, "_rac");
        else {
            WhichError (status);
            sprintf (MsgText, "Skipping %s because PCTECORR is not set to PERFORM", input);
            trlmessage (MsgText);
            continue;
        }

        if (MkName (input, isuffix, osuffix, "", output, SZ_FNAME)) {
            WhichError (status);
            sprintf (MsgText, "Skipping %s", input);
            trlmessage (MsgText);
            continue;
        }

        /* Calibrate the current input file. */
        if (WF3cte (input, output, &ccd_sw, &refnames, printtime, verbose,
                    onecpu)) {
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
