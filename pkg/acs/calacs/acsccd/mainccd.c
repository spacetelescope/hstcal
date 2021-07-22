/* ACSCCD -- basic CCD image reduction */

# include <stdio.h>
# include <stdlib.h>		/* calloc */
# include <time.h>
# include <string.h>

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

/* Standard string buffer for use in messages */
char MsgText[MSG_BUFF_LENGTH]; // Global char auto initialized to '\0'
struct TrlBuf trlbuf = { 0 };

static void printSyntax(void)
{
    printf ("syntax:  acsccd [--help] [-t] [-v] [-q] [--version] [--gitinfo] input output\n");
}
static void printHelp(void)
{
    printSyntax();
}

/* This is the main module for ACSCCD.  It gets the input and output
 file names, calibration switches, and flags, and then calls ACSccd.

 Warren Hack, 1998 June 1:
 Revised for ACS... It will work on 1 ACS image at a time. Based on
 original CALSTIS code by Phil Hodge.

 Pey Lian Lim, 2012 Dec 12:
 Moved FLSHCORR to ACS2D.

 Pey Lian Lim, 2012 Dec 19:
 Commented calibration switches because they are overwritten in
 getacsflags.c by image header anyway.

 Pey Lian Lim, 2013 Aug 9:
 Separated PCTECORR step from ACSCCD.

 */

int main (int argc, char **argv) {

    char *inlist;		/* input file name */
    char *outlist;		/* output blev file name */
    /*int switch_on = 0;*/	/* was any switch specified? */
    int printtime = NO;	/* print time after each step? */
    int verbose = NO;	/* print additional info? */
    int quiet = NO;	/* print additional info? */
    int too_many = 0;	/* too many command-line arguments? */
    int too_long = 0; /* command-line argument too long? */
    int i, j;		/* loop indexes */
    int k;

    IRAFPointer i_imt, o_imt;	/* imt list pointers */
    char *input;		/* name of input science file */
    char *output;		/* name of output file */
    int n_in, n_out;	/* number of files in each list */
    int n;

    /* Input and output suffixes. */
    char isuffix[] = "_raw";
    char osuffix[] = "_blv_tmp";

    /* A structure to pass the calibration switches to ACSCCD */
    CalSwitch ccd_sw;

    /* reference file keywords and names */
    RefFileInfo refnames;

    void InitRefFile (RefFileInfo *);
    void FreeRefFile (RefFileInfo *);
    void initSwitch (CalSwitch *);

    int ACSccd (char *, char *, CalSwitch *, RefFileInfo *, int, int);
    int DefSwitch (char *);
    int MkName (char *, char *, char *, char *, char *, int);
    void WhichError (int);
    int CompareNumbers (int, int, char *);

    /* For image header access */
    Hdr phdr;
    int LoadHdr (char *, Hdr *);
    int GetSwitch (Hdr *, char *, int *);

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

    if (!inlist || !outlist|| !input || !output) {
        printf ("Can't even begin; out of memory.\n");
        freeOnExit(&ptrReg);
        exit (ERROR_RETURN);
    }

    /* Initialize the lists of reference file keywords and names. */
    InitRefFile (&refnames);
    addPtr(&ptrReg, &refnames, &FreeRefFile);

    /* Initial values. */
    initSwitch (&ccd_sw);

    for (i = 1;  i < argc;  i++) {

        /**********
        if (strcmp (argv[i], "-dqi") == 0) {
            ccd_sw.dqicorr = PERFORM;
            switch_on = 1;
        } else if (strcmp (argv[i], "-atod") == 0) {
            ccd_sw.atodcorr = PERFORM;
            switch_on = 1;
        } else if (strcmp (argv[i], "-blev") == 0) {
            ccd_sw.blevcorr = PERFORM;
            switch_on = 1;
        } else if (strcmp (argv[i], "-bias") == 0) {
            ccd_sw.biascorr = PERFORM;
            switch_on = 1;
        } else if (argv[i][0] == '-') {
        **********/
        too_long = strlen(argv[i]) > CHAR_LINE_LENGTH;
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
            for (j = 1;  argv[i][j] != '\0';  j++) {
                if (argv[i][j] == 't') {
                    printtime = YES;
                } else if (argv[i][j] == 'v') {
                    verbose = YES;
                } else if (argv[i][j] == 'q') {
                    quiet = YES;
                } else {
                    printf (MsgText, "Unrecognized option %s\n", argv[i]);
                    printSyntax();
                    freeOnExit(&ptrReg);
                    exit (1);
                }
            }
        } else if (inlist[0] == '\0') {
            strncpy (inlist, argv[i], CHAR_LINE_LENGTH);
        } else if (outlist[0] == '\0') {
            strncpy (outlist, argv[i], CHAR_LINE_LENGTH);
        } else {
            too_many = 1;
        }
    }
    
    if (inlist[0] == '\0' || too_many || too_long) {
        if (too_many) {
            fprintf(stderr, "ERROR: Too many arguments\n");
        }
        
        if (too_long) {
            fprintf(stderr, "ERROR: Input path exceeds maximum supported length (%d)\n",
                    CHAR_LINE_LENGTH);
        }
        
        printSyntax();
        /*
        printf ("  command-line switches:\n");
        printf ("       -dqi -atod -blev -bias\n");
        */
        freeOnExit(&ptrReg);
        exit (ERROR_RETURN);
    }
    /* Initialize the structure for managing trailer file comments */
    InitTrlBuf ();
    addPtr(&ptrReg, &trlbuf, &CloseTrlBuf);
    trlGitInfo();

    /* Copy command-line value for QUIET to structure */
    SetTrlQuietMode(quiet);

    /* Was no calibration switch specified on command line?
       default values (mostly PERFORM) except ATODCORR
    if (!switch_on) {*/
    ccd_sw.dqicorr  = DefSwitch ("dqicorr");
    ccd_sw.atodcorr = DefSwitch ("atodcorr");
    ccd_sw.blevcorr = DefSwitch ("blevcorr");
    ccd_sw.biascorr = DefSwitch ("biascorr");
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
        strcpy(osuffix, "_blv_tmp");

        if (MkName (input, isuffix, osuffix, "", output, CHAR_LINE_LENGTH)) {
            WhichError (status);
            sprintf (MsgText, "Skipping %s", input);
            trlmessage (MsgText);
            continue;
        }

        /* Calibrate the current input file. */
        if (ACSccd (input, output, &ccd_sw, &refnames, printtime, verbose)) {
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
