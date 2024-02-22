/* ACS2d -- basic 2-D image reduction */

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
# include "acsversion.h"
# include "hstcalerr.h"
# include "acscorr.h"		/* calibration switch names for cs1 */
# include "hstcalversion.h"
#include "trlbuf.h"

static void printSyntax(void)
{
    printf("syntax:  acs2d [--help] [-t] [-v] [-q] [--version] [--gitinfo] input output\n");
}
static void printHelp(void)
{
    printSyntax();
}

/* Standard string buffer for use in messages */
char MsgText[MSG_BUFF_LENGTH]; // Global char auto initialized to '\0'
struct TrlBuf trlbuf = { 0 };

/* This is the main module for acs2d.  It gets the input and output
   file names, calibration switches, and flags, and then calls acs2d.

   Warren Hack, 1998 June 9:
   	Initial version for ACS.

   Pey Lian Lim, 2012 Dec 12:
        Moved FLSHCORR from ACSCCD.

   Pey Lian Lim, 2012 Dec 19:
        Added check to see if PCTECORR was performed.

   Pey Lian Lim, 2014 Oct 9:
       Disabled command line keywords. Obtain calibration flags from header.

*/

int main (int argc, char **argv) {

    char *inlist;		/* list of input file names */
    char *outlist;		/* list of output file names */
    /*int switch_on = 0;*/	/* was any switch specified? */
    int printtime = NO;	/* print time after each step? */
    int verbose = NO;	/* print additional info? */
    int quiet = NO;		/* suppress STDOUT messages? */
    int too_many = 0;	/* too many command-line arguments? */
    int i, j;			/* loop indexes */

    IRAFPointer i_imt, o_imt;	/* imt list pointers */
    char *input;		/* name of input science file */
    char *output;		/* optional name of output file */
    int n_in, n_out;	/* number of files in each list */
    int n;

    /* Input and output suffixes. */
    char *isuffix[] = {"_raw", "_blv_tmp", "_blc_tmp", "_crj_tmp", "_crc_tmp"};
    char *osuffix[] = {"_flt", "_flt",     "_flc",     "_crj",     "_crc"};

    int nsuffix = 5;

    /* A structure to pass the calibration switches to acs2d */
    CalSwitch acs2d_sw;

    /* reference file keywords and names */
    RefFileInfo refnames;
    void InitRefFile (RefFileInfo *);
    void FreeRefFile (RefFileInfo *);

    int ACS2d (char *, char *, CalSwitch *, RefFileInfo *, int, int);
    /*int DefSwitch (char *);*/
    int MkOutName (char *, char **, char **, int, char *, int);
    void WhichError (int);
    int CompareNumbers (int, int, char *);
    void initSwitch (CalSwitch *);

    /* For image header access */
    Hdr phdr;
    int pctecorr;
    int LoadHdr (char *, Hdr *);
    int GetSwitch (Hdr *, char *, int *);
    int Get2dSw (CalSwitch *, Hdr *);

    status = 0;
    c_irafinit (argc, argv);

    /* Allocate space for file names. */
    PtrRegister ptrReg;
    initPtrRegister(&ptrReg);
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

    /* Initialize the lists of reference file keywords and names. */
    InitRefFile (&refnames);
    addPtr(&ptrReg, &refnames, &FreeRefFile);

    /* Initial values. */
    initSwitch (&acs2d_sw);

    for (i = 1;  i < argc;  i++) {

        /*********
        if (strcmp (argv[i], "-dqi") == 0) {
            acs2d_sw.dqicorr = PERFORM;
            switch_on = 1;
        } else if (strcmp (argv[i], "-glin") == 0) {
            acs2d_sw.glincorr = PERFORM;
            switch_on = 1;
        } else if (strcmp (argv[i], "-lflg") == 0) {
            acs2d_sw.lflgcorr = PERFORM;
            switch_on = 1;
        } else if (strcmp (argv[i], "-dark") == 0) {
            acs2d_sw.darkcorr = PERFORM;
            switch_on = 1;
        } else if (strcmp (argv[i], "-flash") == 0) {
            acs2d_sw.flashcorr = PERFORM;
            switch_on = 1;
        } else if (strcmp (argv[i], "-flat") == 0) {
            acs2d_sw.flatcorr = PERFORM;
            switch_on = 1;
        } else if (strcmp (argv[i], "-shad") == 0) {
            acs2d_sw.shadcorr = PERFORM;
            switch_on = 1;
        } else if (strcmp (argv[i], "-phot") == 0) {
            acs2d_sw.photcorr = PERFORM;
            switch_on = 1;
        } else if (argv[i][0] == '-') {
        *********/
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
                    printf ("Unrecognized option %s\n", argv[i]);
                    printSyntax();
                    freeOnExit(&ptrReg);
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
        printSyntax();
        /*
        printf ("  command-line switches:\n");
        printf ("       -dqi -glin -lflg -dark\n");
        printf ("       -flash -flat -shad -phot\n");
        */
        freeOnExit(&ptrReg);
        exit (ERROR_RETURN);
    }

    /* Initialize the structure for managing trailer file comments */
    InitTrlBuf ();
    addPtr(&ptrReg, &trlbuf , &CloseTrlBuf);
    trlGitInfo();

    /* Copy command-line value for QUIET to structure */
    SetTrlQuietMode(quiet);

    /* Was no calibration switch specified on command line?
    if (!switch_on) {	default values (mostly PERFORM)
    acs2d_sw.dqicorr   = DefSwitch ("dqicorr");
    acs2d_sw.glincorr  = DefSwitch ("glincorr");
    acs2d_sw.lflgcorr  = DefSwitch ("lflgcorr");
    acs2d_sw.darkcorr  = DefSwitch ("darkcorr");
    acs2d_sw.flashcorr = DefSwitch ("flshcorr");  OMIT
    acs2d_sw.flatcorr  = DefSwitch ("flatcorr");
    acs2d_sw.shadcorr  = DefSwitch ("shadcorr");  OMIT
    acs2d_sw.photcorr  = DefSwitch ("photcorr");
    }*/

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

        j = c_imtgetim (i_imt, input, CHAR_LINE_LENGTH);
        if (n_out > 0)
            j = c_imtgetim (o_imt, output, CHAR_LINE_LENGTH);
        else
            output[0] = '\0';

        /* Open input image in order to read its primary header. */
        if (LoadHdr (input, &phdr)) {
            WhichError (status);
            sprintf (MsgText, "Skipping %s", input);
            trlmessage (MsgText);
            continue;
        }

        /* Get the values for the Calibration Switches from the
           header for processing. */
        if (Get2dSw (&acs2d_sw, &phdr) ) {
            WhichError (status);
            sprintf (MsgText, "Skipping %s", input);
            trlmessage (MsgText);
            continue;
        }

        /* Set PCTECORR flag.
           PERFORM:   Will use DRKCFILE.
           Otherwise: Will use DARKFILE.
        */
        if (GetSwitch (&phdr, "PCTECORR", &pctecorr)) {
            WhichError (status);
            sprintf (MsgText, "Skipping %s", input);
            trlmessage (MsgText);
            continue;
        }
        if (pctecorr == COMPLETE)
            acs2d_sw.pctecorr = PERFORM;
        else
            acs2d_sw.pctecorr = OMIT;

        if (MkOutName (input, isuffix, osuffix, nsuffix, output, CHAR_LINE_LENGTH)) {
            WhichError (status);
            sprintf (MsgText, "Skipping %s", input);
            trlmessage (MsgText);
            continue;
        }

        /* Calibrate the current input file. */
        if (ACS2d (input, output, &acs2d_sw, &refnames, printtime, verbose)) {
            sprintf (MsgText, "Error processing %s.", input);
            trlerror (MsgText);
            WhichError (status);
        }
    }


    freeOnExit(&ptrReg);

    if (status)
        exit (ERROR_RETURN);
    else
        exit (ACS_OK);
}
