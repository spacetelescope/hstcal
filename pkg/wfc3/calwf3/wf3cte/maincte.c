/* WFC3CTE -- CTE loss correction 

This is the routine for running the CTE correction standalone,
as you would wf3ccd or wf32d

MLS 2015


*/

# include <stdio.h>
# include <stdlib.h>		/* calloc */
# include <time.h>
# include <string.h>


# include <c_iraf.h>		/* for c_irafinit */
# include "ximio.h"
# include "hstio.h"

# include "wf3.h"
# include "wf3err.h"
# include "wf3corr.h"		/* calibration switch names for WFC3ccd */
# include "wf3version.h"

# ifdef _OPENMP
#  include <omp.h>
# endif

static void FreeNames (char *, char *, char *, char *);
void FreeRefFile (RefFileInfo *);
void InitRefFile (RefFileInfo *);
int WF3cte (char *, char *, CCD_Switch *, RefFileInfo *, int, int, int);
int MkName (char *, char *, char *, char *, char *, int);
void WhichError (int);
int CompareNumbers (int, int, char *);
int LoadHdr (char *, Hdr *);
int GetSwitch (Hdr *, char *, int *);
void initCCDSwitches (CCD_Switch *);

/* 

This is the main module for WF3CTE.  It gets the input and output
file names, calibration switches, and flags, and then calls WFC3cte.
This is necessary for the task to act standalone.

*/

int main (int argc, char **argv) {

    char *inlist;		/* input blv_tmp file name */
    char *outlist;		/* output blc_tmp file name */
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
    
	/* Initialize status to OK and MsgText to null */
	status     = WF3_OK;

    /* A structure to pass the calibration switches to WFC3CTE */
    CCD_Switch cte_sw;
    RefFileInfo refnames;

    /* For image header access */
    Hdr phdr;


    c_irafinit (argc, argv);
    push_hstioerr(errchk);

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

    
    /*INITIALIZE REFERENCE FILE INFORMATION*/
    InitRefFile (&refnames);
    
    /* Initial values. */
 	initCCDSwitches (&cte_sw);
    
    for (i = 1;  i < argc;  i++) {

        if (argv[i][0] == '-') {
            for (j = 1;  argv[i][j] != '\0';  j++) {
                if (argv[i][j] == 't') {
                    printtime = YES;
                } else if (argv[i][j] == 'v') {
                    verbose = YES;
                } else if (argv[i][j] == '1') {
                    onecpu = YES;
				} else if (argv[i][j] == 'r'){
					printf ("Current version: %s\n", WF3_CAL_VER);
					exit(0);
                } else {
                    printf (MsgText, "Unrecognized option %s\n", argv[i]);
                    FreeNames (inlist, outlist, input, output);
                    exit (ERROR_RETURN);
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
        printf ("syntax:  WF3cte [-v] [-1] input output\n");
        FreeNames (inlist, outlist, input, output);
        exit (ERROR_RETURN);
    }
    /* INITIALIZE THE STRUCTURE FOR MANAGING TRAILER FILE COMMENTS */
    InitTrlBuf ();
    /* COPY COMMAND-LINE VALUE FOR QUIET TO STRUCTURE */
    SetTrlQuietMode(quiet);
           
    /* EXPAND THE TEMPLATES. */
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

    /* LOOP OVER THE LIST OF INPUT FILES. */
    for (n = 0;  n < n_in;  n++) {

        k = c_imtgetim (i_imt, input, SZ_FNAME);

        if (n_out > 0)
            k = c_imtgetim (o_imt, output, SZ_FNAME);
        else
            output[0] = '\0';

        /* OPEN INPUT IMAGE IN ORDER TO READ ITS PRIMARY HEADER. */
        if (LoadHdr (input, &phdr)) {
            WhichError (status);
            sprintf (MsgText, "Skipping %s", input);
            trlmessage (MsgText);
            continue;
        }

        /* GET SWITCHES WE CARE ABOUT. */
        if (GetSwitch (&phdr, "PCTECORR", &cte_sw.pctecorr)) {
            WhichError (status);
        }
        if (GetSwitch (&phdr, "BIASCORR", &cte_sw.biascorr)) {
            WhichError (status);
        }
        if (GetSwitch (&phdr, "BLEVCORR", &cte_sw.blevcorr)) {
            WhichError (status);
        }
        if (GetSwitch (&phdr, "DARKCORR", &cte_sw.darkcorr)){
            WhichError (status);
        }

        /*SIMPLE CHECK*/
        if (cte_sw.biascorr == COMPLETE || cte_sw.blevcorr == COMPLETE || cte_sw.darkcorr == COMPLETE){
            sprintf(MsgText,"An uncalibrated, RAW file must be used as input to CTE corr, skipping %s", input);
            trlmessage(MsgText);
            exit(ERROR_RETURN);
            
        } else if (cte_sw.pctecorr) {

            if (MkName (input, "_raw", "_rac", "", output, SZ_FNAME)) {
                WhichError (status);
                sprintf (MsgText, "Skipping %s, problem making output name", input);
                trlmessage (MsgText);
            }

            /* CALIBRATE THE CURRENT INPUT FILE. */
            if (WF3cte (input, output, &cte_sw, &refnames, printtime, verbose,
                        onecpu)) {
                sprintf (MsgText, "Error processing cte for %s", input);
                trlerror (MsgText);
                WhichError (status);
            }
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



