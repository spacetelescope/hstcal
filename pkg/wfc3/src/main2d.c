/* WF32d -- basic 2-D image reduction */

# include <stdio.h>
# include <stdlib.h>		/* calloc */
# include <time.h>
# include <string.h>

extern int status;

#include "hstcal_memory.h"
#include "hstcal.h"
# include "c_iraf.h"		/* for c_irafinit */
# include "ximio.h"
# include "hstio.h"

# include "wf3.h"
# include "wf3info.h"
# include "hstcalerr.h"
# include "wf3corr.h"		/* calibration switch names for cs1 */
# include "wf3version.h"
# include "hstcalversion.h"
# include "trlbuf.h"

static void printSyntax(void)
{
    printf ("syntax:  wf32d [--help] [-t] [-v] [-q] [-r] [--version] [--gitinfo] input output\n");
    printf ("  command-line switches:\n");
    printf ("       -dqi  -atod\n");
    printf ("       -dark -flat -shad -phot -stat\n");
}
static void printHelp(void)
{
    printSyntax();
}

/* Standard string buffer for use in messages */
char MsgText[MSG_BUFF_LENGTH]; // Global char auto initialized to '\0'

/* This is the main module for wf32d.  It gets the input and output
   file names, calibration switches, and flags, and then calls wf32d.

   Warren Hack, 1998 June 9:
   	Initial version for ACS.
   Howard Bushouse, 2000 Aug 29:
	Initial version for WFC3.
   H.Bushouse, 2002 June 17:
	Remove use of "statcorr" switch - statistics always computed by
	default (in accordance with CALACS changes).
  M Sosey, 2012 December 27:
      Updated to account for a memory leak on linux machines during BuildDth 
      when RPTCORR is off and a new spt is being constructed (#967)       
*/

int main (int argc, char **argv) {

	char *inlist;		/* list of input file names */
	char *outlist;		/* list of output file names */
	int switch_on = 0;	/* was any switch specified? */
	int printtime = NO;	/* print time after each step? */
	int verbose = NO;	/* print additional info? */
	int quiet = NO;		/* suppress STDOUT messages? */
	int too_many = 0;	/* too many command-line arguments? */
	int i, j;		/* loop indexes */
    
	IRAFPointer i_imt, o_imt;	/* imt list pointers */
	char *input;		/* name of input science file */
	char *output;		/* optional name of output file */
	int n_in, n_out;	/* number of files in each list */
	int n;

	/* Input and output suffixes. */
	char *isuffix[] = {"_raw", "_blv_tmp", "_crj_tmp",  "_blc_tmp","_rac_tmp", "_crc_tmp"};
	char *osuffix[] = {"_flt", "_flt",     "_crj",      "_flc",    "_flc", "_crc"};

	int nsuffix = 6;

	/* A structure to pass the calibration switches to wf32d */
	CCD_Switch wf32d_sw;

	/* reference file keywords and names */
	RefFileInfo refnames;
	void InitRefFile (RefFileInfo *);
	void FreeRefFile (RefFileInfo *);

	int WF32d (char *, char *, CCD_Switch *, RefFileInfo *, int, int);
	int DefSwitch (char *);
	int MkOutName (const char *, char **, char **, int, char *, int);
 	void WhichError (int);
	int CompareNumbers (int, int, char *);
	void initCCDSwitches (CCD_Switch *);

/*===========================================================================*/
    status = 0;

	/* Initialize IRAF interface environment */
	c_irafinit (argc, argv);
	
	/* Post HSTIO error handler */
	push_hstioerr (errchk);

    PtrRegister ptrReg;
    initPtrRegister(&ptrReg);
	/* Allocate space for file names. */
	inlist  = calloc (CHAR_LINE_LENGTH+1, sizeof (char));
    addPtr(&ptrReg, inlist, &free);
	outlist = calloc (CHAR_LINE_LENGTH+1, sizeof (char));
    addPtr(&ptrReg, outlist, &free);
	input   = calloc (CHAR_LINE_LENGTH+1, sizeof (char));
    addPtr(&ptrReg, input, &free);
	output  = calloc (CHAR_LINE_LENGTH+1, sizeof (char));
    addPtr(&ptrReg, output, &free);
	if (inlist == NULL || outlist == NULL ||
		input == NULL || output == NULL) {
	    printf ("Can't even begin; out of memory.\n");
	    freeOnExit(&ptrReg);
	    exit (ERROR_RETURN);
	}
        
	/* Initialize the lists of reference file keywords and names. */
	InitRefFile (&refnames);
	addPtr(&ptrReg, &refnames, &FreeRefFile);

	/* Initial values. */
	initCCDSwitches (&wf32d_sw);

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
	    if (strcmp (argv[i], "-dqi") == 0) {	/* turn on */
		wf32d_sw.dqicorr = PERFORM;
		switch_on = 1;
	    } else if (strcmp (argv[i], "-dark") == 0) {
		wf32d_sw.darkcorr = PERFORM;
		switch_on = 1;
	    } else if (strcmp (argv[i], "-flat") == 0) {
		wf32d_sw.flatcorr = PERFORM;
		switch_on = 1;
	    } else if (strcmp (argv[i], "-shad") == 0) {
		wf32d_sw.shadcorr = PERFORM;
		switch_on = 1;
	    } else if (strcmp (argv[i], "-phot") == 0) {
		wf32d_sw.photcorr = PERFORM;
		switch_on = 1;
	    } else if (argv[i][0] == '-') {
		for (j = 1;  argv[i][j] != '\0';  j++) {
		    if (argv[i][j] == 't') {
				printtime = YES;
		    } else if (argv[i][j] == 'v') {
				verbose = YES;
		    } else if (argv[i][j] == 'q') {
				quiet = YES;
            } else if (argv[i][j] == 'r'){
              printf ("Current version: %s\n", WF3_CAL_VER);
              freeOnExit(&ptrReg);
              exit(0);
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
	    freeOnExit(&ptrReg);
	    exit (ERROR_RETURN);
	}
	
	/* Initialize the structure for managing trailer file comments */
	InitTrlBuf ();
	addPtr(&ptrReg, &trlbuf, &CloseTrlBuf);
    trlGitInfo();

	/* Copy command-line value for QUIET to structure */
	SetTrlQuietMode(quiet);

	/* Was no calibration switch specified on command line? */
	if (!switch_on) {	/* default values (mostly PERFORM) */
	    wf32d_sw.dqicorr  = DefSwitch ("dqicorr");
	    wf32d_sw.darkcorr = DefSwitch ("darkcorr");
	    wf32d_sw.flatcorr = DefSwitch ("flatcorr");
	    wf32d_sw.shadcorr = DefSwitch ("shadcorr");
	    wf32d_sw.photcorr = DefSwitch ("photcorr");
	}

	/* Expand the templates. */
	i_imt = c_imtopen (inlist);
	addPtr(&ptrReg, i_imt, &c_imtclose);
	o_imt = c_imtopen (outlist);
	addPtr(&ptrReg, o_imt, &c_imtclose);
	n_in  = c_imtlen (i_imt);
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

	    if (MkOutName (input, isuffix, osuffix, nsuffix, output, CHAR_LINE_LENGTH)) {
		WhichError (status);
		trlmessage("Skipping %s", input);
		continue;
	    }
										
	    /* Calibrate the current input file. */
	    if (WF32d(input, output, &wf32d_sw, &refnames, printtime, verbose)){
		trlerror("Error processing %s.", input);
		WhichError (status);
	    }
	}


	freeOnExit(&ptrReg);

	if (status)
	    exit (ERROR_RETURN);
	else
	    exit (WF3_OK);
}
