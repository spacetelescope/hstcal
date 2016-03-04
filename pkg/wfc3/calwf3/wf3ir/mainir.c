/* WF3IR -- basic IR image reduction */

# include <stdio.h>
# include <stdlib.h>		/* calloc */
# include <time.h>
# include <string.h>

int status = 0;			/* zero is OK */

# include "c_iraf.h"		/* for c_irafinit */
# include "ximio.h"
# include "hstio.h"

# include "wf3.h"
# include "wf3info.h"
# include "wf3err.h"
# include "wf3corr.h"        /* calibration switch names for wf3ir */
# include "wf3version.h"

static void FreeNames (char *, char *, char *, char *);

/* This is the main module for WF3IR.  It gets the input and output
   file names, calibration switches, and flags, and then calls WF3ir.

   Warren Hack, 1998 June 1:
   	Revised for ACS... It will work on 1 ACS image at a time. Based on 
	original CALSTIS code by Phil Hodge.
   Howard Bushouse, 2001 Apr 18:
	Copied from CALWF3 mainccd.c for use in CALWF3 IR branch.
   H.Bushouse, 2002 June 21:
	Removed use of statcorr (always perform doStat).
*/

int main (int argc, char **argv) {

	char *inlist;		/* input file name */
	char *outlist;		/* output blev file name */
	int switch_on = 0;	/* was any switch specified? */
	int printtime = NO;	/* print time after each step? */
	int verbose = NO;	/* print additional info? */
	int quiet = NO;		/* print additional info? */
	int too_many = 0;	/* too many command-line arguments? */
	int i, j;		/* loop indexes */

	IRAFPointer i_imt, o_imt;	/* imt list pointers */
	char *input;		/* name of input science file */
	char *output;		/* name of output file */
	int n_in, n_out;	/* number of files in each list */
	int n;

	/* Input and output suffixes. */
	char isuffix[] = "_raw";
	char osuffix[] = "_flt";

	/* A structure to pass the calibration switches to WF3IR */
	IR_Switch ir_sw;

	/* reference file keywords and names */
	RefFileInfo refnames;
		
	void InitRefFile (RefFileInfo *);
	void FreeRefFile (RefFileInfo *);
	void initIRSwitches (IR_Switch *);

	int WF3ir (char *, char *, IR_Switch *, RefFileInfo *, int, int,
		   float, float, int);
	int DefSwitch (char *);
	int MkName (char *, char *, char *, char *, char *, int);
	void WhichError (int);
	int CompareNumbers (int, int, char *);

/*===========================================================================*/

	/* Initialize IRAF interface environment */
	c_irafinit (argc, argv);

	/* Post HSTIO error handler */
	push_hstioerr (errchk);

	/* Allocate space for file names. */
	inlist  = calloc (SZ_LINE+1, sizeof (char));
	outlist = calloc (SZ_LINE+1, sizeof (char));
	input   = calloc (SZ_LINE+1, sizeof (char));
	output  = calloc (SZ_LINE+1, sizeof (char));

	if (inlist == NULL || outlist == NULL ||
	    input == NULL || output == NULL) {
	    printf ("Can't even begin; out of memory.\n");
	    exit (ERROR_RETURN);
	}
	inlist[0]  = '\0';
	outlist[0] = '\0';
	input[0]   = '\0';
	output[0]  = '\0';
	/* Initialize the lists of reference file keywords and names. */
	InitRefFile (&refnames);

	/* Initial values. */
	initIRSwitches (&ir_sw);

	/* Parse the command-line arguments */
	for (i = 1;  i < argc;  i++) {

	    if (!(strcmp(argv[i],"--version"))) {
		printf("%s\n",WF3_CAL_VER_NUM);
		exit(0);
	    }

	    if (argv[i][0] == '-') {
		for (j = 1;  argv[i][j] != '\0';  j++) {
		    if (argv[i][j] == 't') {
			printtime = YES;
		    } else if (argv[i][j] == 'v') {
			verbose = YES;
		    } else if (argv[i][j] == 'q') {
			quiet = YES;
            } else if (argv[i][j] == 'r'){
                printf ("Current version: %s\n", WF3_CAL_VER);
                exit(0);
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
	    printf ("syntax:  wf3ir [-t] [-v] [-q] [-r] input output\n");
	    printf ("  command-line switches:\n");
	    printf ("       -bseq -pede\n");
	    FreeNames (inlist, outlist, input, output);
	    exit (ERROR_RETURN);
	}

	/* Initialize the structure for managing trailer file comments */
	InitTrlBuf ();
	
	/* Copy command-line value for QUIET to structure */
	SetTrlQuietMode(quiet);

	/* Was no calibration switch specified on command line? */
	if (!switch_on) {	/* default values (mostly PERFORM) */
	    ir_sw.zsigcorr = DefSwitch ("zsigcorr");
	    ir_sw.zoffcorr = DefSwitch ("zoffcorr");
	    ir_sw.dqicorr  = DefSwitch ("dqicorr");
	    ir_sw.blevcorr = DefSwitch ("blevcorr");
	    ir_sw.noiscorr = DefSwitch ("noiscorr");
	    ir_sw.darkcorr = DefSwitch ("darkcorr");
	    ir_sw.nlincorr = DefSwitch ("nlincorr");
	    ir_sw.flatcorr = DefSwitch ("flatcorr");
	    ir_sw.unitcorr = DefSwitch ("unitcorr");
	    ir_sw.photcorr = DefSwitch ("photcorr");
	    ir_sw.crcorr   = DefSwitch ("crcorr");
	}

	/* Expand the templates. */
	i_imt = c_imtopen (inlist);
	o_imt = c_imtopen (outlist);
	n_in  = c_imtlen (i_imt);
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

	    i = c_imtgetim (i_imt, input, SZ_LINE);
		
	    if (n_out > 0)
		i = c_imtgetim (o_imt, output, SZ_LINE);
	    else {
		    output[0] = '\0';

	        if (MkName (input, isuffix, osuffix, "", output, SZ_LINE)) {
		        WhichError (status);
		        sprintf (MsgText, "Skipping %s", input);
		        trlmessage (MsgText);
		        continue;
	        }
        }

	    /* Calibrate the current input file. */
	    if (WF3ir (input, output, &ir_sw, &refnames, printtime, verbose,
		       0.0, 0.0, 0)){
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
