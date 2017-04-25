/* wf3ccd -- basic CCD image reduction

   This file contains:
   WF3ccd
 */

# include <time.h>
# include <string.h>

# include "hstio.h"

# include "wf3.h"
# include "wf3info.h"
# include "hstcalerr.h"
# include "wf3corr.h"		/* calibration switch names */

static int BiasKeywords (WF3Info *);

/* Do basic CCD calibration.

   Warren Hack, 1998 July 28:
   Revised for ACS from CALSTIS1. 
   
   Howard Bushouse, 2000 Aug 29:
   Initial WFC3 version.
   
   H. Bushouse, 2001 May 7:
   Revised to support post-flash processing (adopted equivalent
   calacs changes).
   
   H. Bushouse, 2003 Oct 24:
   Added BiasKeywords routine to populate bias level keywords for
   each amp (CALACS changes).
      
 */

int WF3ccd (char *input, char *output, CCD_Switch *ccd_sw,
        RefFileInfo *refnames, int printtime, int verbose) {

    extern int status;

    WF3Info wf3;	/* calibration switches, reference files, etc */
    int extver;

    Hdr phdr;	/* primary header for input image */

    int DoCCD (WF3Info *, int);
    int FileExists (char *);
    int GetFlags (WF3Info *, Hdr *);
    int GetKeys (WF3Info *, Hdr *);
    void TimeStamp (char *, char *);
    void PrBegin (char *);
    void PrEnd (char *);
    void PrFileName (char *, char *);
    void PrHdrInfo (char *, char *, char *);
    void PrGrpBegin (char *, int);
    void PrGrpEnd (char *, int);
    int LoadHdr (char *, Hdr *);
    void InitCCDTrl (char *, char *);
    void WF3Init (WF3Info *);
    int MkName (char *, char *, char *, char *, char *, int);
    /* ----------------------- Start Code --------------------------------*/

    /* Determine the names of the trailer files based on the input
       and output file names, then initialize the trailer file buffer
       with those names.
     */
    PrBegin ("WF3CCD");

    InitCCDTrl (input, output);
    /* If we had a problem initializing the trailer files, quit... */
    if (status != WF3_OK) 
        return (status);

    if (printtime)
        TimeStamp ("WF3CCD started", "");
        
    /* Initialize structure containing calwf3 information. */
    WF3Init (&wf3);

    /* If we have IR data, do not even proceed here... */
    if (wf3.detector == IR_DETECTOR) {
        trlerror ("Can NOT process IR data with WF3CCD...");
        freeHdr (&phdr);
        return (status = ERROR_RETURN);
    }

    /* Copy command-line arguments into wf3. */
    /* Start by making sure input name is a full filename... 
       
      The input can either be _raw or _rac_tmp or rootname 
        or rootname+unexpected, check for all. This seems
        like rather strange logic, but it works.
    
    */    
    if (strstr(input,"_raw")){
        if (MkName (input, "_raw", "", "", wf3.input, SZ_LINE)) {
            strcpy(wf3.input,input);
            strcat(wf3.input,"_raw.fits");
        } else {
            strcpy(wf3.input,input);
        }
    } else if (strstr(input,"_rac_tmp")){
        if (MkName (input, "_rac_tmp", "", "", wf3.input, SZ_LINE)) {
            strcpy(wf3.input,input);
            strcat(wf3.input,"_rac_tmp.fits");
        } else {
            strcpy(wf3.input,input);
        }
        
    } else {
        strcpy(wf3.input,input);
        strcat(wf3.input,"_raw.fits");
    }    

    /*user specified output*/
    strcpy (wf3.output, output);
    
    wf3.dqicorr  = ccd_sw->dqicorr;
    wf3.atodcorr = ccd_sw->atodcorr;
    wf3.blevcorr = ccd_sw->blevcorr;
    wf3.biascorr = ccd_sw->biascorr;
    wf3.flashcorr = ccd_sw->flashcorr;
    wf3.noiscorr = PERFORM;
    wf3.printtime = printtime;
    wf3.verbose = verbose;
    wf3.refnames = refnames;

    PrFileName ("Input:", wf3.input);
    PrFileName ("Output:", wf3.output);

    /* Check whether the output file already exists. */
    if (FileExists (wf3.output))
        return (status);

    /* Open input image and read its primary header. */
    if (LoadHdr (wf3.input, &phdr) )
        return (status);

    /* Get keyword values from primary header. */
    if (GetKeys (&wf3, &phdr)) {
        freeHdr (&phdr);
        return (status);
    }


    /* PRINT INFORMATION ABOUT THIS IMAGE. */
    PrHdrInfo (wf3.aperture, wf3.filter, wf3.det);

    /* Get reference file names from input image header.  Pedigree is
       checked, and the calibration switch (an internal flag, not the
       header keyword) will be reset from PERFORM to DUMMY if the
       reference file has pedigree = "DUMMY".  Switches that are
       currently set to PERFORM will be reset to OMIT if the value
       in the header is COMPLETE.
     */
    if (GetFlags (&wf3, &phdr)) {
        freeHdr (&phdr);
        return (status);
    }		

    freeHdr (&phdr);

    /* DO BASIC CCD IMAGE REDUCTION. */
    if (wf3.printtime)
        TimeStamp ("Begin processing", wf3.rootname);
 
    /* PROCESS EACH IMSET (CHIP) IN INPUT FILE */
    for (extver = 1;  extver <= wf3.nimsets;  extver++) {
        trlmessage ("\n");
        PrGrpBegin ("imset", extver);

       if (DoCCD (&wf3, extver))
            return (status);
        PrGrpEnd ("imset", extver);
    }

    
    /* UPDATE THE BIASLEVN KEYWORDS IN THE HEADER. THEY MUST BE UPDATED
     ** HERE BECAUSE ONLY SOME ARE COMPUTED FOR EACH SINGLEGROUP AND
     ** HSTIO WILL ONLY ALLOW ONE UPDATE TO THE PRIMARY HEADER WITH
     ** SINGLEGROUP UPDATES.
     **
     ** WHEN THERE IS NO OVERSCAN TO COMPUTE BIAS LEVELS, ALL VALUES WILL
     ** BE ZERO EXCEPT FOR BLEV[AMP] FOR THE AMP USED FOR THE OBSERVATION. */

    trlmessage("Setting Bias Keywords in header");
    BiasKeywords (&wf3);

    trlmessage ("\n");
    PrEnd ("WF3CCD");

    if (wf3.printtime)
        TimeStamp ("WF3CCD completed", wf3.rootname);

    /* WRITE OUT TEMP TRAILER FILE TO FINAL FILE */
    WriteTrlFile ();

    return (status);
}


void InitCCDTrl (char *input, char *output) {

	extern int status;
	int exist;

	char trl_in[SZ_LINE+1]; 	/* trailer filename for input */
	char trl_out[SZ_LINE+1]; 	/* output trailer filename */
	
	int MkOutName (char *, char **, char **, int, char *, int);
	int MkNewExtn (char *, char *);
	void WhichError (int);
	int TrlExists (char *);
	void SetTrlOverwriteMode (int);

	/* Input and output suffixes. */
	char *isuffix[] = {"_raw", "_rac_tmp"};
	char *osuffix[] = {"_blv_tmp", "_blc_tmp"};
	char *trlsuffix[] = {"", ""};

	int nsuffix = 2;
	
	/* Initialize internal variables */
	trl_in[0] = '\0';
	trl_out[0] = '\0';
	exist = EXISTS_UNKNOWN;

	/* Start by stripping off suffix from input/output filenames */
	if (MkOutName (input, isuffix, trlsuffix, nsuffix, trl_in, SZ_LINE)) {
	    WhichError (status);
	    sprintf (MsgText, "Couldn't determine trailer filename for %s",
		     input);
	    trlmessage (MsgText);
	}
	if (MkOutName (output, osuffix, trlsuffix, nsuffix, trl_out, SZ_LINE)) {
	    WhichError (status);
	    sprintf (MsgText, "Couldn't create trailer filename for %s",
		     output);
	    trlmessage (MsgText);
	}

	/* NOW, CONVERT TRAILER FILENAME EXTENSIONS FROM '.FITS' TO '.TRL' */
    
	if (MkNewExtn (trl_in, TRL_EXTN) ) {
	    sprintf (MsgText, "Error with input trailer filename %s", trl_in);
	    trlerror (MsgText);
	    WhichError (status);
	}
	if (MkNewExtn (trl_out, TRL_EXTN) ) {
	    sprintf (MsgText, "Error with output trailer filename %s", trl_out);
	    trlerror (MsgText);
	    WhichError (status);
	}

	/* If we are working with a RAW file, then see if a TRL file 
	   needs to be overwritten after the generic conversion comments.  */
	if (strstr(input, isuffix[0]) != NULL) {
	    /* Test whether the output file already exists */
	    exist = TrlExists(trl_out);            
	    if (exist == EXISTS_YES) {
		/* The output file exists, so we want to add to them 
		** the new trailer comments.  */
		    SetTrlOverwriteMode (NO);	
	    }
	}
    
    

	/* Sets up temp trailer file for output and copies input
		trailer file into it.
	*/
	InitTrlFile (trl_in, trl_out);
	/* Check value of STATUS for errors in calling routine */
}

static int BiasKeywords (WF3Info *wf3) {

    extern int status;

    IODescPtr im;	/* descriptor for output image */
    Hdr phdr;	/* primary header for input image */
    FitsKw key;	/* location of keyword in header */

    int PutKeyFlt (Hdr *, char *, float, char *);

    initHdr (&phdr);

    /* Open input image in order to read its primary header */
    im = openUpdateImage (wf3->output, "", 0, &phdr);
    if (hstio_err()) {
        trlopenerr (wf3->output);
        return (status = OPEN_FAILED);
    }

    key = findKw (&phdr, "BIASLEVA");

    /* Update all bias value keywords for all chips to provide uniform sets
     ** of keywords across all the extensions in multi-chip/multi-extension
     ** data. Only write out non-zero value to avoid overwriting any
     ** previously written BIASLEV keyword values. */

    if (key != NotFound) {

        if (wf3->blev[0] != 0.) {
            if (PutKeyFlt (&phdr, "BIASLEVA", wf3->blev[0],
                        "bias level for amplifier A")) {
                freeHdr (&phdr);
                closeImage (im);
                return (status);
            }
        }

        if (wf3->blev[1] != 0.) {
            if (PutKeyFlt (&phdr, "BIASLEVB", wf3->blev[1],
                        "bias level for amplifier B")) {
                freeHdr (&phdr);
                closeImage (im);
                return (status);
            }
        }

        if (wf3->blev[2] != 0.) {
            if (PutKeyFlt (&phdr, "BIASLEVC", wf3->blev[2],
                        "bias level for amplifier C")) {
                freeHdr (&phdr);
                closeImage (im);
                return (status);
            }
        }

        if (wf3->blev[3] != 0.) {
            if (PutKeyFlt (&phdr, "BIASLEVD", wf3->blev[3],
                        "bias level for amplifier D")) {
                freeHdr (&phdr);
                closeImage (im);
                return (status);
            }
        }

    } else {

        trlmessage ("  No bias level keywords found to be updated.");
        trlmessage ("  Reporting values in trailer file only!");
    }

    /* Write out primary header */
    if (putHeader (im))
        status = HEADER_PROBLEM;
    if (hstio_err() || status) {
        trlreaderr (wf3->output);
        closeImage (im);
        return (status = OPEN_FAILED);
    }

    closeImage (im);
    freeHdr (&phdr);

    return (status);
}


