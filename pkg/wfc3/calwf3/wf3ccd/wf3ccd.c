/* wf3ccd -- basic CCD image reduction

   This file contains:
   WF3ccd
 */

# include <time.h>
# include <string.h>

# include "hstio.h"

# include "wf3.h"
# include "wf3info.h"
# include "wf3err.h"
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
    int SinkDetect (WF3Info *);
    /* ----------------------- Start Code --------------------------------*/

    /* Determine the names of the trailer files based on the input
       and output file names, then initialize the trailer file buffer
       with those names.
     */
    PrBegin ("WF3CCD");


    /* Initialize structure containing calwf3 information. */
    WF3Init (&wf3);

    /* Copy command-line arguments into wf3. */
    /* Start by making sure input name is a full filename... 
      
      The input can either be _raw or _rac or rootname 
        or rootname+unexpected, check for all
    
    */    
    /*if the input doesn't end _rac or _raw add it, this covers rootname user input*/
    strcpy(wf3.input,input);
    
    /*user specified output*/
    strcpy (wf3.output, output);
    
    InitCCDTrl (input, output);


    if (printtime)
        TimeStamp ("WF3CCD started", "");

    /* If we had a problem initializing the trailer files, quit... */
    if (status != WF3_OK) 
        return (status);

    

    wf3.dqicorr  = ccd_sw->dqicorr;
    wf3.atodcorr = ccd_sw->atodcorr;
    wf3.blevcorr = ccd_sw->blevcorr;
    wf3.biascorr = ccd_sw->biascorr;
    wf3.flashcorr = ccd_sw->flashcorr;
    wf3.noiscorr = PERFORM;
    wf3.printtime = printtime;
    wf3.verbose = verbose;

    /* For debugging...
       wf3.dqicorr  = PERFORM;
       wf3.atodcorr = PERFORM;
       wf3.blevcorr = PERFORM;
       wf3.biascorr = PERFORM;
       wf3.noiscorr = PERFORM;
       wf3.printtime = 1;
       wf3.verbose = 1;
     */
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

    /* If we have IR data, do not even proceed here... */
    if (wf3.detector == IR_DETECTOR) {
        trlerror ("Can NOT process IR data with WF3CCD...");
        freeHdr (&phdr);
        return (status = ERROR_RETURN);
    }

    /* Print information about this image. */
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

    /* Do basic CCD image reduction. */
    if (wf3.printtime)
        TimeStamp ("Begin processing", wf3.rootname);
 
    /* Process each imset (chip) in input file */
    for (extver = 1;  extver <= wf3.nimsets;  extver++) {
        trlmessage ("\n");
        PrGrpBegin ("imset", extver);

       if (DoCCD (&wf3, extver))
            return (status);
        PrGrpEnd ("imset", extver);
    }


    /*Update the SINK pixels in the DQ mask of both science image sets
     It's done here with one call to the file because they need to be
     processed in the RAZ format Jay uses
    */
     
     if (wf3.dqicorr) {
        if (SinkDetect(&wf3))
            return(status);          
     }
    
    /* Update the BIASLEVn keywords in the header. They must be updated
     ** here because only some are computed for each SingleGroup and
     ** HSTIO will only allow one update to the Primary header with
     ** SingleGroup updates.
     **
     ** When there is no overscan to compute bias levels, all values will
     ** be zero except for blev[amp] for the amp used for the observation. */

    BiasKeywords (&wf3);

    trlmessage ("\n");
    PrEnd ("WF3CCD");

    if (wf3.printtime)
        TimeStamp ("WF3CCD completed", wf3.rootname);

    /* Write out temp trailer file to final file */
    WriteTrlFile ();

    return (status);
}


void InitCCDTrl (char *input, char *output) {

    extern int status;

    char trl_in[SZ_LINE+1]; 	/* trailer filename for input */
    char trl_out[SZ_LINE+1]; 	/* output trailer filename */
    int exist;

    int MkName (char *, char *, char *, char *, char *, int);
    void WhichError (int);
    int TrlExists (char *);
    void SetTrlOverwriteMode (int);

    /* Initialize internal variables */
    trl_in[0] = '\0';
    trl_out[0] = '\0';
    exist = EXISTS_UNKNOWN;

    /* Start by stripping off suffix from input/output filenames */
    if (strcmp(input,"_raw") == 1){
        if (MkName (input, "_raw", "", TRL_EXTN, trl_in, SZ_LINE)) {
            WhichError (status);
            sprintf (MsgText, "Couldn't determine trailer filename for %s",
                    input);
            trlmessage (MsgText);
        }
    } else {
        if (strcmp(input,"_rac") == 1){
            if (MkName (input, "_rac", "", TRL_EXTN, trl_in, SZ_LINE)) {
                WhichError (status);
                sprintf (MsgText, "Couldn't determine trailer filename for %s",
                        input);
                trlmessage (MsgText);
            }
        }
    }    
    
    
    if (strcmp(input,"_raw") == 1){

        if (MkName (output, "_blv_tmp", "", TRL_EXTN, trl_out, SZ_LINE)) {
            WhichError (status);
            sprintf (MsgText, "Couldn't create trailer filename for %s",
                    output);
            trlmessage (MsgText);
        }
    } else {
        if (strcmp(input,"_rac") == 1){

            if (MkName (output, "_blc_tmp", "", TRL_EXTN, trl_out, SZ_LINE)) {
                WhichError (status);
                sprintf (MsgText, "Couldn't create trailer filename for %s",
                        output);
                trlmessage (MsgText);
            }
        }
    }
    
    
    /* Test whether the output file already exists */
    exist = TrlExists(trl_out);
    if (exist == EXISTS_YES) {
        /* The output file exists, so we want to overwrite them with
         ** the new trailer comments.  */
        SetTrlOverwriteMode (YES);	
    }

    /* Sets up temp trailer file for output and copies input
     ** trailer file into it.  */
    InitTrlFile (trl_in, trl_out);
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


