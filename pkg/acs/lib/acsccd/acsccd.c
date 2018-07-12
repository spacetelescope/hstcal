/* acsccd -- basic CCD image reduction

 This file contains:
 ACSccd
 */

# include <time.h>
# include <string.h>

#include "hstcal.h"
# include "hstio.h"

# include "acs.h"
# include "acsinfo.h"
# include "hstcalerr.h"
# include "acscorr.h"		/* calibration switch names */
# include "trlbuf.h"
# include "getacskeys.h"

static int BiasKeywords (ACSInfo *);
void InitCCDTrl (char *, char *);


/* Do basic CCD calibration.

 Warren Hack, 1998 July 28:
 Revised for ACS...

 Warren Hack, 2000 Sept 11:
 Revised to support post-flash processing.

 Pey Lian Lim, 2012, Dec 12:
 Moved FLSHCORR to ACS2D.

 **	Not much modification from the basic outline
 **		followed in CALSTIS1...
 **
 */
int ACSccd (char *input, char *output, CalSwitch *ccd_sw,
            RefFileInfo *refnames, int printtime, int verbose) {

    extern int status;

    ACSInfo acs;	/* calibration switches, reference files, etc */

    Hdr phdr;		/* primary header for input image */

    int DoCCD (ACSInfo *);
    int FileExists (char *);
    int GetACSFlags (ACSInfo *, Hdr *);
    void TimeStamp (char *, char *);
    void PrBegin (char *);
    void PrEnd (char *);
    void PrFileName (char *, char *);
    void PrHdrInfo (char *, char *, char *, char *);
    void PrGrpBegin (char *, int);
    void PrGrpEnd (char *, int);
    int LoadHdr (char *, Hdr *);
    void ACSInit (ACSInfo *);
    int MkName (char *, char *, char *, char *, char *, int);

    int GetACSGrp (ACSInfo *, Hdr *);
    int GetCCDTab (ACSInfo *, int, int);

    /* ----------------------- Start Code --------------------------------*/

    /* Determine the names of the trailer files based on the input
       and output file names, then initialize the trailer file buffer
       with those names.
    */
    InitCCDTrl (input, output);
    /* If we had a problem initializing the trailer files, quit... */
    if (status != ACS_OK)
        return (status);

    PrBegin ("ACSCCD");

    if (printtime)
        TimeStamp ("ACSCCD started", "");

    /* Initialize structure containing calacs information. */
    ACSInit (&acs);

    /* Copy command-line arguments into acs. */
    /* Start by making sure input name is a full filename... */
    if (MkName (input, "_raw", "_raw", "", acs.input, CHAR_LINE_LENGTH) ) {
        strcpy(acs.input, input);
        strcat (acs.input,"_raw.fits");
    }
    strcpy (acs.output, output);

    acs.dqicorr  = ccd_sw->dqicorr;
    acs.atodcorr = ccd_sw->atodcorr;
    acs.blevcorr = ccd_sw->blevcorr;
    acs.biascorr = ccd_sw->biascorr;
    acs.noisecorr = PERFORM;
    acs.printtime = printtime;
    acs.verbose = verbose;

    /* For debugging...
    acs.dqicorr  = PERFORM;
    acs.atodcorr = PERFORM;
    acs.blevcorr = PERFORM;
    acs.biascorr = PERFORM;
    acs.noisecorr = PERFORM;
    acs.printtime = 1;
    acs.verbose = 1;
    */

    acs.refnames = refnames;

    PrFileName ("input", acs.input);
    PrFileName ("output", acs.output);

    /* Check whether the output file already exists. */
    if (FileExists (acs.output))
        return (status);

    /* Open input image in order to read its primary header. */
    if (LoadHdr (acs.input, &phdr) )
        return (status);

    /* Get keyword values from primary header. */
    if (getAndCheckACSKeys (&acs, &phdr)) {
        freeHdr (&phdr);
        return (status);
    }

    /* If we have MAMA data, do not even proceed here... */
    if (acs.detector == MAMA_DETECTOR) {
        /* Return ACS_OK, since processing can proceed, just with a
           different function */
        trlwarn ("Can NOT process MAMA data with ACSCCD...");
        freeHdr (&phdr);
        return (status);
    }
    /* Print information about this image. */
    PrHdrInfo (acs.aperture, acs.filter1, acs.filter2, acs.det);

    /* Get reference file names from input image header.  Pedigree is
       checked, and the calibration switch (an internal flag, not the
       header keyword) will be reset from PERFORM to DUMMY if the
       reference file has pedigree = "DUMMY".  Switches that are
       currently set to PERFORM will be reset to OMIT if the value
       in the header is COMPLETE.
    */
    if (GetACSFlags (&acs, &phdr)) {
        freeHdr(&phdr);
        return (status);
    }

    freeHdr(&phdr);

     /* Do basic CCD image reduction. */

    if (acs.printtime) {
        TimeStamp("Begin processing", acs.rootname);
    }

    if (DoCCD(&acs)) {
        return (status);
    }

    /*  Commented out until the new keywords can be populated in the headers
        by Generic Conversion. WJH 6 Mar 03.  Enabled: 6-June-03 WJH

        When there is no overscan to compute bias levels, all values will be
        zero except for blev[amp] for the amp used for the observation.

        The BIASLEV keywords need to be updated here since only some are
        computed for each SingleGroup, and HSTIO will only allow one
        update to the Primary header with SingleGroup updates.
    */
    BiasKeywords (&acs);

    trlmessage("\n");
    PrEnd("ACSCCD");

    if (acs.printtime) {
        TimeStamp("ACSCCD completed", acs.rootname);
    }

    /* Write out temp trailer file to final file */
    WriteTrlFile ();

    return (status);
}


void InitCCDTrl (char *input, char *output) {

    extern int status;

    char trl_in[CHAR_LINE_LENGTH+1]; 	/* trailer filename for input */
    char trl_out[CHAR_LINE_LENGTH+1]; 	/* output trailer filename */
    int exist;

    char isuffix[] = "_raw";
    char osuffix[] = "_blv_tmp";
    char trlsuffix[] = "";

    int MkName (char *, char *, char *, char *, char *, int);
    void WhichError (int);
    int TrlExists (char *);

    /* Initialize internal variables */
    trl_in[0] = '\0';
    trl_out[0] = '\0';
    exist = EXISTS_UNKNOWN;

    /* Start by stripping off suffix from input/output filenames */
    if (MkName (input, isuffix, trlsuffix, TRL_EXTN, trl_in, CHAR_LINE_LENGTH)) {
        WhichError (status);
        sprintf (MsgText, "Couldn't determine trailer filename for %s", input);
        trlmessage (MsgText);
    }
    if (MkName (output, osuffix, trlsuffix, TRL_EXTN, trl_out, CHAR_LINE_LENGTH)) {
        WhichError (status);
        sprintf (MsgText, "Couldn't create trailer filename for %s", output);
        trlmessage (MsgText);
    }

    /* Sets up temp trailer file for output and copies input
       trailer file into it.
    */
    InitTrlFile (trl_in, trl_out);
}


static int BiasKeywords (ACSInfo *acs) {

    extern int status;

    IODescPtr im;		/* descriptor for output image */
    int PutKeyFlt (Hdr *, char *, float, char *);
    Hdr phdr;		/* primary header for input image */

    FitsKw key;		/* location of keyword in header */

    initHdr (&phdr);
    /* Open input image in order to read its primary header. */
    im = openUpdateImage (acs->output, "", 0, &phdr);

    if (hstio_err()) {
        trlopenerr (acs->output);
        return (status = OPEN_FAILED);
    }

    key = findKw (&phdr, "BIASLEVA");
    if (key != NotFound) {
        /* Update all bias value keywords for all chips to provide uniform
           sets of keywords acrossl all the extensions in
           multi-chip/multi-extension data.
           Only write out non-zero values to avoid overwriting any previously
           written BIASLEV keyword values. WJH 18-June-2003
        */
        if (acs->blev[0] != 0.) {
            if (PutKeyFlt (&phdr, "BIASLEVA", acs->blev[0],
                           "mean of amp's bias level")) {
                freeHdr (&phdr);
                closeImage (im);
                return (status);
            }
        }
        if (acs->blev[1] != 0.) {
            if (PutKeyFlt (&phdr, "BIASLEVB", acs->blev[1],
                           "mean of amp's bias level")) {
                freeHdr (&phdr);
                closeImage (im);
                return (status);
            }
        }
        if (acs->blev[2] != 0.) {
            if (PutKeyFlt (&phdr, "BIASLEVC", acs->blev[2],
                           "mean of amp's bias level")) {
                freeHdr (&phdr);
                closeImage (im);
                return (status);
            }
        }
        if (acs->blev[3] != 0.) {
            if (PutKeyFlt (&phdr, "BIASLEVD", acs->blev[3],
                           "mean of amp's bias level")) {
                freeHdr (&phdr);
                closeImage (im);
                return (status);
            }
        }
    } else {
        trlmessage ("  NO bias level keywords found to be updated.");
        trlmessage ("  Reporting values in trailer file only!");
    }

    /* write out primary header */
    if (putHeader (im) )
        status = HEADER_PROBLEM;
    if (hstio_err() || status) {
        trlreaderr (acs->output);
        closeImage (im);
        return (status = OPEN_FAILED);
    }

    closeImage (im);
    freeHdr (&phdr);

    return(status);
}
