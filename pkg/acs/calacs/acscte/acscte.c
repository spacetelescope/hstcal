/* acscte -- CTE loss correction

 This file contains:
 ACScte
 */

# include <time.h>
# include <string.h>
# include <stdbool.h>
# include <assert.h>

#include "hstcal.h"
# include "hstio.h"

# include "acs.h"
# include "acsinfo.h"
# include "hstcalerr.h"
# include "acscorr.h"		/* calibration switch names */
# include "trlbuf.h"
# include "getacskeys.h"

void InitCTETrl (char * input, char * output, const char * isuffix, const char * osuffix);


/* Do CTE loss correction.

 Pey Lian Lim, 2013, Aug 9:
 Separated PCTECORR from ACSCCD.

 */
int ACScte (char *input, char *output, CalSwitch *cte_sw,
            RefFileInfo *refnames, int printtime, int verbose,
            const unsigned nThreads, const unsigned cteAlgorithmGen, const char * pcteTabNameFromCmd, const bool forwardModelOnly) {

    extern int status;

    ACSInfo acs;	/* calibration switches, reference files, etc */

    Hdr phdr;		/* primary header for input image */

    int DoCTE (ACSInfo *, const bool forwardModelOnly);
    int FileExists (char *);
    int GetCTEFlags (ACSInfo *, Hdr *);
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

    char * isuffix = "_blv_tmp";
    char osuffix[CHAR_FNAME_LENGTH+1];
    if (forwardModelOnly)
        strcpy(osuffix, "_ctefmod");
    else
        strcpy(osuffix, "_blc_tmp");

    InitCTETrl (input, output, isuffix, osuffix);
    /* If we had a problem initializing the trailer files, quit... */
    if (status != ACS_OK)
        return (status);

    PrBegin ("ACSCTE");

    if (printtime)
        TimeStamp ("ACSCTE started", "");

    /* Initialize structure containing calacs information. */
    ACSInit (&acs);

    /* Copy command-line arguments into acs. */
    strcpy (acs.input, input);
    strcpy (acs.output, output);
    strcpy(acs.pcteTabNameFromCmd, pcteTabNameFromCmd);

    acs.pctecorr = cte_sw->pctecorr;
    acs.printtime = printtime;
    acs.verbose = verbose;
    acs.nThreads = nThreads;
    acs.cteAlgorithmGen = cteAlgorithmGen;

    /* For debugging...
    acs.pctecorr = PERFORM;
    acs.printtime = 1;
    acs.verbose = 1;
    acs.onecpu = 0;
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

    if (forwardModelOnly)
    {
        // Only model full IMSETS defined as ("SCI", "ERR", "DQ") extensions
        if ((status = findTotalNumberOfImsets(acs.input, "SCI", &(acs.nimsets))))
            return status;

        if (acs.nimsets < 1)
        {
            sprintf (MsgText, "N IMSETS found = %d; must be at least %d.", acs.nimsets, 1);
            trlerror (MsgText);
            return INVALID_VALUE;
        }
    }

    /* If we have MAMA data, do not even proceed here... */
    if (acs.detector == MAMA_DETECTOR) {
        /* Return ACS_OK, since processing can proceed, just with a
           different function */
        trlwarn ("Can NOT process MAMA data with ACSCTE...");
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

    if (GetCTEFlags (&acs, &phdr)) {
        freeHdr(&phdr);
        return (status);
    }

    freeHdr(&phdr);

    /* Do basic CCD image reduction. */

    if (acs.printtime) {
        TimeStamp("Begin processing", acs.rootname);
    }

    if (DoCTE(&acs, forwardModelOnly)) {
        return (status);
    }

    trlmessage("\n");
    PrEnd("ACSCTE");

    if (acs.printtime) {
        TimeStamp("ACSCTE completed", acs.rootname);
    }

    /* Write out temp trailer file to final file */
    WriteTrlFile ();

    return (status);
}


void InitCTETrl (char *input, char *output, const char * isuffix, const char * osuffix) {

    extern int status;

    char trl_in[CHAR_LINE_LENGTH+1]; 	/* trailer filename for input */
    char trl_out[CHAR_LINE_LENGTH+1]; 	/* output trailer filename */

    assert(isuffix);
    assert(osuffix);
    char trlsuffix[] = "";

    int MkName (char *, char *, char *, char *, char *, int);
    void WhichError (int);
    int TrlExists (char *);

    /* Initialize internal variables */
    trl_in[0] = '\0';
    trl_out[0] = '\0';

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
