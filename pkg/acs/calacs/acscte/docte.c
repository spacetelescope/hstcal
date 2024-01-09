/* Do CTE loss correction for the current set of extensions.
 The primary header will be updated with history information, and the calibration
 switches in the header will be reset from "PERFORM" to "COMPLETE".

 12-Aug-2012 PLL: Separated PCTECORR from ACSCCD.
 14-Sep-2023 MDD: Updated to accommodate only the "Parallel/Serial PixelCTE 2023"
     (aka Generation 3) correction.

*/
# include <string.h>
# include <stdio.h>
# include <time.h>
# include <assert.h>
# include <stdbool.h>

#include "hstcal_memory.h"
#include "hstcal.h"
# include "hstio.h"

# include "acs.h"
# include "acsinfo.h"
# include "acsversion.h"
# include "hstcalerr.h"
#include "trlbuf.h"

# include "../../../../ctegen2/ctegen2.h"
# include "pcte.h"

static void PCTEMsg (ACSInfo *, int);
static int OscnTrimmed (Hdr*, Hdr *);

/* 
   Typical order of processing is Chip 2 (Amps C and D) and then 
   Chip 1 (Amps A and B). See ctehelpers.c for full description of
   PCTETAB which clarifies the order of the extensions for processing.
*/
/* CTE data as stored in the PCTETAB */ 
int SET_TO_PROCESS[] = {13, 17, 5, 9, 1};  // Starting extension number of a set of qprof/sclbycol/rprof/cprof
char AMP_OF_PROCESS[] = {'C', 'D', 'A', 'B', 'P'}; // Amp corresponding to the set of corrections in SET_TO_PROCESS
# define AMPCALIBORDER "CDAB"

#define SZ_CBUF 24
int getCTE_NAME(char * filename, char * cteName, int cteNameBufferLength);

int DoCTE (ACSInfo *acs_info, const bool forwardModelOnly) {

    /* arguments:
       acs   i: calibration switches and info
    */

    extern int status;

    SingleGroup * x;  /* used for both input and output */
    ACSInfo * acs;    /* hold a copy of the acs_info struct for each extension */
    int option = 0;
    int numCorrection = 0; /* number of corrections, parallel and serial */
    int namp_imsets = 0;   /* number of loops to handle all amp for all chips */
    int ampID;             /* index amp A:0, B:1, etc. */
    int ampIDInCalib;      /* index amp C:0, D:1, etc. */
    char * amploc;         /* pointer to amp character in AMPSORDER */
    char * amplocInCalib;  /* pointer to amp character in calibration file */
    int numamps;

    int i;        /* loop index */
    Bool subarray;
    int CCDHistory (ACSInfo *, Hdr *);
    int doPCTEGen2 (ACSInfo *,  CTEParamsFast * pars, SingleGroup *, const bool forwardModelOnly, const char * corrType, char *ccdamp, char *amploc, int ampID);
    int doPCTESerial (ACSInfo *,  CTEParamsFast * pars, SingleGroup *, const bool forwardModelOnly);
    int pcteHistory (ACSInfo *, Hdr *);
    int GetACSGrp (ACSInfo *, Hdr *);
    int OmitStep (int);
    int PutKeyFlt (Hdr *, char *, float, char *);
    int PutKeyStr (Hdr *, char *, char *, char *);
    void PrSwitch (char *, int);
    void PrRefInfo (char *, char *, char *, char *, char *);
    void TimeStamp (char *, char *);
    void UCalVer (Hdr *);
    void UCteVer(Hdr * hdr, char * cteName, char * cteVersion);
    void UFilename (char *, Hdr *);
    int GetCCDTab (ACSInfo *, int, int);
    int GetKeyBool (Hdr *, char *, int, Bool, Bool *);
    int GetKeyStr (Hdr *, char *, int, char *, char *, int);
    void parseWFCamps (char *acsamps, int chip, char *ccdamp);

    /*========================Start Code here =======================*/

    /* set up an array of SingleGroups to hold all image extensions */
    x = (SingleGroup *) malloc(acs_info->nimsets * sizeof(SingleGroup));
    acs = (ACSInfo *) malloc(acs_info->nimsets * sizeof(ACSInfo));

    for (i = 0; i < acs_info->nimsets; i++) {
        initSingleGroup(&x[i]);
        acs[i] = *acs_info;
    }

    if (acs_info->printtime)
        TimeStamp ("Open SingleGroup now...", "");

    /* Open the input image.
       FOR ACS: Keep this in memory throughout processing.
       Read in reference images line-by-line within the individual
       processing step functions and pass along modified input image
       from one step to the next...
    */
    for (i = 0; i < acs_info->nimsets; i++) {
        getSingleGroup(acs[i].input, i+1, &x[i]);

        if (hstio_err()) {
            freeSingleGroup(&x[i]);
            return (status = OPEN_FAILED);
        }
    }

    if (acs_info->printtime)
        TimeStamp ("Input read into memory", acs_info->rootname);


    /* Get header info that varies from imset to imset necessary
       for reading CCDTAB.

       Probably redundant here but does not hurt.
    */
    for (i = 0; i < acs_info->nimsets; i++) {
        if (GetACSGrp(&acs[i], &x[i].sci.hdr)) {
            freeSingleGroup(&x[i]);
            return (status);
        }

        /* Has the CCD image been overscan trimmed? */
        if (OscnTrimmed (x[i].globalhdr, &x[i].sci.hdr)) {
            freeSingleGroup (&x[i]);
            return (status);
        }

        /* Read in keywords from primary header... */

        if (GetKeyBool(x[i].globalhdr, "SUBARRAY", NO_DEFAULT, 0, &subarray))
            return (status);

        if (subarray)
            acs[i].subarray = YES;
        else
            acs[i].subarray = NO;

        /* For the CCD, update primary header keywords.
           Also reset CRCORR if there's only one image set.
        */
        /* Get values from tables, using same function used in ACSCTE. */
        if (GetCCDTab(&acs[i], x[i].sci.data.nx, x[i].sci.data.ny)) {
            freeSingleGroup(&x[i]);
            return (status);
        }
    }

    if (PutKeyFlt(x[0].globalhdr, "ATODGNA", acs[0].atodgain[0], ""))
        return (status);
    if (PutKeyFlt(x[0].globalhdr, "ATODGNB", acs[0].atodgain[1], ""))
        return (status);
    if (PutKeyFlt(x[0].globalhdr, "ATODGNC", acs[0].atodgain[2], ""))
        return (status);
    if (PutKeyFlt(x[0].globalhdr, "ATODGND", acs[0].atodgain[3], ""))
        return (status);
    if (PutKeyFlt(x[0].globalhdr, "READNSEA", acs[0].readnoise[0], ""))
        return (status);
    if (PutKeyFlt(x[0].globalhdr, "READNSEB", acs[0].readnoise[1], ""))
        return (status);
    if (PutKeyFlt(x[0].globalhdr, "READNSEC", acs[0].readnoise[2], ""))
        return (status);
    if (PutKeyFlt(x[0].globalhdr, "READNSED", acs[0].readnoise[3], ""))
        return (status);

    trlmessage("\n");

    PrRefInfo("ccdtab", acs[0].ccdpar.name, acs[0].ccdpar.pedigree,
              acs[0].ccdpar.descrip, acs[0].ccdpar.descrip2);

    if (CCDHistory(&acs[0], x[0].globalhdr))
        return (status);

    /************************************************************************/
    /* perform CTE correction */
    PCTEMsg(&acs[0], 1);

    if (acs_info->pctecorr == PERFORM)
    {
        PtrRegister ptrReg;//move this up later
        initPtrRegister(&ptrReg);

        /* This variable will hold the serial(x) and then the parallel(y) CTE
           values. */
        CTEParamsFast ctePars;

        //NOTE: The char * below should be const but this would require a massive refactoring.
        char * cteTabFilename = acs->pcte.name;

        // Generation 1 CTE algorithm is obsolete.
        // Generation 3 CTE algorithm is the Generation 2 algorithm applied to both
        // serial and parallel trails.  There is no longer any choice for the
        // generation of the CTE algorithm. It is "generation 3", but it is preferably
        // referenced as the "Serial and Parallel CTE correction".

        // open PCTETAB and check for CTE algorithm version name: CTE_NAME
        char cteName[SZ_CBUF];
        *cteName = '\0';
        int cteAlgorithmGen = 0;
        int tmpStatus = getCTE_NAME(cteTabFilename, cteName, SZ_CBUF);
        if (tmpStatus == KEYWORD_MISSING)
        {
            trlerror("(pctecorr) No CTE_NAME keyword found in the primary header of PCTETAB file.");
            trlerror("(pctecorr) This version of PCTETAB is obsolete - use the 'Parallel and Serial PixelCTE 2023' PCTETAB.");
            freeOnExit(&ptrReg);
            return (status = ERROR_RETURN);
        }
        else if (tmpStatus)
        {
            freeOnExit(&ptrReg);
            return (status = tmpStatus);
        }
        else if (strcmp(cteName, ACS_GEN3_CTE_NAME) == 0)
        {
            cteAlgorithmGen = 3;
            trlmessage("(pctecorr) Parallel and Serial PCTETAB file auto-detected.");
        }
        else if (strcmp(cteName, ACS_GEN2_CTE_NAME) == 0)
        {
            trlerror("(pctecorr) Generation 2 PCTETAB file auto-detected.");
            trlerror("(pctecorr) This version of PCTETAB is obsolete - use the 'Parallel and Serial PixelCTE 2023' PCTETAB.");
            freeOnExit(&ptrReg);
            return (status = ERROR_RETURN);
        }
        else
        {
            trlerror("(pctecorr) Generation 1 PCTETAB file auto-detected.");
            trlerror("(pctecorr) This version of PCTETAB is obsolete - use the 'Parallel and Serial PixelCTE 2023' PCTETAB.");
            freeOnExit(&ptrReg);
            return (status = ERROR_RETURN);
        }
        acs_info->cteAlgorithmGen = cteAlgorithmGen;

        sprintf(MsgText, "(pctecorr) Reading CTE parameters from PCTETAB file: '%s'...", cteTabFilename);
        trlmessage(MsgText);
        if ((status = PutKeyStr(x[0].globalhdr, "PCTETAB", cteTabFilename, "CTE Correction Table")))
        {
            trlerror("(pctecorr) failed to update PCTETAB keyword in image header");
            return status;
        }
    /*
      Loop over the imsets as the CTE is applied per amp
    */
    char ccdamp[strlen(AMPSTR1)+1]; // string to hold amps on current chip
    addPtr(&ptrReg, &ctePars, &freeCTEParamsFast);
    unsigned nScaleTableColumns = N_COLUMNS_FOR_RAZ_CDAB_ALIGNED_IMAGE;
    for (i = 0; i < acs_info->nimsets; i++) {

        sprintf(MsgText, "(docte) Loop nimset: %i", i);
        trlmessage(MsgText);

	    // Determine which amps are on the current chip of the input image 
        ccdamp[0] = '\0'; // "reset" the string for reuse

        // Amps on this chip
        parseWFCamps(acs[i].ccdamp, acs[i].chip, ccdamp);
        numamps = strlen(ccdamp);

        int startOfSetInCalib;
        char corrType[10];
        corrType[0] = '\0';
        for (unsigned nthAmp = 0; nthAmp < numamps; ++nthAmp)
        {
            sprintf(MsgText, "(docte) Loop nthAmp: %i ccdamp: %s", nthAmp, ccdamp);
            trlmessage(MsgText);

            /* Get the amp letter and number where A:0, B:1, etc. as defined in acs.h */
            amploc = strchr(AMPSORDER, ccdamp[nthAmp]);
            ampID = amploc - AMPSORDER; 

            /* 
               Get the amp letter and number where C:0, D:1, etc. as defined in this file.
               This order is based upon the arrangement of the chips in the input science
               image (chip 2 with amps C and D in the first imset, and then chip 1 with
               amps A and B in the second imset. 
            */
            amplocInCalib = strchr(AMPCALIBORDER, ccdamp[nthAmp]); // This is a full string.
            ampIDInCalib = amplocInCalib - AMPCALIBORDER; // This is a number.

            sprintf(MsgText, "(docte) amploc: %s  ampID: %i amplocInCalib: %s ampIDInCalib: %i", amploc, ampID, amplocInCalib, ampIDInCalib);
            trlmessage(MsgText);

            startOfSetInCalib = SET_TO_PROCESS[ampIDInCalib];
            strcpy(corrType, "serial");
            if (startOfSetInCalib == 1) {
                strcpy(corrType, "parallel");
            }

            /* 
               Loop to control the application of the CTE corrections - the serial 
               (Extensions 5-8 [A], 9-12 [B], 13-16 [C], 17-20 [D]) correction for each
               Amp will be done first, then the parallel (Extensions 1-4) correction will
               be done.  The "Extensions" refer to the set of extensions in the input 
               dataset (*raw.fits) pertaining to a particular correction.

               As noted at the top of this file, the typical order of processing is 
               Chip 2 (Amps C and D) and then Chip 1 (Amps A and B). 
            */
            initCTEParamsFast(&ctePars, TRAPS, 0, 0, nScaleTableColumns, acs_info->nThreads);

            ctePars.refAndIamgeBinsIdenticle = True;
            ctePars.verbose = acs->verbose = 0 ? False : True;
            if ((status = allocateCTEParamsFast(&ctePars)))
            {
                freeOnExit(&ptrReg);
                return (status);
            }

            if ((status = loadPCTETAB(cteTabFilename, &ctePars, startOfSetInCalib)))
            {
                freeOnExit(&ptrReg);
                return (status);
            }

            if ((status = getCTEParsFromImageHeader(&x[0], &ctePars)))
            {
                freeOnExit(&ptrReg);
                return (status);
            }

            sprintf(MsgText,"(pctecorr) StartOfSet: %d  CorrType: %s \n", startOfSetInCalib, corrType);
            trlmessage(MsgText);

            ctePars.scale_frac = (acs->expstart - ctePars.cte_date0) / (ctePars.cte_date1 - ctePars.cte_date0);

            if ((status = populateImageFileWithCTEKeywordValues(&x[0], &ctePars)))
            {
                freeOnExit(&ptrReg);
                return (status);
            }

            sprintf(MsgText, "(pctecorr) IGNORING read noise level PCTERNOI from PCTETAB: %f. Using amp dependent values from CCDTAB instead", ctePars.rn_amp);
            trlwarn(MsgText);
            sprintf(MsgText, "(pctecorr) Readout simulation forward modeling iterations PCTENFOR: %i",
                    ctePars.n_forward);
            trlmessage(MsgText);
            sprintf(MsgText, "(pctecorr) Number of iterations used in the parallel transfer PCTENPAR: %i",
                    ctePars.n_par);
            trlmessage(MsgText);
            sprintf(MsgText, "(pctecorr) CTE_FRAC: %f", ctePars.scale_frac);
            trlmessage(MsgText);

            trlmessage("(pctecorr) PCTETAB read");

            /*
               The number of acs_info->nimsets represents the number of chips to process for the
               input file, and each chip can contain two amps. Since the CTE correction is done on 
               an amp basis, and only the correction information for that amp is available, then
               the looping is double the number of nimsets.
            */

            clock_t begin = (double)clock();
            if ((status = doPCTEGen2(&acs[i], &ctePars, &x[i], forwardModelOnly, corrType, &ccdamp[nthAmp], amploc, ampID)))
            {
                return status;
            }
            double time_spent = ((double) clock()- begin +0.0) / CLOCKS_PER_SEC;
            sprintf(MsgText,"(pctecorr) CTE run time for current chip: %.2f(s) with %i procs/threads\n", time_spent/acs_info->nThreads, acs_info->nThreads);
            trlmessage(MsgText);

        }  // End loop to accommodate both the serial and parallel CTE corrections

        freeOnExit(&ptrReg);
    }
    PrSwitch("pctecorr", COMPLETE);
    }

    if (acs_info->printtime)
        TimeStamp ("PCTECORR complete", acs->rootname);

    if (!OmitStep(acs_info->pctecorr)) {
        if (pcteHistory(&acs[0], x[0].globalhdr))
            return (status);

        UCteVer(x[0].globalhdr, ACS_GEN3_CTE_NAME, ACS_GEN3_CTE_VER);

    } else if (acs_info->pctecorr == OMIT) {
        /* this is just to make sure that the PCTECORR flag is set to OMIT in the
         * output image because it's possible that this is the no-CTE correction
         * half of ACSCCD from a raw image that has PCTECORR set to PERFORM.
         * MRD 16 Mar. 2011
         */
        PutKeyStr(x[0].globalhdr, "PCTECORR", "OMIT","");
    }
    /************************************************************************/

    /* Write this imset to the output file.  The
       CAL_VER and FILENAME keywords will be updated, and the primary
       header will be written.
    */
    UCalVer(x[0].globalhdr);
    UFilename(acs_info->output, x[0].globalhdr);

    for (i = 0; i < acs_info->nimsets; i++) {
        if (acs_info->verbose) {
            sprintf(MsgText, "Outputting chip %d", acs[i].chip);
            trlmessage(MsgText);
        }

        putSingleGroup(acs_info->output, i+1, &x[i], option);

        if (hstio_err()) {
            sprintf(MsgText, "Couldn't write imset %d.", i+1);
            trlerror(MsgText);
            return (status = WRITE_FAILED);
        }

        if (acs_info->printtime)
            TimeStamp ("Output written to disk", acs_info->rootname);

        /* Free x, which still has memory allocated. */
        freeSingleGroup (&x[i]);
    }

    free(x);
    free(acs);

    return (status);
}


static void PCTEMsg (ACSInfo *acs, int extver) {
    int OmitStep (int);
    void PrSwitch (char *, int);
    void PrRefInfo (char *, char *, char *, char *, char *);

    trlmessage ("\n");
    PrSwitch ("pctecorr", acs->pctecorr);

    if (extver == 1 && !OmitStep (acs->pctecorr)) {
        PrRefInfo ("pctefile", acs->pcte.name, acs->pcte.pedigree,
                   acs->pcte.descrip, "");
    }
}


/* This function verifies whether the overscan region has been
   trimmed from the image during BLEVCORR.  If not, it tells the user,
   and shuts down gracefully.
   *****************
   Copied from do2d.c
   *****************
*/
static int OscnTrimmed (Hdr *phdr, Hdr *hdr) {
    extern int status;

    double ltv1, ltv2;
    int blevcorr;

    int GetKeyDbl (Hdr *, char *, int, double, double *);
    int GetSwitch (Hdr *, char *, int *);

    if (GetSwitch (phdr, "BLEVCORR", &blevcorr))
        return(status);

    if (GetKeyDbl (hdr, "LTV1", USE_DEFAULT, 0., &ltv1))
        return (status);
    if (GetKeyDbl (hdr, "LTV2", USE_DEFAULT, 0., &ltv2))
        return (status);

    /* If there is any overscan region still specified in the
       image header, then we can not process this image yet...
    */
    if ( blevcorr != COMPLETE || (ltv1 > 0. || ltv2 > 0.)) {
        trlerror ("Overscan region not trimmed.  Please perform BLEVCORR.");
        status = CAL_STEP_NOT_DONE;
    }

    return (status);
}
