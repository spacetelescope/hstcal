/* Do CTE loss correction for the current set of extensions.
 The primary header will be updated with history information, and the calibration
 switches in the header will be reset from "PERFORM" to "COMPLETE".

 12-Aug-2012 PLL: Separated PCTECORR from ACSCCD.

*/
# include <string.h>
# include <stdio.h>
# include <time.h>
# include <assert.h>

# include "hstio.h"

# include "acs.h"
# include "acsinfo.h"
# include "acsversion.h"
# include "hstcalerr.h"

# include "../../../../ctegen2/ctegen2.h"
# include "pcte.h"

static void PCTEMsg (ACSInfo *, int);
static int OscnTrimmed (Hdr*, Hdr *);

#define SZ_CBUF 24
int getCTE_NAME(char * filename, char * cteName, int cteNameBufferLength);

int DoCTE (ACSInfo *acs_info) {

    /* arguments:
       acs   i: calibration switches and info
    */

    extern int status;

    SingleGroup * x;    /* used for both input and output */
    ACSInfo * acs;    /* hold a copy of the acs_info struct for each extension */
    int option = 0;

    int i;        /* loop index */
    Bool subarray;
    int CCDHistory (ACSInfo *, Hdr *);
    int doPCTEGen1 (ACSInfo *, SingleGroup *);
    int doPCTEGen2 (ACSInfo *,  CTEParamsFast * pars, SingleGroup *);
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
        CTEParamsFast ctePars;

        //NOTE: The char * below should be const but this would require a massive refactoring.
        char * cteTabFilename = acs->pcte.name;

        //open PCTETAB and check for CTE algorithm version name: CTE_NAME
        char cteName[SZ_CBUF];
        *cteName = '\0';
        int cteAlgorithmGen = 0;
        int tmpStatus = getCTE_NAME(cteTabFilename, cteName, SZ_CBUF);
        if (tmpStatus == KEYWORD_MISSING)
            cteAlgorithmGen = 1;
        else if (tmpStatus)
        {
            freeOnExit(&ptrReg);
            return (status = tmpStatus);
        }
        else if (strcmp(cteName, ACS_GEN2_CTE_NAME) == 0)
        {
            cteAlgorithmGen = 2;
            trlmessage("(pctecorr) Generation 2 PCTETAB file auto-detected.");
        }
        else
        {
            cteAlgorithmGen = 1;
            trlmessage("(pctecorr) Generation 1 PCTETAB file auto-detected.");
        }

        if (acs_info->cteAlgorithmGen && (acs_info->cteAlgorithmGen != cteAlgorithmGen))
        {
            char msgBuffer[MSG_BUFF_LENGTH];
            sprintf(msgBuffer, "(pctecorr) Cmd line option '--ctegen %d' specified yet gen%d algorithm detected in PCTETAB (CTE_NAME).",
                    acs_info->cteAlgorithmGen, cteAlgorithmGen);
            trlerror(msgBuffer);
            freeOnExit(&ptrReg);
            return (status = CAL_FILE_MISSING);
        }
        else if (!acs_info->cteAlgorithmGen && (cteAlgorithmGen != 2))
        {
            char msgBuffer[MSG_BUFF_LENGTH];
            sprintf(msgBuffer, "(pctecorr) Gen%d algorithm detected in PCTETAB (CTE_NAME). Default is gen2, use '--ctegen %d' to override.",
                    cteAlgorithmGen, cteAlgorithmGen);
            trlerror(msgBuffer);
            freeOnExit(&ptrReg);
            return (status = CAL_FILE_MISSING);
        }
        acs_info->cteAlgorithmGen = cteAlgorithmGen;

        if (acs_info->cteAlgorithmGen == 2)
        {
            sprintf(MsgText, "(pctecorr) Reading CTE parameters from PCTETAB file: '%s'...", cteTabFilename);
            trlmessage(MsgText);
            if ((status = PutKeyStr(x[0].globalhdr, "PCTETAB", cteTabFilename, "CTE Correction Table")))
            {
                trlerror("(pctecorr) failed to update PCTETAB keyword in image header");
                return status;
            }
            //Get parameters from PCTETAB reference file
            addPtr(&ptrReg, &ctePars, &freeCTEParamsFast);
            unsigned nScaleTableColumns = N_COLUMNS_FOR_RAZ_CDAB_ALIGNED_IMAGE;
            initCTEParamsFast(&ctePars, TRAPS, 0, 0, nScaleTableColumns, acs_info->nThreads);
            ctePars.refAndIamgeBinsIdenticle = True;
            ctePars.verbose = acs->verbose = 0 ? False : True;
            if ((status = allocateCTEParamsFast(&ctePars)))
            {
                freeOnExit(&ptrReg);
                return (status);
            }

            if ((status = loadPCTETAB(cteTabFilename, &ctePars)))
            {
                freeOnExit(&ptrReg);
                return (status);
            }

            if ((status = getCTEParsFromImageHeader(&x[0], &ctePars)))
            {
                freeOnExit(&ptrReg);
                return (status);
            }

            ctePars.scale_frac = (acs->expstart - ctePars.cte_date0) / (ctePars.cte_date1 - ctePars.cte_date0);

            if ((status = populateHeaderWithCTEKeywordValues(&x[0], &ctePars)))
            {
                freeOnExit(&ptrReg);
                return (status);
            }

            sprintf(MsgText, "(pctecorr) Read noise level PCTERNCL: %f", ctePars.rn_amp);
            trlmessage(MsgText);
            sprintf(MsgText, "(pctecorr) Readout simulation forward modeling iterations PCTENFOR: %i",
                    ctePars.n_forward);
            trlmessage(MsgText);
            sprintf(MsgText, "(pctecorr) Number of iterations used in the parallel transfer PCTENPAR: %i",
                    ctePars.n_par);
            trlmessage(MsgText);
            sprintf(MsgText, "(pctecorr) CTE_FRAC: %f", ctePars.scale_frac);
            trlmessage(MsgText);

            trlmessage("(pctecorr) PCTETAB read");
        }

        for (i = 0; i < acs_info->nimsets; i++)
        {
            clock_t begin = (double)clock();
            if (acs_info->cteAlgorithmGen == 1)//make explicit as not using bool
            {
                if ((status = doPCTEGen1(&acs[i], &x[i])))
                    return status;
            }
            else
            {
                if ((status = doPCTEGen2(&acs[i], &ctePars, &x[i])))
                    return status;
            }
            double time_spent = ((double) clock()- begin +0.0) / CLOCKS_PER_SEC;
            sprintf(MsgText,"(pctecorr) CTE run time for current chip: %.2f(s) with %i procs/threads\n", time_spent/acs_info->nThreads, acs_info->nThreads);
            trlmessage(MsgText);
        }
        freeOnExit(&ptrReg);

        PrSwitch("pctecorr", COMPLETE);

        if (acs_info->printtime)
            TimeStamp ("PCTECORR complete", acs->rootname);
    }

    if (!OmitStep(acs_info->pctecorr)) {
        if (pcteHistory(&acs[0], x[0].globalhdr))
            return (status);

        if (acs_info->cteAlgorithmGen == 1)
        	UCteVer(x[0].globalhdr, ACS_GEN1_CTE_NAME, ACS_GEN1_CTE_VER);
        else if (acs_info->cteAlgorithmGen == 2)
        	UCteVer(x[0].globalhdr, ACS_GEN2_CTE_NAME, ACS_GEN2_CTE_VER);
        else
        	assert(0);//unimplemented

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
