/* Do basic CCD image reduction for the current set of extensions.
 switches in the header will be reset from "PERFORM" to "COMPLETE".

 01-Jun-1998 WJH: Revised to only perform CCD operations on ACS observations.
 07-Jan-1999 WJH: Revised to only output overscan-trimmed image when BLEVCORR is
     performed successfully. In other cases, full image is written out.
 12-Sep-2000 WJH: Added support for new processing step: post-flash.
 29-Oct-2001 WJH: Updated to use CCDBIAS[A,B,C,D] for default bias levels
     and to use CCDOFST to select CCDBIAS from revised CCDTAB file. Also,
     reports bias levels for each AMP.
 10-Dec-2001 WJH: Added update of SIZAXIS keywords to reflect trimmed size.
 11-Dec-2012 PLL: Moved FLSHCORR stuff to ACS2D.
 12-Aug-2012 PLL: Separated PCTECORR from ACSCCD.
 21-Feb-2017 PLL: Added SINKCORR.
 16-May-2018 MDD: Specific supported subarrays can now be bias shift
     corrected. Clarified the logic that determines if data is processed
     with doBlev or bias_shift/cross_talk/destripe correction.
 14-May-2020: MDD: Modified to apply the full-well saturation flags stored
     as an image to the data instead of in the doDQI step.
 24-May-2021 MDD: Read the post-flashed and unflashed columns from the CCDTAB
     reference file to use for the offset to the DARKTIME FITS keyword value
     under appropropriate date and supported subarray criteria.  The DARKTIME
     keyword value is now the default scaling factor for the superdarks, with
     the offset being an additive correction to DARKTIME under appropriate
     circumstances. The offset values is applicable for WFC and HRC only.
     The code was moved from acs2d/dodark.c which contained the original
     implementation.
 11-Jan-2022 MDD: Only call computeDarktime() once to compute the updated
     DARKTIME primary keyword. Update the keyword in main to avoid passing
     of the reference of the output data unnecessarily. Return type of the
     computeDarktime() is void.  Mods triggered by routine failing under the
     condaforge compilation.
 23-Jul-2025 PLL: Moved some code out of doccd.c in order to lessen
     its complexity .

** This code is a trimmed down version of CALSTIS1 do2d.c.
*/
#include <string.h>
#include "acs.h"
#include "acsinfo.h"
#include "hstcal_memory.h"
#include "hstcalerr.h"
#include "hstio.h"

static void dqiMsg (ACSInfo *, const int);
static void BiasMsg (ACSInfo *, const int);
static void BlevMsg (ACSInfo *, const int);
static void SinkMsg (ACSInfo *, const int);


int DoCCD (ACSInfo *acs_info) {

    /* arguments:
       acs_info   i: calibration switches and info
    */

    extern int status;

    char buff[ACS_FITS_REC+1];

    int to_electrons(ACSInfo *, SingleGroup *);
    void computeDarktime(ACSInfo *, float *);
    int doBias (ACSInfo *, SingleGroup *);
    int biasHistory (ACSInfo *, Hdr *);
    int performBlevCorr(ACSInfo *, ACSInfo *, SingleGroup *, int *, int *, int *, int * , int *);
    int blevHistory (ACSInfo *, Hdr *, int, int);
    int CCDHistory (ACSInfo *, Hdr *);
    int doDQI (ACSInfo *, SingleGroup *);
    int dqiHistory (ACSInfo *, Hdr *);
    int doNoise (ACSInfo *, SingleGroup *, int *);
    int noiseHistory (Hdr *);
    int SinkDetect (ACSInfo *, SingleGroup *);
    int sinkHistory(const ACSInfo *, Hdr *);
    int GetACSGrp (ACSInfo *, Hdr *);
    int OmitStep (int);
    int PutKeyFlt (Hdr *, char *, float, char *);
    int PutKeyInt (Hdr *, char *, int, char *);
    int PutKeyStr (Hdr *, char *, char *, char *);
    void PrSwitch (char *, int);
    void PrRefInfo (char *, char *, char *, char *, char *);
    void UCalVer (Hdr *);
    void UFilename (char *, Hdr *);
    int FindOverscan (ACSInfo *, int, int, int *, int *);
    int GetCCDTab (ACSInfo *, int, int);
    int GetKeyBool (Hdr *, char *, int, Bool, Bool *);
    int doFullWellSat(ACSInfo *, SingleGroup *);

    /*========================Start Code here =======================*/

    PtrRegister ptrReg;
    initPtrRegister(&ptrReg);

    /* set up an array of SingleGroups to hold all image extensions */
    SingleGroup * x = NULL;   /* used for both input and output */
    x = (SingleGroup *) malloc(acs_info->nimsets * sizeof(SingleGroup));
    if (!x)
       return (status = OUT_OF_MEMORY);
    addPtr(&ptrReg, x, &free);

    ACSInfo * acs   = NULL;   /* hold a copy of the acs_info struct for each extension */
    acs = (ACSInfo *) malloc(acs_info->nimsets * sizeof(ACSInfo));
    if (!acs) {
       freeOnExit (&ptrReg);
       return (status = OUT_OF_MEMORY);
    }
    addPtr(&ptrReg, acs, &free);

    {unsigned int i;
    for (i = 0; i < acs_info->nimsets; i++) {
        initSingleGroup(&x[i]);
        acs[i] = *acs_info;
    }}

    if (acs_info->printtime)
        TimeStamp ("Open SingleGroup now...", "");

    /* Open the input image.
       FOR ACS: Keep this in memory throughout processing.
       Read in reference images line-by-line within the individual
       processing step functions and pass along modified input image
       from one step to the next...
    */
    {unsigned int i;
    for (i = 0; i < acs_info->nimsets; i++) {
        getSingleGroup(acs[i].input, i+1, &x[i]);
        if (hstio_err()) {
            freeOnExit (&ptrReg);
            return (status = OPEN_FAILED);
        }
        addPtr(&ptrReg, &x[i], &freeSingleGroup);
    }}

    if (acs_info->printtime)
        TimeStamp ("Input read into memory", acs_info->rootname);

    /* Get header info that varies from imset to imset necessary
       for reading CCDTAB.
    */
    Bool subarray = False;
    {unsigned int i;
    for (i = 0; i < acs_info->nimsets; i++) {
        if (GetACSGrp(&acs[i], &x[i].sci.hdr)) {
            freeOnExit (&ptrReg);
            return status;
        }

        /* Read in keywords from primary header... */

        if (GetKeyBool(x[i].globalhdr, "SUBARRAY", NO_DEFAULT, 0, &subarray)) {
            freeOnExit (&ptrReg);
            return status;
        }

        if (subarray)
            acs[i].subarray = YES;
        else
            acs[i].subarray = NO;

        /* For the CCD, update primary header keywords.
           Also reset CRCORR if there's only one image set.
        */
        /* Get values from tables, using same function used in ACSCCD. */
        if (GetCCDTab(&acs[i], x[i].sci.data.nx, x[i].sci.data.ny)) {
            freeOnExit (&ptrReg);
            return status;
        }
    }}

    int primaryIdx = 0;
    if (PutKeyFlt(x[primaryIdx].globalhdr, "ATODGNA", acs[primaryIdx].atodgain[0], "")) {
        freeOnExit (&ptrReg);
        return status;
    }
    if (PutKeyFlt(x[primaryIdx].globalhdr, "ATODGNB", acs[primaryIdx].atodgain[1], "")) {
        freeOnExit (&ptrReg);
        return status;
    }
    if (PutKeyFlt(x[primaryIdx].globalhdr, "ATODGNC", acs[primaryIdx].atodgain[2], "")) {
        freeOnExit (&ptrReg);
        return status;
    }
    if (PutKeyFlt(x[primaryIdx].globalhdr, "ATODGND", acs[primaryIdx].atodgain[3], "")) {
        freeOnExit (&ptrReg);
        return status;
    }
    if (PutKeyFlt(x[primaryIdx].globalhdr, "READNSEA", acs[primaryIdx].readnoise[0], "")) {
        freeOnExit (&ptrReg);
        return status;
    }
    if (PutKeyFlt(x[primaryIdx].globalhdr, "READNSEB", acs[primaryIdx].readnoise[1], "")) {
        freeOnExit (&ptrReg);
        return status;
    }
    if (PutKeyFlt(x[primaryIdx].globalhdr, "READNSEC", acs[primaryIdx].readnoise[2], "")) {
        freeOnExit (&ptrReg);
        return status;
    }
    if (PutKeyFlt(x[primaryIdx].globalhdr, "READNSED", acs[primaryIdx].readnoise[3], "")) {
        freeOnExit (&ptrReg);
        return status;
    }

    trlmessage("\n");

    PrRefInfo("ccdtab", acs[primaryIdx].ccdpar.name, acs[primaryIdx].ccdpar.pedigree,
              acs[primaryIdx].ccdpar.descrip, acs[primaryIdx].ccdpar.descrip2);

    if (CCDHistory(&acs[primaryIdx], x[primaryIdx].globalhdr)) {
        freeOnExit (&ptrReg);
        return status;
    }

    /* Get overscan region information from OSCNTAB */
    int * overscan = NULL;
    int * virtOverscan = NULL;
    overscan = (int *) malloc(acs_info->nimsets * sizeof(int));
    if (!overscan) {
       freeOnExit (&ptrReg);
       return (status = OUT_OF_MEMORY);
    }
    addPtr (&ptrReg, overscan, &free);

    virtOverscan = (int *) malloc(acs_info->nimsets * sizeof(int));
    if (!virtOverscan) {
       freeOnExit (&ptrReg);
       return (status = OUT_OF_MEMORY);
    }
    addPtr (&ptrReg, virtOverscan, &free);

    {unsigned int i;
    for (i = 0; i < acs_info->nimsets; i++) {
        if (FindOverscan(&acs[i], x[i].sci.data.nx, x[i].sci.data.ny,
                         &overscan[i], &virtOverscan[i])) {
            freeOnExit (&ptrReg);
            return status;
        }
    }}

    /************************************************************************/
    /* Data quality initialization and (for the CCDs) check saturation. */
    const int nextver = 1;
    dqiMsg(&acs[primaryIdx], nextver);

    {unsigned int i;
    for (i = 0; i < acs_info->nimsets; i++) {
        if (acs[i].dqicorr == PERFORM || acs[i].dqicorr == DUMMY) {
            if (doDQI(&acs[i], &x[i])) {
                freeOnExit (&ptrReg);
                return status;
            }
        }
    }}

    PrSwitch("dqicorr", COMPLETE);

    if (acs_info->printtime)
        TimeStamp("DQICORR complete", acs_info->rootname);

    if (!OmitStep(acs[primaryIdx].dqicorr)) {
        if (dqiHistory(&acs[primaryIdx], x[primaryIdx].globalhdr)) {
            freeOnExit (&ptrReg);
            return status;
        }
    }
    /************************************************************************/

    /************************************************************************/
    /* Subtract bias image. */
    BiasMsg(&acs[primaryIdx], nextver);

    {unsigned int i;
    for (i = 0; i < acs_info->nimsets; i++) {
        if (acs[i].biascorr == PERFORM) {
            if (doBias(&acs[i], &x[i])) {
                freeOnExit (&ptrReg);
                return status;
            }
        }
    }}

    PrSwitch("biascorr", COMPLETE);

    if (acs_info->printtime)
        TimeStamp ("BIASCORR complete", acs_info->rootname);

    if (!OmitStep(acs[primaryIdx].biascorr)) {
        if (biasHistory(&acs[primaryIdx], x[primaryIdx].globalhdr)) {
            freeOnExit (&ptrReg);
            return status;
       }
    }
    /************************************************************************/

    /************************************************************************/
    /* convert data to electrons */
    {unsigned int i;
    for (i = 0; i < acs_info->nimsets; i++) {
        if (to_electrons(&acs[i], &x[i])) {
            freeOnExit (&ptrReg);
            return status;
        }

        if (PutKeyStr(&x[i].sci.hdr, "BUNIT", "ELECTRONS", "")) {
            freeOnExit (&ptrReg);
            return status;
        }
        if (PutKeyStr(&x[i].err.hdr, "BUNIT", "ELECTRONS", "")) {
            freeOnExit (&ptrReg);
            return status;
        }
    }}
    /************************************************************************/

    /************************************************************************/
    /* Subtract bias level determined from overscan. */
    BlevMsg(&acs[primaryIdx], nextver);

    int * blevcorr = NULL;
    blevcorr = (int *) malloc(acs_info->nimsets * sizeof(int));
    if (!blevcorr) {
       freeOnExit (&ptrReg);
       return (OUT_OF_MEMORY);
    }
    addPtr (&ptrReg, blevcorr, free);

    {unsigned int i;
    for (i = 0; i < acs_info->nimsets; i++) {
        blevcorr[i] = PERFORM;
    }}

    int done = NO;
    int driftcorr = NO;     /* true means bias level was corrected for drift */

    if (performBlevCorr(acs_info, acs, x, overscan, virtOverscan, blevcorr, &driftcorr, &done)) {
        freeOnExit (&ptrReg);
        return status;
    }

    if (!OmitStep(acs[primaryIdx].blevcorr)) {
        if (blevHistory(&acs[primaryIdx], x[primaryIdx].globalhdr, done, driftcorr)) {
            freeOnExit (&ptrReg);
            return status;
        }
    }
    /************************************************************************/

    /************************************************************************/
    /* Compute the DARKTIME */
    float darktime = 0.0;  /* final darktime value after applicable offset, if any */

    /* Only need to invoke once as the needed keywords are inherited from
       the primary, and the DARKTIME applies to all imsets */
    computeDarktime(&acs[0], &darktime);
    if (PutKeyFlt(x[0].globalhdr, "DARKTIME", darktime, "")) {
        freeOnExit (&ptrReg);
        return status;
    }
    /************************************************************************/

    /************************************************************************/
    /* Apply the saturation image.
       Strictly speaking, the application of the full-well saturation image is
       not a calibration step (i.e., there is no SATCORR), but the application
       of a 2D image to flag pixels versus using a single scalar to flag
       saturated pixels as previously done in dqicorr needs to be done here
       after BIASCORR and BLEVCORR.  This should only be done if both
       BIASCORR and BLEVCORR have been performed, and the data has been
       converted to electrons.  This flagging is only applicable for the
       two CCDs (WFC and HRC). */
    {unsigned int i;
    for (i = 0; i < acs_info->nimsets; i++) {
        if (acs[i].biascorr == PERFORM && acs[i].blevcorr == PERFORM) {
            trlmessage("\nFull-well saturation flagging being performed for imset %d.\n", i+1);
            if (doFullWellSat(&acs[i], &x[i])) {
                freeOnExit (&ptrReg);
                return status;
            }
        } else {
            trlwarn("\nNo Full-well saturation flagging being performed for imset %d.\n", i+1);
        }
    }}
    /************************************************************************/

    /************************************************************************/
    /* Fill in the error array, if it initially contains all zeros. */
    if (acs->noisecorr == PERFORM) {
        {unsigned int i;
        for (i = 0; i < acs_info->nimsets; i++) {
            if (doNoise(&acs[i], &x[i], &done)) {
                freeOnExit (&ptrReg);
                return status;
            }
        }}

        if (done) {
            if (noiseHistory(x[primaryIdx].globalhdr)) {
                freeOnExit (&ptrReg);
                return status;
            }

            trlmessage("\n    Uncertainty array initialized,");
            buff[0] = '\0';

            snprintf(MsgText, sizeof(MsgText), "    readnoise =");
            {unsigned int i;
            for (i=0; i < NAMPS-1; i++) {
                if (acs[primaryIdx].readnoise[i] > 0) {
                    snprintf(buff, sizeof(buff), "%.5g,", acs[primaryIdx].readnoise[i]);
                    strcat(MsgText, buff);
                }
            }}

            if (acs[primaryIdx].readnoise[NAMPS-1] > 0) {
                snprintf(buff, sizeof(buff), "%.5g", acs[primaryIdx].readnoise[NAMPS-1]);
                strcat(MsgText, buff);
            }
            trlmessage(MsgText);

            snprintf(MsgText, sizeof(MsgText), "    gain =");
            {unsigned int i;
            for (i=0; i < NAMPS-1; i++) {
                if (acs[primaryIdx].atodgain[i] > 0) {
                    snprintf(buff, sizeof(buff), "%.5g,", acs[primaryIdx].atodgain[i]);
                    strcat(MsgText, buff);
                }
            }}

            if (acs[primaryIdx].atodgain[NAMPS-1] > 0) {
                snprintf(buff, sizeof(buff), "%.5g", acs[primaryIdx].atodgain[NAMPS-1]);
                strcat(MsgText, buff);
            }
            trlmessage(MsgText);

            snprintf(MsgText, sizeof(MsgText), "   default bias levels =");
            {unsigned int i;
            for (i=0; i < NAMPS-1; i++) {
                if (acs[primaryIdx].ccdbias[i] > 0) {
                    snprintf(buff, sizeof(buff), "%.5g,",
                            acs[primaryIdx].ccdbias[i] * acs[primaryIdx].atodgain[i]);
                    strcat(MsgText, buff);
                }
            }}

            if (acs[primaryIdx].ccdbias[NAMPS-1] > 0) {
                snprintf(buff, sizeof(buff), "%.5g",
                        acs[primaryIdx].ccdbias[NAMPS-1] * acs[primaryIdx].atodgain[NAMPS-1]);
                strcat(MsgText, buff);
            }
            trlmessage("%s", MsgText);

            if (acs_info->printtime)
                TimeStamp("Uncertainty array initialized", acs->rootname);
        }
    }
    /************************************************************************/

    /************************************************************************/
    /* Flag sink pixels.                                                    */
    SinkMsg (&acs[primaryIdx], nextver);

    if (acs->sinkcorr == PERFORM) {
        {unsigned int i;
        for (i = 0; i < acs_info->nimsets; i++) {
            if (SinkDetect(&acs[i], &x[i])) {
                freeOnExit (&ptrReg);
                return status;
            }
        }}

        PrSwitch("sinkcorr", COMPLETE);

        if (acs_info->printtime)
            TimeStamp("SINKCORR complete", acs_info->rootname);
    }

    if ((status = sinkHistory(&acs[primaryIdx], x[primaryIdx].globalhdr))) {
        freeOnExit (&ptrReg);
        return status;
    }
    /************************************************************************/

    /* Write this imset to the output file.  The
       CAL_VER and FILENAME keywords will be updated, and the primary
       header will be written.
    */
    UCalVer(x[primaryIdx].globalhdr);
    UFilename(acs_info->output, x[primaryIdx].globalhdr);

    /* If BLEVCORR was performed, then output trimmed image,
       otherwise output full image...
    */
    int option = 0; /* Write data to disk */
    int sizex = 0, sizey = 0;  /* size of output image */
    {unsigned int i;
    for (i = 0; i < acs_info->nimsets; i++) {
        if (blevcorr[i] == COMPLETE || acs[i].blevcorr == COMPLETE) {
            /* BLEVCORR was completed, so overscan regions can be trimmed... */
            if (acs_info->verbose) {
                trlmessage("Writing out image with trimx = %d,%d, trimy = %d,%d\nOutputting chip %d",
                        acs[i].trimx[0], acs[i].trimx[1],
                        acs[i].trimy[0], acs[i].trimy[1],
                        acs[i].chip);
            }

            /* Output overscan-trimmed, bias-subtracted image
               using the routine 'putSingleGroupSect' from HSTIO V2.0 (or later)
               to write it directly from the full image in memory.
            */
            sizex = x[i].sci.data.nx - (acs[i].trimx[0] + acs[i].trimx[1]);
            sizey = x[i].sci.data.ny - (acs[i].trimy[0] + acs[i].trimy[1]);

            /* Update SIZAXIS keywords to reflect the new trimmed image size. */
            if (PutKeyInt(&x[i].sci.hdr, "SIZAXIS1", sizex,"")) {
                freeOnExit (&ptrReg);
                return status;
            }
            if (PutKeyInt(&x[i].sci.hdr, "SIZAXIS2", sizey,"")) {
                freeOnExit (&ptrReg);
                return status;
            }

            /* Now, write out updated trimmed data to disk... */
            putSingleGroupSect(acs_info->output, i+1, &x[i], acs[i].trimx[0],
                               acs[i].trimy[0], sizex, sizey, option);

        } else {
            /* BLEVCORR was not completed, so keep overscan regions... */
            if (acs_info->verbose) {
                trlmessage("Writing out FULL image with overscan regions\n"
                           "Outputting chip %d", acs[i].chip);
            }

            putSingleGroup(acs_info->output, i+1, &x[i], option);
        }

        /* If putSingleGroup/putSingleGroupSect failed above, then error out. */
        if (hstio_err()) {
            trlerror("Couldn't write imset %d.", i+1);
            freeOnExit (&ptrReg);
            return (status = WRITE_FAILED);
        }

        if (acs_info->printtime)
            TimeStamp ("Output written to disk", acs_info->rootname);

        /* x[i] memory will be freed on exit */
    }}

    freeOnExit (&ptrReg);

    return status;
}


static void dqiMsg (ACSInfo *acs, const int extver) {
    int OmitStep (int);
    void PrSwitch (char *, int);
    void PrRefInfo (char *, char *, char *, char *, char *);

    trlmessage("\n");
    PrSwitch ("dqicorr", acs->dqicorr);

    if (extver == 1 && !OmitStep (acs->dqicorr)) {
        PrRefInfo ("dqitab", acs->bpix.name, acs->bpix.pedigree,
                   acs->bpix.descrip, acs->bpix.descrip2);
    }
}


static void BiasMsg (ACSInfo *acs, const int extver) {
    int OmitStep (int);
    void PrSwitch (char *, int);
    void PrRefInfo (char *, char *, char *, char *, char *);

    trlmessage("\n");
    PrSwitch ("biascorr", acs->biascorr);

    if (extver == 1 && !OmitStep (acs->biascorr)) {
        PrRefInfo ("biasfile", acs->bias.name, acs->bias.pedigree,
                   acs->bias.descrip, "");
    }
}


static void BlevMsg (ACSInfo *acs, const int extver) {
    void PrSwitch (char *, int);
    void PrRefInfo (char *, char *, char *, char *, char *);
    int OmitStep (int);

    trlmessage("\n");
    PrSwitch ("blevcorr", acs->blevcorr);

    if (extver == 1 && !OmitStep (acs->blevcorr)) {
        PrRefInfo ("oscntab", acs->oscn.name, acs->oscn.pedigree,
                   acs->oscn.descrip, acs->oscn.descrip2);
    }
}


static void SinkMsg (ACSInfo *acs, const int extver) {
    void PrSwitch (char *, int);
    void PrRefInfo (char *, char *, char *, char *, char *);
    int OmitStep (int);

    trlmessage("\n");
    PrSwitch ("sinkcorr", acs->sinkcorr);

    if (extver == 1 && !OmitStep (acs->sinkcorr)) {
        PrRefInfo ("snkcfile", acs->sink.name, acs->sink.pedigree,
                   acs->sink.descrip, "");
    }
}
