/* Do basic CCD image reduction for the current set of extensions.
 The primary header will be updated with history information, and the calibration
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

** This code is a trimmed down version of CALSTIS1 do2d.c.
*/
# include <string.h>
# include <stdio.h>

# include "hstio.h"

# include "acs.h"
# include "acsinfo.h"
# include "hstcalerr.h"


static void dqiMsg (ACSInfo *, int);
/*static void AtoDMsg (ACSInfo *, int);*/ /* Not used */
static void BiasMsg (ACSInfo *, int);
static void BlevMsg (ACSInfo *, int);
static void SinkMsg (ACSInfo *, int);


int DoCCD (ACSInfo *acs_info) {

    /* arguments:
       acs   i: calibration switches and info
    */

    extern int status;

    SingleGroup * x;  /* used for both input and output */
    ACSInfo * acs;    /* hold a copy of the acs_info struct for each extension */
    int option = 0;
    float meanblev;   /* mean value of overscan bias (for history) */
    int driftcorr;    /* true means bias level was corrected for drift */
    int done;         /* true means the input SingleGroup has been freed */
    int sizex, sizey;  /* size of output image */

    int i, j;  /* loop index */
    int * overscan;
    int * blevcorr;
    char buff[ACS_FITS_REC+1];
    Bool subarray;

    int to_electrons(ACSInfo *, SingleGroup *);
    int doBias (ACSInfo *, SingleGroup *);
    int biasHistory (ACSInfo *, Hdr *);
    int doBlev (ACSInfo *, SingleGroup *, int, float *, int *, int *);
    int bias_shift_corr(ACSInfo *, SingleGroup *, SingleGroup *);
    void cross_talk_corr(ACSInfo *, SingleGroup *);
    int doDestripe(ACSInfo *, SingleGroup *, SingleGroup *);
    int blevHistory (ACSInfo *, Hdr *, int, int);
    int CCDHistory (ACSInfo *, Hdr *);
    int doDQI (ACSInfo *, SingleGroup *);
    int dqiHistory (ACSInfo *, Hdr *);
    int doNoise (ACSInfo *, SingleGroup *, int *);
    int noiseHistory (Hdr *);
    int SinkDetect (ACSInfo *, SingleGroup *);
    int GetACSGrp (ACSInfo *, Hdr *);
    int OmitStep (int);
    int PutKeyFlt (Hdr *, char *, float, char *);
    int PutKeyInt (Hdr *, char *, int, char *);
    int PutKeyStr (Hdr *, char *, char *, char *);
    void PrSwitch (char *, int);
    void PrRefInfo (char *, char *, char *, char *, char *);
    void TimeStamp (char *, char *);
    void UCalVer (Hdr *);
    void UFilename (char *, Hdr *);
    int FindOverscan (ACSInfo *, int, int, int *);
    int GetCCDTab (ACSInfo *, int, int);
    int GetKeyBool (Hdr *, char *, int, Bool, Bool *);

    /*========================Start Code here =======================*/

    /* set up an array of SingleGroups to hold all image extensions */
    x = (SingleGroup *) malloc(acs_info->nimsets * sizeof(SingleGroup));
    acs = (ACSInfo *) malloc(acs_info->nimsets * sizeof(ACSInfo));
    overscan = (int *) malloc(acs_info->nimsets * sizeof(int));
    blevcorr = (int *) malloc(acs_info->nimsets * sizeof(int));

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
    */
    for (i = 0; i < acs_info->nimsets; i++) {
        if (GetACSGrp(&acs[i], &x[i].sci.hdr)) {
            freeSingleGroup(&x[i]);
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
        /* Get values from tables, using same function used in ACSCCD. */
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

    /* Get overscan region information from OSCNTAB */
    for (i = 0; i < acs_info->nimsets; i++) {
        if (FindOverscan(&acs[i], x[i].sci.data.nx, x[i].sci.data.ny,
                         &overscan[i])) {
            freeSingleGroup(&x[i]);
            return (status);
        }
    }

    /************************************************************************/
    /* Data quality initialization and (for the CCDs) check saturation. */
    dqiMsg(&acs[0], 1);

    for (i = 0; i < acs_info->nimsets; i++) {
        if (acs[i].dqicorr == PERFORM || acs[i].dqicorr == DUMMY) {
            if (doDQI(&acs[i], &x[i]))
                return (status);
        }
    }

    PrSwitch("dqicorr", COMPLETE);

    if (acs_info->printtime)
        TimeStamp("DQICORR complete", acs_info->rootname);

    if (!OmitStep(acs[0].dqicorr))
        if (dqiHistory(&acs[0], x[0].globalhdr))
            return (status);
    /************************************************************************/

    /************************************************************************/
    /* Subtract bias image. */
    BiasMsg(&acs[0], 1);

    for (i = 0; i < acs_info->nimsets; i++) {
        if (acs[i].biascorr == PERFORM) {
            if (doBias(&acs[i], &x[i]))
                return (status);
        }
    }

    PrSwitch("biascorr", COMPLETE);

    if (acs_info->printtime)
        TimeStamp ("BIASCORR complete", acs_info->rootname);

    if (!OmitStep(acs[0].biascorr))
        if (biasHistory(&acs[0], x[0].globalhdr))
            return (status);
    /************************************************************************/

    /************************************************************************/
    /* convert data to electrons */
    for (i = 0; i < acs_info->nimsets; i++) {
        if (to_electrons(&acs[i], &x[i])) {
            freeSingleGroup(&x[i]);
            return (status);
        }

        if (PutKeyStr(&x[i].sci.hdr, "BUNIT", "ELECTRONS", "")) {
            freeSingleGroup(&x[i]);
            return (status);
        }
        if (PutKeyStr(&x[i].err.hdr, "BUNIT", "ELECTRONS", "")) {
            freeSingleGroup(&x[i]);
            return (status);
        }
    }
    /************************************************************************/

    /************************************************************************/
    /* Subtract bias level determined from overscan. */
    BlevMsg(&acs[0], 1);

    for (i = 0; i < acs_info->nimsets; i++) {
        blevcorr[i] = PERFORM;
    }

    /* default blev case, do original bias level subtraction */
    if (acs_info->blevcorr == PERFORM &&
        ((acs_info->expstart < SM4MJD && acs_info->detector == WFC_CCD_DETECTOR) ||
         (acs_info->detector == WFC_CCD_DETECTOR && acs_info->subarray == YES) ||
         (acs_info->detector != WFC_CCD_DETECTOR))) {
        for (i = 0; i < acs_info->nimsets; i++) {
            /* This is set based on results of FindOver */
            done = overscan[i];

            if (doBlev(&acs[i], &x[i], acs[i].chip, &meanblev, &done,
                       &driftcorr))
                return (status);

            if (done) {
                trlmessage("Bias level from overscan has been subtracted;");

                if (PutKeyFlt(&x[i].sci.hdr, "MEANBLEV", meanblev,
                              "mean of bias levels subtracted in electrons")) {
                    return (status);
                }

                sprintf(MsgText, "     mean of bias levels subtracted was %.6g electrons.", meanblev);
                trlmessage(MsgText);

            } else {
                trlmessage("Default bias level from CCDTAB was subtracted.");
            }

            /* Provide immediate feedback to the user on the values computed
               for each AMP, but only for those amps used for the chip being
               processed.
            */
            for (j = 0; j < NAMPS; j++){
                if (acs[i].blev[j] != 0.) {
                    sprintf(MsgText, "     bias level of %.6g electrons was subtracted for AMP %c.", acs[i].blev[j], AMPSORDER[j]);
                    trlmessage(MsgText);

                    acs_info->blev[j] = acs[i].blev[j];
                }
            }

            /* Set this to complete so overscan-trimmed image will be
               written out.
            */
            blevcorr[i] = COMPLETE;
        }

        PrSwitch("blevcorr", COMPLETE);

        if (acs_info->printtime)
            TimeStamp("BLEVCORR complete", acs_info->rootname);

    /* post SM4, full frame WFC case */
    } else if (acs_info->blevcorr == PERFORM && acs_info->expstart > SM4MJD &&
               acs_info->detector == WFC_CCD_DETECTOR &&
               acs_info->subarray == NO) {

        /* This is set based on results of FindOver */
        if (overscan[0] == YES && overscan[1] == YES) {
            done = overscan[0];
        } else {
            done = NO;
        }

        if (done) {
            PrSwitch("blevcorr", PERFORM);

            /* only do bias-shift and cross talk corrections of images taken
               with gain = 2 and in dual-slope integrator mode. */
            if (strcmp(acs_info->jwrotype, "DS_int") == 0 &&
                    acs_info->ccdgain == 2) {
                trlmessage("Performing bias-shift and cross talk corrections.");

                if (bias_shift_corr(acs_info, &x[0], &x[1])) {
                    return status;
                }

                for (i = 0; i < acs_info->nimsets; i++) {
                    cross_talk_corr(&acs[i], &x[i]);
                }
            }

            trlmessage("Performing stripe removal and bias level subtraction.");

            if (doDestripe(acs_info, &x[0], &x[1])) {
                return status;
            }

            for (i = 0; i < acs_info->nimsets; i++) {
                blevcorr[i] = COMPLETE;
            }

            PrSwitch("blevcorr", COMPLETE);

            if (acs_info->printtime) {
                TimeStamp("BLEVCORR complete", acs->rootname);
            }
        } else {
            trlmessage("Overscan missing, no destriping or bias level subtraction possible.");
        }

        driftcorr = NO;
    }

    if (!OmitStep(acs[0].blevcorr))
        if (blevHistory(&acs[0], x[0].globalhdr, done, driftcorr))
            return (status);
    /************************************************************************/

    /************************************************************************/
    /* Fill in the error array, if it initially contains all zeros. */
    if (acs->noisecorr == PERFORM) {
        for (i = 0; i < acs_info->nimsets; i++) {
            if (doNoise(&acs[i], &x[i], &done))
                return (status);
        }

        if (done) {
            if (noiseHistory(x[0].globalhdr))
                return (status);

            trlmessage("\n    Uncertainty array initialized,");
            buff[0] = '\0';

            sprintf(MsgText, "    readnoise =");
            for (i=0; i < NAMPS-1; i++) {
                if (acs[0].readnoise[i] > 0) {
                    sprintf(buff, "%.5g,",acs[0].readnoise[i]);
                    strcat(MsgText, buff);
                }
            }

            if (acs[0].readnoise[NAMPS-1] > 0) {
                sprintf(buff, "%.5g",acs[0].readnoise[NAMPS-1]);
                strcat(MsgText, buff);
            }
            trlmessage(MsgText);

            sprintf(MsgText, "    gain =");
            for (i=0; i < NAMPS-1; i++) {
                if (acs[0].atodgain[i] > 0) {
                    sprintf(buff, "%.5g,",acs[0].atodgain[i]);
                    strcat(MsgText, buff);
                }
            }

            if (acs[0].atodgain[NAMPS-1] > 0) {
                sprintf(buff, "%.5g",acs[0].atodgain[NAMPS-1]);
                strcat(MsgText, buff);
            }
            trlmessage(MsgText);
            sprintf(MsgText, "   default bias levels =");

            for (i=0; i < NAMPS-1; i++) {
                if (acs[0].ccdbias[i] > 0) {
                    sprintf(buff, "%.5g,",
                            acs[0].ccdbias[i] * acs[0].atodgain[i]);
                    strcat(MsgText, buff);
                }
            }

            if (acs[0].ccdbias[NAMPS-1] > 0) {
                sprintf(buff, "%.5g",
                        acs[0].ccdbias[NAMPS-1] * acs[0].atodgain[NAMPS-1]);
                strcat(MsgText, buff);
            }
            trlmessage(MsgText);

            if (acs_info->printtime)
                TimeStamp("Uncertainty array initialized", acs->rootname);
        }
    }
    /************************************************************************/

    /************************************************************************/
    /* Flag sink pixels.                                                    */
    SinkMsg (&acs[0], 1);

    if (acs->sinkcorr == PERFORM) {
        for (i = 0; i < acs_info->nimsets; i++) {
            if (SinkDetect(&acs[i], &x[i]))
                return (status);
        }

        PrSwitch("sinkcorr", COMPLETE);

        if (acs_info->printtime)
            TimeStamp("SINKCORR complete", acs_info->rootname);
    }
    /************************************************************************/

    /* Write this imset to the output file.  The
       CAL_VER and FILENAME keywords will be updated, and the primary
       header will be written.
    */
    UCalVer(x[0].globalhdr);
    UFilename(acs_info->output, x[0].globalhdr);

    /* If BLEVCORR was performed, then output trimmed image,
       otherwise output full image...
    */
    for (i = 0; i < acs_info->nimsets; i++) {
        if (blevcorr[i] == COMPLETE || acs[i].blevcorr == COMPLETE) {
            /* BLEVCORR was completed, so overscan regions can be trimmed... */
            if (acs_info->verbose) {
                sprintf(MsgText,
                        "Writing out image with trimx = %d,%d, trimy = %d,%d",
                        acs[i].trimx[0], acs[i].trimx[1],
                        acs[i].trimy[0],acs[i].trimy[1]);
                trlmessage(MsgText);

                sprintf(MsgText,"Outputting chip %d",acs[i].chip);
                trlmessage(MsgText);
            }

            /* Output overscan-trimmed, bias-subtracted image
               using the routine 'putSingleGroupSect' from HSTIO V2.0 (or later)
               to write it directly from the full image in memory.
            */
            sizex = x[i].sci.data.nx - (acs[i].trimx[0] + acs[i].trimx[1]);
            sizey = x[i].sci.data.ny - (acs[i].trimy[0] + acs[i].trimy[1]);

            /* Update SIZAXIS keywords to reflect the new trimmed image size. */
            if (PutKeyInt(&x[i].sci.hdr, "SIZAXIS1", sizex,""))
                return (status);
            if (PutKeyInt(&x[i].sci.hdr, "SIZAXIS2", sizey,""))
                return (status);

            /* Now, write out updated trimmed data to disk... */
            putSingleGroupSect(acs_info->output, i+1, &x[i], acs[i].trimx[0],
                               acs[i].trimy[0], sizex, sizey, option);

        } else {
            /* BLEVCORR was not completed, so keep overscan regions... */
            if (acs_info->verbose) {
                sprintf(MsgText,"Writing out FULL image with overscan regions");
                trlmessage(MsgText);

                sprintf(MsgText,"Outputting chip %d",acs[i].chip);
                trlmessage(MsgText);
            }

            putSingleGroup(acs_info->output, i+1, &x[i], option);
        }

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
    free(overscan);
    free(blevcorr);

    return (status);
}


static void dqiMsg (ACSInfo *acs, int extver) {
    int OmitStep (int);
    void PrSwitch (char *, int);
    void PrRefInfo (char *, char *, char *, char *, char *);

    trlmessage ("\n");
    PrSwitch ("dqicorr", acs->dqicorr);

    if (extver == 1 && !OmitStep (acs->dqicorr)) {
        PrRefInfo ("dqitab", acs->bpix.name, acs->bpix.pedigree,
                   acs->bpix.descrip, acs->bpix.descrip2);
    }
}


/* Not used */
/*
static void AtoDMsg (ACSInfo *acs, int extver) {
    int OmitStep (int);
    void PrSwitch (char *, int);
    void PrRefInfo (char *, char *, char *, char *, char *);

    trlmessage ("\n");
    PrSwitch ("atodcorr", acs->atodcorr);

    if (extver == 1 && !OmitStep (acs->atodcorr)) {
        PrRefInfo ("atodtab", acs->atod.name, acs->atod.pedigree,
                   acs->atod.descrip, acs->atod.descrip2);
    }
}
*/


static void BiasMsg (ACSInfo *acs, int extver) {
    int OmitStep (int);
    void PrSwitch (char *, int);
    void PrRefInfo (char *, char *, char *, char *, char *);

    trlmessage ("\n");
    PrSwitch ("biascorr", acs->biascorr);

    if (extver == 1 && !OmitStep (acs->biascorr)) {
        PrRefInfo ("biasfile", acs->bias.name, acs->bias.pedigree,
                   acs->bias.descrip, "");
    }
}


static void BlevMsg (ACSInfo *acs, int extver) {
    void PrSwitch (char *, int);
    void PrRefInfo (char *, char *, char *, char *, char *);
    int OmitStep (int);

    trlmessage ("\n");
    PrSwitch ("blevcorr", acs->blevcorr);

    if (extver == 1 && !OmitStep (acs->blevcorr)) {
        PrRefInfo ("oscntab", acs->oscn.name, acs->oscn.pedigree,
                   acs->oscn.descrip, acs->oscn.descrip2);
    }
}


static void SinkMsg (ACSInfo *acs, int extver) {
    void PrSwitch (char *, int);
    void PrRefInfo (char *, char *, char *, char *, char *);
    int OmitStep (int);

    trlmessage ("\n");
    PrSwitch ("sinkcorr", acs->sinkcorr);

    if (extver == 1 && !OmitStep (acs->sinkcorr)) {
        PrRefInfo ("snkcfile", acs->sink.name, acs->sink.pedigree,
                   acs->sink.descrip, "");
    }
}
