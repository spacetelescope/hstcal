/* Do basic CCD image reduction for the current set of extensions.
   If we're processing the first set (extver = 1), then the primary
   header will be updated with history information, and the calibration
   switches in the header will be reset from "PERFORM" to "COMPLETE".

   Warren Hack, 1998 June 1:
   Revised to only perform CCD operations on ACS observations.

   Warren Hack, 1999 Jan 7:
   Revised to only output overscan-trimmed image when BLEVCORR is
   performed successfully.  In other cases, full image is written out.

   Howard Bushouse, 2000 Aug 29:
   Initial revisions for WFC3 use.

   H.Bushouse, 2001 May 7:
   Added support for new processing step: post-flash
   (copy of changes made to calacs).

   H.Bushouse, 2001 Nov 16:
   Updates to track CALACS changes - Changed "postflsh" to "flshcorr"
   everywhere; updated to use CCDBIAS[A,B,C,D] for default bias levels
   and to use CCDOFST to select CCDBIAS from revised CCDTAB file;
   reports bias levels for each amp.

   H.Bushouse, 2002 Jan 28:
   Updated to remove extra WFC3 UVIS CCD serial virtual overscan
   regions from calibrated image data before writing out trimmed
   image.

   H.Bushouse, 2002 Jun 17:
   Added update of SIZAXIS keywords to reflect trimmed size
   (following CALACS changes).

   H.Bushouse, 2003 Oct 23:
   Modified to correctly handle the output of WFC3 subarray images,
   in which there's no serial virtual overscan to remove. Also
   added message info giving mean bias level for each amp (CALACS
   changes).

   M. Sosey, 2013 July 3
   Addded new FLUXCORR step to scale flux units in both chips to
   the same scale using the ratio of PHTFLAM1/PHTFLAM2 ->PHTRATIO

   M. Sosey June 2015
   If DQICORR is performed, then the SINK pixel mask is also propgated
   This is part of the UVIS 2 update. SINK pixels are not currently
   detected for subarrays or binned data.

    C Jones May 2016
    Made changes to the input of the sink-detect so that it will accept
    sub-array data.

    M. sosey July 2016
    Moved the subarray codes to a general function and changed the calculation
    of where the subarray is. Also moved the DQ bit to the singlegroup extension

    M. De La Pena February 2022
    Modified to apply the full-well saturation flags stored as an image
    to the data in doFullWellSat() instead of in the doDQI step.

 */

# include <string.h>
# include <stdio.h>
# include "hstcal.h"
# include "hstio.h"
# include "wf3.h"
# include "wf3info.h"
# include "hstcalerr.h"
# include "doccd.h"

int PutKeyDbl(Hdr *, char *, double , char *);
int PutKeyStr(Hdr *, char *, char *, char *);

int DoCCD (WF3Info *wf3, int extver) {

    /* arguments:
       WF3Info *wf3	i: calibration switches and info
       int extver	i: "imset" number, the current set of extensions
     */

    extern int status;

    SingleGroup x;	/* used for both input and output */
    SingleGroup fullarray;
    int option = 0;
    float meanblev;	/* mean value of overscan bias (for history) */
    float meanflash; /* mean value of post-flash image (for history) */
    int driftcorr;	/* true means bias level was corrected for drift */
    int done;	/* true means the input SingleGroup has been freed */
    int sizex, sizey;	/* size of output image */

    int i, j, x1, dx;		/* loop index */
    int overscan;
    int blevcorr;
    char buff[SZ_FITS_REC+1];
    Bool subarray;


    /*========================Start Code here =======================*/
    initSingleGroup (&x);
    if (wf3->printtime)
        TimeStamp ("Open SingleGroup now...", "");

    /* OPEN THE INPUT IMAGE.
     ** FOR WF3: Keep this in memory throughout processing.
     ** Read in reference images line-by-line within the individual
     ** processing step functions and pass along modified input image
     ** from one step to the next...
     */
    getSingleGroup (wf3->input, extver, &x);
    if (hstio_err())
        return (status = OPEN_FAILED);

    if (wf3->printtime)
        TimeStamp ("Input read into memory", wf3->rootname);

    /* GET HEADER INFO THAT VARIES FROM IMSET TO IMSET NECESSARY
     **	FOR READING CCDTAB.
     */
    if (GetGrp (wf3, &x.sci.hdr)) {
        freeSingleGroup (&x);
        return (status);
    }

    /* READ IN KEYWORDS FROM PRIMARY HEADER... */
    if (GetKeyBool (x.globalhdr, "SUBARRAY", NO_DEFAULT, 0, &subarray))
        return (status);
    if (subarray)
        wf3->subarray = YES;
    else
        wf3->subarray = NO;

    /* FOR THE CCD, UPDATE PRIMARY HEADER KEYWORDS.
     * ALSO RESET CRCORR IF THERE'S ONLY ONE IMAGE SET.
     * GET VALUES FROM TABLES, USING SAME FUNCTION USED IN WF3CCD. */
    if (GetCCDTab (wf3, x.sci.data.nx, x.sci.data.ny)) {
        freeSingleGroup (&x);
        return (status);
    }

    if (PutKeyFlt (x.globalhdr, "ATODGNA", wf3->atodgain[0], ""))
        return (status);
    if (PutKeyFlt (x.globalhdr, "ATODGNB", wf3->atodgain[1], ""))
        return (status);
    if (PutKeyFlt (x.globalhdr, "ATODGNC", wf3->atodgain[2], ""))
        return (status);
    if (PutKeyFlt (x.globalhdr, "ATODGND", wf3->atodgain[3], ""))
        return (status);
    if (PutKeyFlt (x.globalhdr, "READNSEA", wf3->readnoise[0], ""))
        return (status);
    if (PutKeyFlt (x.globalhdr, "READNSEB", wf3->readnoise[1], ""))
        return (status);
    if (PutKeyFlt (x.globalhdr, "READNSEC", wf3->readnoise[2], ""))
        return (status);
    if (PutKeyFlt (x.globalhdr, "READNSED", wf3->readnoise[3], ""))
        return (status);
    trlmessage ("");

    PrRefInfo ("ccdtab", wf3->ccdpar.name, wf3->ccdpar.pedigree,
            wf3->ccdpar.descrip, wf3->ccdpar.descrip2);

    if (extver == 1) {
        if (CCDHistory (wf3, x.globalhdr))
            return (status);
    }

    /* GET OVERSCAN REGION INFORMATION FROM OSCNTAB */
    if (FindOverscan (wf3, x.sci.data.nx, x.sci.data.ny, &overscan))
        return (status);

    /* FILL IN THE ERROR ARRAY, IF IT INITIALLY CONTAINS ALL ZEROS. */
    if (wf3->noiscorr == PERFORM) {
        if (doNoise (wf3, &x, &done))
            return (status);
        if (done) {
            if (extver == 1) {
                if (noiseHistory (x.globalhdr))
                    return (status);
            }
            trlmessage ("    Uncertainty array initialized,");
            buff[0] = '\0';

            sprintf(MsgText, "    readnoise =");
            for (i=0; i < NAMPS-1; i++) {
                if (wf3->readnoise[i] > 0) {
                    sprintf (buff, "%.5g,",wf3->readnoise[i]);
                    strcat (MsgText, buff);
                }
            }
            if (wf3->readnoise[NAMPS-1] > 0) {
                sprintf(buff, "%.5g",wf3->readnoise[NAMPS-1]);
                strcat (MsgText, buff);
            }
            trlmessage (MsgText);

            sprintf(MsgText, "    gain =");
            for (i=0; i < NAMPS-1; i++) {
                if (wf3->atodgain[i] > 0) {
                    sprintf(buff, "%.5g,",wf3->atodgain[i]);
                    strcat (MsgText, buff);
                }
            }
            if (wf3->atodgain[NAMPS-1] > 0) {
                sprintf(buff, "%.5g",wf3->atodgain[NAMPS-1]);
                strcat (MsgText, buff);
            }
            trlmessage (MsgText);

            sprintf (MsgText, "    default bias levels = ");
            for (i=0; i < NAMPS-1; i++) {
                if (wf3->ccdbias[i] > 0) {
                    sprintf (buff, "%.5g,", wf3->ccdbias[i]);
                    strcat (MsgText, buff);
                }
            }
            if (wf3->ccdbias[NAMPS-1] > 0) {
                sprintf (buff, "%.5g", wf3->ccdbias[NAMPS-1]);
                strcat (MsgText, buff);
            }
            trlmessage (MsgText);

            if (wf3->printtime)
                TimeStamp ("Uncertainty array initialized", wf3->rootname);
        }
    }

    /* DATA QUALITY INITIALIZATION AND (FOR THE CCDS) CHECK SATURATION. */
    dqiMsg (wf3, extver);
    if (wf3->dqicorr == PERFORM || wf3->dqicorr == DUMMY) {
        if (doDQI (wf3, &x))
            return (status);
        PrSwitch ("dqicorr", COMPLETE);
        if (wf3->printtime)
            TimeStamp ("DQICORR complete", wf3->rootname);
    }
    if (extver == 1 && !OmitStep (wf3->dqicorr))
        if (dqiHistory (wf3, x.globalhdr))
            return (status);

    /* ANALOG TO DIGITAL CORRECTION. */
    AtoDMsg (wf3, extver);
    if (wf3->atodcorr == PERFORM) {
        if (doAtoD (wf3, &x))
            return (status);
        PrSwitch ("atodcorr", COMPLETE);
        if (wf3->printtime)
            TimeStamp ("ATODCORR complete", wf3->rootname);
    }
    if (extver == 1 && !OmitStep (wf3->atodcorr)){
        if (atodHistory (wf3, x.globalhdr)){
            return (status);
        }
    }

    /* SUBTRACT BIAS FROM OVERSCAN. */
    BlevMsg (wf3, extver);
    blevcorr = PERFORM;
    if (wf3->blevcorr == PERFORM) {
        /* This is set based on results of FindOver */
        done = overscan;

        if (doBlev(wf3, &x, wf3->chip, &meanblev, &done, &driftcorr))
            return (status);

        if (done) {
            trlmessage ("Bias level from overscan has been subtracted;");
            if (PutKeyFlt (&x.sci.hdr, "MEANBLEV", meanblev,
                        "mean of bias levels subtracted")) {
                return (status);
            }
            sprintf (MsgText,
                    "     mean of bias levels subtracted was %.6g.",
                    meanblev);
            trlmessage (MsgText);
        } else {
            trlmessage ("Default bias level from CCDTAB was subtracted.");
        }

        /* Provide immediate feedback to the user on the values computed
         ** for each AMP, but only for those amps used for the chip being
         ** processed. */
        for (i=0; i < 4; i++) {
            if (wf3->blev[i] != 0.) {
                sprintf (MsgText,
                        "     bias level of %.6g was subtracted for AMP %c.",
                        wf3->blev[i], AMPSORDER[i]);
                trlmessage (MsgText);
            }
        }

        PrSwitch ("blevcorr", COMPLETE);

        /* Set this to complete so overscan-trimmed image will be
         ** written out.  */
        blevcorr = COMPLETE;

        if (wf3->printtime)
            TimeStamp ("BLEVCORR complete", wf3->rootname);
    }
    if (extver == 1 && !OmitStep (wf3->blevcorr))
        if (blevHistory (wf3, x.globalhdr, done, driftcorr))
            return (status);

    /* SUBTRACT BIAS IMAGE. */
    BiasMsg (wf3, extver);
    if (wf3->biascorr == PERFORM) {
        if (doBias (wf3, &x))
            return (status);
        PrSwitch ("biascorr", COMPLETE);
        if (wf3->printtime)
            TimeStamp ("BIASCORR complete", wf3->rootname);
    }
    if (extver == 1 && !OmitStep (wf3->biascorr))
        if (biasHistory (wf3, x.globalhdr))
            return (status);

    /* Apply the saturation image.
       Strictly speaking, the application of the full-well saturation image is
       not a calibration step (i.e., there is no SATCORR), but the application
       of a 2D image to flag pixels versus using a single scalar to flag
       saturated pixels as previously done in DQICORR will be done in doFullWellSat()
       after BLEVCORR and BIASCORR.  This correction should only be done if both
       BLEVCORR and BIASCORR have been performed.  This flagging is only applicable
       for the UVIS. */

    if (wf3->biascorr == PERFORM && wf3->blevcorr == PERFORM) {
        SatMsg (wf3, extver);
        sprintf(MsgText, "\nFull-well saturation flagging being performed.");
        trlmessage(MsgText);
        if (doFullWellSat(wf3, &x)) {
            return (status);
        }
    } else {
        sprintf(MsgText, "\nNo Full-well saturation flagging being performed.\n");
        trlwarn(MsgText);
    }

    /*UPDATE THE SINK PIXELS IN THE DQ MASK OF BOTH SCIENCE IMAGE SETS
     IT'S DONE HERE WITH ONE CALL TO THE FILE BECAUSE THEY NEED TO BE
     PROCESSED IN THE RAZ FORMAT JAY USES, THOUGH ONLY ONE CHIP DONE HERE

     IF THE IMAGE IS A SUBARRAY THEN IT MUST BE PLACED INSIDE A FULL FRAME
     ARRAY FOR SINK PIXEL DETECTION

     BINNED IMAGES ARE NOT SUPPORTED AT ALL
    */
    if (wf3->dqicorr == PERFORM) {
        if (checkBinned(&x)){
            sprintf(MsgText,"Binned data not supported for Sink Pixel detection");
            trlmessage(MsgText);
        } else {
            if (wf3->subarray){

                /* CKJ
                 * If the SCI data is a sub array then we are going to
                 * copy the SCI and DQ sub-arrays into a temporary version
                 * of a full-array and pass that to SinkDetect.
                 *
                 * MLS
                 * I created a library function that could be used
                 * by any of the modules
                */

                /* Make temporary full group
                   use the sinkfile as the reference
                   sub2full will reset the pixel values
                */
                initSingleGroup(&fullarray);
                allocSingleGroup(&fullarray,RAZ_COLS/2,RAZ_ROWS, True);
                fullarray.group_num=x.group_num;
                CreateEmptyChip(wf3, &fullarray);
                if (x.group_num == 2){ /*post blev*/
                    if (PutKeyDbl(&fullarray.sci.hdr, "LTV2", 0.0, "offset in Y to light start")) {
                      trlmessage("Error putting LTV1 keyword in header");
                      return (status=HEADER_PROBLEM);
                    }
                }

                if (Sub2Full(wf3, &x, &fullarray, 1, 1, 0))
                  return(status);

                if (SinkDetect(wf3, &fullarray))
                   return(status);

                /*place the marked dq pixels back into x*/
                if (Full2Sub(wf3, &x, &fullarray, 1, 0, 0))
                    return(status);

                /* Free up the group now that we are done with it*/
                freeSingleGroup(&fullarray);

           } else {/* Must be a fullframe image so just do this */
                if (SinkDetect(wf3, &x))
                   return(status);
          }
      }
    }

    /* SUBTRACT POST-FLASH IMAGE. */
    FlashMsg (wf3, extver);
    if (wf3->flashcorr == PERFORM) {

        /* Initialize this to a set value... */
        meanflash = 0.0;

        if (doFlash (wf3, &x, &meanflash))
            return (status);

        /* Report mean of post-flash image subtracted from science image,
         ** if it was performed...*/
        if (meanflash > 0.) {
            sprintf (MsgText, "Mean of post-flash image (MEANFLSH) = %g",
                    meanflash);
            trlmessage(MsgText);

            /* Store mean flash value as keyword */
            if (PutKeyFlt (&x.sci.hdr, "MEANFLSH", meanflash,
                        "mean of post-flash values subtracted"))
                return (status);

            PrSwitch ("flshcorr", COMPLETE);

        } else {
            PrSwitch ("flshcorr", SKIPPED);
        }

        if (wf3->printtime) TimeStamp ("FLSHCORR complete", wf3->rootname);
    }

    if (extver == 1 && !OmitStep (wf3->flashcorr)) {
        if (flashHistory (wf3, x.globalhdr))
            return (status);
    }

    /* WRITE THIS IMSET TO THE OUTPUT FILE.  IF EXTVER IS ONE, THE
       CAL_VER AND FILENAME KEYWORDS WILL BE UPDATED, AND THE PRIMARY
       HEADER WILL BE WRITTEN.  */

    if (extver == 1) {
        UCalVer (x.globalhdr);
        UFilename (wf3->output, x.globalhdr);
    }

    /* IF BLEVCORR WAS PERFORMED, THEN OUTPUT TRIMMED IMAGE,
       OTHERWISE OUTPUT FULL IMAGE...
     */
    if (blevcorr == COMPLETE) {

        /* BLEVCORR WAS COMPLETED, SO OVERSCAN REGIONS CAN BE TRIMMED... */
        if (wf3->verbose) {
            sprintf (MsgText,
                    "Writing out image with trimx = %d,%d,%d,%d, trimy = %d,%d",
                    wf3->trimx[0],wf3->trimx[1],wf3->trimx[2],
                    wf3->trimx[3],wf3->trimy[0],wf3->trimy[1]);
            trlmessage(MsgText);

            sprintf(MsgText,"Outputting chip %d",wf3->chip);
            trlmessage(MsgText);
        }

        /* OUTPUT OVERSCAN-TRIMMED, BIAS-SUBTRACTED IMAGE
         **  USING THE ROUTINE 'PUTSINGLEGROUPSECT' FROM HSTIO
         **  TO WRITE IT DIRECTLY FROM THE FULL IMAGE IN MEMORY.  */
        sizex = x.sci.data.nx - (wf3->trimx[0] + wf3->trimx[1] +
                wf3->trimx[2] + wf3->trimx[3]);
        sizey = x.sci.data.ny - (wf3->trimy[0] + wf3->trimy[1]);

        /* For WFC3 full-frame, 4-amp readouts, we must recopy a subset of
         ** the image in order to remove the serial virtual overscan region
         ** in the middle of the image. */
        if (strlen(wf3->ccdamp) == 4) {
            x1 = x.sci.data.nx/2 + wf3->trimx[2];
            dx = wf3->trimx[2] + wf3->trimx[3];
            for (j=0; j < x.sci.data.ny; j++) {
                for (i=x1; i < x.sci.data.nx; i++) {
                    Pix (x.sci.data, i-dx, j) = Pix (x.sci.data, i, j);
                    Pix (x.err.data, i-dx, j) = Pix (x.err.data, i, j);
                    Pix (x.dq.data,  i-dx, j) = Pix (x.dq.data,  i, j);
                }
            }
        }

        /* UPDATE SIZAXIS KEYWORDS TO REFLECT THE NEW TRIMMED IMAGE SIZE */
        if (PutKeyInt (&x.sci.hdr, "SIZAXIS1", sizex, ""))
            return (status);
        if (PutKeyInt (&x.sci.hdr, "SIZAXIS2", sizey, ""))
            return (status);

        /* NOW WRITE OUT THE APPROPRIATE REMAINING SECTION */
        putSingleGroupSect (wf3->output, extver, &x, wf3->trimx[0],
                wf3->trimy[0], sizex, sizey, option);

    } else {

        /* BLEVCORR WAS NOT COMPLETED, SO KEEP OVERSCAN REGIONS... */
        if (wf3->verbose) {
            sprintf(MsgText,"Writing out FULL image with overscan regions");
            trlmessage(MsgText);

            sprintf(MsgText,"Outputting chip %d",wf3->chip);
            trlmessage(MsgText);
        }

        putSingleGroup (wf3->output, extver, &x, option);

    }

    if (hstio_err()) {
        sprintf (MsgText, "Couldn't write imset %d.", extver);
        trlerror (MsgText);
        return (status = WRITE_FAILED);
    }
    if (wf3->printtime)
        TimeStamp ("Output written to disk", wf3->rootname);

    /* Free x, which still has memory allocated. */
    freeSingleGroup (&x);

    return (status);
}


/*check and see if the image is binned at all, return 1 if binned*/
int checkBinned (SingleGroup *x){
    int rsize=1;
    int bin[2];
    int corner[2];
    int i;

    /*initialize*/
    for(i=0;i<=1;i++){
        bin[i]=0;
        corner[i]=0;
    }
    GetCorner(&x->sci.hdr, rsize, bin, corner);
    if (bin[0] > 1 || bin[1] > 1){
        return 1;
    } else {
        return 0;
    }

}

static void AtoDMsg (WF3Info *wf3, int extver) {

    int OmitStep (int);
    void PrSwitch (char *, int);
    void PrRefInfo (char *, char *, char *, char *, char *);

    trlmessage ("");
    PrSwitch ("atodcorr", wf3->atodcorr);

    if (extver == 1 && !OmitStep (wf3->atodcorr)) {

        PrRefInfo ("atodtab", wf3->atod.name, wf3->atod.pedigree,
                wf3->atod.descrip, wf3->atod.descrip2);
    }
}

static void BiasMsg (WF3Info *wf3, int extver) {

    int OmitStep (int);
    void PrSwitch (char *, int);
    void PrRefInfo (char *, char *, char *, char *, char *);

    trlmessage ("");
    PrSwitch ("biascorr", wf3->biascorr);

    if (extver == 1 && !OmitStep (wf3->biascorr)) {

        PrRefInfo ("biasfile", wf3->bias.name, wf3->bias.pedigree,
                wf3->bias.descrip, "");
    }
}

static void SatMsg (WF3Info *wf3, int extver) {

    int OmitStep (int);
    void PrSwitch (char *, int);
    void PrRefInfo (char *, char *, char *, char *, char *);

    trlmessage ("");
    if (extver == 1 && !OmitStep (wf3->biascorr)) {

        PrRefInfo ("satufile", wf3->satmap.name, wf3->satmap.pedigree,
                wf3->satmap.descrip, "");
    }
}

static void FlashMsg (WF3Info *wf3, int extver) {

    int OmitStep (int);
    void PrSwitch (char *, int);
    void PrRefInfo (char *, char *, char *, char *, char *);

    trlmessage ("");
    PrSwitch ("flshcorr", wf3->flashcorr);

    if (extver == 1 && !OmitStep (wf3->flashcorr)) {

        PrRefInfo ("flshfile", wf3->flash.name, wf3->flash.pedigree,
                wf3->flash.descrip, "");
    }
}

static void BlevMsg (WF3Info *wf3, int extver) {

    void PrSwitch (char *, int);
    void PrRefInfo (char *, char *, char *, char *, char *);
    int OmitStep (int);

    trlmessage ("");
    PrSwitch ("blevcorr", wf3->blevcorr);

    if (extver == 1 && !OmitStep (wf3->blevcorr)) {

        PrRefInfo ("oscntab", wf3->oscn.name, wf3->oscn.pedigree,
                wf3->oscn.descrip, wf3->oscn.descrip2);
    }

}

static void dqiMsg (WF3Info *wf3, int extver) {

    int OmitStep (int);
    void PrSwitch (char *, int);
    void PrRefInfo (char *, char *, char *, char *, char *);

    trlmessage ("");
    PrSwitch ("dqicorr", wf3->dqicorr);

    if (extver == 1 && !OmitStep (wf3->dqicorr)) {

        PrRefInfo ("dqitab", wf3->bpix.name, wf3->bpix.pedigree,
                wf3->bpix.descrip, wf3->bpix.descrip2);
    }
}
