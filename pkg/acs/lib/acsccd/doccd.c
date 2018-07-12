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

** This code is a trimmed down version of CALSTIS1 do2d.c.
*/
# include <string.h>
# include <stdio.h>

#include "hstcal_memory.h"
#include "hstcal.h"
# include "hstio.h"

# include "acs.h"
# include "acsinfo.h"
# include "hstcalerr.h"


static void dqiMsg (ACSInfo *, const int);
static void BiasMsg (ACSInfo *, const int);
static void BlevMsg (ACSInfo *, const int);
static void SinkMsg (ACSInfo *, const int);

/* This variable contains the subarray aperture names for the officially supported
   subarrays since ~Cycle 24.  The polarizer and ramp subarray apertures did NOT 
   change aperture names, even though the dimensions of the data changed to accommodate 
   overscan.  In addition to checking the aperture name, the dimensions will also need 
   to be checked for pol/ramp data to determine if the bias shift correction 
   can be applied.

   The non-pol and non-ramp data adopted new aperture names for easy identification.
   These apertures contain physical overscan.  The new "-2K" subarrays also have virtual 
   overscan so the bias shift correction can be applied to this data. At this time,
   bias shift correction can only be applied to data with both physical and virtual
   overscan so the "fat zero" value can be computed.
*/
static const char *subApertures[] = {"WFC1A-2K", "WFC1B-2K", "WFC2C-2K", "WFC2D-2K",
                                     "WFC1-IRAMPQ", "WFC1-MRAMPQ", "WFC2-MRAMPQ", "WFC2-ORAMPQ",
                                     "WFC1-POL0UV", "WFC1-POL0V", "WFC2-POL0UV", "WFC2-POL0V",
                                     "WFC1-POL60UV", "WFC1-POL60V", "WFC2-POL60UV", "WFC2-POL60V",
                                     "WFC1-POL120UV", "WFC1-POL120V", "WFC2-POL120UV", "WFC2-POL120V"};
static const int numSupportedSubApertures = sizeof(subApertures) / sizeof(subApertures[0]);

static const char *darkSubApertures[] = {"WFC1A-2K", "WFC1B-2K", "WFC2C-2K", "WFC2D-2K",
                                     "WFC1A-1K", "WFC1B-1K", "WFC2C-1K", "WFC2D-1K",
                                     "WFC1A-512", "WFC1B-512", "WFC2C-512", "WFC2D-512",
                                     "WFC1-IRAMPQ", "WFC1-MRAMPQ", "WFC2-ORAMPQ",
                                     "WFC1-POL0UV", "WFC1-POL0V", 
                                     "WFC1-POL60UV", "WFC1-POL60V", 
                                     "WFC1-POL120UV", "WFC1-POL120V", 
                                     "WFC1-SMFL"};

static const int numSupportedDarkSubApertures = sizeof(darkSubApertures) / sizeof(darkSubApertures[0]);

int DoCCD (ACSInfo *acs_info) {

    /* arguments:
       acs   i: calibration switches and info
    */

    extern int status;

    char buff[ACS_FITS_REC+1];

    int to_electrons(ACSInfo *, SingleGroup *);
    void computeDarktime(ACSInfo *, float *);
    int doBias (ACSInfo *, SingleGroup *);
    int biasHistory (ACSInfo *, Hdr *);
    int doBlev (ACSInfo *, SingleGroup *, int, float *, int *, int *);
    int bias_shift_corr(ACSInfo *, int, ...);
    void cross_talk_corr(ACSInfo *, SingleGroup *);
    int doDestripe(ACSInfo *, SingleGroup *, SingleGroup *);
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
    void TimeStamp (char *, char *);
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
            return (status);
        }

        /* Read in keywords from primary header... */

        if (GetKeyBool(x[i].globalhdr, "SUBARRAY", NO_DEFAULT, 0, &subarray)) {
            freeOnExit (&ptrReg);
            return (status);
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
            return (status);
        }
    }}

    int primaryIdx = 0;
    if (PutKeyFlt(x[primaryIdx].globalhdr, "ATODGNA", acs[primaryIdx].atodgain[0], "")) {
        freeOnExit (&ptrReg);
        return (status);
    }
    if (PutKeyFlt(x[primaryIdx].globalhdr, "ATODGNB", acs[primaryIdx].atodgain[1], "")) {
        freeOnExit (&ptrReg);
        return (status);
    }
    if (PutKeyFlt(x[primaryIdx].globalhdr, "ATODGNC", acs[primaryIdx].atodgain[2], "")) {
        freeOnExit (&ptrReg);
        return (status);
    }
    if (PutKeyFlt(x[primaryIdx].globalhdr, "ATODGND", acs[primaryIdx].atodgain[3], "")) {
        freeOnExit (&ptrReg);
        return (status);
    }
    if (PutKeyFlt(x[primaryIdx].globalhdr, "READNSEA", acs[primaryIdx].readnoise[0], "")) {
        freeOnExit (&ptrReg);
        return (status);
    }
    if (PutKeyFlt(x[primaryIdx].globalhdr, "READNSEB", acs[primaryIdx].readnoise[1], "")) {
        freeOnExit (&ptrReg);
        return (status);
    }
    if (PutKeyFlt(x[primaryIdx].globalhdr, "READNSEC", acs[primaryIdx].readnoise[2], "")) {
        freeOnExit (&ptrReg);
        return (status);
    }
    if (PutKeyFlt(x[primaryIdx].globalhdr, "READNSED", acs[primaryIdx].readnoise[3], "")) {
        freeOnExit (&ptrReg);
        return (status);
    }

    trlmessage("\n");

    PrRefInfo("ccdtab", acs[primaryIdx].ccdpar.name, acs[primaryIdx].ccdpar.pedigree,
              acs[primaryIdx].ccdpar.descrip, acs[primaryIdx].ccdpar.descrip2);

    if (CCDHistory(&acs[primaryIdx], x[primaryIdx].globalhdr)) {
        freeOnExit (&ptrReg);
        return (status);
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
            return (status);
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
                return (status);
            }
        }
    }}

    PrSwitch("dqicorr", COMPLETE);

    if (acs_info->printtime)
        TimeStamp("DQICORR complete", acs_info->rootname);

    if (!OmitStep(acs[primaryIdx].dqicorr)) {
        if (dqiHistory(&acs[primaryIdx], x[primaryIdx].globalhdr)) {
            freeOnExit (&ptrReg);
            return (status);
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
                return (status);
            }
        }
    }}

    PrSwitch("biascorr", COMPLETE);

    if (acs_info->printtime)
        TimeStamp ("BIASCORR complete", acs_info->rootname);

    if (!OmitStep(acs[primaryIdx].biascorr)) {
        if (biasHistory(&acs[primaryIdx], x[primaryIdx].globalhdr)) {
            freeOnExit (&ptrReg);
            return (status);
       }
    }
    /************************************************************************/

    /************************************************************************/
    /* convert data to electrons */
    {unsigned int i;
    for (i = 0; i < acs_info->nimsets; i++) {
        if (to_electrons(&acs[i], &x[i])) {
            freeOnExit (&ptrReg);
            return (status);
        }

        if (PutKeyStr(&x[i].sci.hdr, "BUNIT", "ELECTRONS", "")) {
            freeOnExit (&ptrReg);
            return (status);
        }
        if (PutKeyStr(&x[i].err.hdr, "BUNIT", "ELECTRONS", "")) {
            freeOnExit (&ptrReg);
            return (status);
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

    /* The logic here has been re-written (1) to accommodate the new subarrays
       (Ref: ISR ACS 2017-03) which can be bias shift corrected ("*-2K" data only at 
       this time), and (2) the head of the ACS team requested the logic be clarified.  
       This code is verbose and potentially harder to maintain, but the logic is clearer
       for the scientist reading the code.

       NOTE: The variable "done" is overloaded.  When done = overscan[i], this means
       a column was matched in the OSCNTAB for the image and there is physical overscan
       in the image.  Variable "done" is updated in doblev() to indicate the bias has
       been determined from the overscan region rather than using a default bias obtained
       from CCDTAB.
    */

    int done = NO;       
    int driftcorr = NO;     /* true means bias level was corrected for drift */
    Bool isSupportedSubAperture = NO;
    float meanblev = 0.0;   /* mean value of overscan bias (for history) */
    if (acs_info->blevcorr == PERFORM) {

       /* This block handles only the WFC detector */
       if (acs_info->detector == WFC_CCD_DETECTOR) {

          /* If this is subarray data, determine if it is supported for the bias shift 
             correction.  Compare the aperture name against the array of supported
             apertures.  Use the aperture name and not the date as the discriminant 
             as the date is fuzzy - some supported subarrays were acquired before the 
             official name change, and some old style subarrays were acquired after the 
             change.  This scheme works fine for non-polarization and non-ramp data as the
             "*-2K" apertures have both physical and virtual overscan.  These are needed
             to compute the "fat zero" (aka magic square) value needed for the bias
             shift correction.  The "*-1K and *-512" subarrays do not have virtual 
             overscan which is needed at this time for fat zero.

             For polarization and ramp subarrays, the aperture names did not change.  For 
             this data, the dimensions of the image also have to be checked to determine 
             if there are both physical and virtual overscan sections so the bias shift 
             correction can be applied.
          */
          if (acs_info->subarray == YES) {
             {unsigned int i;
             for (i = 0; i < numSupportedSubApertures; i++) {
                 if (strcmp(acs_info->aperture, subApertures[i]) == 0) {

                    /* If this is polarization or ramp data, then the size of the science
                       array must also be checked to make sure there is physical and
                       virtual overscan data.  This does not need to be checked explicitly
                       for non-polarization or non-ramp data, but the comparison works for 
                       this data too.  Only need to check one virtOverscan value in the array. 

                       NOTE:Just because the aperture is supported does not mean the bias
                       shift correction will be applied.  The DS_int and gain criteria also
                       have to be satisfied.
                    */
                    if (virtOverscan[0]) {
                       isSupportedSubAperture = YES;
                       trlmessage("This subarray data is supported for the bias shift correction.");
                    }
                 }
             }}
          }

          /* Only Pre-SM4 WFC data - this is the original bias level subtraction */
          if (acs_info->expstart < SM4MJD) {

             {unsigned int i;
             for (i = 0; i < acs_info->nimsets; i++) {

                 /* The variable done is set to the results of FindOver indicating there is
                    physical overscan. */
                 done = overscan[i];

                 if (doBlev(&acs[i], &x[i], acs[i].chip, &meanblev, &done, &driftcorr)) {
                    freeOnExit (&ptrReg);
                    return (status);
                 }

                 /* Variable done is overwritten in doBlev if bias is determined from overscan region */
                 if (done) {
                    trlmessage("Bias level from overscan has been subtracted;");

                    if (PutKeyFlt(&x[i].sci.hdr, "MEANBLEV", meanblev,
                                   "mean of bias levels subtracted in electrons")) {
                       freeOnExit (&ptrReg);
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
                 {unsigned int j;
                 for (j = 0; j < NAMPS; j++) {
                     if (acs[i].blev[j] != 0.) {
                         sprintf(MsgText, "     bias level of %.6g electrons was subtracted for AMP %c.", acs[i].blev[j], AMPSORDER[j]);
                         trlmessage(MsgText);

                         acs_info->blev[j] = acs[i].blev[j];
                     }
                 }}

                 /* Set this to complete so overscan-trimmed image will be written out. */
                 blevcorr[i] = COMPLETE;

             }} /* End loop over imsets */
  
             PrSwitch("blevcorr", COMPLETE);
  
             if (acs_info->printtime)
                TimeStamp("BLEVCORR complete", acs_info->rootname);

          /* End Pre-SM4 WFC data */
          } else {

             /* Post-SM4 WFC data - Fullframe data with DS_int and gain of 2 will 
                have bias shift, cross talk, and destripe corrections applied.  Supported
                subarray data with DS_int and gain of 2 will have the bias shift applied as
                cross talk and destripe do not apply to subarray data.
             */

             /* The variable done is set based on results of FindOver */
             done = NO;
             if ((overscan[0] == YES) && (overscan[1] == YES)) 
                done = YES;
           
             /* Process full frame data */
             if (acs_info->subarray == NO) {

                if (done) {
                   PrSwitch("blevcorr", PERFORM);

                   /* Only do bias-shift and cross talk corrections of images taken
                      with gain = 2 and in dual-slope integrator mode.  These corrections
                      are done on a per chip basis (one chip per imset).
                   */
                   if (strcmp(acs_info->jwrotype, "DS_INT") == 0 &&
                           acs_info->ccdgain == 2) {
                      trlmessage("Performing bias-shift and cross talk corrections for full frame data.");

                      if (bias_shift_corr(acs_info, acs_info->nimsets, &x[0], &x[1])) {
                         freeOnExit (&ptrReg);
                         return status;
                      }
        
                      {unsigned int i;
                      for (i = 0; i < acs_info->nimsets; i++) {
                          cross_talk_corr(&acs[i], &x[i]);
                      }}
                   } else {
                      trlmessage("WFC readout type/gain not set as needed, no bias shift nor cross talk correction done for full frame data.");
                   }

                   trlmessage("Performing stripe removal and bias level subtraction for full frame data.");

                   if (doDestripe(acs_info, &x[0], &x[1])) {
                      freeOnExit (&ptrReg);
                      return status;
                   }

                   {unsigned int i;
                   for (i = 0; i < acs_info->nimsets; i++) {
                       blevcorr[i] = COMPLETE;
                   }}

                   PrSwitch("blevcorr", COMPLETE);
   
                   if (acs_info->printtime) {
                       TimeStamp("BLEVCORR complete", acs->rootname);
                   }
                } else {
                   trlmessage("Overscan missing, no destriping or bias level subtraction possible for full frame data.");
                }

                driftcorr = NO;

             /* Process subarray data */
             } else {

                /* Supported subarray data with physical and virtual overscan can be
                   processed in the same manner as full frame data.  The data must
                   also be DS_int with a gain of 2.

                   No need to check "done" here as the isSupportedSubAperture has been set based upon
                   tha array names of supported subarrays AND the determination they have virtual overscan 
                   according to the image size (e.g., WFC1A-2K with size 2072x2068).
                */
                if ((isSupportedSubAperture == YES) && (strcmp(acs_info->jwrotype, "DS_INT") == 0) &&
                              (acs_info->ccdgain == 2)) {
                   PrSwitch("blevcorr", PERFORM);

                   trlmessage("Performing bias-shift correction for subarray data.");
                      
                   /* Only bias shift correction is done for subarray data. For a supported
                      subarray, there is only one amp on a single chip in use. 
                   */
                   if (bias_shift_corr(acs_info, 1, &x[0])) {
                      freeOnExit (&ptrReg);
                      return status;
                   }
                     
                   /* There is only one imset for supported subarrays */
                   blevcorr[primaryIdx] = COMPLETE;
                   PrSwitch("blevcorr", COMPLETE);
   
                   if (acs_info->printtime) {
                      TimeStamp("BLEVCORR complete", acs->rootname);
                   }

                   driftcorr = NO;

                   /* Unsupported subarray apertures and other oddities.  These
                      are processed using the original/old bias correction algorithm.
                   */ 
                } else {
                   trlmessage("WFC readout is not a supported subarray or type/gain not set as needed for new bias level algorithm for subarray data.");
                   trlmessage("Data to be processed with original bias level algorithm.");

                   {unsigned int i;
                   for (i = 0; i < acs_info->nimsets; i++) {

                       /* This is set based on results of FindOver */
                       done = overscan[i];

                       if (doBlev(&acs[i], &x[i], acs[i].chip, &meanblev, &done, &driftcorr)) {
                          freeOnExit (&ptrReg);
                          return (status);
                       }

                       /* Variable done is overwritten in doBlev if bias is determined from overscan region */
                       if (done) {
                          trlmessage("Bias level from overscan has been subtracted;");
 
                          if (PutKeyFlt(&x[i].sci.hdr, "MEANBLEV", meanblev,
                                         "mean of bias levels subtracted in electrons")) {
                             freeOnExit (&ptrReg);
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
                       {unsigned int j;
                       for (j = 0; j < NAMPS; j++) {
                           if (acs[i].blev[j] != 0.) {
                               sprintf(MsgText, "     bias level of %.6g electrons was subtracted for AMP %c.", acs[i].blev[j], AMPSORDER[j]);
                               trlmessage(MsgText);

                               acs_info->blev[j] = acs[i].blev[j];
                           }
                       }}

                       /* Set this to complete so overscan-trimmed image will be written out. */
                       blevcorr[i] = COMPLETE;
                   }}
  
                   PrSwitch("blevcorr", COMPLETE);

                   if (acs_info->printtime)
                      TimeStamp("BLEVCORR complete", acs_info->rootname);
 
                } /* End subarray processing */

             } /* End of full frame/subarray data */

          } /* End of Pre-/Post-SM4 data */

       } /* End of only WFC detector */

       /* All HRC data - use the original bias level subtraction correction */
       else {
            {unsigned int i;
            for (i = 0; i < acs_info->nimsets; i++) {
                /* This is set based on results of FindOver */
                done = overscan[i];

                if (doBlev(&acs[i], &x[i], acs[i].chip, &meanblev, &done, &driftcorr)) {
                   freeOnExit (&ptrReg);
                   return (status);
                }
  
                /* Variable done is overwritten in doBlev if bias is determined from overscan region */
                if (done) {
                   trlmessage("Bias level from overscan has been subtracted;");

                   if (PutKeyFlt(&x[i].sci.hdr, "MEANBLEV", meanblev,
                                  "mean of bias levels subtracted in electrons")) {
                      freeOnExit (&ptrReg);
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
                {unsigned int j;
                for (j = 0; j < NAMPS; j++) {
                    if (acs[i].blev[j] != 0.) {
                        sprintf(MsgText, "     bias level of %.6g electrons was subtracted for AMP %c.", acs[i].blev[j], AMPSORDER[j]);
                        trlmessage(MsgText);

                        acs_info->blev[j] = acs[i].blev[j];
                    }
                }}

                /* Set this to complete so overscan-trimmed image will be written out. */
                blevcorr[i] = COMPLETE;
            }}

            PrSwitch("blevcorr", COMPLETE);
 
            if (acs_info->printtime)
               TimeStamp("BLEVCORR complete", acs_info->rootname);

       } /* End for HRC data */
 
    } /* End of BLEVCORR = PERFORM */

    if (!OmitStep(acs[primaryIdx].blevcorr)) {
        if (blevHistory(&acs[primaryIdx], x[primaryIdx].globalhdr, done, driftcorr)) {
            freeOnExit (&ptrReg);
            return (status);
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
        return (status);
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
            sprintf(MsgText, "\nFull-well saturation flagging being performed for imset %d.\n", i+1);
            trlmessage(MsgText);
            if (doFullWellSat(&acs[i], &x[i])) {
                freeOnExit (&ptrReg);
                return (status);
            }
        } else {
            sprintf(MsgText, "\nNo Full-well saturation flagging being performed for imset %d.\n", i+1);
            trlwarn(MsgText);
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
                return (status);
            }
        }}

        if (done) {
            if (noiseHistory(x[primaryIdx].globalhdr)) {
                freeOnExit (&ptrReg);
                return (status);
            }

            trlmessage("\n    Uncertainty array initialized,");
            buff[0] = '\0';

            sprintf(MsgText, "    readnoise =");
            {unsigned int i;
            for (i=0; i < NAMPS-1; i++) {
                if (acs[primaryIdx].readnoise[i] > 0) {
                    sprintf(buff, "%.5g,",acs[primaryIdx].readnoise[i]);
                    strcat(MsgText, buff);
                }
            }}

            if (acs[primaryIdx].readnoise[NAMPS-1] > 0) {
                sprintf(buff, "%.5g",acs[primaryIdx].readnoise[NAMPS-1]);
                strcat(MsgText, buff);
            }
            trlmessage(MsgText);

            sprintf(MsgText, "    gain =");
            {unsigned int i;
            for (i=0; i < NAMPS-1; i++) {
                if (acs[primaryIdx].atodgain[i] > 0) {
                    sprintf(buff, "%.5g,",acs[primaryIdx].atodgain[i]);
                    strcat(MsgText, buff);
                }
            }}

            if (acs[primaryIdx].atodgain[NAMPS-1] > 0) {
                sprintf(buff, "%.5g",acs[primaryIdx].atodgain[NAMPS-1]);
                strcat(MsgText, buff);
            }
            trlmessage(MsgText);
            sprintf(MsgText, "   default bias levels =");

            {unsigned int i;
            for (i=0; i < NAMPS-1; i++) {
                if (acs[primaryIdx].ccdbias[i] > 0) {
                    sprintf(buff, "%.5g,",
                            acs[primaryIdx].ccdbias[i] * acs[primaryIdx].atodgain[i]);
                    strcat(MsgText, buff);
                }
            }}

            if (acs[primaryIdx].ccdbias[NAMPS-1] > 0) {
                sprintf(buff, "%.5g",
                        acs[primaryIdx].ccdbias[NAMPS-1] * acs[primaryIdx].atodgain[NAMPS-1]);
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
    SinkMsg (&acs[primaryIdx], nextver);

    if (acs->sinkcorr == PERFORM) {
        {unsigned int i;
        for (i = 0; i < acs_info->nimsets; i++) {
            if (SinkDetect(&acs[i], &x[i])) {
                freeOnExit (&ptrReg);
                return (status);
            }
        }}

        PrSwitch("sinkcorr", COMPLETE);

        if (acs_info->printtime)
            TimeStamp("SINKCORR complete", acs_info->rootname);
    }

    if ((status = sinkHistory(&acs[primaryIdx], x[primaryIdx].globalhdr))) {
        freeOnExit (&ptrReg);
        return (status);
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
            if (PutKeyInt(&x[i].sci.hdr, "SIZAXIS1", sizex,"")) {
                freeOnExit (&ptrReg);
                return (status);
            }
            if (PutKeyInt(&x[i].sci.hdr, "SIZAXIS2", sizey,"")) {
                freeOnExit (&ptrReg);
                return (status);
            }

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

        /* If putSingleGroup/putSingleGroupSect failed above, then error out. */
        if (hstio_err()) {
            sprintf(MsgText, "Couldn't write imset %d.", i+1);
            trlerror(MsgText);
            freeOnExit (&ptrReg);
            return (status = WRITE_FAILED);
        }

        if (acs_info->printtime)
            TimeStamp ("Output written to disk", acs_info->rootname);

        /* x[i] memory will be freed on exit */
    }}

    freeOnExit (&ptrReg);

    return (status);
}

/*
 * This routine updates the DARKTIME value such that every BLV_TMP file will
 * have the correct DARKTIME keyword.
 */
void computeDarktime(ACSInfo *acs, float *darktime) {

    float darktimeFromHeader;  /* base darktime value from FITS keyword DARKTIME in image header */
    float darktimeOffset;      /* overhead offset time (s) based upon full-frame/subarray and post-flashed/unflashed */

    /* Get the DARKTIME FITS keyword value stored in the ACSInfo structure */
    darktimeFromHeader = (float)acs->darktime;

    /*
       The overhead offset time is a function of full-frame vs subarray and post-flash vs unflashed 
       and has been determined empirically for CCD data.  Both the unflashed and post-flash
       overhead values for full-frame or subarray were extracted from the calibration file 
       during the table read.  Now it is just necessary to determine which offset actually
       applies.
    */

    /* Unflashed observation */
    darktimeOffset = acs->overhead_unflashed;
 
    /* Post-flashed observation */
    if ((acs->flashdur > 0.0) && (strcmp(acs->flashstatus, "SUCCESSFUL") == 0)) {
          darktimeOffset = acs->overhead_postflashed;
    }

    /* 
       Compute the final darktime based upon the date of full-frame or subarray data.
       The full-frame overhead offset is applicable to all data post-SM4.  The subarray 
       overhead offset is applicable to all data post-CYCLE24 and ONLY for supported
       subarrays. 
  
       Effectively the additive factor only applies to ACS/WFC as the HRC was no longer
       operational by SM4MJD or CYCLE24, and SBC is a MAMA detector.
    */
    *darktime = darktimeFromHeader;  /* Default */

    /* Full-frame data */
    if (acs->subarray == NO) {
      if (acs->expstart >= SM4MJD) {
        *darktime = darktimeFromHeader + darktimeOffset;
      }
      sprintf(MsgText, "Full Frame adjusted Darktime: %f\n", *darktime);
      trlmessage(MsgText);

    /* Subarray data */
    } else {
      if (acs->expstart >= CYCLE24) {
        for (unsigned int i = 0; i < numSupportedDarkSubApertures; i++) {
            if (strcmp(acs->aperture, darkSubApertures[i]) == 0) {
                *darktime = darktimeFromHeader + darktimeOffset;
                sprintf(MsgText, "Supported Subarray adjusted Darktime: %f for aperture: %s\n", *darktime, darkSubApertures[i]);
                trlmessage(MsgText);
                break;
            }
        }
      }
    }
    sprintf(MsgText, "DARKTIME from SCI header: %f  Offset from CCDTAB: %f  Final DARKTIME: %f", darktimeFromHeader, darktimeOffset, *darktime);
    trlmessage(MsgText);
}


static void dqiMsg (ACSInfo *acs, const int extver) {
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


static void BiasMsg (ACSInfo *acs, const int extver) {
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


static void BlevMsg (ACSInfo *acs, const int extver) {
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


static void SinkMsg (ACSInfo *acs, const int extver) {
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
