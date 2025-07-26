#include <string.h>
#include "acs.h"
#include "acsinfo.h"


/*
 23-Jul-2025 PLL: Created performBlevCorr() in order to move code out of
     doccd.c to simplify the complexity of doccd.c
*/
int performBlevCorr(ACSInfo *acs_info, ACSInfo *acs, SingleGroup *x,
                    int *overscan, int *virtOverscan,
                    int *blevcorr, int *driftcorr, int *done) {
    extern int status;
    const int primaryIdx = 0;

    int doBlev (ACSInfo *, SingleGroup *, int, int *, int *);
    void blevSubTrlMessage(ACSInfo *, ACSInfo *);
    int bias_shift_corr(ACSInfo *, int, ...);
    void cross_talk_corr(ACSInfo *, SingleGroup *);
    int doDestripe(ACSInfo *, SingleGroup *, SingleGroup *);
    int isValidBiasShiftSubArrWithVirtOscn(int, char *, int);
    void PrSwitch (char *, int);

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
    if (acs_info->blevcorr == PERFORM) {

       /* This block handles only the WFC detector */
       if (acs_info->detector == WFC_CCD_DETECTOR) {

          /* Only Pre-SM4 WFC data - this is the original bias level subtraction */
          if (acs_info->expstart < SM4MJD) {

             {int i;
             for (i = 0; i < acs_info->nimsets; i++) {

                 /* The variable done is set to the results of FindOver indicating there is
                    physical overscan. */
                 *done = overscan[i];

                 if (doBlev(&acs[i], &x[i], acs[i].chip, done, driftcorr)) {
                    return status;
                 }

                 /* Variable done is overwritten in doBlev if bias is determined from overscan region */
                 if (!*done) {
                     trlmessage("Default bias level from CCDTAB was subtracted.");
                 }

                 blevSubTrlMessage(acs_info, &acs[i]);

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
             *done = NO;

             /* Process full frame data */
             if (acs_info->subarray == NO) {
                if ((overscan[0] == YES) && (overscan[1] == YES))
                   *done = YES;

                if (*done) {
                   PrSwitch("blevcorr", PERFORM);

                   /* Only do bias-shift and cross talk corrections of images taken
                      with gain = 2 and in dual-slope integrator mode.  These corrections
                      are done on a per chip basis (one chip per imset).
                   */
                   if (strcmp(acs_info->jwrotype, "DS_INT") == 0 &&
                           acs_info->ccdgain == 2) {
                      trlmessage("Performing bias-shift and cross talk corrections for full frame data.");

                      if (bias_shift_corr(acs_info, acs_info->nimsets, &x[0], &x[1])) {
                         return status;
                      }

                      {int i;
                      for (i = 0; i < acs_info->nimsets; i++) {
                          cross_talk_corr(&acs[i], &x[i]);
                      }}
                   } else {
                      trlmessage("WFC readout type/gain not set as needed, no bias shift nor cross talk correction done for full frame data.");
                   }

                   trlmessage("Performing stripe removal and bias level subtraction for full frame data.");

                   if (doDestripe(acs_info, &x[0], &x[1])) {
                      return status;
                   }

                   {int i;
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

                   No need to check "done" here as we will decide based upon
                   tha array names of supported subarrays AND the determination they have virtual overscan
                   according to the image size (e.g., WFC1A-2K with size 2072x2068).
                   See comments in dobiasshift.c fo this function for more info.
                */
                if ((isValidBiasShiftSubArrWithVirtOscn(acs_info->subarray, acs_info->aperture, virtOverscan[0]) == YES) &&
                        (strcmp(acs_info->jwrotype, "DS_INT") == 0) &&
                        (acs_info->ccdgain == 2)) {
                   PrSwitch("blevcorr", PERFORM);

                   trlmessage("Performing bias-shift correction for subarray data.");

                   if (overscan[0] == YES)  /* Not used as BLEVCORR switch, just informational */
                      *done = YES;

                   /* Only bias shift correction is done for subarray data. For a supported
                      subarray, there is only one amp on a single chip in use.
                   */
                   if (bias_shift_corr(acs_info, 1, &x[0])) {
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

                   {int i;
                   for (i = 0; i < acs_info->nimsets; i++) {

                       /* This is set based on results of FindOver */
                       *done = overscan[i];

                       if (doBlev(&acs[i], &x[i], acs[i].chip, done, driftcorr)) {
                          return status;
                       }

                       /* Variable done is overwritten in doBlev if bias is determined from overscan region */
                       if (!*done) {
                          trlmessage("Default bias level from CCDTAB was subtracted.");
                       }

                       blevSubTrlMessage(acs_info, &acs[i]);

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
            {int i;
            for (i = 0; i < acs_info->nimsets; i++) {
                /* This is set based on results of FindOver */
                *done = overscan[i];

                if (doBlev(&acs[i], &x[i], acs[i].chip, done, driftcorr)) {
                   return status;
                }

                /* Variable done is overwritten in doBlev if bias is determined from overscan region */
                if (!*done) {
                   trlmessage("Default bias level from CCDTAB was subtracted.");
                }

                blevSubTrlMessage(acs_info, &acs[i]);

                /* Set this to complete so overscan-trimmed image will be written out. */
                blevcorr[i] = COMPLETE;
            }}

            PrSwitch("blevcorr", COMPLETE);

            if (acs_info->printtime)
               TimeStamp("BLEVCORR complete", acs_info->rootname);

       } /* End for HRC data */

    } /* End of BLEVCORR = PERFORM */

    return status;
}
