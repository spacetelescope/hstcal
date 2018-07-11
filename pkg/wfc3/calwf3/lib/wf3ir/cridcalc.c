# include <math.h>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>

#include "hstcal.h"
# include "hstio.h"    /* defines HST I/O functions */
# include "wf3.h"
# include "wf3info.h"
# include "wf3dq.h"
# include "trlbuf.h"
# include "wf3rej.h"
# include "rej.h"

# define max_CRs         4
# define equal_weight    0
# define optimum_weight  1

# define DEBUG 0
# define DEBUG2 0
# define X1     605-1
# define Y1     108-1
# define X2     398-1
# define Y2     15-1

static short DQIGNORE = SATPIXEL;

/** Function Instantiation **/

static int RejSpikes (float *, float *, short *, float *, short, float, int *);
static int RejFirstRead (short *, float *, short, float);
static int RejCRs (short *, float *, short, float, int *);
static void EstimateDarkandGlow(const short nsamp,float *time,float,float *tot_ADUs);

extern int status;

static int  crrej (WF3Info *, MultiNicmosGroup *, SingleNicmosGroup *);
static void fitsamps (const short, float *, float *, short *,float *, float *, float,
        float *, float *, short *, float *, short, short, float, float);
static void linfit (short, float *, float *, float *, float *, short, float, 
        float, float *, float *, float *, float *, int, int);


/** THE DRIVER ROUTINE **/

/* CRIDCALC: Identify cosmic ray hits in MultiAccum images.
 **
 ** Revision history:
 ** H.Bushouse  Oct. 2000   Initial CALNICA to CALWF3 port.  
 ** M.Sosey July 2008   Extensive update to use improvements in CALNICA
 **                              made in 2007 to the way CR's are detected.
 */

int cridcalc (WF3Info *wf3, MultiNicmosGroup *input, SingleNicmosGroup *crimage) {

    /* Arguments:
     ** wf3  i: WFC3 info structure
     ** input   io: input image stack
     ** crimage  o: CR reject output image
     **     flat    i: input flatfield image -- REMOVED
     */

    /* Function definitions */
    void PrSwitch (char *, int);

    if (wf3->crcorr == PERFORM) {

        if (crrej (wf3, input, crimage))
            return (status);

        PrSwitch ("crcorr", COMPLETE);
    }

    /* Successful return */
    return (status = 0);
}


/* CRREJ: Identify cosmic ray hits in a stack of MultiACCUM samples.
 ** Detected hits are flagged in the DQ arrays of the individual samples,
 ** but the SCI images in the individual samples are NOT changed.
 ** A final single-group image is produced using only unflagged data
 ** from the individual samples. The accumlated counts vs. time are fitted
 ** with a first-order (linear) polynomial using standard linear regression.
 ** Outliers are identified by their distance from the fit and rejected.
 **
 ** Revision history:
 ** H.Bushouse  10-Apr-2002 Modified to skip reference pixels by setting
 **             loop limits based on trim values from OSCNTAB.
 ** H.Bushouse  08-May-2002 Switched DQ flag macros over from calnica
 **             names to calwf3 names (SATURATED->SATPIXEL,
 **             CR_HIT->DATAREJECT, BADPIX->SPIKE).
 **             Removed code involving "GROT" flags (not used
 **             for wfc3). Added use of wf3dq.h.
 **
 ** M. Sosey    Updated to use improvements in CALNICA pipeline.
 ** July-2008   This includes adding the new CR detection and fitting routines 
 **     using the Fixsen et al standards. Adding the HIGH_CURVATURE
 **     DQ flag for pixels with higher than max cr hits (16384, which
 **     was previously reserved. fitsamps and  linfit are significantly
 **     changed, so was crrej, and RejFirstRead was added. 
 **
 ** M. Sosey    OPUS 2008.3 testing revealed a problem with 
 ** 25-July-2008    the estimation of the gain that was fed to linfit later on. If 
 **     FLATCORR was set to perform then the FLATFIELD file was 
 **     referenced directly using the Pix function to return indiv. 
 **     pixel values. However, if the reference flatfield image has a 
 **     null ERR extension it would cause a memory error. Also, WFC3 
 **     can use up to 3 flatfields, which need to be combined to get
 **     the true value. After examining the code and realizing that 
 **     this step was being performed for NICMOS to back out the 
 **     flatfield processing so that the errors made sense, we decided 
 **     to take out the loop all together and just return the value of 
 **     the gain in all cases. This is reasonable because WFC3 performs
 **     flatfielding AFTER the CRIDCALC stage, so the step is 
 **     unneccesary. It would also have mandated a bunch of changes to 
 **     the code.
 ** 
 ** M. Sosey    Corrected last few lines of code in linfit, which contained an
 ** 01-Aug-2008 if statement that wasn't supposed to be there. This was 
 **     discovered by Dahlen and Reagan in the NICMOS cridcalc code, 
 **     per ticket #233 in stsci_python.
 **
 ** H.Bushouse  Reinstated code that had been inadvertantly removed from the
 ** 10-Dec-2008 calnica code ported to calwf3, which propogates CR DQ flags to
 **     all samples following a hit.
 **
 ** H.Bushouse  Enhanced to get CR rejection threshold from CRREJTAB, rather
 ** 13-Mar-2009  than having it hardwired in the code. This involves a call to
 **     the new crrpar_in routine, which is located in refdata.c.
 **     Decreased max_CRs from 6 to 4, now that we have a good
 **     non-linearity correction for WFC3.
 **     Reinstated old code that adjusted loop limits so that ramp
 **     fitting is not applied to reference pixels. The old code had
 **     gotten removed when porting the calnica algorithms into calwf3.
 **     Fixed bug in IF syntax that checks for pixels saturated already
 **     in first read.
 **     Don't set HIGH_CURVATURE value in output DQ arrays; use
 **     UNSTABLE instead, and change all messages to say UNSTABLE.
 **     Don't set ZEROSIG value in output crimage (flt file) DQ array,
 **     because these pixels are still good (assuming they have no
 **     other flags also set).
 **     Removed unnecessary call to EstimateDarkandGlow at end of
 **     processing.
 **     Fixed many bugs in fitsamps that came with the calnica port:
 **     include time of zero-to-first read interval in output TIME;
 **     include zero read in ouput SAMP; don't double count samples
 **     at ends of multiple intervals in output SAMP; fixed bug in 
 **     syntax of IF statement that starts block of output computations
 **     for pixels that have only 1 good sample (was causing outputs to
**      to be zero, instead of using the 1 good sample); fixed bug in
**      "firstgood", "lastgood" assignments that was causing output TIME
**      to sometimes be negative.
**
** H.Bushouse   Modified to use badinpdq value from crrejtab to set DQ mask
** 13-Oct-2009  that is used to reject samples and zero-out output pixels,
    **      rather than just hardwiring the DQIGNORE constant to a value.
    **
    ** H.Bushouse   Modified the handling of output DQ values to match the new
    ** 26-Oct-2009  handling of sample rejection implemented in previous update,
    **      such that DQ values from any sample are carried through to
    **      the output if a pixel has no good samples, rather than the
    **      previous method of only carry those DQ bits that appear in
    **      ALL samples. And don't reset output DQ to SATPIXEL if all
    **      samples are rejected.
    **
    ** H.Bushouse   Added check for pixels already saturated in zeroth read (set
            ** Feb-Mar 2010 all outputs to zero). Some general cleanup and removal of
    **      unused variables. Switch from using commanded ccdgain to
    **      actual mean_gain. Modified linfit to include readnoise in
    **      sample weights and Poisson noise from source in final fit
    **      uncertainty. Added SPIKE_THRESH in RejSpikes to use as
    **      different threshold than CR thresh. Updated hardwired dark
**      and readnoise values to use SMOV results. (calwf3 2.0)
    **
    ** H.Bushouse   Changed crrej to always call EstimateDarkandGlow, regardless of
    ** 21-Oct-2010  darkcorr switch setting, because for WFC3 we use a static dark
    **      value and therefore don't need access to the darkfile.
    **
    ** H.Bushouse   Updated crrej to free memory for tot_ADUs before returning.
    ** 01-Aug-2011  (PR 68993; Trac #748)
    **
    ** H.Bushouse   Fixed fitsamps routine to correctly accumulate int_time when
    ** 31-Aug-2011  first valid interval contains a single sample or starts in
    **      something other than the first read. (PR 69230; Trac #770)
    **
    ** M.Sosey      Added code to detect negative cr hits, and set the detected SPIKES to 1024
    ** May 2012     during detection so that the 4 value isn't overloaded and ignored during the fit
    **              its saved to the output as 4 though. Also redid the logic for handling
    **              saturation in zeroth and first reads; they are now treated the same so that
    **              output pixels are never zeroed out. There's probably a lot of diffs from the previous
    **              version because I re-indented the entire file to help figure out the logic
    **                
    */

    int crrej (WF3Info *wf3, MultiNicmosGroup *input, SingleNicmosGroup *crimage) {

        /* Arguments:
         ** wf3  i: WFC3 info structure
         ** input   io: input image stack
         ** crimage  o: CR reject output image (single-group)
         **      flat     i: input flatfield image - removed
         */

        /* Local variables */
        short i, j, k, l;       /* pixel and loop indexes */
        short ibeg, iend, jbeg, jend;   /* loop limits */
        short nsamp;            /* number of samples for pixel */
        float current_sci;      /* current sci value */
        float current_err;      /* current err value */
        short current_dq;       /* current dq  value */
        float *sci;         /* list of sci  values for pixel */
        float *err;         /* list of err  values for pixel */
        short *dq;          /* list of dq   values for pixel */
        float *time;            /* list of time values for pixel */
        float *tot_ADUs;        /* list of total ADUs from dark current and 
                           amp glow for pixel */
        float out_sci;      /* output sci  value */
        float out_err;      /* output err  value */
        short out_dq;           /* output dq   value */
        short out_samp;     /* output samp value */
        float out_time;     /* output time value */
        float flat_value;             /* value to convert from flat fielded ADUs to 
                         electrons */
        float flat_uncertainty;       /* unitless rms flat field uncertainty */
        int ncurved;            /* Number of pixels with high curvature */
        int   niter = 0;        /* number of rejection iterations */
        float sigma[MAX_ITER];  /* list of sigma values for rejection */
        clpar par;          /* parameters used */
        int   newpar[MAX_PAR+1];    /* user specifiable parameters */

        /* Function definitions */
        int copyGroup (SingleNicmosGroup *, SingleNicmosGroup *);
        int crrpar_in (clpar *, int [], int, float, int *, float []);
        void rej_reset (clpar *, int []);

        /* Allocate memory for local arrays */
        sci    = NULL;
        err    = NULL;
        dq     = NULL;
        time   = NULL;
        sci    = (float *) calloc(wf3->ngroups, sizeof(float));
        err    = (float *) calloc(wf3->ngroups, sizeof(float));
        dq     = (short *) calloc(wf3->ngroups, sizeof(short));
        time   = (float *) calloc(wf3->ngroups, sizeof(float));
        tot_ADUs   = (float *) calloc(wf3->ngroups, sizeof(float));

        /* Read in the parameters from the crrejtab */
        rej_reset (&par, newpar);
        strcpy (par.tbname, wf3->crrej.name);
        strcpy (par.expname, "EXPTIME");
        par.verbose = wf3->verbose;
        if (crrpar_in (&par, newpar, wf3->ngroups, (float)wf3->exptime[0], &niter,
                    sigma))
            return (status);
        DQIGNORE = par.badinpdq;

        /* Initialize the output crimage by copying the first input
         ** group to the crimage */
        if (copyGroup (crimage, &(input->group[0])))
            return (status);

        /* Set CR reject threshold; use default if no user input */
        if (wf3->crthresh == 0)
            wf3->crthresh = sigma[0];

        /* Print info to processing log */
        if (wf3->samp_rej > 0) {
            sprintf (MsgText, "CRIDCALC using %g sigma rejection threshold",
                    wf3->crthresh);
            trlmessage (MsgText);
            sprintf (MsgText, "               %d bad DQ mask", DQIGNORE);
            trlmessage (MsgText);
            sprintf (MsgText, "               %d max CRs for UNSTABLE",max_CRs);
            trlmessage (MsgText);
            sprintf (MsgText,
                    "         and rejecting first %d samples for all pixels",
                    wf3->samp_rej);
            trlmessage (MsgText);
        } else {
            sprintf (MsgText, "CRIDCALC using %g sigma rejection threshold",
                    wf3->crthresh);
            trlmessage (MsgText);
            sprintf (MsgText, "               %d bad DQ mask", DQIGNORE);
            trlmessage (MsgText);
            sprintf (MsgText, "               %d max CRs for UNSTABLE",max_CRs);
            trlmessage (MsgText);
        }
        ncurved = 0;

        /* Loop over image array, computing mean countrate at each pixel;
         ** the loop limits are set so that reference pixels are skipped;
         ** this also handles subarrays, because limits are based on size
         ** of input science image. */
        ibeg = wf3->trimx[0];
        iend = input->group[0].sci.data.nx - wf3->trimx[1];
        jbeg = wf3->trimy[0];
        jend = input->group[0].sci.data.ny - wf3->trimy[1];

        for (j=jbeg; j<jend; j++) {
            for (i=ibeg; i<iend; i++) {

                /* Get the DQ value in the zeroth-read */
                k = wf3->ngroups-1;
                current_dq = DQPix(input->group[k].dq.data,i,j);

                /* Back out the ZEROSIG  dq flags */
                if (current_dq & ZEROSIG)
                    current_dq -= ZEROSIG;

                /* Check for pixels that are saturated already in first read */
                if ((wf3->zsigcorr == PERFORM || wf3->zsigcorr == COMPLETE) &&
                    (DQPix(input->group[k-1].dq.data,i,j) & SATPIXEL)) {

                    /* For these pixels, just set the output values equal
                     ** to what's in the input zeroth-read image, regardless
                     ** of whether the zeroth-read is saturated or not. If it
                     ** is saturated in the zeroth-read, the DQ flag will get
                     ** carried over to the output to indicate it's bad. */
                    Pix(crimage->sci.data,i,j) = Pix(input->group[k].sci.data,i,j);
                    Pix(crimage->err.data,i,j) = Pix(input->group[k].err.data,i,j);
                    DQSetPix(crimage->dq.data,i,j,current_dq);
                    Pix(crimage->smpl.data,i,j) = 1;
                    Pix(crimage->intg.data,i,j) = wf3->sampzero;

                    continue;
                }

                /* Initialize the output image DQ value */
                out_dq = 0;

                /* Create list of samples for this pixel */
                nsamp = 0;
                for (k = wf3->ngroups-1; k >= 0; k--) {

                    /* Retrieve SCI, ERR, and DQ values for this sample */
                    current_sci = Pix(input->group[k].sci.data,i,j);
                    current_err = Pix(input->group[k].err.data,i,j);
                    current_dq  = DQPix(input->group[k].dq.data,i,j);

                    /* Temporarily convert countrates back to counts */
                    if (wf3->bunit[0] == COUNTRATE) {
                        current_sci *= Pix(input->group[k].intg.data,i,j);
                        current_err *= Pix(input->group[k].intg.data,i,j);
                    }

                    /* Remove ZEROSIG bits from DQ, if present,
                     ** because those samples are OK to use here (Vsn 3.2) */
                    if (current_dq & ZEROSIG)
                        current_dq -= ZEROSIG;

                    /* Propagate DQ values to output DQ */
                    out_dq = out_dq | current_dq;

                    /* Temporarily flag first samp_rej samples to exclude
                     ** them from the fit (Version 3.3) */
                    if (wf3->samp_rej > 0 && nsamp > 0 && nsamp <= wf3->samp_rej)
                        current_dq = current_dq | RESERVED2;

                    /* Copy all values to the fitting arrays */
                    sci[nsamp]  = current_sci;
                    err[nsamp]  = current_err;
                    dq[nsamp]   = current_dq;
                    time[nsamp] = Pix(input->group[k].intg.data,i,j);
                    nsamp++;
                }

                /* Get dark and amp glow values for each sample of this pixel */
                /*if (wf3->darkcorr == PERFORM)*/
                EstimateDarkandGlow (nsamp, time, wf3->mean_gain, tot_ADUs);

                if (wf3->flatcorr == PERFORM) {
                    flat_value = wf3->mean_gain;
                    flat_uncertainty = 0.0;
                } else {
                    flat_value = wf3->mean_gain;
                    flat_uncertainty = 0.0;
                }

                /* Do iterative rejection and computation of slope */
                fitsamps (nsamp, sci, err, dq, time, tot_ADUs, wf3->crthresh,
                        &out_sci, &out_err, &out_samp, &out_time, i, j, flat_value,
                        flat_uncertainty);

                /* Propagate all DQ flags to output EXCEPT for SATPIXEL
                 ** and DATAREJECT  (Version 4.2) */
                out_dq = out_dq - (out_dq & SATPIXEL);

                /* Convert results back to counts if necessary */
                if (wf3->bunit[0] == COUNTS) {
                    out_sci *= out_time;
                    out_err *= out_time;
                }

                /* Set DATAREJECT DQ for all samples following a hit
                   This is done so that people looking at the imas in the
                   future know that the absolute value of the pixel is wrong
                   after the first hit, but it smears the location of any hits
                   which occurred in addition to the first one. */
                for (k = 0; k < wf3->ngroups-1; k++) {
                    if (dq[k] & DATAREJECT) {
                        dq[k+1] = dq[k+1] | DATAREJECT;
                    }
                }  

                /* If the HIGH-CURVATURE bit is set anywhere, set it for all
                   groups and unset the DATAREJECT bit */
                for (k = 0; k < wf3->ngroups-1; k++) {
                    if (dq[k] & HIGH_CURVATURE) {
                        for (l=0; l < wf3->ngroups; l++) {
                            dq[l] = dq[l] | HIGH_CURVATURE;
                            dq[l] = dq[l] - (dq[l] & DATAREJECT);
                        }
                        break;
                    }
                }

                /* Add UNSTABLE bit to output crimage if the pixel is marked as
                 ** HIGH_CURVATURE and unset the DATAREJECT bit. */
                if (dq[nsamp-1] & HIGH_CURVATURE) {
                    out_dq = out_dq | UNSTABLE;
                    out_dq = out_dq - (out_dq & DATAREJECT);
                    /* And increment the ncurved counter */
                    ncurved++;
                }

                /* Store final values in output crimage */
                Pix(crimage->sci.data,i,j)  = out_sci;
                Pix(crimage->err.data,i,j)  = out_err;
                DQSetPix(crimage->dq.data,i,j,out_dq);
                Pix(crimage->smpl.data,i,j) = out_samp;
                Pix(crimage->intg.data,i,j) = out_time;

                /* Update input DQ values for detected outliers */
                for (k = wf3->ngroups-1; k >= 0; k--) {

                    if (dq[wf3->ngroups-1-k] & DATAREJECT)
                        DQSetPix (input->group[k].dq.data,i,j,
                                DQPix (input->group[k].dq.data,i,j) | DATAREJECT);

                    if (dq[wf3->ngroups-1-k] & SPIKE)
                        DQSetPix (input->group[k].dq.data,i,j,
                                DQPix (input->group[k].dq.data,i,j) | DETECTORPROB);

                    if (dq[wf3->ngroups-1-k] & HIGH_CURVATURE)
                        DQSetPix (input->group[k].dq.data,i,j,
                                DQPix (input->group[k].dq.data,i,j) | UNSTABLE);
                }

            } /* end of loop over nx */
        } /* end of loop over ny */

        if (ncurved > 0) {
            sprintf (MsgText, "%d pixels detected as unstable", ncurved);
            trlmessage (MsgText);
        }

        /* Free memory allocated locally */
        free(sci);
        free(err);
        free(dq);
        free(time);
        free(tot_ADUs);

        /* Successful return */
        return (status = 0);
    }


/* FITSAMPS: Fit accumulating counts vs. time to compute mean countrate,
 ** iteratively rejecting CR hits and refitting until no new samples are
 ** rejected. The rejection test is based on the input error value for
 ** each sample.
 */

static void fitsamps (const short nsamp, float *sci, float *err, short *dq,
        float *time, float *darkandglow, float thresh, 
        float *out_sci, float *out_err, 
        short *out_samp, float *out_time,
        short i, short j, float gain, float flat_uncertainty) {
    /* Local variables */
    int k,idx,l;            /* loop index */
    int nrej, rej;      /* number & index of rejected samples */
    float (*tsci)[nsamp];       /* list of cleaned fluxes */
    float (*terr)[nsamp];       /* list of cleaned errors */
    float (*ttime)[nsamp];  /* list of cleaned times */
    short (*tdq)[nsamp];        /* list of cleaned dqa */
    float (*tdarkandglow)[nsamp]; /* poisson counts for this interval */
    int   (*tlookup)[nsamp];      /* poisson counts for this interval */
    short   interval;             /* number of intervals that have to be fit */
    int   sample;       
    int   number_CRs;             /* total number of CRs found in the set of reads */
    float a[max_CRs+1], siga[max_CRs+1];    /* linear fit zeropoint and error */
    float b[max_CRs+1], sigb[max_CRs+1];    /* linear fit slope and error */
    short tcount[max_CRs+1];        /* length in each interval */
    float *diff;                /* sample differences relative to fit */
    float sumwts,sum;             /* used to combine the slopes of the CR
                     created intervals */
    float int_time;               /* the sum of all the times without CR hits */
    float sigsquared;
    short iter_count;             /* for debugging*/
    int found_error;              /* Flag for too many cosmic rays */
    int no_samples;     /* Flag for no valid samples */
    int firstgood;                /* First pixel with OK DQ */
    int lastgood;                 /* Last pixel with OK DQ */
    int ngood;
    int nsflag;                   /* Flag set when we go through no_samples=1 */
    int naccum;

    /* Initialize the output values */
    *out_sci  = 0;
    *out_err  = 0;
    *out_samp = 0;
    *out_time = 0;

    iter_count = 0;
    found_error = 0;
    no_samples = 0;
    nsflag = 0;

    /* Allocate memory for local arrays */
    tsci   = calloc(nsamp*(max_CRs+1), sizeof(float));
    terr   = calloc(nsamp*(max_CRs+1), sizeof(float));
    tdq   =  calloc(nsamp*(max_CRs+1), sizeof(short));
    ttime  = calloc(nsamp*(max_CRs+1), sizeof(float));
    tdarkandglow   = calloc(nsamp*(max_CRs+1), sizeof(float));
    tlookup   = calloc(nsamp*(max_CRs+1), sizeof(int));
    diff   = (float *) calloc(nsamp, sizeof(float));


    /* Fit and test samples for this pixel until no new samples
     ** are rejected:

     start_index and end_index are paired variables defining the start
     and end of intervals along the samples. If no CRs were found, then 
     start_index[0] would be 1 (skip the 0th read) and
     end_index[0] would be nsamp-1.

     Finding a CR ends one interval and begins the next. The CR must
     be included in the next interval.

     This is an iterative process:
     - first, intervals are defined based on existing CRs
     - then, each interval is fitted separately
     - then, each interval is inspected for spikes
     - then, each interval is inspected for CRs

     If any spikes or CRs are found, the entire procedure repeats: new
     intervals are defined, etc.

     After the iteration ends because no new spikes or CRs are found,
     each interval is fitted separately once again, with optimum weighting,
     and then the results for each interval are combined to obtain the
     final solution for the pixel.

     */
    number_CRs=0;
    nrej = 1; 
    while (nrej > 0) {
        nrej = 0;  /* RejSpikes and RejCRs may increment nrej */
        iter_count++;

        /*---------------------------------------------------------------
          IDENTIFY THE INTERVALS
          ----------------------------------------------------------------*/
        interval = 0;
        sample = 0;
        for (k=1; k <nsamp; k++) {
            if (dq[k] == 0) { /*accumulate the sample*/
                tsci[interval][sample]     =     sci[k];
                ttime[interval][sample]    =     time[k];
                terr[interval][sample]     =     err[k];
                tdq[interval][sample]      =     dq[k];
                tdarkandglow[interval][sample] = darkandglow[k];
                tlookup[interval][sample]  =     k;
                sample++;
            } else if (dq[k] == DATAREJECT) { /* When a CR is found, */
                tcount[interval]=sample; /*save the num. of points in this interval*/
                interval++;               /*start a new interval*/
                sample=0;                 /*resetting its sample count*/
                tsci[interval][sample]     =     sci[k];
                ttime[interval][sample]    =     time[k];
                terr[interval][sample]     =     err[k];
                tdq[interval][sample]      =     dq[k];
                tdarkandglow[interval][sample] = darkandglow[k];
                tlookup[interval][sample]  =     k;
                sample++;
            }
        }

        /*If the last sample we kept has a CR, then the last interval
          has nothing in it but one CR: it can't be used, so ignore it*/
        if (sample >= 1) {
            if (tdq[interval][sample-1] == DATAREJECT) {
                interval--;
            } else {  /*Otherwise close the last interval*/
                tcount[interval]=sample;
            }
        } else {
            /*There are no valid samples in this pixel. Bail out.*/
            no_samples = 1;
            nrej = 0;
            break;
        }

        /******  Print statements that can be turned on for debugging *****/
        if (DEBUG) {
            /* Set the pixel you want to examine up top */

            if ((i==X1 && j==Y1) || (i==X2 && j==Y2)) {
                printf("========== Debug info: iteration %d ====================\n",
                        iter_count);
                printf("Pixel i,j=%d,%d:    nsamp %d\n",i,j,nsamp);
                printf("Sci: ");
                for (k=0; k<nsamp; k++){
                    printf("% 10.6f",sci[k]);
                }
                printf("\nDQ: ");
                for (k=0; k<nsamp; k++){
                    printf(" %5d",dq[k]);
                }
                printf("\nTIME: ");
                for (k=0; k<nsamp; k++){
                    printf(" %10.6f",time[k]);
                }

                printf("\n\n Intervals: %d  \n", interval+1);
                for (idx=0; idx<=interval; idx++){
                    printf("Interval %d, tcount %d\n sci: ",idx+1,tcount[idx]);
                    for (k=0; k<tcount[idx]; k++){
                        printf(" %10.6f",tsci[idx][k]);
                    }
                    printf("\n  dq:");
                    for (k=0; k<tcount[idx]; k++){
                        printf(" %5d",tdq[idx][k]);
                    }
                    printf("\ntime:");
                    for (k=0; k<tcount[idx]; k++){
                        printf(" %10.6f",ttime[idx][k]);
                    }
                    printf("\n...........................\n");
                }  /* End for loop over idx  */
            }    /* End if DEBUG pixel */
        }      /* End if DEBUG */

        /*****************end debugging***********************/


        /*---------------------------------------------------------------
          FIT EACH (non-zero-length) INTERVAL SEPARATELY
          ----------------------------------------------------------------*/

        /* We need a unique solution for each interval. */
        for (idx=0;(idx<=interval)&&nrej==0;idx++){
            if (tcount[idx] > 1){
                /* Compute mean countrate using linear fit and equal weighting 
                   to best find SPIKES and CRs
                   
                   
                   using a straight difference between samples to detect hits before
                   calculating the slope might be more mathmatically sound, see
                   Karls paper
                   */
                linfit (equal_weight, ttime[idx], tsci[idx], terr[idx], 
                        tdarkandglow[idx], tcount[idx], gain, flat_uncertainty, 
                        &a[idx], &b[idx], &siga[idx], &sigb[idx],i,j);

                /*---------------------------------------------------------------
                  EXAMINE EACH INTERVAL FOR SPIKES AND CRS
                  ----------------------------------------------------------------*/

                /* Compute difference of each sample from the local fit 
                   Note that the fit is done with the CR included, this might
                   throw off the actual detection. 
                */
                for (k = 0; k <tcount[idx]; k++) {
                    if (terr[idx][k] != 0 && tdq[idx][k] == 0){
                        diff[k] =(tsci[idx][k]-(a[idx]+b[idx]*ttime[idx][k]))
                            / terr[idx][k];
                        if (DEBUG && ((i==X1 && j==Y1) || (i==X2 && j==Y2))){
                            printf (" samp %d, tsci=%g, terr=%g, a=%g, b=%g, fit=%g, diff=%g\n",
                                    k,tsci[idx][k],terr[idx][k],a[idx],b[idx],
                                    (a[idx]+b[idx]*ttime[idx][k]),diff[k]);
                        }
                    } else { /* ignore diffs for already-identified CRs & spikes */
                        diff[k] = 0;
                    }
                 }

                /* Look for and flag new electronic noise spikes and CRs*/
                nrej=RejSpikes (tsci[idx], terr[idx], tdq[idx], diff, tcount[idx],
                        thresh, &rej);
                if (nrej > 0) {
                    dq[tlookup[idx][rej]] = dq[tlookup[idx][rej]] | SPIKE;
                } else if (nrej == 0) {
                    /* Only look for new CR hits if no spike was found */
                    nrej = RejCRs (tdq[idx], diff, tcount[idx], thresh, &rej);
                    if (nrej > 0) {
                        dq[tlookup[idx][rej]] = dq[tlookup[idx][rej]] | DATAREJECT;
                        number_CRs +=1;
                        if (number_CRs >= max_CRs) {
                            if (DEBUG && ((i==X1 && j==Y1) || (i==X2 && j==Y2)))
                                printf (" hit max_CRs limit\n");
                            found_error = 1;
                            dq[tlookup[idx][rej]] = dq[tlookup[idx][rej]] - DATAREJECT;
                            dq[tlookup[idx][rej]] = dq[tlookup[idx][rej]] | HIGH_CURVATURE;
                            nrej = 0; /*To break out of the while loop*/
                            break;
                        }
                    } else if (nrej == 0) {
                        /* Only check the first read if nothing else has been found */
                        nrej = RejFirstRead(tdq[idx], diff, tcount[idx], thresh);
                        if (nrej > 0){
                            dq[tlookup[idx][0]] = dq[tlookup[idx][0]] | SPIKE;
                        }
                    }
                }
            } else if (tcount[idx] == 1) {
                if (tdq[idx][0] != DATAREJECT) { /*ignore CR hits if they're alone*/
                    b[idx]=tsci[idx][0]/ttime[idx][0];
                    sigb[idx]=terr[idx][0]/ttime[idx][0];
                }
            }   /*endif tcount for this interval*/
        }     /*for each interval*/
    }       /*while new rejections occur*/

    /*--------------------------------------------------------------------
      If there are no good pixels, repeat this process for all
      pixels that don't have bits in the DQIGNORE pattern set
      --------------------------------------------------------------------*/

    if (no_samples == 1) {
        no_samples = 0;      /* Reset because we're starting again   */
        nsflag = 1;
        if (DEBUG && ((i==X1 && j==Y1) || (i==X2 && j==Y2))) 
            printf ("no_samples = 1 at i = %d, j = %d\n", i, j);
        number_CRs=0;
        nrej = 1;
        iter_count = 0;
        while (nrej > 0) {
            nrej = 0;  /* RejSpikes and RejCRs may increment nrej */
            iter_count++;

            /*---------------------------------------------------------------
              IDENTIFY THE INTERVALS
              ----------------------------------------------------------------*/
            interval = 0;
            sample = 0;
            for (k=1; k <nsamp; k++) {
                if ((dq[k] & (DQIGNORE | DATAREJECT)) == 0) { /*accumulate the sample*/
                    tsci[interval][sample]     =     sci[k];
                    ttime[interval][sample]    =     time[k];
                    terr[interval][sample]     =     err[k];
                    tdq[interval][sample]      =     dq[k];
                    tdarkandglow[interval][sample] = darkandglow[k];
                    tlookup[interval][sample]  =     k;
                    sample++;
                } else if ((dq[k] & DATAREJECT) != 0) { /* When a CR is found, */
                    tcount[interval]=sample; /*save the num. of points in this interval*/
                    interval++;               /*start a new interval*/
                    sample=0;                 /*resetting its sample count*/
                    tsci[interval][sample]     =     sci[k];
                    ttime[interval][sample]    =     time[k];
                    terr[interval][sample]     =     err[k];
                    tdq[interval][sample]      =     dq[k];
                    tdarkandglow[interval][sample] = darkandglow[k];
                    tlookup[interval][sample]  =     k;
                    sample++;
                }
            }

            /*If the last sample we kept has a CR, then the last interval
              has nothing in it but one CR: it can't be used, so ignore it*/
            if (sample >= 1) {
                if (tdq[interval][sample-1] == DATAREJECT) {
                    interval--;
                } else {  /*Otherwise close the last interval*/
                    tcount[interval]=sample;
                }
            } else {
                /*There are no valid samples in this pixel. Bail out.*/
                no_samples = 1;
                nrej = 0;
                if (DEBUG2) printf ("Still no valid samples\n");
                break;
            }
            if (DEBUG2) {
                printf ("Iter count = %d. %d intervals\n", iter_count, interval+1);
                for (k=0; k <= interval; k++) {
                    for (l=0; l < tcount[k]; l++) {
                        printf (" %d", tdq[k][l]);
                    }
                    printf ("\n");
                    for (l=0; l < tcount[k]; l++) {
                        printf (" %7.2f", tsci[k][l]);
                    }
                    printf ("\n");
                }
            }
            /******  Print statements that can be turned on for debugging *****/
            if (DEBUG) {
                /* Set the pixel you want to examine up top */

                if ((i==X1 && j==Y1) || (i==X2 && j==Y2)){
                    printf("========== Debug info: iteration %d ====================\n",
                            iter_count);
                    printf("Pixel i,j=%d,%d:    nsamp %d\n",i,j,nsamp);
                    printf("Sci: ");
                    for (k=0; k<nsamp; k++){
                        printf(" %10.6f",sci[k]);
                    }
                    printf("\nDQ: ");
                    for (k=0; k<nsamp; k++){
                        printf(" %5d",dq[k]);
                    }
                    printf("\nTIME: ");
                    for (k=0; k<nsamp; k++){
                        printf(" %10.6f",time[k]);
                    }

                    printf("\n\n Intervals: %d  \n", interval+1);
                    for (idx=0; idx<=interval; idx++){
                        printf("Interval %d, tcount %d\n sci: ",idx+1,tcount[idx]);
                        for (k=0; k<tcount[idx]; k++){
                            printf(" %10.6f",tsci[idx][k]);
                        }
                        printf("\n  dq:");
                        for (k=0; k<tcount[idx]; k++){
                            printf(" %5d",tdq[idx][k]);
                        }
                        printf("\ntime:");
                        for (k=0; k<tcount[idx]; k++){
                            printf(" %10.6f",ttime[idx][k]);
                        }
                        printf("\n...........................\n");
                    }  /* End for loop over idx  */
                }    /* End if DEBUG pixel */
            }      /* End if DEBUG */

            /*****************end debugging***********************/


            /*---------------------------------------------------------------
              FIT EACH (non-zero-length) INTERVAL SEPARATELY
              ----------------------------------------------------------------*/

            /* We need a unique solution for each interval. */
            for (idx=0;(idx<=interval)&&nrej==0;idx++){
                if (tcount[idx] > 1){
                    /* Compute mean countrate using linear fit and equal weighting 
                       to best find spikes and CRs*/
                    linfit (equal_weight, ttime[idx], tsci[idx], terr[idx], 
                            tdarkandglow[idx], tcount[idx], gain, flat_uncertainty, 
                            &a[idx], &b[idx], &siga[idx], &sigb[idx],i,j);
                    if (DEBUG2) {
                        printf ("Slope = %8.2f, intercept = %8.2f ", b[idx], a[idx]);
                        printf ("(equal weight)\n");
                    }

                    /*---------------------------------------------------------------
                      EXAMINE EACH INTERVAL FOR SPIKES AND CRS
                      ----------------------------------------------------------------*/

                    /* Compute difference of each sample from the local fit */
                    naccum = 0;
                    for (k = 0; k <tcount[idx]; k++) {
                        if ((terr[idx][k] != 0) && ((tdq[idx][k] & DQIGNORE) == 0)){
                            diff[k] =(tsci[idx][k]-(a[idx]+b[idx]*ttime[idx][k]))
                                / terr[idx][k];
                            naccum +=1;
                            if (DEBUG && ((i==X1 && j==Y1) || (i==X2 && j==Y2)))
                                printf (" diff[%d]=%g\n",k,diff[k]);
                        }
                        else{ /* ignore diffs for already-identified CRs & spikes*/
                            diff[k] = 0;
                        }
                    }
                    if (DEBUG2) printf ("Accumulated %d diffs\n", naccum);

                    /* Look for and flag new electronic noise spikes and CRs*/
                    nrej=RejSpikes (tsci[idx], terr[idx], tdq[idx], diff, tcount[idx],
                            thresh, &rej);
                    if (nrej > 0) {
                        dq[tlookup[idx][rej]] = dq[tlookup[idx][rej]] | SPIKE;
                        if (DEBUG2) printf ("Rejected %d points\n", nrej);
                    } else if (nrej == 0) {
                        /* Only look for new CR hits if no spike was found */
                        nrej = RejCRs (tdq[idx], diff, tcount[idx], thresh, &rej);
                        if (nrej > 0) {
                            dq[tlookup[idx][rej]] = dq[tlookup[idx][rej]] | DATAREJECT;
                            if (DEBUG2) printf ("Rejected %d cosmic rays\n", nrej);
                            number_CRs +=1;
                            if (number_CRs >= max_CRs) {
                                if (DEBUG && ((i==X1 && j==Y1) || (i==X2 && j==Y2)))
                                    printf (" hit max_CRs limit\n");
                                found_error = 1;
                                dq[tlookup[idx][rej]] = dq[tlookup[idx][rej]] - DATAREJECT;
                                dq[tlookup[idx][rej]] = dq[tlookup[idx][rej]] | HIGH_CURVATURE;
                                nrej = 0; /*To break out of the while loop*/
                                break;
                            }
                        } else if (nrej == 0) {
                            /* Only check the first read if nothing else has been found */
                            nrej = RejFirstRead(tdq[idx], diff, tcount[idx], thresh);
                            if (nrej > 0){
                                dq[tlookup[idx][0]] = dq[tlookup[idx][0]] | SPIKE;
                                if (DEBUG2) printf ("Rejected first read\n");
                            }
                        }
                    }
                } else if (tcount[idx] == 1) {
                    if ((tdq[idx][0] & DATAREJECT) != 0) { /*ignore CR hits if alone*/
                        if (ttime[idx][0] == 0) {
                            if (DEBUG2) printf ("Exposure time is zero\n");
                        } else {
                            b[idx]=tsci[idx][0]/ttime[idx][0];
                            sigb[idx]=terr[idx][0]/ttime[idx][0];
                        }
                    }
                }   /*endif tcount for this interval*/
            }     /*for each interval*/
        }       /*while new rejections occur*/
    }         /*if no_samples == 1 */

    /*--------------------------------------------------------------------
      ITERATION COMPLETE: PERFORM FINAL FIT OF EACH INTERVAL
      ---------------------------------------------------------------------*/

    /* Now go back through the intervals and find the best slopes using 
       the optimum weighting */
    int_time = 0;
    sumwts = 0;
    sum = 0;
    
    /* has more that one sample and doesn't contain too many cosmic ray hits */
    if (found_error == 0 && no_samples == 0) { 

        /* Loop over all intervals */
        for (idx=0; idx<=interval; idx++) {

            /* If there's more than 1 data point in the interval, do a slope fit */
            if (tcount[idx] > 1) {

                /* Initialize integration time counter with delta time
                   of first valid sample in first interval */
                if (idx==0) {
                    k = tlookup[idx][0];
                    int_time += time[k] - time[k-1];
                }

                /* Compute mean countrate using linear fit, with optimum weighting 
                   this time to get best estimate of the slope */
                linfit (optimum_weight, ttime[idx], tsci[idx], terr[idx], 
                        tdarkandglow[idx],tcount[idx], gain, flat_uncertainty, 
                        &a[idx], &b[idx], &siga[idx], &sigb[idx],i,j);
                if (nsflag == 1) {
                    if (DEBUG2) {
                        printf ("Slope = %8.2f, intercept = %8.2f", b[idx], a[idx]);
                        printf (" (optimum weight)\n");
                    }
                }

                /* also accumulate the integration time and sums  */
                int_time += ttime[idx][tcount[idx]-1] - ttime[idx][0];
                sigsquared = sigb[idx]*sigb[idx];
                if (sigsquared == 0.0) {
                    if (DEBUG2) printf ("sigb = 0\n");
                } else {
                    sumwts += 1.0/sigsquared;
                    sum    += b[idx]/sigsquared;
                }

                /* If there's only 1 data point, compute countrate directly */
            } else if (tcount[idx] == 1) {

                if ((tdq[idx][0] & DATAREJECT) == 0) {
                    if (ttime[idx][0] == 0.0) {
                        if (DEBUG2) printf ("Exposure time is zero\n");
                    } else {
                        b[idx]=tsci[idx][0]/ttime[idx][0];
                        sigb[idx]=terr[idx][0]/ttime[idx][0];
                    }

                    /* accumulate the integration time and sums */
                    k = tlookup[idx][0];
                    int_time += time[k] - time[k-1];
                    sigsquared = sigb[idx]*sigb[idx];
                    if (sigsquared == 0.0) {
                        if (DEBUG2) printf ("sigb[idx] = 0 and tcount[idx] = 1\n");
                    } else {
                        sumwts += 1/sigsquared;
                        sum    += b[idx]/sigsquared;
                    }
                }
            }
        }  /* end of loop over intervals */
    }    /* end test against found_error */

    /*----------------------------------------------------------------
      COMBINE ALL THE FITS TO GET THE FINAL SOLUTION FOR THE PIXEL
      ---------------------------------------------------------------*/

    /* Set output SCI and ERR values */
    if (no_samples == 0) {

        /* Good sums to work with ... */
        if ((found_error == 0) && (sumwts > 0.0)) {
            *(out_sci) = sum/sumwts;
            *(out_err) = sqrt(1/sumwts);
            *(out_time) = int_time;
            if (nsflag == 1) {
                if (DEBUG2) printf ("Output = %8.2f\n", *(out_sci));
            }
            /* Compute output SAMP  */
            for (idx = 0; idx <= interval; idx++) {
                if (idx==0) {
                    /* Add 1 to the number of samps in the first interval, because
                       it doesn't take into account the zeroth read sample. */
                    (*out_samp) += tcount[idx] + 1;
                } else {
                    /* Subtract 1 from the number of samps in all remaining intervals
                       to avoid double counting the end points that are shared. */
                    (*out_samp) += tcount[idx] - 1;
                }
            }

            /* No good sums to work with ... */
        } else {
            /*Set the flux and errors to the last good read*/
            firstgood = 0;
            lastgood = 0;
            ngood = 0;
            for (k=1; k < nsamp; k++) {
                if ((dq[k] & DQIGNORE) == 0) {
                    firstgood = k;
                    lastgood = k;
                    ngood++;
                    break;
                }
            }
            for (k=firstgood; k < nsamp; k++) {
                if (dq[k] == 0) {
                    lastgood = k;
                    ngood++;
                }
            }
            if (lastgood != firstgood) {
                *(out_sci) = (sci[lastgood] - sci[firstgood]) /
                    (time[lastgood] - time[firstgood]);
                *(out_err) = (err[lastgood] - err[firstgood]) /
                    (time[lastgood] - time[firstgood]);
                *(out_time) = time[lastgood] - time[firstgood];
                *(out_samp) = ngood;
            } else {
                *(out_sci) = 0;
                *(out_err) = 0;
                *(out_time) = 0;
                *(out_samp) = 0;
            }
        }
    } else {
        *(out_err) = 1e5;
    }

    if (DEBUG && ((i==X1 && j==Y1) || (i==X2 && j==Y2))){
        printf("Output values SCI/ERR/SAMP/TIME %f / %f / %d / %f \n",
                *out_sci,*out_err,*out_samp,*out_time);
        printf("========== End debug info ====================\n");
    }
    /*--------------------------------------------------------------
      Clean up and go home
      ----------------------------------------------------------------*/

    free (*tsci);
    free (*terr);
    free (*ttime);
    free (*tdq);
    free (*tdarkandglow);
    free (*tlookup);
    free (diff);
}


/* LINFIT: Compute a linear fit of the form y = a + bx to a set of
 ** data points x, y with individual standard deviations sig. Returned
 ** are a, b and their respective probable uncertainties siga and sigb.
 ** Taken from Numerical Recipes (with modifications to handle masked data).
 **
 ** Modified to remove mask: only data to be used should be passed in to
 ** this routine. Optimal weighting (SNR-dependent) may now be used in
 ** the fit. Calculation of uncertainty now includes contributions
 ** from Poisson noise (source and dark) and read noise.
 */

static void linfit (short weight_type, float *x, float *y, float *sig, 
        float *darkandglow, short ndata, 
        float gain, float flat_uncert,
        float *a, float *b, float *siga, float *sigb,
        int i, int j) {

    /* Local variables */
    int k;
    float  S, Sx, Sy, Sxx, Sxy, wt;
    float  snr,power,fit_uncert,dx,dy,ddg,denom;
    float  terma,termb,termc,errterms;
    float  rdns, invrdns2;

    /* WFC3 SMOV mean CDS read noise is ~21 e/pixel */
    rdns = 21.0 / gain; /* DN/pixel/CDS */
    invrdns2 = 1. / (rdns*rdns);

    /* Initialize accumulators and results */
    Sx=0; Sy=0; Sxx=0; Sxy=0; S=0;
    *a=0; *b=0; *siga=0; *sigb=0;

    /* Check for trivial cases */
    if (ndata == 0) {
        (*b)=0.0;
        (*sigb)=1e10;
        return;
    } else if (ndata == 1) {
        if (x[0] == 0.0) {
            (*b)=0.0;
            (*sigb)=1e10;
        } else {
            (*b)    = y[0] / x[0];
            (*sigb) = sig[0] / x[0];
        }
        return;
    }

    /*equal weight is a boolean set up top
      when the weights are calculated, the rejection isn't
      taken into account. Only the good data should be passed
      into this routine!!
    */
    if (weight_type == equal_weight){ 
        power=0;
        snr=0;
    }else{
        /* Determine the correct power term for the optimum weighting.
         * First determine the SNR over the whole interval */
        /* The sig term of the final read is the total rms in the final read
           (Poisson and read) */
        if (sig[ndata-1] > 0.0 ){
            snr=(y[ndata-1]-y[0])/sig[ndata-1];
        }else{
            snr=0.0;
        }
        /* Then use that SNR to select the exponent for the weighting. 
         * These numbers come from a paper by Fixsen (ref. TBA)*/
        if (snr > 100) {
            power =10.0;
        }else{
            if (snr > 50){
                power =6.0;
            }else{
                if(snr > 20){
                    power=3.0;
                }else{
                    if(snr>10){
                        power=1.0;
                    }else{
                        if(snr>5){
                            power=0.4;
                        }else{
                            power=0;
                        }
                    }
                }
            }
        }
    }

    /* Accumulate sums */
    for (k = 0; k < ndata; k++) {
        if (sig[k] != 0) {

            /*wt = 1.0 / (sig[k]*sig[k]);*/

            /* Compute dimensionless weight (ranges from 0 to 1) */
            wt=fabs(pow(fabs(k-((ndata-1)/2.))/((ndata-1)/2.),power)); 
            /* Renormalize dimensionless weight by inverse readnoise squared */
            wt*=invrdns2;

            S += wt;
            Sy += y[k] * wt;
            Sx += x[k] * wt;
            Sxx += x[k] * x[k] * wt;
            Sxy += x[k] * y[k] * wt;
        }
    }

    if (S>0) {
        dx = x[ndata-1]-x[0]; /* total time span of interval being fit */
        dy = y[ndata-1]-y[0]; /* total signal span of interval being fit */
        ddg = darkandglow[ndata-1]-darkandglow[0]; /* total dark of interval */
        denom = (S*Sxx - (Sx*Sx));  

        /* Avoid dividing by zero or other problematic numerical issues */
        if (denom < 1e-6)
            denom = 1e-6;

        if (ndata == 2) { 
            *b = dy/dx;
            *a = y[0] - (*b)*x[0];
            *siga = 1.0; /*not strictly correct, but not really used*/
            fit_uncert = 0;
            *sigb = sqrt((sig[0]/x[0])*(sig[0]/x[0]) + (sig[1]/x[1])*(sig[1]/x[1]));

        } else { /* we can do a proper calculation */
            *b = (S*Sxy - Sx*Sy)/denom;
            *a = (Sxx*Sy - Sx*Sxy)/denom;
            *siga = sqrt(Sxx/denom);
            fit_uncert = sqrt(S/denom);
            *sigb = fit_uncert;
        }

        /* Contributions to the slope uncertainty are:
           - the fit uncertainty, which is uncorrelated read noise
           - the Poisson noise from dark
           - the Poisson noise from source
         */
        terma = (gain*fit_uncert*dx)*(gain*fit_uncert*dx);
        termb =  gain*ddg;
        termc =  gain*dy;

        errterms = terma + termb + termc;

        if (errterms > 0)
            *sigb = (sqrt(errterms)/dx)/gain;   

    } else { /* No points had any weight */
        *a=0.0;
        *b=0.0;
        *sigb=1e9;
        *siga=1e9;
    }
}


/* Find and flag electronic noise spikes. These show up as
 ** large positive or negative deviations for a single sample relative to
 ** the samples on either side of it.
 **
 ** July 2008: Updated to conform to the current method used
 **            in the NICMOS calibration routine. MLS
 ** Mar 2010: Modified to use SPIKE_THRESH, separate from CR_THRESH. HAB
 */

# define SPIKE_THRESH 6.0   /* sigma threshold for spike rejection */ 

static int RejSpikes (float *tsci, float*terr, short *dq, float *diff,
        short nsamp, float thresh, int *max_samp) {

    /* Local variables */
    int   k, nrej;  /* loop index and rejection counter */
    float max_diff;     /* max sample deviation */

    max_diff = 0.0; (*max_samp) = 0;
    nrej = 0;

    /* Loop over samples skipping the first and last samples */
    /* The first sample is treated separately and a spike in the last sample */
    /* is treated as a cosmic ray. */
    for (k = 2; k < nsamp; k++) {
        if (((dq[k-1] & (DQIGNORE | DATAREJECT | SPIKE)) == 0) &&
                ((dq[k]   & (DQIGNORE | DATAREJECT | SPIKE)) == 0) &&
                ((dq[k+1] & (DQIGNORE | DATAREJECT | SPIKE)) == 0)) {

            /* Check for positive or negative spikes */
            /* Here we look for a large difference that has small differences */
            /* on both sides. This spike is just ignored in the fit and */
            /* doesn't affect further reads. */
            if (diff[k-1] - diff[k-2] > SPIKE_THRESH &&
                    diff[k-1] - diff[k]   > SPIKE_THRESH &&
                    diff[k-2] < 0 && diff[k] < 0) {
                if (fabs(diff[k-1]) > max_diff) {
                    max_diff = diff[k-1];
                    (*max_samp) = k-1;
                    nrej++;
                }
            }

        }
    }

    /* If a spike was found, tell the caller*/
    if (nrej > 0)
        return 1;
    else
        return 0;

}


/* REJCRS: Find and flag Cosmic Ray hits. This is done by looking for
 ** large negative-to-positive going changes in the sample values from one
 ** sample to the next.
 **
 ** July 2008: updated to conform to current NICMOS methods  MLS
 ** April 2012: Updated to try and catch negative cr trends as well ML
 */

static int RejCRs (short *dq, float *diff, short nsamp, float thresh,
        int *cr_index) {

    /* Local variables */
    int k;      /* loop index */
    float prev_diff;    /* last useable difference value */
    float max_diff; /*the largest jump */
    float current_diff; /*temp to hold value */

    prev_diff = 0.;
    max_diff = 0.;
    current_diff=0.;

    /* Loop over samples */
    for (k = 1; k < nsamp; k++) {
        if ((dq[k-1] & (DQIGNORE | DATAREJECT )) == 0)
            prev_diff =  diff[k-1];
        if ((dq[k] & (DQIGNORE | DATAREJECT )) == 0) {
            /* A large absolute increase in the difference is used as */
            /* a sign of a CR. */
            current_diff =  fabs(diff[k] - prev_diff);
            if ( current_diff > thresh) {
                if ( current_diff> max_diff) {
                    max_diff = fabs(diff[k] -prev_diff);
                    (*cr_index) = k;
                } 
                prev_diff = diff[k];
            }
        }
    }

    /* If a CR was found, tell the caller */
    if (max_diff != 0)
        return 1;
    else
        return 0;

}


/* REJFIRSTREAD: Check for spikes in the first read of the interval only.
 ** This was formerly part of RejSpikes. 
 */
static int RejFirstRead(short *dq, float *diff, short nsamp, float thresh) {

    /* Local variables */
    int   nrej;     /* return value */

    nrej = 0;

    if (((dq[0] & (DQIGNORE | DATAREJECT | SPIKE)) == 0) &&
            ((dq[1] & (DQIGNORE | DATAREJECT | SPIKE)) == 0)){
        /* first read of interval is not flagged already */
        if (fabs(diff[0]-diff[1]) > thresh){ /* a positive or negative spike */
            nrej++;
        }
    }
    /* If a spike was found, flag the sample that has the max deviation */
    if (nrej > 0)
        return 1;
    else
        return 0;

}


/* This subroutine is used to estimate the dark and amp glow component of the
   signal for a pixel.  It determines the average rate of each pixel and 
   stores that in the rate structure. When called it will return the tot_ADUs 
   for each read of a given pixel. 

   At the moment a static value for the dark and ampglow are used. This
   could be upgraded in the future to use pixel- and read-dependent values
   from an actual dark reference image.
 */  

static void EstimateDarkandGlow (const short nsamp, float *time, float gain,
        float *tot_ADUs) {

    int   i;
    float lineardark, ampglow;

    /* WFC3 SMOV values: average dark = 0.04 e/s */
    lineardark = 0.036 / gain; /** in dn/s **/
    ampglow = 0.0;

    /* create the rate array */
    for (i=0; i <nsamp; i++)
        tot_ADUs[i] = time[i]*lineardark + i*ampglow;
}

