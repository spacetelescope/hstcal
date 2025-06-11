/* This file contains:
 doFlash
 */

# include <stdio.h>
# include <stdlib.h>  /* calloc */
# include <math.h>    /* fabs */
# include <stdbool.h> /* true, false */

# include "hstcal.h"
# include "hstio.h"
# include "acs.h"
# include "acsinfo.h"
# include "hstcalerr.h"

/* This routine subtracts the post-flash image from x (in-place).
 For CCD data, the post-flash image is multiplied by flash time
 before subtracting.  The flash time is given by the keyword
 FLASHDUR, and may represent an interrupted exposure indicated
 by the keyword FLASHSTA.

 MAMA detector should not be using this.

 Reference image should have been selected to have
 the same binning factor as the science image, so
 assume ratio of bin factors to be 1.

 The value of MEANFLSH is calculated based on the weighted average
 of each lines' post-flash value.  The weighting is based on the percent of
 good pixels in each line, so only pixels not flagged BAD (in some way)
 will contribute to the average, and each line will contribute only
 as much as the line has good pixels.

 Warren Hack, 2000 Sept 12:
 Initial ACS Version.
 Warren Hack, 2000 Nov 10:
 Revised to fully support multi-amp configurations.
 Warren Hack, 2001 April 16:
 Revised to scale FLSHFILE by FLASHDUR instead of EXPTIME. Required
 modification to calling sequence to 'multgn1d'.
 Pey Lian Lim, 2012 Dec 12:
 Moved FLSHCORR from ACSCCD to ACS2D.
 Modified algorithm to be similar to dodark.c.
 Pey Lian Lim, 2025 May 8:
 Applied scaling for linear decrease over time of flash rate.
 */

int doFlash (ACSInfo *acs2d, SingleGroup *x, float *meanflash) {

    /* arguments:
       ACSInfo *acs2d     i: calibration switches, etc
       SingleGroup *x    io: image to be calibrated; written to in-place
       float *meanflash   o: mean of post-flash image values subtracted
    */

    extern int status;

    Hdr hdr_ptr;
    SingleGroupLine y, z;  /* y and z are scratch space */
    int extver = 1;        /* get this imset from post-flash image */
    int rx, ry;            /* for binning post-flash down to size of x */
    int x0, y0;            /* offsets of sci image */
    int same_size;         /* true if no binning of ref image required */
    int avg = 0;           /* bin2d should sum values within each bin */
    int scilines;          /* number of lines in science image */
    int i, j;
    float mean, flash;
    float weight, wflash;  /* weights for line averages */
    double slope;          /* m of linear decrease over time */
    double t0;             /* t_0 of linear decrease over time */
    double avglastflash;   /* F of linear decrease over time */
    float s_lastflash = 1.0; /* s of linear decrease over time */
    bool is_lastflash = true; /* flag to apply s or not */
    int update;

    int LoadHdr (char *, Hdr *);
    int GetKeyDbl(Hdr *, char *, int, double, double *);
    int FindLine (SingleGroup *, SingleGroupLine *, int *, int *, int *, int *, int *);
    int sub1d (SingleGroup *, int, SingleGroupLine *);
    int trim1d (SingleGroupLine *, int, int, int, int, int, SingleGroupLine *);
    int DetCCDChip (char *, int, int, int *);
    void AvgSciValLine (SingleGroupLine *, short, float *, float *);
    int multk1d (SingleGroupLine *a, float k);
    int streq_ic (char *, char *);
    void WhichError (int);

    /* Check to see whether we need to do any processing at all */
    if (acs2d->flashdur <= 0.) {
        snprintf(MsgText, sizeof(MsgText), "Post-flash exposure was 0 seconds. FLSHCORR not performed.");
        trlwarn(MsgText);
        addHistoryKw (x->globalhdr, MsgText);
        acs2d->flashcorr = IGNORED;

        /* This is not an error condition, so continue with the remainder
           of the calibrations... */
        return (status);
    }

    /* Flag an aborted Post-Flash exposure in the trailer file comments. */
    /* Add this message to the image header as well. */
    if (streq_ic(acs2d->flashstatus, "ABORTED")) {
        snprintf (MsgText, sizeof(MsgText), "Post-flash STATUS was ABORTED. Post-flash may be compromised.");
        trlwarn (MsgText);
        addHistoryKw (x->globalhdr, MsgText);
    }

    /* Start with the actual post-flash subtraction now. */

    initSingleGroupLine (&y);

    scilines = x->sci.data.ny;

    /* Compute correct extension version number to extract from
       reference image to correspond to CHIP in science data.

       Also load reference image primary header for scaling that is
       needed because no new reference image is generated going forward
       and we know that the flash rate decreases linearly over time.
    */
    if (acs2d->pctecorr == PERFORM) {
        if (DetCCDChip (acs2d->flashcte.name, acs2d->chip, acs2d->nimsets,
                        &extver))
            return (status);
        if (LoadHdr(acs2d->flashcte.name, &hdr_ptr)) {
            WhichError (status);
            return (status);
        }
    } else {
        if (DetCCDChip (acs2d->flash.name, acs2d->chip, acs2d->nimsets,
                        &extver))
            return (status);
        if (LoadHdr(acs2d->flash.name, &hdr_ptr)) {
            WhichError (status);
            return (status);
        }
    }

    /* Get linear decrease parameters from reference image primary header.

       The scaling factor (s) by which to multiply the SCI and ERR arrays
       of the last available post-flash reference image prior to
       flash correction in CALACS:

       s = (m * (t_obs - t_0)) / F + 1

       where:

       m = SLOPE in reference image primary header
       t_obs = DATE-OBS of input image (MJD), approx. by acs2d->expstart
       t_0 = T0 in reference image primary header (MJD)
       F = AVGFLASH in reference image primary header (e/s)

       Example reference image header content:

       SLOPE   =         -6.52415E-05
       T0      =                60126
       AVGFLASH=          11.76780605
    */
    if (GetKeyDbl(&hdr_ptr, "SLOPE", NO_DEFAULT, 0.0, &slope) ||
            GetKeyDbl(&hdr_ptr, "T0", NO_DEFAULT, 0.0, &t0) ||
            GetKeyDbl(&hdr_ptr, "AVGFLASH", NO_DEFAULT, -1.0, &avglastflash)) {
        status = 0;  /* Keep going */
        is_lastflash = false;
    } else {
        s_lastflash = (float)((slope * (acs2d->expstart - t0)) / avglastflash + 1.0);
    }

    if (acs2d->verbose) {
        trlmessage("Performing post-flash subtraction on chip %d in imset %d", acs2d->chip, extver);
        if (is_lastflash) {
            trlmessage("Linear decrease over time scaling for post-flash: %e", s_lastflash);
        } else {
            trlmessage("Linear decrease over time scaling for post-flash: SKIPPED");
        }
    }

    /* Get the post-flash image data. */
    if (acs2d->pctecorr == PERFORM) {
        openSingleGroupLine (acs2d->flashcte.name, extver, &y);
    } else {
        openSingleGroupLine (acs2d->flash.name, extver, &y);
    }
    if (hstio_err())
        return (status = OPEN_FAILED);

    /* Compare binning of science image and reference image;
       get same_size and high_res flags, and get info about
       binning and offset for use by bin2d.
    */
    if (FindLine (x, &y, &same_size, &rx, &ry, &x0, &y0))
        return (status);

    if (rx != 1 || ry != 1) {
        trlwarn("Reference image and input are not binned to the same pixel size!");
    }

    if (acs2d->verbose) {
        trlmessage("Image has an offset of %d,%d", x0, y0);
    }

    mean = 0.0;
    weight = 0.0;

    /* Bin the post-flash image down to the actual size of x. */

    initSingleGroupLine(&z);
    allocSingleGroupLine(&z, x->sci.data.nx);
    for (i=0, j=y0; i < scilines; i++,j++) {

        /* We are working with a sub-array and need to apply the
           proper section from the reference image to the science image.
        */
        if (acs2d->pctecorr == PERFORM) {
            getSingleGroupLine (acs2d->flashcte.name, j, &y);
	} else {
            getSingleGroupLine (acs2d->flash.name, j, &y);
	}

        update = NO;

        if (trim1d(&y, x0, y0, rx, avg, update, &z)) {
            trlerror("(flshcorr) size mismatch.");
            return (status);
        }

        if (is_lastflash) {
            multk1d(&z, acs2d->flashdur * s_lastflash);
        } else {
            multk1d(&z, acs2d->flashdur);
        }

        AvgSciValLine(&z, acs2d->sdqflags, &flash, &wflash);

        /* Sum the contribution from each line */
        mean += flash * wflash;
        weight += wflash;

        status = sub1d(x, i, &z);
        if (status)
            return (status);
    }
    freeSingleGroupLine(&z);  /* done with z */

    closeSingleGroupLine(&y);
    freeSingleGroupLine(&y);

    /* Compute the mean for the entire image */
    if (scilines > 0)
        *meanflash = mean / weight;
    else
        *meanflash = 0.;

    return (status);
}
