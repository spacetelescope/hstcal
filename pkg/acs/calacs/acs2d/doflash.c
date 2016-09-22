/* This file contains:
     doFlash
*/
# include <stdio.h>
# include <stdlib.h>  /* calloc */
# include <math.h>    /* exp, pow, sinh */

# include "hstio.h"
# include "acs.h"
# include "acsinfo.h"
# include "acserr.h"


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
   Pey Lian Lim, 2016 Sep 22:
     Added correction factor based on AMP and FLASHDUR.
 */
int doFlash (ACSInfo *acs2d, SingleGroup *x, float *meanflash) {
    /* arguments:
       ACSInfo *acs2d     i: calibration switches, etc
       SingleGroup *x    io: image to be calibrated; written to in-place
       float *meanflash   o: mean of post-flash image values subtracted
    */

    extern int status;

    SingleGroupLine y, z;  /* y and z are scratch space */
    int extver = 1;        /* get this imset from post-flash image */
    int rx, ry;            /* for binning post-flash down to size of x */
    int x0, y0;            /* offsets of sci image */
    int same_size;         /* true if no binning of ref image required */
    int avg = 0;           /* bin2d should sum values within each bin */
    int dimx, dimy;        /* number of cols and lines in science image */
    int offsetx, offsety;
    int ampx;        /* border column for 2-amp readout regions */
    int ampy;        /* Boundary values corrected for trim regions */
    int i, j;
    float corrfac[NAMPS]; /* Correction factor to be multiplied to postflash */
    float mean, flash;
    float weight, wflash;  /* weights for line averages */
    const int update = NO;

    int FindLine (SingleGroup *, SingleGroupLine *, int *, int *, int *, int *,
                  int *);
    int sub1d (SingleGroup *, int, SingleGroupLine *);
    int trim1d (SingleGroupLine *, int, int, int, int, int, SingleGroupLine *);
    int DetCCDChip (char *, int, int, int *);
    void AvgSciValLine (SingleGroupLine *, short, float *, float *);
    int multk1d (SingleGroupLine *a, float k);
    int streq_ic (char *, char *);

    /* Check to see whether we need to do any processing at all */
    if (acs2d->flashdur <= 0.) {
        sprintf(MsgText,
            "Post-flash exposure was 0 seconds. FLSHCORR not performed.");
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
        sprintf (MsgText,
            "Post-flash STATUS was ABORTED. Post-flash may be compromised.");
        trlwarn (MsgText);
        addHistoryKw (x->globalhdr, MsgText);
    }

    /* Start with the actual post-flash subtraction now. */

    initSingleGroupLine (&y);

    offsetx = (int)(acs2d->offsetx > 0) ? acs2d->offsetx : 0;
    offsety = (int)(acs2d->offsety > 0) ? acs2d->offsety : 0;

    dimx = x->sci.data.nx;
    dimy = x->sci.data.ny;

    /* Correct the AMP readout boundaries for this offset */
    ampx = ((acs2d->ampx == 0) ? 0 : (int)(acs2d->ampx + offsetx));
    ampy = ((acs2d->ampy == 0) ? 0 : (int)(acs2d->ampy + offsety));

    /* Bounds checking to make sure we don't try to apply flash
       correction outside the bounds of the image.
       This would apply if using only 1 AMP on each WFC chip when
       ampx is given in CCDTAB as something greater image size.

       We need to make sure that if the ampx value extends into the
       overscan at the end of the line, ampx gets automatically
       moved to cover the whole line. This allows all AMPX and AMPY
       values to be specified in CCDTAB in trimmed coordinates.
    */
    if (ampx >= (dimx - acs2d->offsetx) || ampx > dimx)
        ampx = dimx;
    if (ampy >= (dimy - acs2d->offsety) || ampy > dimy)
        ampy = dimy;

    /* Set up correction factor that depends on AMP and FLASHDUR.
       https://github.com/spacetelescope/hstcal/issues/12 */
    corrfac[AMP_A] = 1.0 - (-0.27143973 *
                            sinh(0.89499787 * acs2d->flashdur - 1.0) /
                            exp(1.64930466 * acs2d->flashdur - 1.0));
    corrfac[AMP_B] = 1.0 - (-0.38018589 *
                            sinh(1.63426073 * acs2d->flashdur - 1.0) /
                            exp(2.19362160 * acs2d->flashdur - 1.0));
    corrfac[AMP_C] = 1.0 + (0.04802005 /
                            (pow(acs2d->flashdur, 2.12968265) + 0.01922919));
    corrfac[AMP_D] = 1.0 - (-0.32215275 *
                            sinh(1.19033580 * acs2d->flashdur - 1.0) /
                            exp(1.59900158 * acs2d->flashdur - 1.0));

    /* Compute correct extension version number to extract from
       reference image to correspond to CHIP in science data. */
    if (acs2d->pctecorr == PERFORM) {
        if (DetCCDChip (acs2d->flashcte.name, acs2d->chip, acs2d->nimsets,
                        &extver))
            return (status);
    } else {
        if (DetCCDChip (acs2d->flash.name, acs2d->chip, acs2d->nimsets,
                        &extver))
            return (status);
    }

    if (acs2d->verbose) {
        sprintf(MsgText,
                "Performing post-flash subtraction on chip %d in imset %d",
                acs2d->chip, extver);
        trlmessage(MsgText);
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
       binning and offset for use by bin2d. */
    if (FindLine (x, &y, &same_size, &rx, &ry, &x0, &y0))
        return (status);

    if (rx != 1 || ry != 1) {
        sprintf(MsgText,
            "Reference image and input are not binned to the same pixel size!");
        trlmessage(MsgText);
    }

    if (acs2d->verbose) {
        sprintf(MsgText, "Image has an offset of %d,%d", x0, y0);
        trlmessage(MsgText);
    }

    mean = 0.0;
    weight = 0.0;

    /* Bin the post-flash image down to the actual size of x. */

    initSingleGroupLine(&z);
    allocSingleGroupLine(&z, dimx);

    for (j = 0; j < ampy; j++) {
        /* We are working with a sub-array and need to apply the
           proper section from the reference image to the science image. */
        if (acs2d->pctecorr == PERFORM) {
            getSingleGroupLine (acs2d->flashcte.name, y0 + j, &y);
	} else {
            getSingleGroupLine (acs2d->flash.name, y0 + j, &y);
	}

        if (trim1d(&y, x0, y0, rx, avg, update, &z)) {
            trlerror("(flshcorr) size mismatch.");
            return (status);
        }

        /* This region corresponds to AMP_C,
           if it is even used for this observation. */
        for (i = 0; i < ampx; i++) {
            /* Apply correction factor to reference file data. */
            z->sci.line[i] = corrfac[AMP_C] * z->sci.line[i];
            z->err.line[i] = corrfac[AMP_C] * z->err.line[i];
        }

        /* This region corresponds to AMP_D,
           if it is even used for this observation. */
        for (i = ampx; i < dimx; i++) {
            /* Apply correction factor to reference file data. */
            z->sci.line[i] = corrfac[AMP_D] * z->sci.line[i];
            z->err.line[i] = corrfac[AMP_D] * z->err.line[i];
        }

        multk1d(&z, acs2d->flashdur);

        AvgSciValLine(&z, acs2d->sdqflags, &flash, &wflash);

        /* Sum the contribution from each line */
        mean += flash * wflash;
        weight += wflash;

        status = sub1d(x, j, &z);
        if (status)
            return (status);
    }

    for (j = ampy; j < dimy; j++) {
        /* We are working with a sub-array and need to apply the
           proper section from the reference image to the science image. */
        if (acs2d->pctecorr == PERFORM) {
            getSingleGroupLine (acs2d->flashcte.name, y0 + j, &y);
	} else {
            getSingleGroupLine (acs2d->flash.name, y0 + j, &y);
	}

        if (trim1d(&y, x0, y0, rx, avg, update, &z)) {
            trlerror("(flshcorr) size mismatch.");
            return (status);
        }

        /* This region corresponds to AMP_A,
           if it is even used for this observation. */
        for (i = 0;  i < ampx;  i++) {
            /* Apply correction factor to reference file data. */
            z->sci.line[i] = corrfac[AMP_A] * z->sci.line[i];
            z->err.line[i] = corrfac[AMP_A] * z->err.line[i];
        }

        /* This region corresponds to AMP_B,
           if it is even used for this observation. */
        for (i = ampx;  i < dimx;  i++) {
            /* Apply correction factor to reference file data. */
            z->sci.line[i] = corrfac[AMP_B] * z->sci.line[i];
            z->err.line[i] = corrfac[AMP_B] * z->err.line[i];
        }

        multk1d(&z, acs2d->flashdur);

        AvgSciValLine(&z, acs2d->sdqflags, &flash, &wflash);

        /* Sum the contribution from each line */
        mean += flash * wflash;
        weight += wflash;

        status = sub1d(x, j, &z);
        if (status)
            return (status);
    }

    freeSingleGroupLine(&z);  /* done with z */

    closeSingleGroupLine(&y);
    freeSingleGroupLine(&y);

    /* Compute the mean for the entire image */
    if (dimy > 0)
        *meanflash = mean / weight;
    else
        *meanflash = 0.;

    return (status);
}
