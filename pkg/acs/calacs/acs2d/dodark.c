# include <stdio.h>
# include <stdlib.h>  /* calloc */
# include <math.h>    /* fabs */

# include "hstio.h"
# include "acs.h"
# include "acsinfo.h"
# include "acserr.h"

static float calc_darktime(ACSInfo *);


/* This routine subtracts the dark image from x (in-place).
 For CCD data, the dark image is multiplied by the dark time
 before subtracting.  The dark time is the DARKTIME value, if available.
 Else, it is exposure time + FLASHDUR;
 for post-SM4 non-BIAS WFC, it also INCLUDES the idle time since the last
 flushing of the chip or the readout time.

 For MAMA data, the dark image is just multiplied by the exposure time
 before subtracting.

 Reference image should have been selected to have
 the same binning factor as the science image, so
 assume ratio of bin factors to be 1.

 The value of MEANDARK is calculated based on the weighted average
 of each lines' dark value.  The weighting is based on the percent of
 good pixels in each line, so only pixels not flagged BAD (in some way)
 will contribute to the average, and each line will contribute only
 as much as the line has good pixels.

 Date       Author      Description
 ----       ------      -----------
 1998-06-11 W.J. Hack   Initial ACS Version.
 1999-07-27 W.J. Hack   Added weighting to calculation of MEANDARK.
                        Use same code for all bin cases.
 2000-04-14 W.J. Hack   Fixed the code which applies the gain to each line
                        so that if there is only 1 AMP used per line,
                        it doesn't try to apply a second amp with gain = 0.
 2001-04-16 W.J. Hack   Revised calling sequence for 'multgn1d' to be more
                        general  and accept arbitrary scaling factor instead
                        of being hardwired to 'exptime'.
 2001-12-04 W.J. Hack   Compute darktime from exposure time keywords
                        EXPSTART, EXPEND and use it instead of EXPTIME for
                        scaling dark image.
 2002-04-24 W.J. Hack   Reverted back to using EXPTIME for scaling dark image.
 2012-12-11 P.L. Lim    Now use DARKTIME variable for scaling dark image.
 2016-01-15 P.L. Lim    Now use DARKTIME keyword, if available.
*/
int doDark (ACSInfo *acs2d, SingleGroup *x, float *meandark) {
    /* Parameters:

       ACSInfo     *acs       i: Calibration switches, etc.
       SingleGroup *x        io: Image to be calibrated; written to in-place.
       float       *meandark  o: Mean of dark image values subtracted.
    */
    extern int status;

    SingleGroupLine y, z; /* y and z are scratch space */
    int extver = 1;       /* get this imset from dark image */
    int rx, ry;           /* for binning dark down to size of x */
    int x0, y0;           /* offsets of sci image */
    int same_size;        /* true if no binning of ref image required */
    int avg = 0;          /* bin2d should sum values within each bin */
    int scilines;         /* number of lines in science image */
    int i, j;
    float mean, dark;
    float weight, wdark;  /* weights for line averages */
    int update;
    float darktime;

    int FindLine (SingleGroup *, SingleGroupLine *, int *, int *, int *, int *,
                  int *);
    int sub1d (SingleGroup *, int, SingleGroupLine *);
    int trim1d (SingleGroupLine *, int, int, int, int, int, SingleGroupLine *);
    int DetCCDChip (char *, int, int, int *);
    void AvgSciValLine (SingleGroupLine *, short, float *, float *);
    int multk1d (SingleGroupLine *, float);

    initSingleGroupLine (&y);
    scilines = x->sci.data.ny;

    /* Use DARKTIME from header if available. It is negative if unavailable,
       which then we calculate it instead.
       See https://trac.stsci.edu/ssb/stsci_python/ticket/1167#comment:3 */
    if (acs2d->darktime < 0) {
        darktime = calc_darktime(acs2d);
    } else {
        darktime = acs2d->darktime;
    }


    /* Compute correct extension version number to extract from
       reference image to correspond to CHIP in science data. */
    if (acs2d->pctecorr == PERFORM) {
        if (DetCCDChip (acs2d->darkcte.name, acs2d->chip, acs2d->nimsets,
                        &extver) )
            return (status);
    } else {
        if (DetCCDChip (acs2d->dark.name, acs2d->chip, acs2d->nimsets,
                        &extver) )
            return (status);
    }

    /* Extra info, useful for debugging. */
    if (acs2d->verbose) {
        sprintf(MsgText,
                "Performing dark subtraction on chip %d in imset %d",
                acs2d->chip, extver);
        trlmessage(MsgText);

        if (acs2d->darktime < 0) {
            sprintf(MsgText, "DARKTIME=%f (calculated)", darktime);
        } else {
            sprintf(MsgText, "DARKTIME=%f (direct keyword)", darktime);
        }
        trlmessage(MsgText);
    }

    /* Get the dark image data. */
    if (acs2d->pctecorr == PERFORM) {
        openSingleGroupLine (acs2d->darkcte.name, extver, &y);
    } else {
        openSingleGroupLine (acs2d->dark.name, extver, &y);
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
    if (acs2d->verbose){
        sprintf(MsgText, "Image has an offset of %d,%d", x0, y0);
        trlmessage(MsgText);
    }

    mean = 0.0;
    weight = 0.0;

    /* Bin the dark image down to the actual size of x. */

    initSingleGroupLine (&z);
    allocSingleGroupLine (&z, x->sci.data.nx);

    for (i=0, j=y0; i < scilines; i++,j++) {
        /* We are working with a sub-array and need to apply the
           proper section from the reference image to the science image. */
        if (acs2d->pctecorr == PERFORM) {
            getSingleGroupLine (acs2d->darkcte.name, j, &y);
        } else {
            getSingleGroupLine (acs2d->dark.name, j, &y);
        }

        update = NO;

        if (trim1d (&y, x0, y0, rx, avg, update, &z)) {
            trlerror ("(darkcorr) size mismatch.");
            return (status);
        }

        multk1d(&z, darktime);
        AvgSciValLine (&z, acs2d->sdqflags, &dark, &wdark);

        /* Sum the contribution from each line */
        mean += dark * wdark;
        weight += wdark;

        status = sub1d (x, i, &z);
        if (status)
            return (status);
    }

    freeSingleGroupLine (&z);  /* done with z */
    closeSingleGroupLine (&y);
    freeSingleGroupLine (&y);

    /* Compute the mean for the entire image */
    if (scilines > 0)
        *meandark = mean / weight;
    else
        *meandark = 0.;

    return (status);
}


/* Calculate dark time if there is no direct DARKTIME keyword. */
static float calc_darktime(ACSInfo *acs2d) {
    const float darkscaling = 3.0;  /* Extra idle time (seconds) */
    float darktime;

    /* SBC does not have FLASH keywords. */
    if (acs2d->detector == MAMA_DETECTOR) {
        darktime = acs2d->exptime;
    } else {
        darktime = acs2d->exptime + acs2d->flashdur;

        /* Post-SM4 non-BIAS WFC only.
           TARGNAME unavailable, assume EXPTIME=0 means BIAS */
        if (acs2d->detector == WFC_CCD_DETECTOR && acs2d->expstart > SM4MJD &&
                acs2d->exptime > 0) {
            darktime += darkscaling;
        }
    }

    return (darktime);
}
