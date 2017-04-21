/* ACS -- Detect and mark SINK pixels in the DQ array. */
# include <stdio.h>

# include "hstio.h"
# include "acs.h"
# include "acsinfo.h"
# include "err.h"

/* local mask values */
# define SINKPIXEL (short)1024

static int DetSinkChip (char *, int, int *);


/* See Section 6 of ACS ISR by Ryon & Grogin (2017).
   SINKCORR is based on WFC3/UVIS but not identical.
   We do not mark flags by dates nor use crazy RAZ format.
   This should be the last step in ACSCDD, before ACSCTE.
   Also see https://github.com/spacetelescope/hstcal/issues/44

   Relevant info from Section 6:

   * Pix value from SNKCFILE = sink pixel depth from INS calib
   * Downstream pixel (-1, closer to amp) = -1.0 (optional)
   * Upstream pixel (+1) = 999.0 (always)
   * Upstream pixels (+2 to +12) depends on sink pixel depth.

   Algorithm from Section 6 along a column:

   1. Look for pixel < -1 in SNKCFILE, set DQ to 1024.
   2. If downstream (-1) pixel == -1, set DQ to 1024.
   3. If upstream (+1 to +n) pixel > SCI, set DQ to 1024 until
      pixel <= SCI or pixel == 0 or pixel == another sink pixel.
*/
int SinkDetect(ACSInfo *acs, SingleGroup *x) {
    /* arguments:
       ACSInfo *acs       i: calibration switches, etc
       SingleGroup *x    io: image to be calibrated; written to in-place
    */
    extern int status;

    int i, j, jdone, jend, jstep, n;  /* counters */
    short *keep_going, dqval;
    int extver;      /* get this imset from sink image */
    int dimx, refx, dimy;
    int rx, ry;      /* for binning sink image down to size of x */
    int x0, y0;      /* offsets of sci image */
    int same_size;   /* true if no binning of ref image required */
    float cur_sinkpix, *cur_sci;
    FloatHdrData sinkref;  /* array to store sink image */

    int FindLineHdr(Hdr *, Hdr *, int, int, int *, int *, int *, int *, int *);

    if (acs->sinkcorr != PERFORM)
        return (status);

    /* Initialize local variables */
    extver = 1;
    rx = 1;
    ry = 1;
    x0 = 0;
    y0 = 0;
    same_size = 1;
    jdone = 0;
    n = 0;

    /* Compute correct extension version number to extract from
       reference image to correspond to CHIP in science data. */
    if (DetSinkChip(acs->sink.name, acs->chip, &extver))
        return (status);

    /* Get sink reference image. */
    initFloatHdrData(&sinkref);
    getFloatHD(acs->sink.name, "SCI", extver, &sinkref);

    /* Extract relevant portion from reference image. */
    dimx = x->sci.data.tot_nx;
    dimy = x->sci.data.tot_ny;
    refx = sinkref.data.tot_nx;

    if (FindLineHdr (&x->sci.hdr, &sinkref.hdr, dimx, refx,
                     &same_size, &rx, &ry, &x0, &y0))
        return (status);

    if (acs->verbose) {
        sprintf(MsgText, "Ratio of (%d,%d) with offset =(%d,%d)",
                rx, ry, x0, y0);
        trlmessage(MsgText);
        if (same_size) {
            sprintf(MsgText, "SINK image and input are the same size ");
        } else {
            sprintf(MsgText, "SINK image and input are NOT the same size ");
        }
        trlmessage(MsgText);
    }

    if ((rx != 1) || (ry != 1)) {
        trlerror ("(sinkcorr) binned data is not supported.");
        return (status = INVALID_VALUE);
    }

    /* Always traverse from sink pixel head to tail. */
    if (extver == 1) {
        j = 0;
        jend = dimy - 1;
        jstep = 1;
    } else {  /* extver == 2 */
        j = dimy - 1;
        jend = 0;
        jstep = -1;
    }

    /* Allocate data array */
    keep_going = (short *) calloc (dimx, sizeof(short));
    if (keep_going == NULL) {
        trlerror ("Couldn't allocate memory for scratch array in SINKCORR.");
        return (status = OUT_OF_MEMORY);
    }
    cur_sci = (float *) calloc (dimx, sizeof(float));
    if (cur_sci == NULL) {
        trlerror ("Couldn't allocate memory for scratch array in SINKCORR.");
        return (status = OUT_OF_MEMORY);
    }

    /* Flag sink pixels in input DQ. */
    while (jdone == 0) {
        for (i = 0; i < dimx; i++) {

            /* Sink image is always untrimmed fullframe.
               While, since image at this point is still untrimmed. */
            cur_sinkpix = Pix(sinkref.data, i + x0, j + y0);

            /* This is either sink or downstream pixel, flag it.
               It comes with a tail upstream.
               Does not matter if it sits on another tail. */
            if (cur_sinkpix < 0) {
                if (keep_going[i] == 0) {
                    keep_going[i] = 1;
                }
                if (cur_sinkpix < -1) {  /* The actual sink pixel */
                    cur_sci[i] = Pix(x->sci.data, i, j);
                }
                dqval = SINKPIXEL | DQPix(x->dq.data, i, j);
                DQSetPix(x->dq.data, i, j, dqval);
                n++;

            /* This is upstream pixel. */
            } else if (cur_sinkpix > 0) {
                /* Tail is still continuous. */
                if (keep_going[i] == 1) {
                    /* This upstream pixel > SCI at sink pixel, flag it. */
                    if (cur_sinkpix > cur_sci[i]) {
                        dqval = SINKPIXEL | DQPix(x->dq.data, i, j);
                        DQSetPix(x->dq.data, i, j, dqval);
                        n++;
                    /* Tail ended, stop flagging. */
                    } else {
                        keep_going[i] = 0;
                    }
                }

            /* No sink pixel, stop flagging. */
            } else if (keep_going[i] == 1) {
                keep_going[i] = 0;
            }
        } /* end i loop */

        j += jstep;
        if (jstep >= 0) {
            if (j > jend)
                jdone = 1;
        } else {
            if (j < jend)
                jdone = 1;
        }
    } /* end j loop */

    if (acs->verbose) {
        sprintf(MsgText, "Sink pixels flagged = %d", n);
        trlmessage(MsgText);
    }

    free(keep_going);
    free(cur_sci);
    freeFloatHdrData(&sinkref);

    return (status);
}


/* Find the corresponding EXT from SNKCFILE. Adapted from DetCCDChip. */
static int DetSinkChip (char *fname, int chip, int *extver) {
    /* parameters:
       char *fname   i: name of file
       int chip      i: CHIP ID to be found
       int nimsets;  i: number of IMSETS in image
       int *extver   o: extension (IMSET) from file corresponding
       to input image chip ID
    */
    extern int status;

    Hdr scihdr, prihdr;
    IODescPtr ip;
    int ccdchip;        /* CHIP id from reference header */
    int n, foundit;
    int nextend;

    int GetKeyInt (Hdr *, char *, int, int, int *);

    initHdr (&scihdr);
    initHdr (&prihdr);
    *extver = 0;
    ip = NULL;
    foundit = NO;
    ip = openInputImage (fname, "", 0);
    getHeader (ip, &prihdr);
    closeImage (ip);

    /* Find out how many extensions there are in this file. */
    if (GetKeyInt (&prihdr, "NEXTEND", USE_DEFAULT, 1, &nextend))
        return (status);

    /* Loop over all the extensions in the reference file
       to search for the extension which corresponds to the desired
       CCDCHIP id of the exposure. */
    for (n = 1; n <= nextend ; n++) {
        ip = openInputImage (fname, "SCI", n);

        getHeader (ip, &scihdr);

        if (ip != NULL)
            closeImage (ip);

        /* Get CCD-specific parameters. */
        if (GetKeyInt (&scihdr, "CCDCHIP", USE_DEFAULT, 1, &ccdchip))
            return (status);

        if (ccdchip == chip) {
            /* Found it! */
            *extver = n;
            foundit = YES;
            break;
        } else {
            /* Check next extension for CHIP id */
            ccdchip = 0;
        }
    }

    freeHdr(&scihdr);
    freeHdr(&prihdr);

    if (foundit == NO) {
        sprintf (MsgText, "No Reference Data found for chip %d", chip);
        trlerror (MsgText);
        return (status = NO_CHIP_FOUND);
    }
    return (status);
}
