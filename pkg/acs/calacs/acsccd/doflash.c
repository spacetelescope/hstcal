/* This file contains:
	doDark
	AvgSciVal

*/

# include <stdio.h>
# include <stdlib.h>		/* calloc */
# include <math.h>		/* fabs */

# include "hstio.h"
# include "acs.h"
# include "acsinfo.h"
# include "acserr.h"

/* This routine subtracts the post-flash image from x (in-place).
   For CCD data, the post-flash image is multiplied by the exposure time and
   divided by the gain before subtracting.  The flash time is given by
   the keyword FLASHDUR, and may represent an interrupted exposure indicated
   by the keyword FLASHSTA.

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

*/

int doFlash (ACSInfo *acsccd, SingleGroup *x, float *meanflash) {

/* arguments:
ACSInfo *acs     	i: calibration switches, etc
SingleGroup *x     io: image to be calibrated; written to in-place
float *meanflash    o: mean of post-flash image values subtracted
*/

    extern int status;

    SingleGroupLine y, z;	/* y and z are scratch space */
    int extver = 1;		/* get this imset from post-flash image */
    int rx, ry;		/* for binning post-flash down to size of x */
    int x0, y0;		/* offsets of sci image */
    int same_size;		/* true if no binning of ref image required */
    int avg = 0;		/* bin2d should sum values within each bin */
    int scilines; 			/* number of lines in science image */
    int i, j;
    float mean, flash;
    float weight, wflash;    /* weights for line averages */
    int update;
    float gain[NAMPS];
    float rn2[NAMPS];       /* only need this to call get_nsegn */
	int ampx;		/* border column for 2amp readout regions */
	int ampy;		/* Boundary values corrected for trim regions */
    int dimx, dimy;
    int offsetx, offsety;

    int FindLine (SingleGroup *, SingleGroupLine *, int *, int *, int *, int *, int *);
    int sub1d (SingleGroup *, int, SingleGroupLine *);
    int trim1d (SingleGroupLine *, int, int, int, int, int, SingleGroupLine *);
    int DetCCDChip (char *, int, int, int *);
    void get_nsegn (int, int, int, int, float *, float*, float *, float *);
	void AvgSciValLine (SingleGroupLine *, short, float *, float *);
	void multgn1d (SingleGroupLine *, int, int, int, float *, float);
        int streq_ic (char *, char *);

	/* Check to see whether we need to do any processing
		at all...  */
	if (acsccd->flashdur <= 0.) {
		sprintf(MsgText,"Post-flash exposure was 0 seconds. FLSHCORR not performed.");
		trlwarn(MsgText);
		addHistoryKw (x->globalhdr, MsgText);
		acsccd->flashcorr = IGNORED;

		/* This is not an error condition, so continue with the remainder
			of the calibrations... */
		return (status);
	}

	/* Flag an aborted Post-Flash exposure in the trailer file comments. */
	if (streq_ic(acsccd->flashstatus,"ABORTED")){
		sprintf (MsgText,"Post-flash STATUS was ABORTED. Post-flash may be compromised.");
		trlwarn (MsgText);
		/* Add this message to the image header as well... */
		addHistoryKw (x->globalhdr, MsgText);
	}

	/* Start with the actual post-flash subtraction now... */
	initSingleGroupLine (&y);

	scilines = x->sci.data.ny;

	/* Compute correct extension version number to extract from
		reference image to correspond to CHIP in science data.
	*/
	if (DetCCDChip (acsccd->flash.name, acsccd->chip, acsccd->nimsets, &extver) )
		return (status);

	if (acsccd->verbose) {
		sprintf(MsgText,"Performing post-flash subtraction on chip %d in imset %d",acsccd->chip, extver);
		trlmessage(MsgText);
	}

	/* Get the post-flash image data. */
	openSingleGroupLine (acsccd->flash.name, extver, &y);
	if (hstio_err())
	    return (status = OPEN_FAILED);

	/* Compare binning of science image and reference image;
	   get same_size and high_res flags, and get info about
	   binning and offset for use by bin2d.
	*/
	if (FindLine (x, &y, &same_size, &rx, &ry, &x0, &y0))
	    return (status);

    if (rx != 1 || ry != 1) {
		sprintf(MsgText,"Reference image and input are not binned to the same pixel size!");
		trlmessage(MsgText);
    }
	if (acsccd->verbose){
		sprintf(MsgText,"Image has an offset of %d,%d",x0,y0);
		trlmessage(MsgText);
	}

    /*
		AMPX,AMPY initialization
	*/
    dimx = x->sci.data.nx;
    dimy = x->sci.data.ny;

    offsetx = (int)(acsccd->offsetx > 0) ? acsccd->offsetx : 0;
    offsety = (int)(acsccd->offsety > 0) ? acsccd->offsety : 0;

	/* Correct the AMP readout boundaries for this offset */
	ampx = ((acsccd->ampx == 0) ? 0 : (int)(acsccd->ampx + offsetx) );
	ampy = ((acsccd->ampy == 0) ? 0 : (int)(acsccd->ampy + offsety) );
	/* Bounds checking to make sure we don't try to apply gain
		and noise outside the bounds of the image.
	*/
	if (ampx >= (dimx - acsccd->offsetx) || ampx > dimx ) ampx = dimx;
	if (ampy >= (dimy - acsccd->offsety) || ampy > dimy ) ampy = dimy;

	acsccd->ampx = ampx;
	acsccd->ampy = ampy;
	mean = 0.0;
    weight = 0.0;

	/* Multiply the post-flash image by the exposure time and divide by the
	   atodgain, and subtract it from x.
	*/

    for (i = 0; i < NAMPS; i++) {
        gain[i] = 0.;
        rn2[i] = 0.;
    }
    get_nsegn (acsccd->detector, acsccd->chip, acsccd->ampx, acsccd->ampy, acsccd->atodgain, acsccd->readnoise, gain, rn2);

	/* Bin or sub-sample the post-flash image down to the actual size of x. */

	initSingleGroupLine (&z);
	allocSingleGroupLine (&z, x->sci.data.nx);
	for (i=0, j=y0; i < scilines; i++,j++) {

	/* We are working with a sub-array and need to apply the
		proper section from the reference image to the science image.
	*/
		getSingleGroupLine (acsccd->flash.name, j, &y);

        /*
		rx = 1;
		*/
        update = NO;

	    if (trim1d (&y, x0, y0, rx, avg, update, &z)) {
			trlerror ("(flshcorr) size mismatch.");
			return (status);
	    }

	    multgn1d (&z, j, acsccd->ampx, acsccd->ampy, gain, acsccd->flashdur);

	    AvgSciValLine (&z, acsccd->sdqflags, &flash, &wflash);
		/* Sum the contribution from each line */
		mean += flash * wflash;
		weight += wflash;

	    status = sub1d (x, i, &z);
		if (status)
			return (status);
	}
	freeSingleGroupLine (&z);			/* done with z */

/*	} */

	closeSingleGroupLine (&y);
	freeSingleGroupLine (&y);

	/* Compute the mean for the entire image */
	if (scilines > 0)
		*meanflash = mean / weight;
	else
		*meanflash = 0.;

	return (status);
}
