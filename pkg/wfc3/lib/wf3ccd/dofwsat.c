# include <stdio.h>
# include <string.h>

# include "hstcal.h"
# include "hstio.h"
# include "wf3.h"
# include "wf3dq.h"
# include "wf3info.h"
# include "hstcalerr.h"

/*
  This routine compares the full-well saturation image pixels against
  those of the SingleGroup x science array pixels.  If the science array
  pixel value >= the value in the saturation image, the corresponding
  dq array pixel of the SingleGroup x is ORed with the full-well saturation
  flag, SATPIXEL, and the dq array is modified in-place.
  Note the saturation image data is only valid in the region of a
  trimmed full-frame image.  In order to apply the saturation image to the
  science data (full-frame or subarray), the starting X and Y pixels, as well
  as the X and Y sizes need to be obtained, and the saturation image properly
  overlaid on the science data.  Note that formerly the full-well saturation 
  flags were applied during doDQI.

  Michele De La Pena, 2022 March 14
  Initial implementation of the full-well saturation flagging done
  via an image.  Based on a similar routine written for ACS, as well as
  WFC3 doflash to accommodate binned or subarray data.

  Michele De La Pena, 2023 April 17
  Replace hardcoded values with variables isolated to this module.

  Michele De La Pena, 
  The ComputeLimits function has been moved to its own file in /lib.
  
 */

/* SIZE_SV_OVERSCAN is the size of the serial virtual overscan region 
   (in pixels) between amps on the same chip.  
   END_PIX_AC_AMP is the last pixel value (0-based system) of amp A 
   or C (CCD area plus serial phyical overscan */
# define SIZE_SV_OVERSCAN 60
# define END_PIX_AC_AMP 2072 

int doFullWellSat(WF3Info *wf3, SingleGroup *x) {

    /*
      WF3Info *wf3     i: calibration switches, etc.
      SingleGroup *x   io: image to be modified; dq array is updated in-place
     */

    extern int status;

    SingleGroupLine y;		/* scratch space */
    SingleGroup satimage;	/* storage for entire saturation image */
    int extver = 1;			/* get this imset from bias image */
    int rx, ry;				/* for binning bias image down to size of x */
    int x0, y0;				/* offsets of sci image */
    int same_size;			/* true if no binning of ref image required */
    int xdim;				/* number of columns in science image */
    int ydim;				/* number of lines in science image */
    short sum_dq;			/* total of the DQ bits for each pixel */
    int is_subarray = 0;	/* identification of data as a subarray */
    int straddle = 0;		/* subarray starts in A or C and straddles the virtual overscan in the reference image */
    int overstart = -1;		/* location where the overscan starts in the cut science image */
    int xbeg[2], ybeg[2];
    int xend[2], yend[2];

    int FindLine (SingleGroup *, SingleGroupLine *, int *, int *, int *, int *, int *);
    int DetCCDChip (char *, int, int, int *);
    void ComputeLimits(WF3Info *, int, int, int *, int *, int *, int *);

    {unsigned int i;
    for (i = 0; i < 2; i++) {
        xbeg[i] = -1;
        xend[i] = -1;
        ybeg[i] = -1;
        yend[i] = -1;
    }}

    /* Initialize local variables */
    rx = 1;
    ry = 1;
    x0 = 0;
    y0 = 0;

    /* 
       Compute correct extension version number to extract from
       reference image to correspond to CHIP in science data.
       Note: wf3->nimsets is no longer used in DetCCDChip() routine.
    */
    if (DetCCDChip (wf3->satmap.name, wf3->chip, wf3->nimsets, &extver))
        return (status);

    /* Get the first line of saturation image data */
    initSingleGroupLine (&y);
    openSingleGroupLine (wf3->satmap.name, extver, &y);
    if (hstio_err())
        return (status = OPEN_FAILED);

    /* 
       Reference image should already be selected to have the
       same binning factor as the science image.  All we need to
       make sure of is whether the science array is a subarray.

       x0,y0 is the location of the start of the 
       subimage in the reference image.
    */
    if (FindLine (x, &y, &same_size, &rx, &ry, &x0, &y0))
        return (status);
    sprintf(MsgText,"Image has starting location of %d,%d in the reference image", x0, y0);
    trlmessage(MsgText);

    /* Clean up the SingleGroupLine object */
    closeSingleGroupLine (&y);
    freeSingleGroupLine (&y);

    /* If reference data not binned same as input, return with an error */
    if (rx != 1 || ry != 1) {
        sprintf (MsgText,
                "Saturation image and input are not binned to the same pixel size!");
        trlerror (MsgText);
        return (status = SIZE_MISMATCH);
    }

    /* If the science and reference images are not the same size, check to
       see if we are calibrating a subarray image against the untrimmed,
       full-frame, saturation image. The full-frame saturation image will contain
       serial virtual overscan columns in the middle of it, which must be accounted
       for in the value of x0 if the science image subarray is located to the
       right (higher column numbers) of the virtual overscan.  This will be the
       case for science image subarrays located in the amp B and D quadrants.

       For subarrays which START in B or D we can just move the x0 over 60 pixels
       to avoid the serial virtual overscan in the reference image, which is always
       an untrimmed, full-frame image. 

       Otherwise, only part of the subarray overlaps the serial virtual overscan
       and special measures must be taken to avoid it. 
    */

    /* Science image dimensions */
    xdim = x->sci.data.nx;
    ydim = x->sci.data.ny;

    /* Full-frame */
    if (same_size) {
        sprintf (MsgText,"Saturation image and input are the same size.");
        trlmessage (MsgText);
    /* Subarray */
    } else {
        is_subarray = 1;
        sprintf (MsgText,"Saturation image and input are NOT the same size - SUBARRAY found, amp %s", wf3->ccdamp);
        trlmessage (MsgText);

        /*
           ONLY 1 AMP is used to read subarrays, so ampx and ampy should be 
           set to the size of the image for all cases
        */ 
        wf3->ampx = xdim;
        wf3->ampy = ydim;

        /* The image starts in the B or D regions so we can just shift the starting pixel */
        if (x0 > END_PIX_AC_AMP) {
            if (wf3->verbose) {
                sprintf (MsgText,"Subarray starts in B or D region, moved from (%d,%d) to ", x0, y0);
                trlmessage (MsgText);
            }
            x0 += SIZE_SV_OVERSCAN;
            if (wf3->verbose) {
                sprintf (MsgText,"(%d,%d) to avoid virtual overscan in reference image", x0, y0);
                trlmessage (MsgText);
            }

        /* The subarray starts somewhere in A or C and might straddle the virtual overscan region */
        } else {
            if ((x0 + xdim) > END_PIX_AC_AMP) {
                straddle = 1;
                overstart = (END_PIX_AC_AMP + 1) - x0;
            }
        }
    }

    if (wf3->verbose) {
        sprintf(MsgText,"ccdamp = %s, straddle = %d, (xdim,ydim) = (%d,%d), (ampx,ampy) = (%d,%d), (x0,y0) = (%d,%d)",
                wf3->ccdamp, straddle, xdim, ydim, wf3->ampx, wf3->ampy, x0, y0);
        trlmessage(MsgText);
    }

    /* Get the full saturation image */
    initSingleGroup(&satimage);
    getSingleGroup(wf3->satmap.name, extver, &satimage);
    if (hstio_err()) {
        freeSingleGroup(&satimage);
        return (status = OPEN_FAILED);
    }

    /* If the science image is binned, it will be assumed to have
       same size as the reference image, since reading subarrays
       of a binned chip is not supported in current flight software.

       Need to ignore the columns which represent the serial virtual
       overscan between Amps A and B and Amps C and D as there is no
       saturation information for this data.  The values in the
       saturation image are set to zero in this location.
       Bin size: 1x1 ==> overscan = 60 pixels, pixels (2073 - 2132)
       Bin size: 2x2 ==> overscan = 28 pixels, pixels (1037 - 1064)
       Bin size: 3x3 ==> overscan = 20 pixels, pixels (691 - 710)
       The pixel ranges documented here (begin - end) are zero-based.
    */

    /* Determine the limits of real data to process - note the end point
       in the loops is "<" and NOT "<=".
    */
    ComputeLimits(wf3, xdim, ydim, xbeg, ybeg, xend, yend);
    if (wf3->verbose) {
        sprintf(MsgText,"xbeg[0]: %d xend[0]: %d  xbeg[1]: %d  xend[1]: %d", xbeg[0], xend[0], xbeg[1], xend[1]);
        trlmessage(MsgText);
        sprintf(MsgText,"ybeg[0]: %d yend[0]: %d  ybeg[1]: %d  yend[1]: %d", ybeg[0], yend[0], ybeg[1], yend[1]);
        trlmessage(MsgText);
    }

    /* The saturation image already has the gain applied, but the science data
       at this stage in the WFC3 pipeline is still in counts.  It is 
       necessary to divide out the gain from the saturation data before
       any comparison is done. 
    */

    if (wf3->verbose) {
        sprintf(MsgText, "Mean gain: %f", wf3->mean_gain);
        trlmessage(MsgText);
    }
    
    /* Full-frame */
    if (same_size) {

        /* Loop over the lines in the science image */
        {unsigned int  j;
        for (j=ybeg[0]; j < yend[0]; j++) {

            /* Loop over the indices in the line in the science image */
            {unsigned int  i;
            for (i = xbeg[0];  i < xend[0];  i++) {
                /* Flag full-well saturated pixels with 256 dq bit*/             
                if (Pix(x->sci.data, i, j) > (Pix(satimage.sci.data, i, j) / wf3->mean_gain)) {
                    sum_dq = DQPix(x->dq.data, i, j) | SATPIXEL;
			        DQSetPix(x->dq.data, i, j, sum_dq);
                }
            }}

            /* If there is a second Amp in play, complete the processing of the line */ 
            {unsigned int  i;
            for (i = xbeg[1];  i < xend[1];  i++) {
                /* Flag full-well saturated pixels with 256 dq bit*/             
                if (Pix(x->sci.data, i, j) > (Pix(satimage.sci.data, i, j) / wf3->mean_gain)) {
                    sum_dq = DQPix(x->dq.data, i, j) | SATPIXEL;
			        DQSetPix(x->dq.data, i, j, sum_dq);
                }
            }}
        }}
        sprintf(MsgText, "Full-frame full-well saturation image flagging step done.\n");
        trlmessage(MsgText);

    /* Subarray */
    } else {

        /* Loop over all the indices and lines in the science array, and
           match the data to the appropriate indices and lines in the reference image.

           i - index in science image
           k - index in reference image
           j - line in the science image
           l - line in reference image
        */

        {unsigned int j, l;
        for (j = 0, l = y0; j < ydim; j++, l++) {

            /* Working with a subarray so need to apply the proper 
               section from the reference image to the science image.
            */

            {unsigned int i, k;
            for (i = 0, k = x0; i < xdim; i++, k++) {

                /* Increase the value of l to jump over the virtual overscan */
                if (i == overstart)
                    k += SIZE_SV_OVERSCAN;

                /* Flag full-well saturated pixels with 256 dq bit*/             
		        if (Pix(x->sci.data, i, j) > (Pix(satimage.sci.data, k, l) / wf3->mean_gain)) {
			        sum_dq = DQPix(x->dq.data, i, j) | SATPIXEL;
			        DQSetPix(x->dq.data, i, j, sum_dq);
		        }
            }}
        }}
        sprintf(MsgText, "Subarray full-well saturation image flagging step done.\n");
        trlmessage(MsgText);
    }

    freeSingleGroup(&satimage);

    return (status);
}
