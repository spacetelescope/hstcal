# include <stdio.h>

# include "hstcal.h"
# include "hstio.h"
# include "acs.h"
# include "acsdq.h"
# include "acsinfo.h"
# include "hstcalerr.h"

/*
  This routine compares the full-well saturation image pixels against
  those of the SingleGroup x science array pixels.  If the science array
  pixel value >= the value in the saturation image, the corresponding
  SingleGroup x dq array pixel is ORed with the full-well saturation
  flag, SATPIXEL, and the dq array is modified in-place.
  The ACS science data must have already been corrected for bias (BIASCORR)
  and bias level (BLEVCORR), as well as converted to electrons.
  Note the saturation image data is only valid in the region of a
  trimmed full-frame image.  In order to apply the saturation image to the
  bias/blev corrected science data (full-frame or subarray), the starting
  X and Y pixels, as well as the X and Y sizes need to be obtained, and
  the saturation image properly overlaid on the science data.
  Note that formerly the full-well saturation flags were applied during
  doDQI.
  *** The process of applying a saturation image is an awkward retro-fit. ***

  Michele De La Pena, 2020 May 14
  Initial implementation of the full-well saturation flagging done
  via an image.
 */

int doFullWellSat(ACSInfo *acs, SingleGroup *x) {

    /*
      ACSInfo *acs     i: calibration switches, etc.
      SingleGroup *x   io: image to be modified; dq array is updated in-place
     */

    extern int status;

    SingleGroupLine y, z;	/* y and z are scratch space */
    int extver = 1;			/* get this imset from bias image */
    int rx, ry;				/* for binning bias image down to size of x */
    int x0, y0;				/* offsets of sci image */
    int same_size;			/* true if no binning of ref image required */
    int avg = 0;			/* bin2d should sum within each bin */
    int xdim;
    int ydim;				/* number of lines in science image */
    int i, j, k;
	short sum_dq;
    int xbeg, ybeg;			/* Beginning pixels for saturation image overlay */
    int xend, yend;			/* Beginning pixels for saturation image overlay */
    int rsize = 1;
    int sci_bin[2];
    int sci_corner[2];
    int ref_bin[2];
    int ref_corner[2];

    int FindLine (SingleGroup *, SingleGroupLine *, int *, int *,int *, int *, int *);
    int trim1d (SingleGroupLine *, int, int, int, int, int, SingleGroupLine *);
    int DetCCDChip (char *, int, int, int *);
	int GetCorner (Hdr *, int, int *, int *);

    initSingleGroupLine (&y);  /* A single line of the saturation image data */

    /* Initialize local variables */
    rx = 1;
    ry = 1;
    x0 = 0;
    y0 = 0;
    same_size = 1;

    /* 
       Compute correct extension version number to extract from
       reference image to correspond to CHIP in science data.
       Note: acs->nimsets is no longer used in DetCCDChip() routine.
    */
    if (DetCCDChip (acs->satmap.name, acs->chip, acs->nimsets, &extver) )
        return (status);

    /* Get the first line of saturation image data */
    openSingleGroupLine (acs->satmap.name, extver, &y);
    if (hstio_err())
        return (status = OPEN_FAILED);

    /*
      It is necessary to determine if the science array is full-frame
      or a subarray.

      x0,y0 is the location of the start of the subimage in the 
      reference image.

    **
    **	FindLine is a modified version of FindBin routine from CALSTIS.
    **
    */
    if (FindLine (x, &y, &same_size, &rx, &ry, &x0, &y0))
        return (status);

    /* Report information regarding the size of the input data and calibration image */
    if (acs->verbose) {
        sprintf(MsgText,"Ratio of (%d,%d) with offset =(%d,%d)", rx, ry, x0, y0);
        trlmessage(MsgText);
        if (same_size) {
            sprintf(MsgText,"Full-well Saturation image and input are the same size ");
        } else {
            sprintf(MsgText,"Full-well Saturation image and input are NOT the same size ");
        }
        trlmessage(MsgText);
    }

    /* Get the bin factor and "corner" (xpixel and ypixel) where science data actually 
       begins (versus overscan) in the science and the saturation image data */
	if (GetCorner (&x->sci.hdr, rsize, sci_bin, sci_corner))
	    return (status);
	if (GetCorner (&y.sci.hdr, rsize, ref_bin, ref_corner))
	    return (status);

    /* Compute the output size */
    xdim = x->sci.data.nx - (acs->trimx[0] + acs->trimx[1]);
    ydim = x->sci.data.ny - (acs->trimy[0] + acs->trimy[1]);

    if (same_size) {
        sprintf(MsgText, "Full-frame full-well saturation image flagging step.\n");
        trlmessage(MsgText);

        /* Define the beginning and ending pixels */
        xbeg = abs(sci_corner[0]);	
        ybeg = abs(sci_corner[1]);
        xend = xbeg + xdim;
        yend = ybeg + ydim;

        /* Loop over the lines in the science image, excluding the overscan lines */
        {unsigned int  j;
        for (j = ybeg; j < yend; j++) {
            {unsigned int  i;
            for (i = xbeg;  i < xend;  i++) {
                status = getSingleGroupLine (acs->satmap.name, j, &y);
                if (status) {
                    sprintf(MsgText,"Could not read line %d from satmap image.", j+1);
                    trlerror(MsgText);
                }

                /* Flag full-well saturated pixels with 256 dq bit*/             
		        if (Pix (x->sci.data, i, j) > y.sci.line[i]) {
			        sum_dq = DQPix (x->dq.data, i, j) | SATPIXEL;
			        DQSetPix (x->dq.data, i, j, sum_dq);
                }
            }}
        }}
    } else {  /* subarray */
        sprintf(MsgText, "Subarray full-well saturation image flagging step.\n");
        trlmessage(MsgText);

        /* Compute the end point of the X and Y dimension loops */
        xend = x0 + xdim;
        yend = y0 + ydim;

        /* Loop over all the lines in the science array, and
           match them to the appropriate line in the reference image...

           j - index for line in science image
           k - index for line in reference image
           y0 - line in reference image corresponding to line in input image
        */
        {unsigned int j, k;
        for (j = 0, k = y0; j < ydim; j++, k++) {

            /* 
               Working with a sub-array so need to apply the proper 
               section from the reference image to the science image.
            */
            status = getSingleGroupLine (acs->satmap.name, k, &y);
            if (status) {
                sprintf(MsgText,"Could not read line %d from Full-well Saturation image.", k+1);
                trlerror(MsgText);
            }

            {unsigned int i, l;
            for (i = 0, l = x0; i < xdim; i++, l++) {
                /* Flag full-well saturated pixels with 256 dq bit*/             
		        if (Pix (x->sci.data, i, j) > y.sci.line[l]) {
			        sum_dq = DQPix (x->dq.data, i, j) | SATPIXEL;
			        DQSetPix (x->dq.data, i, j, sum_dq);
		        }
            }}
        }}
    }

    closeSingleGroupLine (&y);
    freeSingleGroupLine (&y);

    return (status);
}
