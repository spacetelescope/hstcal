# include <stdio.h>
# include <string.h>

# include "hstio.h"
# include "wf3.h"
# include "wf3info.h"
# include "err.h"

/* This routine subtracts the bias image from x (in-place).
   For WF3 science data, it will normally be the
   case that this correction will be performed before the data has
   been combined for cosmic-ray rejection.
   Note that the bias image is assumed to have already been scaled by
   the gain, and binned the same as the input image.

   In order to accomodate the larger WFC images, the bias images will be
   applied line-by-line to the input image.  

   Warren Hack, 1998 June 3:
   Initial version for ACS. 
   
   Howard Bushouse, 2000 Aug 29:
   Initial version for WFC3.
   
   H.Bushouse, 2001 Dec 4:
   Added error message and return after FindLine if reference data
   are found to be binned differently than science image.
   
   H.Bushouse, 2009 Mar 9:
   Updated to compute correct x-offset values for amp B and D subarray
   science images, which needs to account for the columns of serial 
   virtual overscan that are in the middle of a full-frame bias reference
   image.
   
   M. Sosey, 2013 Feb 7:
   Updated to take into account the bias reference files who have values for their
   CCDAMP of SINGLE_AMP or SINGLE_OR_ALL. This was not being tested against and their
   are now reference files in CDBS with this designation when any single amp is used. 
   
 */

int doBias (WF3Info *wf3, SingleGroup *x) {

    /* arguments:
       WF3Info *wf3     i: calibration switches, etc
       SingleGroup *x	io: image to be calibrated; written to in-place
     */

    extern int status;

    SingleGroupLine y, z;	/* y and z are scratch space */
    int extver = 1;		/* get this imset from bias image */
    int rx, ry;		/* for binning bias image down to size of x */
    int x0, y0;		/* offsets of sci image */
    int same_size;		/* true if no binning of ref image required */
    int avg = 0;		/* bin2d should sum within each bin */
    int scilines; 		/* number of lines in science image */
    int i, j;
    int update;
    char biasamp[SZ_CBUF+1];	/* CCDAMP string for bias ref file */

    int FindLine (SingleGroup *, SingleGroupLine *, int *, int *,int *,
            int *, int *);
    int sub1d (SingleGroup *, int, SingleGroupLine *);
    int trim1d (SingleGroupLine *, int, int, int, int, int,
            SingleGroupLine *);
    int DetCCDChip (char *, int, int, int *);
    int GetKeyStr (Hdr *, char *, int, char *, char *, int);
    int streq_ic (char *, char *); /* case insensitive string equal */

    if (wf3->biascorr != PERFORM)
        return (status);

    initSingleGroupLine (&y);

    scilines = x->sci.data.ny;

    /* Initialize local variables */
    rx = 1;
    ry = 1;
    x0 = 0;
    y0 = 0;
    same_size = 1;

    /* Compute correct extension version number to extract from
     ** reference image to correspond to CHIP in science data. */
    if (DetCCDChip (wf3->bias.name, wf3->chip, wf3->nimsets, &extver))
        return (status);		

    /* Get the first line of bias image data. */
    openSingleGroupLine (wf3->bias.name, extver, &y);
    if (hstio_err())
        return (status = OPEN_FAILED);

    /* 
       Reference image should already be selected to have the
       same binning factor as the science image.  All we need to
       make sure of is whether the science array is a sub-array of
       the bias image.  

       x0,y0 is the location of the start of the 
       subimage in the reference image.
     */
    if (FindLine (x, &y, &same_size, &rx, &ry, &x0, &y0))
        return (status);

    /* If the science and reference images are not the same size, check to
       see if we're calibrating a subarray image with a full-frame, 4-amp
       bias image. The full-frame bias will contain serial virtual overscan
       columns in the middle of it, which must be accounted for in the value
       of x0 if the science image subarray is located to the right (higher
       column numbers) of the virtual overscan, which will be the case for
       science image subarrays located in the amp B and D quadrants.
     */
    if (!same_size) {
        sprintf(MsgText,"Bias ref and science image not the same size");
        trlmessage(MsgText);
        /* Retrieve CCDAMP value from bias image header */
        if(GetKeyStr(y.globalhdr, "CCDAMP", NO_DEFAULT, "", biasamp, SZ_CBUF))
            return (status);

        /* If the bias image uses all 4 amps and the science image uses
           either amp B or D, then add the virtual overscan columns to x0 */
        if ((streq_ic(biasamp,"ABCD") || streq_ic(biasamp,"ANY") \
                    || streq_ic(biasamp,"SINGLE_AMP") ||streq_ic(biasamp,"SINGLE_OR_ALL") ) &&
                (streq_ic(wf3->ccdamp,"B") || streq_ic(wf3->ccdamp,"D"))) {
            x0 += 60;

        }
    }

    if (wf3->verbose) {
        sprintf (MsgText, "Ratio of (%d,%d) with offset =(%d,%d)",
                rx,ry,x0,y0);
        trlmessage(MsgText);
        if (same_size) {
            sprintf(MsgText,"BIAS image and input are the same size ");
        } else {
            sprintf(MsgText,"BIAS image and input are NOT the same size ");
        }
        trlmessage(MsgText);
    }

    /* Return with error if reference data not binned same as input */
    if (rx != 1 || ry != 1) {
        closeSingleGroupLine (&y);
        freeSingleGroupLine (&y);
        sprintf (MsgText,
                "BIAS image and input are not binned to the same pixel size!");
        trlerror (MsgText);
        return (status = SIZE_MISMATCH);
    }

    /* Subtract the bias image from x. */

    /* If the science image is binned, it will be assumed to have
       same size as the reference image, since reading subarrays
       of a binned chip is not supported in current flight software.
     */
    if (same_size) {

        /* Loop over all the lines in the science image */
        for (i=0; i < scilines; i++) {
            status = getSingleGroupLine (wf3->bias.name, i, &y);
            if (status) {
                sprintf(MsgText,"Could not read line %d from bias image.",
                        i+1);
                trlerror(MsgText);
            }

            /* No trimming required. */
            status = sub1d(x, i, &y);
            if (status) {
                trlerror ("(biascorr) size mismatch.");
                return (status);
            }
        }

    } else {

        /* Loop over all the lines in the science array, and
         ** match them to the appropriate line in the reference image. */
        /* 
           i - index for line in science image
           j - index for line in reference image
           y0 - line in reference image corresponding to
           line in input image
         */
        initSingleGroupLine (&z);
        allocSingleGroupLine (&z, x->sci.data.nx);
        for (i=0, j=y0; i < scilines; i++,j++) { 

            /* We are working with a sub-array and need to apply the
             ** proper section from the reference image to the science
             ** image.  */
            status = getSingleGroupLine (wf3->bias.name, j, &y);
            if (status) {
                sprintf (MsgText,"Could not read line %d from bias image.",
                        j+1);
                trlerror(MsgText);
            }			

            update = NO;

            if (trim1d (&y, x0, y0, rx, avg, update, &z)) {
                trlerror ("(biascorr) size mismatch.");
                return (status);
            }

            status = sub1d (x, i, &z);
            if (status)
                return (status);
        }
        freeSingleGroupLine (&z);			/* done with z */
    }

    closeSingleGroupLine (&y);
    freeSingleGroupLine (&y);

    return (status);
}
