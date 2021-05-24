# include <stdio.h>

#include "hstcal.h"
# include "hstio.h"
# include "acs.h"
# include "acsinfo.h"
# include "hstcalerr.h"

/* This routine subtracts the bias image from x (in-place).
 For ACS science data, it will normally be the
 case that this correction will be performed before the data has
 been combined for cosmic-ray
 rejection.
 Note that the bias image is assumed to have already been scaled by
 the gain, and binned the same as the input image.
 
 In order to accomodate the larger WFC images, the bias images will be
 applied line-by-line to the input image.  
 
 Warren Hack, 1998 June 3:
 Initial version for ACS... 
 
 ** Revised to read Bias image in line-by-line.  All code written around
 **	SingleGroupLine routines are new to ACS...
 */

int doBias (ACSInfo *acs, SingleGroup *x) {
  
  /* arguments:
   ACSInfo *acs     i: calibration switches, etc
   SingleGroup *x    io: image to be calibrated; written to in-place
   */
  
	extern int status;
  
	SingleGroupLine y, z;	/* y and z are scratch space */
	int extver = 1;			/* get this imset from bias image */
	int rx, ry;					/* for binning bias image down to size of x */
	int x0, y0;				/* offsets of sci image */
	int same_size;			/* true if no binning of ref image required */
	int avg = 0;			/* bin2d should sum within each bin */
	int scilines; 			/* number of lines in science image */
	int i, j;
	int update;
  
	int FindLine (SingleGroup *, SingleGroupLine *, int *, int *,int *, int *, int *);
	int sub1d (SingleGroup *, int, SingleGroupLine *);
	int trim1d (SingleGroupLine *, int, int, int, int, int, SingleGroupLine *);
	int DetCCDChip (char *, int, int, int *);
  
	if (acs->biascorr != PERFORM)
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
   reference image to correspond to CHIP in science data.
   */
	if (DetCCDChip (acs->bias.name, acs->chip, acs->nimsets, &extver) )
		return (status);		
  
	/* Get the first line of bias image data. */
	openSingleGroupLine (acs->bias.name, extver, &y);
	if (hstio_err())
    return (status = OPEN_FAILED);
  
	/* 
   Reference image should already be selected to have the
   same binning factor as the science image.  All we need to
   make sure of is whether the science array is a sub-array of
   the bias image.  
   
   x0,y0 is the location of the start of the 
   subimage in the reference image.
   
   **   
   **	FindLine is a modified version of FindBin routine from CALSTIS.
   **
   */
	if (FindLine (x, &y, &same_size, &rx, &ry, &x0, &y0))
    return (status);
  
	if (acs->verbose) {
		sprintf(MsgText,"Ratio of (%d,%d) with offset =(%d,%d)",rx,ry,x0,y0);
		trlmessage(MsgText);
		if (same_size) {
			sprintf(MsgText,"BIAS image and input are the same size ");
		} else {
			sprintf(MsgText,"BIAS image and input are NOT the same size ");
		}
    trlmessage(MsgText);
	}
	
	/* Subtract the bias image from x. */
  
	/* If the science image is binned, it will be assumed to have
   same size as the reference image, since reading subarrays
   of a binned chip is not supported in current flight software.
   */
	if (same_size) {
    
		/* Loop over all the lines in the science image */
		
		for (i=0; i < scilines; i++) {
			status = getSingleGroupLine (acs->bias.name, i, &y);
			if (status) {
				sprintf(MsgText,"Could not read line %d from bias image.",i+1);
				trlerror(MsgText);
			}
			
      /* No binning required. */
			status = sub1d(x, i, &y);
      if (status) {
				trlerror ("(biascorr) size mismatch.");
				return (status);
      }
		}
	} else {
    
        /* Loop over all the lines in the science array, and
           match them to the appropriate line in the reference image... 
           i - index for line in science image
           j - index for line in reference image
           y0 - line in reference image corresponding to line in input image
        */
        initSingleGroupLine (&z);
        allocSingleGroupLine (&z, x->sci.data.nx);
        for (i=0, j=y0; i < scilines; i++,j++) { 
            /* We are working with a sub-array and need to apply the
               proper section from the reference image to the science image.
            */
            status = getSingleGroupLine (acs->bias.name, j, &y);
            if (status) {
                sprintf(MsgText,"Could not read line %d from bias image.",j+1);
                trlerror(MsgText);
            }			

            update = NO;
     
            /* rx = 1; */
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
