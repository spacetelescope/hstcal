/* This file contains:
	doDark
*/

# include <stdio.h>
# include <stdlib.h>		/* calloc */
# include <math.h>		/* fabs */

#include "hstcal.h"
# include "hstio.h"

# include "wf3.h"
# include "wf3info.h"
# include "hstcalerr.h"

/* This routine subtracts the dark image from x (in-place).
   For CCD data, the dark image is multiplied by the exposure time and
   divided by the gain before subtracting.  The dark time is just the
   exposure time; it DOES NOT INCLUDE the idle time since the last
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
    
   Warren Hack, 1998 June 11:
   	Initial ACS Version.
   Howard Bushouse, 2000 Aug 28:
	Initial WFC3 Version.
   H. Bushouse, 2001 May 7:
        Updated to track ACS changes - moved AvgSciValLine and multgn1d
        routines to lib/multk1d.c module.
   H. Bushouse, 2001 Nov 15:
	Updated to track CALACS changes - revised calling sequence for
	'multgn1d' to be more general and accept arbitrary scaling factor
	instead of being hardwired to 'exptime'.
*/

int doDark (WF3Info *wf32d, SingleGroup *x, float *meandark) {

/* arguments:
WF3Info *wf3       i: calibration switches, etc
SingleGroup *x    io: image to be calibrated; written to in-place
float *meandark	   o: mean of dark image values subtracted
*/

    extern int status;

    SingleGroupLine y, z;	/* y and z are scratch space */
    int extver = 1;		/* get this imset from dark image */
    int rx, ry;			/* for binning dark down to size of x */
    int x0, y0;			/* offsets of sci image */
    int same_size;		/* true if no binning of ref image required */
    int avg = 0;		/* bin2d should sum values within each bin */
    int scilines; 		/* number of lines in science image */
    int i, j;
    float mean, dark;
    float weight, wdark;    	/* weights for line averages */
    int update;
    float gain[NAMPS];
    float rn2[NAMPS];       	/* only need this to call get_nsegn */

    int FindLine (SingleGroup *, SingleGroupLine *, int *, int *, int *,
		  int *, int *);
    int sub1d (SingleGroup *, int, SingleGroupLine *);
    int trim1d (SingleGroupLine *, int, int, int, int, int, SingleGroupLine *);
    int DetCCDChip (char *, int, int, int *);
    void get_nsegn (int, int, int, int, float *, float*, float *, float *);
    void AvgSciValLine (SingleGroupLine *, short, float *, float *);
    void multgn1d (SingleGroupLine *, int, int, int, float *, float);


	initSingleGroupLine (&y);
	
	scilines = x->sci.data.ny;

	/* Compute correct extension version number to extract from
	   reference image to correspond to CHIP in science data.  */
	if (DetCCDChip(wf32d->dark.name, wf32d->chip, wf32d->nimsets, &extver))
	    return (status);	
	
	if (wf32d->verbose) {
	    sprintf (MsgText,
		     "Performing dark subtraction on chip %d in imset %d",
		     wf32d->chip, extver);
	    trlmessage(MsgText);
	}

	/* Get the dark image data. */
	openSingleGroupLine (wf32d->dark.name, extver, &y);
	if (hstio_err())
	    return (status = OPEN_FAILED);

	/* Compare binning of science image and reference image;
	   get same_size flag, and get info about binning and offset
	   for use by bin2d.
	*/
	if (FindLine (x, &y, &same_size, &rx, &ry, &x0, &y0))
	    return (status);
    
	/* Return with error if reference data not binned same as input */
	if (rx != 1 || ry != 1) {
	    closeSingleGroupLine (&y);
	    freeSingleGroupLine (&y);
	    sprintf (MsgText,
	    "DARK image and input are not binned to the same pixel size!");
	    trlerror(MsgText);
	    return (status = SIZE_MISMATCH);
	}
	if (wf32d->verbose){
	    sprintf(MsgText,"Image has an offset of %d,%d",x0,y0);
	    trlmessage(MsgText);
	}

	mean = 0.0;
	weight = 0.0;
    
	/* Multiply the dark image by the exposure time and divide by the
	   atodgain (or just by exposure time for the MAMAs), and
	   subtract it from x.
	*/
    
	for (i = 0; i < NAMPS; i++) {
	     gain[i] = 0.;
	     rn2[i] = 0.;
	}
	get_nsegn (wf32d->detector, wf32d->chip, wf32d->ampx, wf32d->ampy,
		   wf32d->atodgain, wf32d->readnoise, gain, rn2);

	initSingleGroupLine (&z);
	allocSingleGroupLine (&z, x->sci.data.nx);
	for (i=0, j=y0; i < scilines; i++,j++) { 

	     /* We are working with a sub-array and need to apply the
		proper section from the reference image to the science image.
	     */
	     getSingleGroupLine (wf32d->dark.name, j, &y);

             update = NO;

	     if (trim1d (&y, x0, y0, rx, avg, update, &z)) {
		 trlerror ("(darkcorr) size mismatch.");
		 return (status);
	     }

	     multgn1d(&z, j, wf32d->ampx, wf32d->ampy, gain, wf32d->exptime[0]);

	     AvgSciValLine (&z, wf32d->sdqflags, &dark, &wdark);

	     /* Sum the contribution from each line */			
	     mean += dark * wdark;
	     weight += wdark;

	     status = sub1d (x, i, &z);
	     if (status)
		 return (status);
	}
	freeSingleGroupLine (&z);			/* done with z */

	closeSingleGroupLine (&y);
	freeSingleGroupLine (&y);

	/* Compute the mean for the entire image */	
	if (scilines > 0) 
	    *meandark = mean / weight; 
	else 
	    *meandark = 0.;
	
	return (status);
}

