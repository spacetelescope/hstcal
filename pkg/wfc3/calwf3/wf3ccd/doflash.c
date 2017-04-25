/* This file contains:
	doFlash

*/

# include <stdio.h>
# include <stdlib.h>		/* calloc */
# include <math.h>		/* fabs */

# include "hstio.h"
# include "wf3.h"
# include "wf3info.h"
# include "hstcalerr.h"

/* This routine subtracts the post-flash image from x (in-place).
   For CCD data, the post-flash image is multiplied by the postflash exposure time and
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
   Howard Bushouse, 2001 May 7:
	Initial WFC3 Version.
   H.Bushouse, 2001 Nov 16:
	Updates to track CALACS changes - Revised to scale FLSHFILE by 
	FLASHDUR instead of EXPTIME. This also required modification to
	calling sequence to 'multgn1d'.
   H.Bushouse, 2001 Dec 4:
	Added error return if reference image not binned same as the
	science image.

   M. Sosey 2013 Septamber 09:
    Updated to correctly deal with both canned and user specified sub arrays
    adding 2 new functions to the multk1d and sub1d to deal with the virtual 
    overscan and applying the correct gain value for subarrays
    
*/

int doFlash (WF3Info *wf3ccd, SingleGroup *x, float *meanflash) {

/* arguments:
WF3Info *wf3     	i: calibration switches, etc
SingleGroup *x     io: image to be calibrated; written to in-place
float *meanflash    o: mean of post-flash image values subtracted

    Logic used in other parts of the pipeline

        IF AMPY > 0., then more than one AMP is used in the Y direction
        If AMPX > 0., then more than one AMP is used in the X direction
        
        So, if line being processed > AMPY then,  
            For all pixels up to AMPX, use AMP_C values 
                If AMPX is ZERO, don't use AMP_C values for any pixel
            For remaining pixels in line (or all), Use AMP_D values 
            
        However, if AMPY is ZERO or the line number > AMPY (and AMPY > 0) then,
            For all pixels up to AMPX, Use AMP_A values
                If AMPX is ZERO, don't use AMP_A values for any pixel
            Then, for remaining pixels in line (or all of them), 
                    Use AMP_B values
                    
        AMPX=AMPY=0 for all subarrays since only 1 amp is allowed according to the above logic
        which messes up the scaling call to multk1d
*/

    extern int status;

    SingleGroupLine y, z;	/* y and z are scratch space */
    int extver = 1;		/* get this imset from post-flash image */
    int rx, ry;			/* for binning post-flash down to size of x */
    int x0, y0;			/* offsets of sci image relative to reference image */
    int same_size;		/* true if no binning of ref image required */
    int avg = 0;		/* bin2d should sum values within each bin */
    int scilines; 		/* number of lines in science image */
    int i, j;
    float mean, flash;
    float weight, wflash;	/* weights for line averages */
    int update;
    float gain[NAMPS];
    float rn2[NAMPS];		/* only need this to call get_nsegn */
    int ampx;			/* border column for 2amp readout regions, set to size of image in ccdtab */
    int ampy;			/* Boundary values corrected for trim regions, set to size of image in ccdtab */
    int dimx, dimy;     /*dimensions of science image */
    int offsetx, offsety;

    int FindLine (SingleGroup *, SingleGroupLine *, int *, int *, int *, int *,
		  int *);
    int sub1d (SingleGroup *, int, SingleGroupLine *);
    int sub1dreform (SingleGroup *, int, int, SingleGroupLine *);
    int trim1d (SingleGroupLine *, int, int, int, int, int, SingleGroupLine *);
    int DetCCDChip (char *, int, int, int *);
    void get_nsegn (int, int, int, int, float *, float*, float *, float *);
    void AvgSciValLine (SingleGroupLine *, short, float *, float *);
    void multgn1d (SingleGroupLine *, int, int, int, float *, float);
    void multgn1dsub(SingleGroupLine *a, int , float *, float , char *);
    
    int streq_ic (char *, char *);
    int subarray;
    int straddle; /* the subarray starts in A or C and straddles the virtual overscan in the reference image */
    int overstart; /*where the overscan starts in the cut science image */

    /*init variables*/
    offsetx=0;
    offsety=0;
    ampx=0;
    ampy=0;
    dimx=0;
    dimy=0;
    mean=0.;
    flash=0.;
	weight = 0.0;
    wflash=0.;
    x0=0.;
    y0=0.;
    subarray=0;
    straddle=0;
    overstart=0;

    /*This will be set down below */
	for (i = 0; i < NAMPS; i++) {
	     gain[i] = 0.;
	     rn2[i] = 0.;
	}
 
	/* Check to see whether we need to do any processing at all...  */
	if (wf3ccd->flashdur <= 0.) {
	    sprintf(MsgText,
		  "Post-flash exposure was 0 seconds. FLSHCORR not performed.");
	    trlwarn(MsgText);
	    addHistoryKw (x->globalhdr, MsgText);
	    wf3ccd->flashcorr = IGNORED;
	
	    /* This is not an error condition, so continue with the remainder
	    ** of the calibrations... */
	    return (status);
	} 
	
	/* Flag an aborted Post-Flash exposure in the trailer file comments. */
	if (streq_ic(wf3ccd->flashstatus,"ABORTED")){
	    sprintf (MsgText,
	       "Post-flash STATUS was ABORTED. Post-flash may be compromised.");
	    trlwarn (MsgText);
	    /* Add this message to the image header as well... */
	    addHistoryKw (x->globalhdr, MsgText);
	}
    
    /*return an array of valid gain and readnoise values 
      This returns 2 amps regardless of subarray */
	get_nsegn (wf3ccd->detector, wf3ccd->chip, wf3ccd->ampx, wf3ccd->ampy,
		   wf3ccd->atodgain, wf3ccd->readnoise, gain, rn2);
    
    if (wf3ccd->verbose){
     sprintf(MsgText,"**gain,flashdur** = ([%f,%f,%f,%f],%f)",gain[0],gain[1],gain[2],gain[3],wf3ccd->flashdur);
     trlmessage(MsgText);
    }        
    
	
	/* Start with the actual post-flash subtraction now... */
	
	scilines = x->sci.data.ny;
    
	/* Compute correct extension version number to extract from
	** reference image to correspond to CHIP in science data.  */
	if (DetCCDChip (wf3ccd->flash.name, wf3ccd->chip, wf3ccd->nimsets, &extver) )
	    return (status);	
	
	if (wf3ccd->verbose) {
	    sprintf (MsgText,
		     "Performing post-flash subtraction on chip %d in imset %d",
		     wf3ccd->chip, extver);
	    trlmessage(MsgText);
	}

	/* Get the post-flash image data. */
	initSingleGroupLine (&y);
	openSingleGroupLine (wf3ccd->flash.name, extver, &y);

	if (hstio_err())
	    return (status = OPEN_FAILED);

	/* Compare binning of science image and reference image;
	   get same_size and high_res flags, and get info about
	   binning and offset for use by bin2d.
       
       x0 and y0 are the starting location of the subarrray in the reference image
	*/
	if (FindLine (x, &y, &same_size, &rx, &ry, &x0, &y0))
	    return (status);
    
	/* Return with error if reference image not binned same as input */
	if (rx != 1 || ry != 1) {
	    closeSingleGroupLine (&y);
	    freeSingleGroupLine (&y);
	    sprintf (MsgText,
	        "FLASH image and input are not binned to the same pixel size!");
	    trlerror (MsgText);
	    return (status = SIZE_MISMATCH);
	}

	if (wf3ccd->verbose){
	    sprintf(MsgText,"Image has starting location of %d,%d in the reference image",x0,y0);
	    trlmessage(MsgText);
	}

	/* AMPX,AMPY initialization, needed for the gain+flashdur step
       Different gains may be applied depending on what amp region the science data is in
    */
	dimx = x->sci.data.nx;
	dimy = x->sci.data.ny;
    

    /*For subarrays which START in B or D we can just move the x0 over 60 pixels to avoid the overscan in the reference image,
      which is always the full frame 4 amp image 

    Otherwise, only part of the subarray overlaps the overscan and special measures must be taken to avoid it. 
    
    ONLY 1 AMP is used to read subarrays, ever, so ampx and ampy should be set to the size of the image for all cases
    
    2073=2048+25 for physical overscan in reference image
    
    */
     if (!same_size){ 
        subarray=1;
        if (wf3ccd->verbose){
	        sprintf(MsgText,"SUBARRAY FOUND, amp=%s",wf3ccd->ccdamp);
            trlmessage(MsgText);
        }
        
        /*now we need to figure out where the subarray starts and if it straddles the virtual overscan in the middle
        of the reference postflash image 
        
        Since this is a subarray and only 1 amp is ever used, the value of ampx and ampy is set to the size of the image
        */
        
                
        wf3ccd->ampx=dimx;
        wf3ccd->ampy=dimy;
             
             
        if (x0 > 2072){ /*image starts in B or D regions and we can just shift the starting pixel*/
            if (wf3ccd->verbose){
	            sprintf(MsgText,"Subarray starts in B or D region, moved from (%d,%d) to ",x0,y0);
                trlmessage(MsgText);
	        }            
                x0 += 60;
            if (wf3ccd->verbose){
	            sprintf(MsgText,"(%d,%d) to avoid virtual overscan in reference",x0,y0);
                trlmessage(MsgText);
	        }            
        } else { /*the subarray starts somewhere in A or C and might straddle the virtual overscan region */
         
            if ( (x0 + dimx) > 2072){
                straddle=1;
                overstart=2073-x0;
            }
            
        }
    }
	
    if (wf3ccd->verbose){
	    sprintf(MsgText,"ccdamp=%s, straddle=%d, offset=(%d,%d),ampx,ampy=(%d,%d),x0,y0=(%d,%d)",wf3ccd->ccdamp,straddle,offsetx,offsety,wf3ccd->ampx,wf3ccd->ampy,x0,y0);
        trlmessage(MsgText);
    }
    
	initSingleGroupLine (&z);
    if (straddle){
    	allocSingleGroupLine (&z, x->sci.data.nx+60);
    } else {
        allocSingleGroupLine (&z, x->sci.data.nx);
    }
    
    
    for (i=0, j=y0; i < scilines; i++,j++) { 

	/* IF we are working with a sub-array and need to apply the
		proper section from the reference image to the science image.
        
        y is the reference image
        z is the reference image after it has been trimmed
        x is the science image
	*/
	    getSingleGroupLine (wf3ccd->flash.name, j, &y);

	    update = NO;
        
	    if (trim1d (&y, x0, j, rx, avg, update, &z)) {
			trlerror ("(flshcorr)reference file size mismatch.");
			return (status);
	    }
        
        if(subarray){
            multgn1dsub (&z, j, gain, wf3ccd->flashdur, wf3ccd->ccdamp);
        } else {
    	    multgn1d (&z, j, wf3ccd->ampx, wf3ccd->ampy, gain, wf3ccd->flashdur);
        }    

        /*currently this also uses the overscan pixels in the average because it's not flagged as bad
         But at this point, the reference image has already been trimmed to the science image, so
         it makes computing where the relative overscan pixels are that much harder. They are not flagged
         as bad in the DQ image, which is what the routine depends on to throw out bad pixels */
         
	    AvgSciValLine (&z, wf3ccd->sdqflags, &flash, &wflash);


		/* Sum the contribution from each line */			
		mean += flash * wflash;
		weight += wflash;

        if (straddle) {
    	    status = sub1dreform (x, i, overstart,&z);
	    } else {
            status = sub1d (x, i, &z);
        }
        if (status)
		return (status);
	}
    
	freeSingleGroupLine (&z);			/* done with z */

	closeSingleGroupLine (&y);
	freeSingleGroupLine (&y);

	/* Compute the mean for the entire image for save to the header*/	
	if (scilines > 0) 
	    *meanflash = mean / weight; 
	else 
	    *meanflash = 0.;
	
	return (status);
}
