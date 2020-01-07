# include <stdio.h>

#include "hstcal.h"
# include "hstio.h"
# include "acs.h"
# include "acsinfo.h"
# include "hstcalerr.h"

/* This routine applies the shutter shading correction:

	corrected = uncalibrated * EXPTIME / (EXPTIME + SHADFILE)

   The value of EXPTIME will (locally) be divided by the number of
   exposures that were combined by cosmic-ray rejection to make the
   current image.
   
    Warren J. Hack, 1999 Nov 10:
        Corrected problem with handling sub-arrays.  It now expands a 
        shadfile section to full size always, and then pulls out appropriate
        sub-section.
*/

int doShad (ACSInfo *acs2d, int extver, SingleGroup *x) {

/* arguments:
ACSInfo *acs2d     i: calibration switches, etc
SingleGroup *x    io: image to be calibrated; written to in-place
*/

	extern int status;

	SingleGroupLine y, z, zl;	/* z and zl are scratch space */
	ACSsect ysect, zsect;	/* scratch space to use for expansion */
	float exptime;		/* exposure time / ncombine */
	int rx, ry;		/* for binning dark down to size of x */
	int x0, y0;		/* offsets of sci image */
	int same_size;		/* true if no binning of ref image required */
	int avg = 1;		/* bin2d should average within each bin */
	int i,j, zline;
	int chipext;
	int update = NO;
    int zsecty, zsectx;

	int FindLine (SingleGroup *, SingleGroupLine *, int *, 
					int *, int *, int *, int *);
	int addk1d (SingleGroupLine *, float);
	int multk1d (SingleGroupLine *, float);
	int div1d (SingleGroup *, int, SingleGroupLine *);
	int trim1d (SingleGroupLine *, int, int, int, int, int, SingleGroupLine *);
	int allocACSsect (ACSsect *, int, int);
	void initACSsect (ACSsect *);
	void freeACSsect (ACSsect *);
	void copySectLine (ACSsect *, int, SingleGroupLine *);
	void getACSsect (char *, SingleGroupLine *, int, int, ACSsect *);
	int unbinsect (ACSsect *, int, ACSsect *);
	int DetCCDChip (char *, int, int, int *);
	
	if (acs2d->shadcorr != PERFORM)
	    return (status);

	if (acs2d->exptime <= 0.) {
	    trlwarn ("EXPTIME must be positive if SHADCORR = PERFORM.");
	    acs2d->shadcorr = IGNORED;
	    return (status);
	}
	
	initSingleGroupLine (&y);
		
	/* Correct for the number of exposures that were combined. */
	exptime = acs2d->exptime / acs2d->ncombine;

	/* Compute correct extension version number to extract from
		reference image to correspond to CHIP in science data.
	*/
	chipext = extver;
	if (DetCCDChip (acs2d->shad.name, acs2d->chip, acs2d->nimsets, &chipext) )
		return (status);	

	/* Get the shutter shading image data. */
	openSingleGroupLine (acs2d->shad.name, chipext, &y);
	if (hstio_err())
	    return (status = OPEN_FAILED);


	/* Compare binning of science image and reference image;
	   get the same_size flag, and get info about binning and offset
	   for use by bin2d.
	*/
	if (FindLine (x, &y, &same_size, &rx, &ry, &x0, &y0))
	    return (status);

	/* Apply the shutter shading image to the input image. */

	if (rx == 1. && ry == 1) {
		/* Loop over all the lines in the science image.
            The shading image may be binned the same, but the SCI
            image may be a sub-array, so it will need to be trimmed.
        */
        
	    initSingleGroupLine (&z);
	    allocSingleGroupLine (&z, x->sci.data.nx);
		for (i=0,j=y0; i < x->sci.data.ny; i++,j++) {
		
			getSingleGroupLine (acs2d->shad.name, j, &y);
				
	    	if (trim1d (&y, x0, y0, rx, avg, update, &z)) {
				trlerror ("(doShad) size mismatch.");
				return (status);
	    	}

	    	if (multk1d (&z, 1. / exptime))
				return (status);
	    	if (addk1d (&z, 1.))
				return (status);
	    	if (div1d (x, i, &z))
				return (status);
		
		}
        freeSingleGroupLine (&z);			/* done with z */
	} else {
	/* 
		Shading correction image is binned down more than image, and
			needs to be expanded to match the image size
	*/
		/* Set up scratch spaces to expand shadfile		*/
		/* Input file section */
		initACSsect (&ysect);
		if (allocACSsect (&ysect, y.sci.tot_nx, SECTLINES) ) {
			trlerror ("(doShad) Out of memory.");
			return (status = OUT_OF_MEMORY);
		}
		/* Output (expanded) file section - fully expanded section
            based on LTV/LTM values.  Sub-sections are taken out later.
        */
		initACSsect (&zsect);
        
        zsecty = SECTLINES * ry;
        zsectx = y.sci.tot_nx * rx;
		if (allocACSsect (&zsect, zsectx, zsecty) ){
			trlerror ("(doShad) Out of memory. ");
			return (status = OUT_OF_MEMORY);
		}
		
		/* Set up individual shadfile line to be applied to image */
		initSingleGroupLine (&zl);
	    allocSingleGroupLine (&zl, zsectx);
		initSingleGroupLine (&z);
	    allocSingleGroupLine (&z, x->sci.data.nx);
				
		if (acs2d->verbose) {
		    sprintf(MsgText, "Shad file will be expanded to %d pixels.", y.sci.tot_nx*rx);
		    trlmessage(MsgText);
		}

		/* Initialize row counter for input image 	*/
		i = 0;
		
		/* Initialize row counter for shadfile, based on offset
			calculated by FindLine.					*/
		j = y0;
		
		/* Read in shadfile and apply to input image here */
		while (i < x->sci.data.ny) {
			
			getACSsect (acs2d->shad.name, &y, j, SECTLINES, &ysect);
 
			/* Increment row counter for reference image */
			j += SECTLINES;			

			/* Expand binned reference data */
			unbinsect (&ysect, update, &zsect);
			
			/* For each line in expanded section, copy out a 
				SingleGroupLine and apply it to the science image. */
			for (zline=0; zline < zsecty; zline++) {
				
				/* Copy out individual expanded lines from reference data */
				copySectLine (&zsect, zline, &zl); 
				
	    		if (trim1d (&zl, x0, y0, 1, avg, update, &z)) {
					trlerror ("(doShad) size mismatch.");
					return (status);
	    		}
				
	    		if (multk1d (&z, 1. / exptime))
					return (status);
	    		if (addk1d (&z, 1.))
					return (status);
	    		if (div1d (x, i, &z))
					return (status);
					
				/* Increment row counters for science image */
				i++;
                /* In case we are processing a sub-array SCI image... */
                if (i >= x->sci.data.ny) break;                
			}
		}
		
		/* Clean up scratch space... */
		freeSingleGroupLine (&z);		
		freeSingleGroupLine (&zl);
		freeACSsect (&zsect);
		freeACSsect (&ysect);
	}
	
	closeSingleGroupLine (&y);
	freeSingleGroupLine (&y);
	
	return (status);
}
