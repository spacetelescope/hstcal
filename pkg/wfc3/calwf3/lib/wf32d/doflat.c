# include <stdio.h>
# include <stdlib.h>		/* for calloc */
# include <string.h>		/* for strncmp */

#include "hstcal.h"
# include "hstio.h"

# include "wf3.h"
# include "wf3info.h"
# include "hstcalerr.h"

static int divFlat (SingleGroup *, char *, WF3Info *, int);

/* This routine divides x in-place by the flat fields.
   There are up to three flat fields.  Each flat field is loaded
   and divided into the science image line-by-line, one flat at
   a time.

   The low-order flat, which is usually stored at a lower resolution
   than the other flats, is rebinned to match the science image before
   being applied.

   Warren Hack, 1998 June 12:
   	Initial ACS version.
   Howard Bushouse, 2000 Aug 29:
	Initial WFC3 version.
   H. Bushouse, 2001 Nov 15:
	Updated to track changes in CALACS - Fixed a major problem with
	applying reference data to sub-arrays; the indexing for the loop
	in divFlat was made consistent with the code in doDark.
   H. Bushouse, 2002 June 17:
	Revised divFlat to divide the PFLTFILE by gain before dividing
	the data by the flat-field. This converts the science data to
	units of ELECTRONS from DN (in accordance with CALACS changes).
   H. Bushouse, 2008 Oct 9:
	Added 'applygain' switch to divFlat to turn on/off the gain
	correction so that the gain will only be used to correct one ref
	file and not both, otherwise the gain will be applied twice to
	the science data.
   H. Bushouse, 2009 Jan 27:
	Fixed bugs that were causing crashes when expanding binned L-flat
	data and added capability to apply an L-flat that is the same size as
	the science image. Also inserted a forced error return if the L-flat
	needs to be interpolated, until the interpolation methods are improved.
   H. Bushouse, 2009 Oct 21:
	Updated divFlat to use the mean gain of all amps when applying
	the gain correction, except in the case of grism images, which still
	use the amp-dependent gain values.
   H. Bushouse, 2009 Nov 24:
	Updated divFlat to use the mean gain for all images, including grism.
	(calwf3 v2.0)

*/

int doFlat (WF3Info *wf32d, int extver, SingleGroup *x) {

/* arguments:
WF3Info *wf32d     i: calibration switches, etc
int extver	   i: extension/imset to be flat-fielded
SingleGroup *x    io: image to be calibrated; written to in-place
*/

	extern int status;

	SingleGroupLine w, zl, ztrim;	/* scratch space */
	WF3sect lfsect, elfsect;
	int rx, ry;		/* for binning dark down to size of x */
	int x0, y0;		/* offsets of sci image */
	int same_size;		/* true if no binning of ref image required */
	int xline;		/* counter for science image lines */
	int chipext;		/* Reference file IMSET corresponding to
				** CCD chip id for science image */
	int lf, zline;
	Hdr phdr;
	int update = NO;	/* Flag to determine whether hdr info needs
				** to be updated*/
	int scilines;
	int applygain;		/* Flag to determine whether to apply the gain
				** to ref file */
    
	int FindLine (SingleGroup *, SingleGroupLine *, int *, int *, int *,
		      int *, int *);
	int div1d (SingleGroup *, int, SingleGroupLine *);
	int allocWF3sect (WF3sect *, int, int);
	void initWF3sect (WF3sect *);
	void freeWF3sect (WF3sect *);
	void closeWF3sect (WF3sect *);
	void copySectLine (WF3sect *, int, SingleGroupLine *);
	void getWF3sect (char *, SingleGroupLine *, int, int, WF3sect *);
	int unbinsect (WF3sect *, int, WF3sect *);
	int DetCCDChip (char *, int, int, int *);
	
	/* Initialize applygain so that correction gets applied */
	applygain = 1;

	/* Apply pixel-to-pixel flat, if set to PERFORM */
	if (wf32d->pfltcorr == PERFORM) {
	    if (divFlat (x, wf32d->pflt.name, wf32d, applygain)) {
		sprintf (MsgText, "Problem applying PFLTFILE %s... ",
			 wf32d->pflt.name);
		trlerror(MsgText);
		return(status);
	    }
	    /* Turn off applygain so that it doesn't get applied again */
	    applygain = 0;
	}

	/* Apply delta flat, if set to PERFORM */
	if (wf32d->dfltcorr == PERFORM) {
	    if (divFlat (x, wf32d->dflt.name, wf32d, applygain)) {
		sprintf (MsgText, "Problem applying DFLTFILE %s... ",
			 wf32d->dflt.name);
		trlerror(MsgText);
		return (status);
	    }
	    /* Turn off applygain so that it doesn't get applied again */
	    applygain = 0;
	}
	
	initHdr (&phdr);

	/* low-order flat */
	if (wf32d->lfltcorr == PERFORM) {
	    if (wf32d->verbose) {
		sprintf(MsgText, "Reading in LFLAT file...");
		trlmessage (MsgText);
	    }

	    /* Compute correct extension version number to extract from
	    ** reference image to correspond to CHIP in science data.  */
	    chipext = extver;
	    if (DetCCDChip (wf32d->lflt.name, wf32d->chip, wf32d->nimsets,
			    &chipext))
		return (status);

	    /* Get the low-order flat field image data. */
	    initSingleGroupLine (&w);
	    openSingleGroupLine (wf32d->lflt.name, chipext, &w);
	    if (hstio_err())
    		return (status = OPEN_FAILED);

	    /* Compare binning of science image and reference image;
	    ** get the same_size flag, and get info about binning and offset
	    ** for use by bin2d.  */
	    if (FindLine (x, &w, &same_size, &rx, &ry, &x0, &y0))
		return (status);

	    /* If the L-flat is the same size as the science data, then just do
	    ** a straight division of the two images. */
	    if (same_size) {
		if (divFlat (x, wf32d->lflt.name, wf32d, applygain)) {
		    sprintf (MsgText, "Problem applying LFLTFILE %s... ",
			     wf32d->lflt.name);
		    trlerror(MsgText);
		    return (status);
		}

	    /* If the L-flat is binned, then we have to interpolate it to
	    ** match the science image before doing the division. */
	    } else {

	    /* Now we know what IMSET corresponds to the appropriate chip,
	    ** we can read in the flat-field file, expand it, combine it
	    ** with the other flats, then apply it to the science data... */

	    /* The current L-flat interpolation methods are known to be
	    ** non-optimal and we're not going to allow them to be used until
	    ** they're improved. So put in a temporary stub here that forces
	    ** a return with a warning if the L-flat is binned. */
	    sprintf (MsgText, "LFLTFILE %s size does not match science data.",
		     wf32d->lflt.name);
	    trlerror (MsgText);
	    sprintf (MsgText, 
	     "LFLTFILE interpolation methods are not available at this time.");
	    trlerror (MsgText);
	    sprintf (MsgText,
	     "Please use an LFLTFILE that matches size of science image.");
	    trlerror (MsgText);
	    closeSingleGroupLine (&w);
	    freeSingleGroupLine (&w);
	    return (status = SIZE_MISMATCH);

	    /* Low-order flat is binned down more than science image, and
	    ** needs to be expanded to match the image size */
	    /* Set up scratch spaces to expand low-order flat */
	    /* Input file section */
	    initWF3sect (&lfsect);
	    if (allocWF3sect (&lfsect, w.sci.tot_nx, SECTLINES) ) {
		trlerror ("(doFlat) Out of memory. ");
		return (status = OUT_OF_MEMORY);
	    }
	    /* Output (expanded) file section.  */
	    initWF3sect (&elfsect);
	    if (allocWF3sect (&elfsect, w.sci.tot_nx * rx, SECTLINES * ry) ){
		trlerror ("(doFlat) Out of memory. ");
		return (status = OUT_OF_MEMORY);
	    }

	    /* Set up individual lfltfile line to be applied to image */
	    initSingleGroupLine (&zl);
	    allocSingleGroupLine (&zl, w.sci.tot_nx * rx);

	    initSingleGroupLine (&ztrim);
	    allocSingleGroupLine (&ztrim, x->sci.data.nx);

	    xline = 0;	lf = 0;
	    scilines = x->sci.data.ny;

	    /* Loop over the science image lines, reading in and expanding the
	    ** appropriate L-flat image lines to apply to it. */
	    while (xline < scilines) {

		   /* Read in lfltfile and apply to input image here */
                   if (scilines-xline <= ry) lf=lf-1;
		   getWF3sect (wf32d->lflt.name, &w, lf, SECTLINES, &lfsect);

		   /* Increment row counter for reference image */
		   lf += 1;			

		   /* Expand binned reference data accounting for
		   ** offsets of subarrays */
		   unbinsect (&lfsect, update, &elfsect);

		   /* For each line in expanded section, copy out a 
		   SingleGroupLine and apply it to the science image. */
		   for (zline=0; (zline < ry) && (xline<scilines); zline++) {

			/* Copy out individual averaged, expanded lines from 
        		** reference data. */
			if (scilines-xline <= ry) {
			    copySectLine (&elfsect, zline+ry, &zl); 
			} else {
			    copySectLine (&elfsect, zline, &zl); 
			}

			/* We now have 1 expanded low-order flat line to apply 
			** Let's check to see if we have any other
			** flat-fields to apply...  */

			/* Now, apply flat-field */
			div1d (x, xline, &zl);
			xline++;

		   } /* End loop over expanded lflt lines, zline loop */
	    } /* End loop over input image lines, xline loop */

	    /* Clean up scratch areas that were used... */
	    freeWF3sect (&lfsect);
	    freeWF3sect (&elfsect);	
	    freeSingleGroupLine (&zl);
	    freeSingleGroupLine (&ztrim);
	    }
	    closeSingleGroupLine (&w);
	    freeSingleGroupLine (&w);
						
	} /* End if (lfltcorr) */

	return (status);
}

static int divFlat (SingleGroup *x, char *flatname, WF3Info *wf32d,
		    int applygain) {

	extern int status;

	int pchipext;
	SingleGroupLine y, ytrim; /* scratch space */
	int i, line;		  /* counters for science and ref image lines */
	int y_rx, y_ry;		  /* for binning dark down to size of x */
	int y_x0, y_y0;		  /* offsets of sci image */
	int ysame_size;		/* true if no binning of ref image required */
	int avg = 1;		/* bin2d should average within each bin */
	int update = NO;	/* Flag to determine whether hdr info needs
				** to be updated*/

	int scilines;           /* Number of lines in 'x' */
	float gain[NAMPS];
	float rn2[NAMPS];	/* only need this to call get_nsegn */
	float gnscale;
    
	int FindLine (SingleGroup *, SingleGroupLine *, int *, int *, int *,
		      int *, int *);
	int DetCCDChip (char *, int, int, int *);
	int trim1d (SingleGroupLine *, int, int, int, int, int,
		    SingleGroupLine *);
	int div1d (SingleGroup *, int, SingleGroupLine *);
	void get_nsegn (int, int, int, int, float *, float *, float *, float *);
	void multgn1d (SingleGroupLine *, int, int, int, float *, float);

	initSingleGroupLine (&y);

	/* Compute correct extension version number to extract from
	** reference image to correspond to CHIP in science data.  */
	if (DetCCDChip (flatname, wf32d->chip, wf32d->nimsets, &pchipext) )
	    return (status);	

	openSingleGroupLine (flatname, pchipext, &y);
	if (hstio_err())
	    return (status = OPEN_FAILED);

	if (FindLine (x, &y, &ysame_size, &y_rx, &y_ry, &y_x0, &y_y0))
	    return (status);
	if (wf32d->verbose) {
	    sprintf (MsgText,"ratio of flat/input = %d,%d offset by %d,%d",
		     y_rx,y_ry,y_x0,y_y0);
	    trlmessage(MsgText);
	}

	/* Return with error if reference data not binned same as input */
	if (y_rx != 1 || y_ry != 1) {
	    closeSingleGroupLine (&y);
	    freeSingleGroupLine (&y);
	    sprintf (MsgText,
		"FLAT image and input are not binned to the same pixel size!");
	    trlerror (MsgText);
	    return (status = SIZE_MISMATCH);
	}

	/* Perform one pointer operation here, instead of one for each line
	** in the image...  */
	scilines = x->sci.data.ny;

	/* Divide the reference image by the calibrated gain value atodgain,
	** and divide the flat into x. */

	for (i = 0; i < NAMPS; i++) {
	     gain[i] = 0.;
	     rn2[i]  = 0.;
	}
	get_nsegn (wf32d->detector, wf32d->chip, wf32d->ampx, wf32d->ampy,
		   wf32d->atodgain, wf32d->readnoise, gain, rn2);
	gnscale = 1.0;

	/* Load the average gain value for use in all quadrants. */
	for (i=0; i<NAMPS; i++)
	     gain[i] = wf32d->mean_gain;

	/* For the sake of run-time speed, the loop over lines
	** is performed differently and separately depending on
	** whether it is the same size or not...  */
	if (ysame_size) {

	    for (line = 0; line < scilines; line++) {
		 getSingleGroupLine (flatname, line, &y);

		 /* Divide flat by gain, if requested */
		 if (applygain)
		   multgn1d (&y, line, wf32d->ampx, wf32d->ampy, gain, gnscale);

		 /* Now, apply flat-field */
		 div1d (x, line, &y);

	    } /* End loop over input image lines, xline loop */

	} else {    

	    /* We are working with a sub-array image and need to trim
	    ** down the flat field lines.
	    ** So, we need a buffer for these flat-field lines... */
	    initSingleGroupLine (&ytrim);
	    allocSingleGroupLine (&ytrim, x->sci.data.nx);

	    for (i=0,line=y_y0; i < scilines; i++,line++) {

		   getSingleGroupLine (flatname, line, &y);

		   /* Make sure it is the same length as science image */
		   trim1d (&y, y_x0, y_y0, y_rx, avg, update, &ytrim);

		   /* Divide flat by gain, if requested */
		   if (applygain) {
		       multgn1d (&ytrim, line, wf32d->ampx, wf32d->ampy, gain,
				 gnscale);
		   }

		   /* Now, apply flat-field */
		   div1d (x, i, &ytrim);

	    } /* End loop over input image lines, xline loop */

	    /* Clean up buffers... */
	    freeSingleGroupLine (&ytrim);    
	}
    
	closeSingleGroupLine (&y);
	freeSingleGroupLine (&y);

	return (status);

}

