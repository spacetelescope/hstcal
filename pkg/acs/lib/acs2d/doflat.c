# include <stdio.h>
# include <stdlib.h>		/* for calloc */
# include <time.h>

#include "hstcal.h"
# include "hstio.h"
# include "acs.h"
# include "acsinfo.h"
# include "hstcalerr.h"

static int divFlat (SingleGroup *, char *, ACSInfo *);

/* This routine divides x in-place by the flat fields.
 There are up to three flat fields.  They are read into SingleGroups,
 multiplied together (leaving the result each time in y), and the
 product of the flat fields is divided into the input x, with the
 final quotient left in x.
 
 Warren Hack, 1998 June 12:
 Initial ACS version.
 WJH, 2001 Oct 17: Fixed a major problem with applying reference
 data to sub-arrays; the indexing for the loop in divFlat was
 made consistent with the code in doDark.
 WJH, 2002 Jan 29: Revised divFlat to divide the PFLTFILE by gain
 before dividing the data by the flat-field.  This converts
 the science data to units of ELECTRONS from DN.
 WJH, 2008 Oct 9: Added 'applygain' switch to divFlat to turn on/off
 the gain correction so that the gain will only be used to correct
 one ref file and not both, otherwise the gain will be applied twice
 to the science data.
 */

int doFlat (ACSInfo *acs2d, int extver, SingleGroup *x) {
  
  /* arguments:
   ACSInfo *acs2d     i: calibration switches, etc
   int extver		   i: extension/imset to be flat-fielded
   SingleGroup *x    io: image to be calibrated; written to in-place
   */
  
	extern int status;
  
	SingleGroupLine w, zl, ztrim;	/* scratch space */
	ACSsect lfsect, elfsect;
	int rx, ry;		/* for binning dark down to size of x */
	int x0, y0;		/* offsets of sci image */
	int same_size;		/* true if no binning of ref image required */
	int xline;		/* counter for science image lines */
	int chipext;	/* Reference file IMSET corresponding to
                 CCD chip id for science image */
	int lf, zline;
	Hdr phdr;
	int update = NO;	/* Flag to determine whether hdr info needs to be updated*/
	int scilines;
  
	int FindLine (SingleGroup *, SingleGroupLine *, int *,
                int *, int *, int *, int *);
	int div1d (SingleGroup *, int, SingleGroupLine *);
	int allocACSsect (ACSsect *, int, int);
	void initACSsect (ACSsect *);
	void freeACSsect (ACSsect *);
	void closeACSsect (ACSsect *);
	void copySectLine (ACSsect *, int, SingleGroupLine *);
	void getACSsect (char *, SingleGroupLine *, int, int, ACSsect *);
	int unbinsect (ACSsect *, int, ACSsect *);
	int DetCCDChip (char *, int, int, int *);
  
	/* apply pixel-to-pixel flat, if set to PERFORM */
	if (acs2d->pfltcorr == PERFORM) {
    if (divFlat (x, acs2d->pflt.name, acs2d)){
      sprintf (MsgText, "Problem applying PFLTFILE %s... ", acs2d->pflt.name);
      trlerror(MsgText);
      return(status);
    }
	}
  
  
	/* apply delta flat, if set to PERFORM */
	if (acs2d->dfltcorr == PERFORM) {
    if (divFlat (x, acs2d->dflt.name, acs2d) ){
      sprintf (MsgText, "Problem applying DFLTFILE %s... ", acs2d->dflt.name);
      trlerror(MsgText);
      return (status);
    }
	}
  
	lf = 0;
	initHdr (&phdr);
  
	/* low-order flat */
	if (acs2d->lfltcorr == PERFORM) {
		if (acs2d->verbose) {
			sprintf(MsgText, "Reading in LFLAT file...");
			trlmessage (MsgText);
		}
		initSingleGroupLine (&w);
    
		/* Compute correct extension version number to extract from
     reference image to correspond to CHIP in science data.
     */
		chipext = extver;
		if (DetCCDChip (acs2d->lflt.name, acs2d->chip, acs2d->nimsets, &chipext) )
			return (status);
    
		/* Get the low-order flat field image data. */
		openSingleGroupLine (acs2d->lflt.name, chipext, &w);
		if (hstio_err())
      return (status = OPEN_FAILED);
    
		/* Compare binning of science image and reference image;
     get the same_size flag, and get info about binning and offset
     for use by bin2d.
     */
		if (FindLine (x, &w, &same_size, &rx, &ry, &x0, &y0))
      return (status);
    
		/* Now we know what IMSET corresponds to the appropriate chip,
     we can read in the flat-field file, expand it, combine it
     with the other flats, then apply it to the science data... */
    
		/* Low-order flat is binned down more than science image, and
     needs to be expanded to match the image size
     */
		/* Set up scratch spaces to expand low-order flat				*/
		/* Input file section */
		initACSsect (&lfsect);
		if (allocACSsect (&lfsect, w.sci.tot_nx, SECTLINES) ) {
			trlerror ("(doFlat) Out of memory. ");
			return (status = OUT_OF_MEMORY);
		}
		/* Output (expanded) file section.  */
		initACSsect (&elfsect);
		if (allocACSsect (&elfsect, w.sci.tot_nx * rx, SECTLINES * ry) ){
			trlerror ("(doFlat) Out of memory. ");
			return (status = OUT_OF_MEMORY);
		}
    
		/* Set up individual lfltfile line to be applied to image */
		initSingleGroupLine (&zl);
    allocSingleGroupLine (&zl, w.sci.tot_nx * rx);
    
		initSingleGroupLine (&ztrim);
    allocSingleGroupLine (&ztrim, x->sci.data.nx);
    
		xline = 0;
    scilines = x->sci.data.ny;
    
		while (xline <= scilines) {
      
			/* Read in lfltfile and apply to input image here */
			getACSsect (acs2d->lflt.name, &w, lf, SECTLINES, &lfsect);
      
			/* Increment row counter for reference image */
			lf += 1;
      
			/* Expand binned reference data
       accounting for offsets of subarrays */
      
			unbinsect (&lfsect, update, &elfsect);
      
			/* For each line in expanded section, copy out a
       SingleGroupLine and apply it to the science image. */
			for (zline=0; zline < ry; zline++) {
        
				/* Copy out individual averaged, expanded lines from
         reference data. */
				copySectLine (&elfsect, zline, &zl);
        
				/* We now have 1 expanded low-order flat line to apply
         Let's check to see if we have any other flat-fields to
         apply...
         */
        
				/* Now, apply flat-field */
				div1d (x, xline, &ztrim);
				xline++;
        
			} /* End loop over expanded lflt lines, zline loop */
      
		} /* End loop over input image lines, xline loop */
		/* Clean up scratch areas that were used... */
    closeACSsect (&lfsect);
    freeACSsect (&lfsect);
    /*closeACSsect (&elfsect); */
    freeACSsect (&elfsect);
    closeSingleGroupLine (&w);
    freeSingleGroupLine (&w);
    freeSingleGroupLine (&zl);
    freeSingleGroupLine (&ztrim);
    
    /* End if (lfltcorr) */
  }
  
	return (status);
}

static int divFlat (SingleGroup *x, char *flatname, ACSInfo *acs2d) {
  
  extern int status;
  
  int pchipext;
  SingleGroupLine y, ytrim;                  /* scratch space */
	int i,line;		            /* counters for science and ref image lines */
	int y_rx, y_ry;		        /* for binning dark down to size of x */
	int y_x0, y_y0;		        /* offsets of sci image */
	int ysame_size;		/* true if no binning of ref image required */
	int avg = 1;		/* bin2d should average within each bin */
	int update = NO;	/* Flag to determine whether hdr info needs to be updated*/
  
  int scilines;           /* Number of lines in 'x' */
  SingleGroup inspot,outspot;
  SingleGroupLine spotline, strim;
  float shiftx, shifty;
  time_t date;
  
	int FindLine (SingleGroup *, SingleGroupLine *, int *,
                int *, int *, int *, int *);
	int DetCCDChip (char *, int, int, int *);
	int trim1d (SingleGroupLine *, int, int, int, int, int, SingleGroupLine *);
	int div1d (SingleGroup *, int, SingleGroupLine *);
  void get_nsegn (int, int, int, int, float *, float*, float *, float *);
	void multgn1d (SingleGroupLine *, int, int, int, float *, float);
  
  int GetSpotTab(char *, time_t, float *, float *);
  int shiftSpot(SingleGroup *, float, float, SingleGroup *);
  int readSpotImage(char *, SingleGroup *, SingleGroupLine *);
  void copySpotLine(SingleGroup *, int, SingleGroupLine *);
  int parseObsDate (Hdr *, time_t *);
  int multlines (SingleGroupLine *, SingleGroupLine *);
  
	initSingleGroupLine (&y);
  
  /* Compute correct extension version number to extract from
   reference image to correspond to CHIP in science data.
   */
  if (DetCCDChip (flatname, acs2d->chip, acs2d->nimsets, &pchipext) )
    return (status);
  
  openSingleGroupLine (flatname, pchipext, &y);
  if (hstio_err())
    return (status = OPEN_FAILED);
  
  if (FindLine (x, &y, &ysame_size, &y_rx, &y_ry, &y_x0, &y_y0))
    return (status);
  if (acs2d->verbose){
    sprintf(MsgText,"ratio of flat/input = %d,%d offset by %d,%d",y_rx,y_ry,y_x0,y_y0);
    trlmessage(MsgText);
  }
  
  /* Perform one pointer operation here, instead of one for each line
   in the image...
   */
  scilines = x->sci.data.ny;
  
  /* Check for existence of SPOTFLAT file... */
	if (acs2d->cfltcorr == PERFORM) {
    
    status = readSpotImage(acs2d->cflt.name, &inspot, &spotline);
    
    /* Initialize shifted spot arrays */
    initSingleGroup(&outspot);
    allocSingleGroup(&outspot, inspot.sci.data.nx, inspot.sci.data.ny, True);
    
    parseObsDate(x->globalhdr, &date);
    status = GetSpotTab(acs2d->spot.name, date, &shiftx, &shifty);
    
    sprintf(MsgText,"SPOTTAB:  Using spot shift of: %0.2g  %0.2g",shiftx,shifty);
    trlmessage(MsgText);
    
    /* Shift input spot flat
     Multiply PFLT with this shifted spot, outspot, when
     applying it to the data (line-by-line)...
     */
    status = shiftSpot(&inspot, shiftx, shifty, &outspot);
    
  }
  
  /* For the sake of run-time speed, the loop over lines
   is performed differently and separately depending on
   whether it is the same size or not...
   */
  if (ysame_size) {
    for (line = 0; line < scilines; line++) {
      
      getSingleGroupLine (flatname, line, &y);
      
      /* Apply spot flat if one was specified... */
      if (acs2d->cfltcorr == PERFORM) {
        copySpotLine(&outspot, line, &spotline);
        multlines(&y, &spotline);
      }
      
      /* Now, apply flat-field */
      div1d (x, line, &y);
      
    } /* End loop over input image lines, xline loop */
    
  } else {
    /* We are working with a sub-array image and need to trim
     down the flat field lines
     
     So, we need a buffer for these flat-field lines...
     */
    initSingleGroupLine (&ytrim);
    allocSingleGroupLine (&ytrim, x->sci.data.nx);
    if (acs2d->cfltcorr == PERFORM) {
      initSingleGroupLine (&strim);
      allocSingleGroupLine (&strim, x->sci.data.nx);
    }
    
    for (i=0,line=y_y0; i < scilines; i++,line++) {
      
      getSingleGroupLine (flatname, line, &y);
      
      /* Make sure it is the same length as science image */
      trim1d (&y, y_x0, y_y0, y_rx, avg, update, &ytrim);
      
      if (acs2d->cfltcorr == PERFORM) {
        copySpotLine(&outspot, line, &spotline);
        trim1d(&spotline, y_x0, y_y0, y_rx, avg, update, &strim);
        multlines(&ytrim, &strim);
      }
      
      /* Now, apply flat-field */
      div1d (x, i, &ytrim);
      
    } /* End loop over input image lines, xline loop */
    
    /* Clean up buffers... */
    freeSingleGroupLine (&ytrim);
    if (acs2d->cfltcorr == PERFORM) {
      freeSingleGroupLine (&ytrim);
    }
  }
  
	if (acs2d->cfltcorr == PERFORM) {
    freeSingleGroup(&inspot);
    freeSingleGroupLine(&spotline);
    freeSingleGroup(&outspot);
  }
  
	closeSingleGroupLine (&y);
	freeSingleGroupLine (&y);
  
	return (status);
  
}
