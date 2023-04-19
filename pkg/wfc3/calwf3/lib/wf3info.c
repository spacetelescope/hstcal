# include <string.h>			/* for strncmp, strlen */

#include "hstcal.h"
# include "hstio.h"
# include "wf3.h"
# include "wf3info.h"
# include "wf3corr.h"		/* calibration switch names */
# include "hstcalerr.h"			/* defines error codes */

/* These functions initialize values for the WF3Info structure.
	The values set here are NOT default values, just initialization...
   Warren Hack, 1998 June 10:
   	Initial ACS version.
	
   Howard Bushouse, 2000 Aug 24:
	Initial WFC3 version.

   Howard Bushouse, 2001 Apr 17:
	Changed exptime from scalar to array for WF3 IR channel

   H.Bushouse, 2001 May 8:
	Added support for post-flash processing step;
	Revised FILTER1,FILTER2 usage to just FILTER for WFC3.

   H.Bushouse, 2001 Nov 16: (CALACS changes)
	Removed 'filtcorr', 'apertab', and 'phottab'. Changed to use of
	'graphtab' and 'comptab' for new PHOTCORR scheme. Upgraded to
	use individual bias values for each CCD amp.

   H.Bushouse, 2002 Jan 28: 
	Added biassectc, biassectd, and increased dimesions of trimx
	to accomodate extra serial virtual overscan regions in WFC3
	UVIS CCD's.

   H.Bushouse, 2002 Mar 1:
	Increased dimensions of vx and vy from 2 to 4 to allow for two
	separate parallel overscan regions to be used - one for each amp
	on the UVIS CCD.

   H.Bushouse, 2002 Jun 20:
	Added expend; removed use of statcorr (following CALACS changes).

   H.Bushouse, 2003 Oct 27:
	Added initialization of wf3->blev array (following CALACS).

   H.Bushouse, 2009 Jan 08:
	Added initialization of new 'type' parameter in InitRefImg and 
	InitRefTab. Also added new CheckImgType, CheckTabType, CheckFilter,
	CheckDetector, and CheckGain routines.

   H.Bushouse, 2009 Oct 21:
	Added initialization of new 'mean_gain' parameter in WF3Init.

   H.Bushouse, 2011 Jul 29:
	Fixed the logic in the CheckGain routine so that the ref image gets
	closed before returning when keyval=-1. (PR68983; Trac 745).

   H.Bushouse, 2011 Sep 7:
	Modified WF3Init to initialize new phot tab instead of graph and
	comp tabs.
    
  M. Sosey 2013 December 4
    Added fluxcorr switches for new uvis photometry scaling of chips
    
  M. sosey 2014
    Added pctecorr switches for new CTE correction

  M. Sosey 2015
    Added chip keyowrds to fully support all phot keys in all headers

  M. De La Pena 2020 March:
    Added reading of PCTERNOI keyword from the raw primary header for 
    possible use in the CTE reduction

  M. De La Pena 2022 February:
    Added variable, satmap - the reference image for full-well
    saturation. Use of this image rendered wf3->saturate variable
    obsolete.  Removed "wf3->saturate" as part of the cleanup operation.
*/

void WF3Init (WF3Info *wf3) {

	int i;

	void InitRefTab (RefTab *tab);
	void InitRefImg (RefImage *img);
	
	/* Assign default values. */

	wf3->input[0] = '\0';
	wf3->output[0] = '\0';
	wf3->rootname[0] = '\0';
	wf3->printtime = 0;
	wf3->verbose = 0;
	
	wf3->det[0] = '\0';
	wf3->aperture[0] = '\0';
	wf3->filter[0]  = '\0';
	wf3->obstype[0] = '\0';
	wf3->detector = UNKNOWN_DETECTOR;
	wf3->chip = 1;
    wf3->chip1_flam=1.;
    wf3->chip2_flam=1.;
	wf3->ncombine = 1;
	wf3->nimsets = 1;
	wf3->members = 1;
	wf3->mtype[0] = '\0';

	wf3->exptime[0] = 1.;
	wf3->expstart = 0.;
	wf3->expend = 0.;

	wf3->sdqflags = MAX_DQ;			/* 16 bits set */

	wf3->subarray = 0;
	wf3->bin[0] = 1;
	wf3->bin[1] = 1;
	wf3->offsetx = 0.;
	wf3->offsety = 0.;

	/* CCD-specific info */
	wf3->ccdamp[0] = '\0';
	wf3->ccdgain = 1;
	wf3->binaxis[0] = 1;
	wf3->binaxis[1] = 1;
	wf3->mean_gain = 0;
	for (i=0; i < NAMPS; i++) {
	     wf3->ccdoffset[i] = 0;
	     wf3->ccdbias[i] = 0.;
	     wf3->atodgain[i] = 0.;
	     wf3->readnoise[i] = 0.;
	     wf3->blev[i] = 0.;
	}
	wf3->ampx = 0;
	wf3->ampy = 0;
	wf3->trimx[0] = 0;
	wf3->trimx[1] = 0;
	wf3->trimx[2] = 0;
	wf3->trimx[3] = 0;
	wf3->trimy[0] = 0;
	wf3->trimy[1] = 0;
	wf3->vx[0] = 0;
	wf3->vx[1] = 0;
	wf3->vx[2] = 0;
	wf3->vx[3] = 0;
	wf3->vy[0] = 0;
	wf3->vy[1] = 0;
	wf3->vy[2] = 0;
	wf3->vy[3] = 0;
	wf3->biassecta[0] = 0;
	wf3->biassecta[1] = 0;
	wf3->biassectb[0] = 0;
	wf3->biassectb[1] = 0;
	wf3->biassectc[0] = 0;
	wf3->biassectc[1] = 0;
	wf3->biassectd[0] = 0;
	wf3->biassectd[1] = 0;
	wf3->flashdur = 0;
	wf3->flashstatus[0] = '\0';
    wf3->pcternoi = 0.;

	/* Initialize Calibration switches */
	wf3->dqicorr  = OMIT;
	wf3->atodcorr = OMIT;
	wf3->blevcorr = OMIT;
	wf3->biascorr = OMIT;
	wf3->noiscorr = OMIT;
	wf3->darkcorr = OMIT;
	wf3->nlincorr = OMIT;
	wf3->flatcorr = OMIT;
	wf3->flashcorr = OMIT;
    wf3->fluxcorr = OMIT;
	wf3->pfltcorr = OMIT;
	wf3->dfltcorr = OMIT;
	wf3->lfltcorr = OMIT;
	wf3->shadcorr = OMIT;
	wf3->photcorr = OMIT;
	wf3->zoffcorr = OMIT;
	wf3->zsigcorr = OMIT;
    wf3->pctecorr = OMIT;

	/* Initialize reference images and tables for WF3CCD */
	InitRefImg (&(wf3->bias));
    InitRefImg (&(wf3->sink));
    InitRefImg (&(wf3->biac));
	InitRefImg (&(wf3->flash));
	InitRefTab (&(wf3->bpix));
	InitRefTab (&(wf3->ccdpar));
	InitRefTab (&(wf3->oscn));
	InitRefTab (&(wf3->atod));
    InitRefTab (&(wf3->pctetab));
    InitRefImg (&(wf3->satmap));

	/* Initialize reference images and tables for WF32D */
	InitRefImg (&(wf3->dark));
	InitRefImg (&(wf3->pflt));
	InitRefImg (&(wf3->dflt));
	InitRefImg (&(wf3->lflt));
	InitRefImg (&(wf3->shad));
	InitRefTab (&(wf3->phot));
	
	/* Initialize reference images and tables for WF3IR */
	InitRefImg (&(wf3->nlin));
}

void InitRefTab (RefTab *tab) {

	tab->name[0] = '\0';
	tab->type[0] = '\0';
	tab->pedigree[0] = '\0';
	tab->descrip[0] = '\0';
	tab->descrip2[0] = '\0';
	tab->exists = EXISTS_UNKNOWN;
	tab->goodPedigree = PEDIGREE_UNKNOWN;

}


void InitRefImg (RefImage *img) {

	img->name[0] = '\0';
	img->type[0] = '\0';
	img->pedigree[0] = '\0';
	img->descrip[0] = '\0';
	img->exists = EXISTS_UNKNOWN;
	img->goodPedigree = PEDIGREE_UNKNOWN;

}


/* This routine gets the name of a reference image, checks whether it
   exists, and gets filetype, pedigree and descrip if they are present.  The 
   image name will be null if the keyword is not present in the header; this is
   not an error.
*/

int GetImageRef (RefFileInfo *refnames, Hdr *phdr,
		char *keyword, RefImage *image, int *calswitch) {

	extern int status;

	int GetRefName (RefFileInfo *, Hdr *, char *, char *);
	int ImgPedigree (RefImage *);

	/* Get the reference image name. */
	if (GetRefName (refnames, phdr, keyword, image->name))
	    return (status);

	/* ImgPedigree opens the image to verify that it exists, and if so,
	   gets filetype, pedigree & descrip.
	*/
	if (ImgPedigree (image))
	    return (status);
	if (image->exists == EXISTS_YES) {
	    if (image->goodPedigree != GOOD_PEDIGREE)
		*calswitch = DUMMY;
	}

	return (status);
}



/* This routine gets the name of a reference table, checks whether it
   exists, and gets filetype, pedigree and descrip if they are present.  The 
   table name will be null if the keyword is not present in the header; this is
   not an error.
*/

int GetTabRef (RefFileInfo *refnames, Hdr *phdr,
		char *keyword, RefTab *table, int *calswitch) {

	extern int status;

	int GetRefName (RefFileInfo *, Hdr *, char *, char *);
	int TabPedigree (RefTab *);

	/* Get the reference table name. */    
	if (GetRefName (refnames, phdr, keyword, table->name))
	    return (status);

	/* TabPedigree opens the table to verify that it exists, and if so,
	   gets filetype, pedigree & descrip.
	*/
	if (TabPedigree (table))
	    return (status);
	if (table->exists == EXISTS_YES) {
	    if (table->goodPedigree != GOOD_PEDIGREE)
		*calswitch = DUMMY;
	}

	return (status);
}


void MissingFile (char *keyword, char *filename, int *missing) {

	sprintf (MsgText, "%s `%s' not found or can't open.", keyword, filename);
	trlerror (MsgText);
	(*missing)++;
}


void CheckImgType (RefImage *image, char *filetype, char *keyword, int *badtype) {

	if (strncmp (image->type, filetype, strlen(filetype)) != 0) {
	    sprintf (MsgText, "%s FILETYPE='%s' is invalid for %s", 
		     image->name, image->type, keyword);
	    trlerror (MsgText);
	    (*badtype)++;
	}
}


void CheckTabType (RefTab *table, char *filetype, char *keyword, int *badtype) {

	if (strncmp (table->type, filetype, strlen(filetype)) != 0) {
	    sprintf (MsgText, "%s FILETYPE='%s' is invalid for %s", 
		     table->name, table->type, keyword);
	    trlerror (MsgText);
	    (*badtype)++;
	}
}


int CheckFilter (char *image, char *filter, char *keyword, int *badtype) {

	extern int status;

	FitsKw key;		/* location of keyword in header */
	IODescPtr im;		/* descriptor for primary header unit */
	Hdr phdr;		/* primary header */

	char keyval[SZ_FITS_REC+1];

	keyval[0] = '\0';
	initHdr (&phdr);

	/* Open the primary header of the reference file. */
	im = openInputImage (image, "", 0);
	getHeader (im, &phdr);
	if (hstio_err())
	    return (status = HEADER_PROBLEM);

	/* Get the FILTER keyword. */
	key = findKw (&phdr, keyword);
	if (key == NotFound) {
	    trlkwerr (keyword, image);
	    return (status = KEYWORD_MISSING);
	} else {
	    getStringKw (key, keyval, SZ_FITS_REC);
	    if (hstio_err()) {
		trlkwerr (keyword, image);
		return (status = KEYWORD_MISSING);
	    }
	}

	/* Does the ref file FILTER value match the science image?
	   Note: make a special allowance for a value of 'ANY'.
	*/
	if ((strncmp (keyval, "ANY", strlen(keyval)) != 0) &&
	    (strncmp (keyval, filter, strlen(keyval)) != 0)) {
	    sprintf (MsgText, "%s %s='%s' does not match science data",
		     image, keyword, keyval);
	    trlerror (MsgText);
	    (*badtype)++;
	}

	/* Close the reference file. */
	closeImage (im);
	freeHdr (&phdr);

	return (status);
}


int CheckDetector (char *image, int detector, char *keyword, int *badtype) {

	extern int status;

	FitsKw key;		/* location of keyword in header */
	IODescPtr im;		/* descriptor for primary header unit */
	Hdr phdr;		/* primary header */

	char keyval[SZ_FITS_REC+1];

	keyval[0] = '\0';
	initHdr (&phdr);

	/* Open the primary header of the reference file. */
	im = openInputImage (image, "", 0);
	getHeader (im, &phdr);
	if (hstio_err())
	    return (status = HEADER_PROBLEM);

	/* Get the DETECTOR keyword. */
	key = findKw (&phdr, keyword);
	if (key == NotFound) {
	    trlkwerr (keyword, image);
	    return (status = KEYWORD_MISSING);
	} else {
	    getStringKw (key, keyval, SZ_FITS_REC);
	    if (hstio_err()) {
		trlkwerr (keyword, image);
		return (status = KEYWORD_MISSING);
	    }
	}

	/* Does the ref file DETECTOR value match the science image? */
	if (detector == IR_DETECTOR) {
	    if (strncmp (keyval, "IR", strlen(keyval)) != 0) {
		sprintf (MsgText, "%s %s='%s' does not match science data",
			 image, keyword, keyval);
		trlerror (MsgText);
		(*badtype)++;
	    }
	} else {
	    if (strncmp (keyval, "UVIS", strlen(keyval)) != 0) {
		sprintf (MsgText, "%s %s='%s' does not match science data",
			 image, keyword, keyval);
		trlerror (MsgText);
		(*badtype)++;
	    }
	}

	/* Close the reference file. */
	closeImage (im);
	freeHdr (&phdr);

	return (status);
}


int CheckGain (char *image, float gain, char *keyword, int *badtype) {

	extern int status;

	FitsKw key;		/* location of keyword in header */
	IODescPtr im;   /* descriptor for primary header unit */
	Hdr phdr;		/* primary header */

	float keyval;

	initHdr (&phdr);

	/* Open the primary header of the reference file. */
	im = openInputImage (image, "", 0);
	getHeader (im, &phdr);
	if (hstio_err())
	    return (status = HEADER_PROBLEM);

	/* Get the CCDGAIN keyword. */
	key = findKw (&phdr, keyword);
	if (key == NotFound) {
	    trlkwerr (keyword, image);
	    return (status = KEYWORD_MISSING);
	} else {
	    keyval = getFloatKw (key);
	    if (hstio_err()) {
		trlkwerr (keyword, image);
		return (status = KEYWORD_MISSING);
	    }
	}

	/* Does the ref file CCDGAIN value match the science image? */
	/* A value of -1 is considered to be OK */
	if ((keyval != -1) && (gain != keyval)) {
	    sprintf (MsgText, "%s %s=%g does not match science data",
		     image, keyword, keyval);
	    trlerror (MsgText);
	    (*badtype)++;
	}

	/* Close the reference file. */
	closeImage (im);
	freeHdr (&phdr);

	return (status);
}
void initCCDSwitches (CCD_Switch *sw) {

	sw->atodcorr = OMIT;
	sw->biascorr = OMIT;
	sw->blevcorr = OMIT;
	sw->crcorr   = OMIT;
	sw->darkcorr = OMIT;
	sw->dqicorr  = OMIT;
	sw->flatcorr = OMIT;
	sw->flashcorr = OMIT;
    sw->fluxcorr = OMIT;
	sw->photcorr = OMIT;
	sw->rptcorr  = OMIT;
	sw->shadcorr = OMIT;
    sw->pctecorr = OMIT;
}

void initIRSwitches (IR_Switch *sw) {

	sw->zsigcorr = OMIT;
	sw->zoffcorr = OMIT;
	sw->dqicorr  = OMIT;
	sw->blevcorr = OMIT;
	sw->noiscorr = OMIT;
	sw->darkcorr = OMIT;
	sw->nlincorr = OMIT;
	sw->flatcorr = OMIT;
	sw->unitcorr = OMIT;
	sw->photcorr = OMIT;
	sw->crcorr   = OMIT;
	sw->rptcorr  = OMIT;

}

