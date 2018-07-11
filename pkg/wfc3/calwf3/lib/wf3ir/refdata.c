/* REFDATA: Contains various routines for reading reference data,
** as well as initializing and freeing the reference data structures.
**
** Revision history:
** H.Bushouse	Oct. 2000	Initial CALNICA to CALWF3 port.
** H.Bushouse	Mar. 2002	Cleaned out remaining unecessary NICMOS
**				procedures.
** H.Bushouse	20-Mar-2002	Created routine RebinRef.
** H.Bushouse	28-Mar-2002	Created routine getFlats.
** H.Bushouse	08-May-2002	Modified to use trlkwerr, trlopenerr, and
				trlreaderr.
** H.Bushouse	14-Aug-2003	Modified to use only 1 node array for nlinfile,
				and to combine pflt, dflt, lflt data into a
				master flat (same as UVIS).
** H.Bushouse	14-Feb-2007	Reduced ALLOWDIFF from 0.1 to 0.01 for use
				with IR subarray exptimes. Added check for
				SUBTYPE in getDarkInfo.
** H.Bushouse	29-Aug-2008	Modified getNlinData and freeNlinData to use
				ncoeff and nerr members of NlinData struct.
** H.Bushouse	09-Jan-2009	Removed FILTER check from getFlatImage because
				that's now handled by checkFlat in getIRFlags.
** H.Bushouse	25-Feb-2009	Added crrpar_in to read parameters for CRCORR.
*/

# include <math.h>
# include <ctype.h>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>

#include "hstcal.h"
# include "hstio.h"	/* defines HST I/O functions */
# include "wf3.h"
# include "wf3info.h"
# include "hstcalerr.h"
# include "trlbuf.h"

# define  FATAL   1
# define  WARNING 2

extern int status;

/* GETNLINDATA: Read linearity coefficients and their associated
** data quality flags and node values from the NLINFILE reference file.
** The coefficient, DQ, and node arrays are returned in a single NlinData 
** structure. The super zero read SCI and ERR arrays are also read and loaded
** into the NlinData structure. Also check to make sure that this reference
** file is the appropriate one to use with the science data being processed.
*/

int getNlinData (WF3Info *wf3, NlinData *nlin) {

/* Arguments:
**	wf3 	 i: WFC3 info structure
**	nlin	 o: data structure containing linearity coeff's and DQ flags
*/

	/* Local variables */
	int i;		/* loop index */
	IODescPtr im;	/* image pointer */

	/* Initialize non-linearity data structure */
	nlin->ncoeff	= 0;
	nlin->nerr	= 0;
	nlin->globalhdr = NULL;
	nlin->coeff     = NULL;
	nlin->error     = NULL;
	nlin->dqual     = NULL;
	nlin->nodes     = NULL;
	nlin->zsci	= NULL;
	nlin->zerr	= NULL;

	/* Open the nonlinearity reference file */
	im = openInputImage (wf3->nlin.name, "", 0);
	if (hstio_err()) {
	    trlopenerr (wf3->nlin.name);
	    return (status = 0);	/* don't abort yet */
	}

	/* Allocate memory for the file's primary header */
	nlin->globalhdr = (Hdr *)calloc(1,sizeof(Hdr));
	if (nlin->globalhdr == NULL) {
	    sprintf (MsgText, "Can't allocate memory for NLIN ref file header");
	    trlerror (MsgText);
	    return (status = 1);
	}

	/* Read the primary header */
	getHeader (im, nlin->globalhdr);
	if (hstio_err()) {
	    trlreaderr (wf3->nlin.name);
	    closeImage (im);
	    freeHdr (nlin->globalhdr);
	    return (status = 1);
	}

	/* Close the image */
	closeImage (im);

	/* Read the NCOEF keyword from the NLINFILE */
        if (getKeyI (nlin->globalhdr, "NCOEF", &nlin->ncoeff)) {
            trlkwerr ("NCOEF", wf3->nlin.name);
            return (status = 1);
        }

	/* Read the NERR keyword from the NLINFILE */
        if (getKeyI (nlin->globalhdr, "NERR", &nlin->nerr)) {
            trlkwerr ("NERR", wf3->nlin.name);
            return (status = 1);
        }

	/* Allocate memory for the nlin data structure */
	nlin->coeff = (FloatHdrData *)calloc(nlin->ncoeff,sizeof(FloatHdrData));
	nlin->error = (FloatHdrData *)calloc(nlin->nerr,sizeof(FloatHdrData));
	nlin->dqual = (ShortHdrData *)calloc(1,sizeof(ShortHdrData));
	nlin->nodes = (FloatHdrData *)calloc(1,sizeof(FloatHdrData));
	nlin->zsci  = (FloatHdrData *)calloc(1,sizeof(FloatHdrData));
	nlin->zerr  = (FloatHdrData *)calloc(1,sizeof(FloatHdrData));
	if (nlin->coeff==NULL || nlin->error==NULL || nlin->dqual==NULL ||
	    nlin->nodes==NULL || nlin->zsci ==NULL || nlin->zerr ==NULL) {
	    sprintf (MsgText, "Can't allocate memory for NLIN ref data");
	    trlerror (MsgText);
	    return (status = 1);
	}

	/* Get the coefficient images from the NLINFILE */
	for (i=0; i < nlin->ncoeff; i++) { 
	     initFloatHdrData (&nlin->coeff[i]);
	     if (getFloatHD (wf3->nlin.name, "COEF", i+1, &nlin->coeff[i]))
		 return (status=1);
	}

	/* Get the error images from the NLINFILE */
	for (i=0; i < nlin->nerr; i++) {
	     initFloatHdrData (&nlin->error[i]);
	     if (getFloatHD (wf3->nlin.name, "ERR", i+1, &nlin->error[i]))
		 return (status=1);
	}

	/* Get the DQ array from the NLINFILE */
	initShortHdrData (&nlin->dqual[0]);
	if (getShortHD (wf3->nlin.name, "DQ", 1, &nlin->dqual[0]))
	    return (status=1);

	/* Get the saturation node array from the NLINFILE */
	for (i=0; i<1; i++) {
	     initFloatHdrData (&nlin->nodes[i]);
	     if (getFloatHD (wf3->nlin.name, "NODE", i+1, &nlin->nodes[i]))
		 return (status=1);
	}

	/* Get the super zero-read SCI array from the NLINFILE */
	initFloatHdrData (&nlin->zsci[0]);
	if (getFloatHD (wf3->nlin.name, "ZSCI", 1, &nlin->zsci[0]))
	    return (status=1);

	/* Get the super zero-read ERR array from the NLINFILE */
	initFloatHdrData (&nlin->zerr[0]);
	if (getFloatHD (wf3->nlin.name, "ZERR", 1, &nlin->zerr[0]))
	    return (status=1);

	/* Successful return */
	return (status = 0);
}

/* GETFLATS: Load each of the 3 possible flat field reference images
** that exist and combine them into a single image to be used when
** applying the flat field correction to the science data.
**
** The low-order flat, which is usually stored at a lower resolution
** than the other 2 flats, is rebinned to match the others before
** being combined.
*/

int getFlats (WF3Info *wf3, SingleNicmosGroup *in, SingleNicmosGroup *flat) {

	/* Local variables */
	int dummy;
	int rx, ry;
	SingleNicmosGroup dflt;
	SingleNicmosGroup lflt;

	/* Function definitions */
	int getFlatImage (WF3Info *, RefImage *, SingleNicmosGroup *);
	int FindBinIR (SingleNicmosGroup *, SingleNicmosGroup *, int *, int *,
		       int *, int *, int *);
	void amul (SingleNicmosGroup *, SingleNicmosGroup *);
	int  unbin2d_ir (SingleNicmosGroup *, SingleNicmosGroup *);

	initSingleNicmosGroup (&dflt);
	initSingleNicmosGroup (&lflt);

	/* Load the pixel-to-pixel flat */
	if (wf3->pfltcorr == PERFORM) {

	    if (getFlatImage (wf3, &wf3->pflt, flat))
		return (status);
	}

	/* Load the delta flat */
	if (wf3->dfltcorr == PERFORM) {

	    /* If pflt also exists, load dflt into separate temporary image */
	    if (wf3->pfltcorr == PERFORM) {

		if (getFlatImage (wf3, &wf3->dflt, &dflt))
		    return (status);

		/* Are the pflt and dflt the same size? */
		if (flat->sci.data.nx != dflt.sci.data.nx ||
		    flat->sci.data.ny != dflt.sci.data.ny) {
		    sprintf (MsgText,
		    "Pixel-to-pixel flat and delta flat are no the same size.");
		    trlerror (MsgText);
		    return (status = SIZE_MISMATCH);
		}

		/* Multiply the pflt and dflt together, leaving
		** the result in flat */
		amul (flat, &dflt);

		/* Free the dflat temporary image */
		freeSingleNicmosGroup (&dflt);

	    /* Otherwise, just load dflt by itself */
	    } else {

		if (getFlatImage (wf3, &wf3->dflt, flat))
		    return (status);
	    }
	}
		
	/* Load the low-order flat */
	if (wf3->lfltcorr == PERFORM) {

	    /* If either pflt or dlft also exist, load lflt into separate
	    ** temporary image */
	    if (wf3->pfltcorr == PERFORM || wf3->dfltcorr == PERFORM) {

		/* This is the normal case; we already have a product
		** in flat */
		if (getFlatImage (wf3, &wf3->lflt, &lflt))
		    return (status);

		/* Allocate an image the same size as flat */
		allocSingleNicmosGroup (&dflt, flat->sci.data.nx,
					flat->sci.data.ny);
		if (hstio_err()) {
		    sprintf (MsgText, "Can't allocate memory for FLAT image");
		    trlerror (MsgText);
		    return (status = 1);
		}
		
		/* Resample lflt into dflt by linear interpolation */
		if ( (status = unbin2d_ir (&lflt, &dflt)))
		    return (status);
		freeSingleNicmosGroup (&lflt);

		/* Combine the resampled lflt with others */
		amul (flat, &dflt);

		freeSingleNicmosGroup (&dflt);

	    /* Otherwise low-order flat is only flat we have */
	    } else {

		if (getFlatImage (wf3, &wf3->lflt, &lflt))
		    return (status);

		/* Figure out how much to expand the low-order flat
		** to match the science data binning. */
		if (FindBinIR (in, &lflt, &dummy, &rx, &ry, &dummy, &dummy)) {
		    if (status == REF_TOO_SMALL)
			status = 0;
		    else
			return (status);
		}

		/* Create an expanded version of the low-order flat,
		** to match binning of science data */
		allocSingleNicmosGroup (flat, rx*lflt.sci.data.nx,
					      ry*lflt.sci.data.ny);
		if (hstio_err()) {
		    sprintf (MsgText, "Can't allocate memory for FLAT image");
		    trlerror (MsgText);
		    return (status = 1);
		}

		if ( (status = unbin2d_ir (&lflt, flat)))
		    return (status);
		freeSingleNicmosGroup (&lflt);

	    }
	}

	/* Now flat contains the product of up to 3 flats. */

	return (status = 0);
}

/* GETFLATIMAGE: Load the flat field reference file image.
*/

int getFlatImage (WF3Info *wf3, RefImage *ref, SingleNicmosGroup *flat) {

/* Arguments:
**	nic	i: NICMOS info structure
**	flat	o: flat field image
*/

	/* Local variables */

	/* Function definitions */
	int getRefImage (RefImage *, int, SingleNicmosGroup *);

	/* Read the reference image */
	if (getRefImage (ref, 1, flat))
	    return (status);

	/* Successful return */
	return (status = 0);
}

/* GETDARKINFO: Load information from the DARKFILE reference file.
*/

int getDarkInfo (WF3Info *wf3) {

/* Arguments:
**	nic	 i: NICMOS info structure
*/

	/* Local variables */
	int i;			/* loop index */
	char kword[8+1];	/* keyword name */

	SingleNicmosGroup dark;

	/* Function definitions */
	int  getRefImage (RefImage *, int, SingleNicmosGroup *);
	void checkKeyS (char *, Hdr *, char *, char *, int);

	/* Read the first group of the dark reference file */
	if (getRefImage (&wf3->dark, 1, &dark))
	    return (status);

	/* Check the DARK file SAMP_SEQ value */
	checkKeyS (wf3->dark.name, dark.globalhdr, "SAMP_SEQ",
		   wf3->sampseq, FATAL);
	if (status)
	    return (status);

	/* Check the DARK file SUBTYPE value */
	checkKeyS (wf3->dark.name, dark.globalhdr, "SUBTYPE",
		   wf3->subtype, FATAL);
	if (status)
	    return (status);

	/* Find out how many dark images are in this file */
	wf3->ndarks = 0;
	if (getKeyI (dark.globalhdr, "NUMEXPOS", &wf3->ndarks)) {
	    trlkwerr ("NUMEXPOS", wf3->dark.name);
	    return (status = 1);
	}

	/* Read the list of exposure times */
	for (i = 0; i < wf3->ndarks; i++) {
	     sprintf (kword, "EXPOS_%d", i+1);
	     wf3->dtimes[i] = 0;
	     if (getKeyD (dark.globalhdr, kword, &(wf3->dtimes[i]))) {
		 trlkwerr (kword, wf3->dark.name);
		 return (status = 1);
	     }
	}

	freeSingleNicmosGroup (&dark);

	/* Successful return */
	return (status = 0);

}

# define ALLOWDIFF 0.01	/* Max allowed exptime difference */

int getDarkImage (WF3Info *wf3, SingleNicmosGroup *dark1, int ngroup) {

/* Arguments:
**	nic	 i: NICMOS info structure
**	dark1	 o: dark current image
**	ngroup	 i: group number
*/

	/* Local variables */
	int i;					/* loop index */
	double etime_lower, etime_upper;	/* dark exposure times */
	SingleNicmosGroup dark2;		/* temporary dark image */

	/* Function definitions */
	int  getRefImage (RefImage *, int, SingleNicmosGroup *);
	void dark_interp (SingleNicmosGroup *, SingleNicmosGroup *, 
			  double, double, double);
	void dark_extrap (SingleNicmosGroup *, double, double);

    etime_lower=0.0f;
    etime_upper=0.0f;
    
	/* Initialize frame interpolation information */
	wf3->DarkType   = 0;
	wf3->darkframe1 = 0;
	wf3->darkframe2 = 0;

	/* Find the appropriate match between science and dark exp time */
	if (strncmp (wf3->sampseq, "NONE", 4) != 0) {

	    /* If the science data are from a MultiAccum standard exposure
	    ** time sequence, then find the ref file group that matches the
	    ** exposure time of the science group. */
	    for (i = 0; i < wf3->ndarks; i++) {
		 if (fabs(wf3->dtimes[i]-wf3->exptime[ngroup-1]) <= ALLOWDIFF){
		     wf3->DarkType = MATCH;
		     wf3->darkframe1 = i+1;
		     break;
		 }
	    }

	    /* Make sure we found a match */
	    if (wf3->darkframe1 == 0) {
		sprintf (MsgText,
			 "Can't find matching dark time in %s for group %d\n",
			 wf3->dark.name, ngroup);
		trlerror (MsgText);
		return (status = 1);
	    }

	} else {

	    /* If the science data aren't from a MultiAccum standard exposure
	    ** time sequence, then find an exposure time that matches the
	    ** science data, or two times that bracket it */

	    for (i = 0; i < wf3->ndarks; i++) {
		 if (fabs(wf3->dtimes[i]-wf3->exptime[ngroup-1]) <= ALLOWDIFF){
		     wf3->DarkType = MATCH;
		     wf3->darkframe1 = i+1;
		     break;
		 } else if (wf3->dtimes[i] < wf3->exptime[ngroup-1]) {
		     wf3->DarkType = INTERP;
		     wf3->darkframe1 = i+1;
		     etime_lower = wf3->dtimes[i];
		 } else if (wf3->dtimes[i] > wf3->exptime[ngroup-1] &&
		     wf3->darkframe2 == 0) {
		     wf3->DarkType = INTERP;
		     wf3->darkframe2 = i+1;
		     etime_upper = wf3->dtimes[i];
		 }
	    }
	}

	/* If there's a dark image with a matching exposure time, load it. */
	if (wf3->DarkType == MATCH) {
	    if (getRefImage (&wf3->dark, wf3->darkframe1, dark1))
		return (status);

	/* Otherwise, load two images and interpolate between them */
	} else if (wf3->darkframe1 != 0  &&  wf3->darkframe2 != 0) {
	    if (getRefImage (&wf3->dark, wf3->darkframe1, dark1))
		return (status);

	    if (getRefImage (&wf3->dark, wf3->darkframe2, &dark2))
		return (status);

	    dark_interp (dark1, &dark2, etime_lower, etime_upper,
			 wf3->exptime[ngroup-1]);

	    freeSingleNicmosGroup (&dark2);

	/* Otherwise, extrapolate outside the available range */
	} else if (wf3->darkframe1 !=0) {
	    wf3->DarkType = EXTRAP;
	    if (getRefImage (&wf3->dark, wf3->darkframe1, dark1))
		return (status);

	    dark_extrap (dark1, etime_lower, wf3->exptime[ngroup-1]);

	} else if (wf3->darkframe2 !=0) {
	    wf3->DarkType = EXTRAP;
	    if (getRefImage (&wf3->dark, wf3->darkframe2, dark1))
		return (status);

	    dark_extrap (dark1, etime_upper, wf3->exptime[ngroup-1]);
	}

	/* Successful return */
	return (status = 0);

}

/* DARK_INTERP: Do a linear interpolation of two dark images */

void dark_interp (SingleNicmosGroup *d1, SingleNicmosGroup *d2, 
	     double etime_d1, double etime_d2, double exptime) {

/* Arguments:
**	d1		io: dark image #1; overwritten by interpolated image
**	d2	 	 i: dark image #2
**	etime_d1   	 i: exposure time of dark image #1
**	etime_d2   	 i: exposure time of dark image #2
**	exptime		 i: desired exposure time of interpolated image
*/
	/* Local variables */
	float frac;

	/* Function definitions */
	void aadd  (SingleNicmosGroup *, SingleNicmosGroup *);
	void asub  (SingleNicmosGroup *, SingleNicmosGroup *);
	void amulk (SingleNicmosGroup *, float);

	/* Compute the fractional exposure time */
	frac = (exptime - etime_d1) / (etime_d2 - etime_d1);

	/* Interpolate the images */
	asub  (d2, d1);
	amulk (d2, frac);
	aadd  (d1, d2);
}

/* DARK_EXTRAP: Do a linear extrapolation of a dark image */

void dark_extrap (SingleNicmosGroup *d1, double etime, double exptime) {

/* Arguments:
**	d1		io: dark image; overwritten by extrapolated image
**	etime		 i: exposure time of dark image
**	exptime		 i: desired exposure time of extrapolated image
*/
	/* Local variables */
	float frac;

	/* Function definitions */
	void amulk (SingleNicmosGroup *, float);

	/* Compute the fractional exposure time */
	frac = exptime / etime;

	/* Extrapolate the image */
	amulk (d1, frac);
}

/* GETREFIMAGE: Load the data from a reference image. */

int getRefImage (RefImage *step, int extver, SingleNicmosGroup *image) {

/* Arguments:
**	step	io: calibration step info structure
**	extver	 i: extension version number to load
**	image	 o: reference image data
*/

	/* Function definitions */

	/* Initialize the image data structure */
	initSingleNicmosGroup (image);

	/* Read the data */
	if (getSingleNicmosGroup (step->name, extver, image))
	    status = 1;
	if (hstio_err() || status) {
	    trlreaderr (step->name);
	    freeSingleNicmosGroup (image);
	    return (status = 1);
	}

	/* Successful return */
	return (status = 0);
}

void freeNlinData (NlinData *nlin) {

	/* Local variables */
	int i;				/* loop index */

	/* Free non-linearity data structure */
	if (nlin->globalhdr != NULL)
	    free (nlin->globalhdr);
	if (nlin->coeff != NULL) {
	    for (i=0; i<nlin->ncoeff; i++)
		 freeFloatHdrData (&nlin->coeff[i]);
	    free (nlin->coeff);
	}
	if (nlin->error != NULL) {
	    for (i=0; i<nlin->nerr; i++)
		 freeFloatHdrData (&nlin->error[i]);
	    free (nlin->error);
	}
	if (nlin->dqual != NULL) {
	    freeShortHdrData   (&nlin->dqual[0]);
	    free (nlin->dqual);
	}
	if (nlin->nodes != NULL) {
	    freeFloatHdrData (&nlin->nodes[0]);
	    /*freeFloatHdrData (&nlin->nodes[1]);*/
	    free (nlin->nodes);
	}
	if (nlin->zsci  != NULL) {
	    freeFloatHdrData (&nlin->zsci[0]);
	    free (nlin->zsci);
	}
	if (nlin->zerr  != NULL) {
	    freeFloatHdrData (&nlin->zerr[0]);
	    free (nlin->zerr);
	}
}

void checkKeyI (char *filename, Hdr *header, char *keyword, int scival,
	       int severity) {

	/* Local variables */
	int keyval;

	keyval = 0;

	/* Read the keyword from the ref file */
	if (getKeyI (header, keyword, &keyval)) {
	    trlkwerr (keyword, filename);
	    status = 1;
	    return;
	}

	/* Does it match the science image value? */
	if (keyval != scival) {
	    sprintf (MsgText, "%s %s=%d does not match science data",
		     filename, keyword, keyval);
	    if (severity == FATAL) {
		trlerror (MsgText);
	    } else
		trlwarn  (MsgText);
	}
}

void checkKeyS (char *filename, Hdr *header, char *keyword, char *scival,
	       int severity) {

	/* Local variables */
	char keyval[SZ_FITS_VAL+1];

	keyval[0] = '\0';

	/* Read the keyword from the ref file */
	if (getKeyS (header, keyword, keyval)) {
	    trlkwerr (keyword, filename);
	    status = 1;
	    return;
	}

	/* Does it match the science image value? */
	if (strncmp (keyval, scival, strlen(keyval)) != 0) {
	    sprintf (MsgText, "%s %s=\"%s\" does not match science data",
		     filename, keyword, keyval);
	    if (severity == FATAL) {
		trlerror (MsgText);
		status = 1;
	    } else
		trlwarn  (MsgText);
	}
}

int RebinRef (SingleNicmosGroup *in, SingleNicmosGroup *ref, int avg) {

        /* Local variables */
        int rx, ry;             /* bin ratio of dark vs. science image */
        int x0, y0;             /* offset of dark vs. science image */
        int same_size;          /* true if no binning of dark image required */
        SingleNicmosGroup z;    /* scratch array for dark image */

        /* Function definitions */
        int FindBinIR (SingleNicmosGroup *, SingleNicmosGroup *, int *, int *,
                       int *, int *, int *);
        int bin2d_ir (SingleNicmosGroup *, int, int, int, int, int,
                      SingleNicmosGroup *);
	int copyGroup (SingleNicmosGroup *, SingleNicmosGroup *);

	/* Compare binning and dimensions of science and reference image;
	** get same_size flag, and get info about binning and offset
	** for use by bin2d. */
	if ( (status = FindBinIR (in, ref, &same_size, &rx, &ry, &x0, &y0)))
	    return (status);

	/* Rebin or extract subarray of reference image, if necessary */
	if (!same_size) {
	    initSingleNicmosGroup (&z);
	    allocSingleNicmosGroup (&z, in->sci.data.nx, in->sci.data.ny);
	    if ( (status = bin2d_ir (ref, x0, y0, rx, ry, avg, &z))) {
		trlerror ("Reference image size mismatch.");
		return (status);
	    }
	    freeSingleNicmosGroup (ref);
	    if (copyGroup (ref, &z))
		return (status);
	    freeSingleNicmosGroup (&z);
	}

	return (status = 0);
}

# include   "xtables.h"

# include   "wf3rej.h"
# include   "wf3dq.h"
# include   "rej.h"

static int strtor (char *, float []);

/*  crrpar_in -- Read parameters either from user input or table.

  Description:
  ------------
  Reads CL parameters and does necessary checkings
  
  Input parameters from crrej reference table:
  -------------------------------------------
  Col. Name     Parameter
  "skysub"      sky         Sky levels subtraction scheme
  "crsigmas"    sigmas      Rejection thresholds
  "crradius"    radius      Radius (in pixels) to propagate the cosmic ray
  "crthresh"    thresh      Propagation factor
  "initgues"    initgues    Scheme of computing initial-guess image
  "scalense"    scalense    multiplicative noise in percents
  "badinpdq"    badinpdq    Data quality pset
  "crmask"      mask        flag CR-rejected pixels in input files?

  Input parameters from input image primary header:
  ------------------------------------------------


  Date          Author          Description
  ----          ------          -----------
  24-Sep-1998   W.J. Hack       Initial ACS Version	
  29-Aug-2000	H.A. Bushouse	Initial WFC3 Version
  25-Feb-2009   H. Bushouse	Modified version of wf3rej/rejpar_in to use with
				IR CRCORR.
  06-Mar-2009	H. Bushouse	Added logic to use new IRRAMP column info in
				CRREJTAB to identify rows that apply to IR 
				CRCORR.
  10-Apr-2009	H. Bushouse	Fixed bug in calls to read irramp in each row
				(variable "row" was used before being set).
  13-Oct-2009	H. Bushouse	Report badinpdq value, now that it's used in
				the cridcalc routine.
  02-Mar-2010	H. Bushouse	Fixed initialization of maxcrsplit.
				(calwf3 v2.0)
*/

int crrpar_in (clpar *par, int newpar[], int nimgs, float exptot, int *niter,
	       float sigma[]) {

    extern int status;

    IRAFPointer     tp;
    IRAFPointer     colptr, colptr1, colptr2;
    int             i, nrows, nmatch, row;
    int             crsplit_in, crsplit, maxcrsplit, irramp;
    float           exp_in, meanexp, mindiff, diff;

    void    PrRefInfo (char *, char *, char *, char *, char *);
    void    WhichError (int);

/* -------------------------------- begin ---------------------------------- */
    row=0;
    mindiff=0.0f;
    crsplit_in = nimgs;

    exp_in = exptot;
    par->meanexp = exp_in;

    /* if all parameters are specified by the user, no need to open the 
	    reference CRREJ table */
    if (newpar[0] < MAX_PAR) {

        tp = c_tbtopn (par->tbname, IRAF_READ_ONLY, 0);
        if (c_iraferr() != 0) {
            sprintf (MsgText,"CRREJTAB table '%s' does not exist", par->tbname);
            trlerror (MsgText);
            return (status = TABLE_ERROR);
        }
        nrows = c_tbpsta (tp, TBL_NROWS);

        /* read the columns IRRAMP, CRSPLIT and MEANEXP */
        c_tbcfnd1 (tp, "crsplit", &colptr);
        if (colptr == 0) {
            trlerror ("column CRSPLIT does not exist in CRREJTAB");
            return (status = COLUMN_NOT_FOUND);
        }
        c_tbcfnd1 (tp, "meanexp", &colptr1);
        if (colptr1 == 0) {
            trlerror ("column MEANEXP does not exist in CRREJTAB\n");
            return (status = COLUMN_NOT_FOUND);
        }
        c_tbcfnd1 (tp, "irramp", &colptr2);
        if (colptr2 == 0) {
            trlerror ("column IRRAMP does not exist in CRREJTAB");
            return (status = COLUMN_NOT_FOUND);
        }
        nmatch = 0;

        /* find the largest value in the CRSPLIT column for rows in
           in which IRRAMP=yes */
	maxcrsplit = 0;
	for (i = 1; i <= nrows; i++) {
	    c_tbegti (tp, colptr2, i, &irramp);
	    if (!irramp) continue;
	    c_tbegti (tp, colptr, i, &crsplit);
	    if (crsplit > maxcrsplit) maxcrsplit = crsplit;
	}
        if (crsplit_in >= maxcrsplit) crsplit_in = maxcrsplit;

        /* find the matching row in CRREJTAB to use */
        for (i = 1; i <= nrows; i++) {
            c_tbegti (tp, colptr2, i, &irramp);
	    if (!irramp) continue;
            c_tbegti (tp, colptr, i, &crsplit);
            c_tbegtr (tp, colptr1, i, &meanexp);
            diff = meanexp - exp_in;
            if (crsplit_in == crsplit && diff >= 0.) {
                nmatch++;
                if (nmatch == 1) mindiff = diff;
                if (diff <= mindiff) {
                    row = i;
                    mindiff = diff;
                }
            }
        }
        if (nmatch == 0) {
            trlerror (" No matching CRSPLIT and MEANEXP in CRREJTAB");
            return (status = ROW_NOT_FOUND);
        }

        /* read the sigmas parameter */ 
        if (newpar[CRSIGMAS] == 0) {
            c_tbcfnd1 (tp, "crsigmas", &colptr);
            if (colptr == 0) {
                trlerror ("column CRSIGMAS does not exist in CRREJTAB");
                return (status = COLUMN_NOT_FOUND);
            }
            c_tbegtt (tp, colptr, row, par->sigmas, CHAR_LINE_LENGTH);
        }

        /* read other parameters */
        if (newpar[SKYSUB] == 0) {
            c_tbcfnd1 (tp, "skysub", &colptr);
            if (colptr == 0) {
                trlerror ("column SKYSUB does not exist in CRREJTAB");
                return (status = COLUMN_NOT_FOUND);
            }
            c_tbegtt (tp, colptr, row, par->sky, SZ_FITS_REC);
        }

        /* CR propagation parameters */
        if (newpar[CRRADIUS] == 0) {
            c_tbcfnd1 (tp, "crradius", &colptr);
            if (colptr == 0) {
                trlerror ("column CRRADIUS does not exist in CRREJTAB");
                return (status = COLUMN_NOT_FOUND);
            }
            c_tbegtr (tp, colptr, row, &par->radius);
        }
        if (newpar[CRTHRESH] == 0) {
            c_tbcfnd1 (tp, "crthresh", &colptr);
            if (colptr == 0) {
                trlerror ("column CRTHRESH does not exist in CRREJTAB");
                return (status = COLUMN_NOT_FOUND);
            }
            c_tbegtr (tp, colptr, row, &par->thresh);
        }

            /* figure out how to do initial comparison image */
        if (newpar[INITGUES] == 0) {
            c_tbcfnd1 (tp, "initgues", &colptr);
            if (colptr == 0) {
                trlerror ("column INITGUES does not exist in CRREJTAB");
                return (status = COLUMN_NOT_FOUND);
            }
            c_tbegtt (tp, colptr, row, par->initgues, SZ_FITS_REC);
        }

        /* read the noise model */
        if (newpar[SCALENSE] == 0) {
            c_tbcfnd1 (tp, "scalense", &colptr);
            if (colptr == 0) {
                trlerror ("column SCALENSE does not exist in CRREJTAB");
                return (status = COLUMN_NOT_FOUND);
            }
            c_tbegtr (tp, colptr, row, &par->scalense);
        }

        if (newpar[BADINPDQ] == 0) {
            c_tbcfnd1 (tp, "badinpdq", &colptr);
            if (colptr == 0) {
                trlerror ("column BADINPDQ does not exist in CRREJTAB");
                return (status = COLUMN_NOT_FOUND);
            }
            c_tbegts (tp, colptr, row, &par->badinpdq);
        }

        if (newpar[CRMASK] == 0) {
            c_tbcfnd1 (tp, "crmask", &colptr);
            if (colptr == 0) {
                trlerror ("column CRMASK does not exist in CRREJTAB");
                return (status = COLUMN_NOT_FOUND);
            }
            c_tbegti (tp, colptr, row, &par->mask);
        }

        c_tbtclo (tp);
    }

    PrRefInfo ("crrejtab", par->tbname, "", "", "");

    /* parse the sigmas string into numbers */
    *niter = strtor (par->sigmas, sigma);
    if (status != WF3_OK) {
        WhichError (status);
        return (status);
    }
    if (*niter > MAX_ITER) {
        sprintf (MsgText,"No more than %d iterations permitted.", MAX_ITER);
        trlerror (MsgText);
        return (status = ERROR_RETURN);
    }
    if (*niter <= 0) {
        trlerror ("Number of iterations is ZERO.");
        return (status = ERROR_RETURN);
    }

    /* other fixed (for now) parameters */
    par->crval = (short) DATAREJECT;
    par->fillval = 0.;

    /* print out which parameter are used */
    if (par->verbose) {
        sprintf (MsgText,"\n number of samples = %d", nimgs);
        trlmessage (MsgText);
        sprintf (MsgText," CRREJ ref table used: %s", par->tbname);
        trlmessage (MsgText);
/*
        sprintf (MsgText," initial guess method: %s", par->initgues);
        trlmessage (MsgText);
*/
        sprintf (MsgText," total exposure time = %0.1f", exptot);
        trlmessage (MsgText);
        sprintf (MsgText," sigmas used: %s", par->sigmas);
        trlmessage (MsgText);
/*
        sprintf (MsgText," sky subtraction used: %s", par->sky);
        trlmessage (MsgText);
        sprintf (MsgText," rejection radius = %0.1f", par->radius);
        trlmessage (MsgText);
        sprintf (MsgText," propagation threshold = %0.1f", par->thresh);
        trlmessage (MsgText);
        sprintf (MsgText," scale noise = %0.1f%%", par->scalense);
        trlmessage (MsgText);
*/
        sprintf (MsgText," input bad bits value = %d", par->badinpdq);
        trlmessage (MsgText);
    }
    return (status);
}

/*  strtor -- convert a string of real numbers into a real array 

    Description:
    ------------
    If the input string is blank(s), this routine will return 0.  If there are
    characters other than digits, decimal point, comma, semi-colon, colon, or
    slash in the input string, this routine will issue an error message.
    
    NOTE: This function sets 'status' upon error, but returns a 
        different variable.

    Date            Author          Description
    ----            ------          -----------
    09-May-1996     J.-C. Hsu       Adapt from the SPP code strtor.x
*/

static int strtor (char *str, float arr[]) {

    extern int status;

    int	    ip;         /* index of the string to be searched */
    int	    n, i, ipx;
    double  rval;
    char    tmp[100];

    n = 0;
    ip = 0;
    ipx = 0;

    /* Initialize value to allow for proper error-checking
        after this function returns to the calling routine */
    status = WF3_OK;

    while (1) {
        if (str[ip] == ',' || str[ip] == ';' || str[ip] == '/' || 
            str[ip] == '\0') {
            for (i = 0; i < (ip-ipx); ++i)
                tmp[i] = str[ipx+i];
            
            tmp[ip-ipx] = '\0';
            rval = strtod(tmp, (char **)NULL);
            
            if (!rval && (ip-ipx) != 0) {
                sprintf (MsgText, "illegal input string '%s'", str);
                trlerror (MsgText);
                status = INVALID_VALUE;
                return (0);
            }
            arr[n] = (float) rval;
            ++n;
            ipx = ip + 1;
        }
        if (str[ip] == '\0') return (n);
        ++ip;
    }
}
