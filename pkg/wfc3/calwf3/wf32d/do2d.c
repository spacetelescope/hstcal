/* Do basic 2-D image reduction for the current set of extensions.
   If we're processing the first set (extver = 1), then the primary
   header will be updated with history information, and the calibration
   switches in the header will be reset from "PERFORM" to "COMPLETE".

   Warren Hack, 1998 June 10:
	Initial ACS version.
   Warren Hack, 1999 Jan 7:
   	Revised to not proceed unless it detects a overscan-trimmed image.
	The keywords LTV1, LTV2 should be zero for an overscan-trimmed
	image, otherwise it will return without processing.
   Howard Bushouse, 2000 Aug 28:
	Initial revisions for WFC3 version.
   H. Bushouse, 2001 May 21:
	Modified for WFC3 so that doPhot gets called for both UVIS imsets.
   H. Bushouse, 2001 Nov 15:
        Updated for use of ccdbias array and new PHOTCORR methods
	(in accordance with changes in CALACS).
   H. Bushouse, 2002 May 21:
	Updated PhotMode calling sequence to pass only image header instead
	of entire image.
   H. Bushouse, 2002 June 17:
	Added update to BUNIT keyword after flat-field correction to reflect
	change in datatype for science data after calibration. Removed all
	references to "statcorr" and STATFLAG - always perform "doStat"
	(in accordance with CALACS changes).
   H. Bushouse, 2003 Oct 24:
	Added update to EXPSCORR switch header keyword in output file
	(CALACS change).
   H. Bushouse, 2005 Mar 2:
	Modified LTV-based test in OscnTrimmed routine to make it compatible
	with WFC3 binned images.
   H. Bushouse, 2007 Feb 16:
	Modified call to PhotMode to use science extension header, rather
	than primary header, because that's where phot keywords are.
   H. Bushouse, 2011 Sep 7:
	Modified PhotMsg to use new phot table instead of graph and comp tabs.
   M. Sosey, 2012 May 7:
    changed "electrons" to "ELECTRONS" to be consistent, also there are places
    in the code that are doing a case sensative check. I'm going to try and make
    those case insensitive as well
   M. Sosey, 2013 July 3
   Added code to implement new FLUXCORR step, see #1011

*/


# include <string.h>
# include <stdio.h>

# include "hstio.h"

# include "wf3.h"
# include "wf3info.h"
# include "err.h"

static void DarkMsg (WF3Info *, int);
static void dqiMsg (WF3Info *, int);
static void PhotMsg (WF3Info *);
static void ShadMsg (WF3Info *, int);
static void FlatMsg (WF3Info *, int);
static int OscnTrimmed (Hdr*, Hdr *);

Hdr phdr; /*primary header for input image, all output information saved here*/
Hdr scihdr; /*science header in case of subarray image to detect chip*/
IODescPtr ip = NULL;

int Do2D (WF3Info *wf32d, int extver) {

/* arguments:
WF3Info *wf32d   i: calibration switches and info
int extver       i: "imset" number, the current set of extensions
*/

	extern int status;

	SingleGroup x;	/* used for both input and output */
	int option = 0;
	float meandark;		/* mean value of dark (for history) */
	int done;	/* true means error array was initialized in donoise */
	int i;

	int logit;
	char buff[SZ_FITS_REC+1];
	Bool subarray;

	int CCDHistory (WF3Info *, Hdr *);
	int doDark (WF3Info *, SingleGroup *, float *);
	int darkHistory (WF3Info *, Hdr *);
	int doDQI (WF3Info *, SingleGroup *);
	int dqiHistory (WF3Info *, Hdr *);
	int doFlat (WF3Info *, int, SingleGroup *);
	int flatHistory (WF3Info *, Hdr *);
	int doNoise (WF3Info *, SingleGroup *, int *);
	int noiseHistory (Hdr *);
	int doPhot (WF3Info *, SingleGroup *);
	int PhotMode (WF3Info *, Hdr *);
	int photHistory (WF3Info *, Hdr *);
    int fluxHistory (WF3Info *, Hdr *);
	int doShad (WF3Info *, int, SingleGroup *);
	int shadHistory (WF3Info *, Hdr *);
	int doStat (SingleGroup *, short);
	int GetGrp (WF3Info *, Hdr *);
	int OmitStep (int);
	int PutKeyFlt (Hdr *, char *, float, char *);
	int PutKeyStr (Hdr *, char *, char *, char *);
	void PrSwitch (char *, int);
	void PrRefInfo (char *, char *, char *, char *, char *);
	void TimeStamp (char *, char *);
	void UCalVer (Hdr *);
	void UFilename (char *, Hdr *);
	int UpdateSwitch (char *, int, Hdr *, int *);
	int GetCCDTab (WF3Info *, int, int);
	int GetKeyBool (Hdr *, char *, int, Bool, Bool *);

	initSingleGroup (&x);

	/* Open the input image. */
	getSingleGroup (wf32d->input, extver, &x);
	if (hstio_err())
	    return (status = OPEN_FAILED);
	if (wf32d->printtime)
	    TimeStamp ("Input read into memory", wf32d->rootname);

	/* Get header info that varies from imset to imset. */
	if (GetGrp (wf32d, &x.sci.hdr)) {
	    freeSingleGroup (&x);
	    return (status);
	}


	/* For the CCD, update key primary header keywords.	*/
	if (wf32d->detector != IR_DETECTOR) {

	    /* Has the CCD image been overscan trimmed? */
	    if (OscnTrimmed (x.globalhdr, &x.sci.hdr)) {
		    freeSingleGroup (&x);
		    return (status);
	    }

            /* Read in keywords from primary header... */
	    if (GetKeyBool (x.globalhdr, "SUBARRAY", NO_DEFAULT, 0, &subarray))
		return (status);

	    if (subarray)
		wf32d->subarray = YES;
	    else
		wf32d->subarray = NO;

	    /* Get values from tables, using same function used in WF3CCD. */
	    if (GetCCDTab (wf32d, x.sci.data.nx, x.sci.data.ny)) {
		freeSingleGroup (&x);
		return (status);
            }

	    if (PutKeyFlt (x.globalhdr, "ATODGNA", wf32d->atodgain[0], ""))
		return (status);
	    if (PutKeyFlt (x.globalhdr, "ATODGNB", wf32d->atodgain[1], ""))
		return (status);
	    if (PutKeyFlt (x.globalhdr, "ATODGNC", wf32d->atodgain[2], ""))
		return (status);
	    if (PutKeyFlt (x.globalhdr, "ATODGND", wf32d->atodgain[3], ""))
		return (status);
	    if (PutKeyFlt (x.globalhdr, "READNSEA", wf32d->readnoise[0], ""))
		return (status);
	    if (PutKeyFlt (x.globalhdr, "READNSEB", wf32d->readnoise[1], ""))
		return (status);
	    if (PutKeyFlt (x.globalhdr, "READNSEC", wf32d->readnoise[2], ""))
		return (status);
	    if (PutKeyFlt (x.globalhdr, "READNSED", wf32d->readnoise[3], ""))
		return (status);
	    trlmessage ("");

	    PrRefInfo ("ccdtab", wf32d->ccdpar.name, wf32d->ccdpar.pedigree,
		wf32d->ccdpar.descrip, wf32d->ccdpar.descrip2);

	    if (extver == 1 ) {
		if (CCDHistory (wf32d, x.globalhdr))
		    return (status);
	    }

	}

	/* Fill in the error array, if it initially contains all zeros. */
	if (wf32d->noiscorr == PERFORM) {
	    if (doNoise (wf32d, &x, &done))
		return (status);
	    if (done) {
		if (extver == 1) {
		    if (noiseHistory (x.globalhdr))
			return (status);
		}
		trlmessage ("         Uncertainty array initialized,");
		buff[0] = '\0';

		if (wf32d->detector != IR_DETECTOR) {
		    sprintf (MsgText, "    readnoise =");
		    for (i=0; i < NAMPS-1; i++) {
			 if (wf32d->readnoise[i] > 0) {
			     sprintf (buff, "%.5g,",wf32d->readnoise[i]);
			     strcat (MsgText, buff);
			 }
		    }
		    if (wf32d->readnoise[NAMPS-1] > 0) {
			sprintf (buff, "%.5g",wf32d->readnoise[NAMPS-1]);
			strcat (MsgText, buff);
		    }
		    trlmessage (MsgText);

		    sprintf(MsgText, "    gain =");
		    for (i=0; i < NAMPS-1; i++) {
			 if (wf32d->atodgain[i] > 0) {
			     sprintf (buff, "%.5g,",wf32d->atodgain[i]);
			     strcat (MsgText, buff);
			 }
		    }
		    if (wf32d->atodgain[NAMPS-1] > 0) {
			sprintf(buff, "%.5g",wf32d->atodgain[NAMPS-1]);
			strcat (MsgText, buff);
		    }
		    trlmessage (MsgText);
		}
		if (wf32d->printtime)
		TimeStamp ("Uncertainty array initialized", wf32d->rootname);
	    }
	}

	/* Data quality initialization and (for the CCD) check saturation. */
	dqiMsg (wf32d, extver);
	if (wf32d->dqicorr == PERFORM ||
	    (wf32d->dqicorr == DUMMY && wf32d->detector != IR_DETECTOR)) {
	    if (doDQI (wf32d, &x))
		return (status);
	    PrSwitch ("dqicorr", COMPLETE);
	    if (wf32d->printtime)
		TimeStamp ("DQICORR complete", wf32d->rootname);
	}
	if (extver == 1 && !OmitStep (wf32d->dqicorr))
	    if (dqiHistory (wf32d, x.globalhdr))
		return (status);


	/* Subtract dark image. */
	DarkMsg (wf32d, extver);
	if (wf32d->darkcorr == PERFORM) {
	    if (doDark (wf32d, &x, &meandark))
		return (status);

	    sprintf(MsgText,"Mean of dark image (MEANDARK) = %g",meandark);
	    trlmessage(MsgText);

	    if (PutKeyFlt (&x.sci.hdr, "MEANDARK", meandark,
			   "mean of dark values subtracted"))
		return (status);
	    PrSwitch ("darkcorr", COMPLETE);

	    if (wf32d->printtime)
		TimeStamp ("DARKCORR complete", wf32d->rootname);
	}
	if (extver == 1 && !OmitStep (wf32d->darkcorr))
	    if (darkHistory (wf32d, x.globalhdr))
		return (status);

	/*  These corrections are not yet ready to be run...  They require
		the use of Sections...

		Multiply by flat field(s).
	*/

	FlatMsg (wf32d, extver);
	if (wf32d->flatcorr == PERFORM) {
	    if (doFlat (wf32d, extver, &x))
		return (status);
	    PrSwitch ("flatcorr", COMPLETE);
	    if (wf32d->printtime)
		TimeStamp ("FLATCORR complete", wf32d->rootname);
	    if (PutKeyStr (&x.sci.hdr, "BUNIT", "ELECTRONS", ""))
		return (status);
	    if (PutKeyStr (&x.err.hdr, "BUNIT", "ELECTRONS", ""))
		return (status);
	}
	if (extver == 1 && !OmitStep (wf32d->flatcorr))
	    if (flatHistory (wf32d, x.globalhdr))
		return (status);

	/*
		Apply shutter shading correction.
	*/
	ShadMsg (wf32d, extver);
	if (wf32d->shadcorr == PERFORM) {
	    if (doShad (wf32d, extver, &x))
		return (status);
	    PrSwitch ("shadcorr", COMPLETE);
	    if (wf32d->printtime)
		TimeStamp ("SHADCORR complete", wf32d->rootname);
	}
	if (extver == 1 && !OmitStep (wf32d->shadcorr))
	    if (shadHistory (wf32d, x.globalhdr))
		return (status);

	/* Assign values to photometry keywords. */
	trlmessage ("");
	PrSwitch ("photcorr", wf32d->photcorr);
	trlmessage ("");
	PrSwitch ("fluxcorr", wf32d->fluxcorr);

	/* Update the PHOTMODE keyword regardless of PHOTCORR. */
	if (PhotMode (wf32d, &x.sci.hdr))
	    return (status);

	/* Calculate the photometry keyword values if set to PERFORM... */
	if (wf32d->photcorr == PERFORM) {
	    if (doPhot (wf32d, &x))
		    return (status);
	    PhotMsg (wf32d);
	    PrSwitch ("photcorr", COMPLETE);
	    if (wf32d->printtime)
		    TimeStamp ("PHOTCORR complete",wf32d->rootname);

	}

	if (extver == 1 && !OmitStep (wf32d->photcorr))
	    if (photHistory (wf32d, x.globalhdr))
    		return (status);

    if (extver == 1 && !OmitStep (wf32d->fluxcorr))
        if (fluxHistory (wf32d, x.globalhdr))
            return(status);

	/* Compute min, max, mean, etc. of good science data. */
	if (doStat (&x, wf32d->sdqflags))
	    return (status);
	if (wf32d->printtime)
	    TimeStamp ("Image statistics computed", wf32d->rootname);

	/* Write this imset to the output file.  If extver is one, the
	   CAL_VER and FILENAME keywords will be updated, and the primary
	   header will be written.
	*/
	if (extver == 1) {
	    UCalVer (x.globalhdr);
	    UFilename (wf32d->output, x.globalhdr);
	}

	/* Update EXPSCORR switch upon completion */
	logit = 0;
	if (UpdateSwitch ("EXPSCORR", wf32d->expscorr, x.globalhdr, &logit))
	    return (status);

	putSingleGroup (wf32d->output, extver, &x, option);
	if (hstio_err()) {
	    sprintf (MsgText, "Couldn't write imset %d.", extver);
	    trlerror (MsgText);
	    return (status = 1001);
	}
	if (wf32d->printtime)
	    TimeStamp ("Output written to disk", wf32d->rootname);

	freeSingleGroup (&x);

	return (status);
}

static void DarkMsg (WF3Info *wf32d, int extver) {

	int OmitStep (int);
	void PrSwitch (char *, int);
	void PrRefInfo (char *, char *, char *, char *, char *);

	trlmessage ("");
	PrSwitch ("darkcorr", wf32d->darkcorr);

	if (extver == 1 && !OmitStep (wf32d->darkcorr)) {

	    PrRefInfo ("darkfile", wf32d->dark.name, wf32d->dark.pedigree,
			wf32d->dark.descrip, "");
	}
}

static void dqiMsg (WF3Info *wf32d, int extver) {

	int OmitStep (int);
	void PrSwitch (char *, int);
	void PrRefInfo (char *, char *, char *, char *, char *);

	trlmessage ("");
	PrSwitch ("dqicorr", wf32d->dqicorr);

	if (extver == 1 && !OmitStep (wf32d->dqicorr)) {

	    PrRefInfo ("dqitab", wf32d->bpix.name, wf32d->bpix.pedigree,
			wf32d->bpix.descrip, wf32d->bpix.descrip2);
	}
}

static void FlatMsg (WF3Info *wf32d, int extver) {

	int GotFileName (char *);
	int OmitStep (int);
	void PrSwitch (char *, int);
	void PrRefInfo (char *, char *, char *, char *, char *);

	trlmessage ("");
	PrSwitch ("flatcorr", wf32d->flatcorr);

	if (extver == 1 && !OmitStep (wf32d->flatcorr)) {

	    if (GotFileName (wf32d->pflt.name)) {	/* pixel-to-pixel */
		PrRefInfo ("pfltfile", wf32d->pflt.name,
			wf32d->pflt.pedigree, wf32d->pflt.descrip, "");
	    }
	    if (GotFileName (wf32d->dflt.name)) {		/* delta flat */
		PrRefInfo ("dfltfile", wf32d->dflt.name,
			wf32d->dflt.pedigree, wf32d->dflt.descrip, "");
	    }
	    if (GotFileName (wf32d->lflt.name)) {	/* low-order flat */
		PrRefInfo ("lfltfile", wf32d->lflt.name,
			wf32d->lflt.pedigree, wf32d->lflt.descrip, "");
	    }
	}
}

static void PhotMsg (WF3Info *wf32d) {

	int OmitStep (int);
	void PrRefInfo (char *, char *, char *, char *, char *);

	if (!OmitStep (wf32d->photcorr)) {

	    PrRefInfo ("imphttab", wf32d->phot.name, wf32d->phot.pedigree,
			wf32d->phot.descrip, wf32d->phot.descrip2);
	}
}


static void ShadMsg (WF3Info *wf32d, int extver) {

	int OmitStep (int);
	void PrSwitch (char *, int);
	void PrRefInfo (char *, char *, char *, char *, char *);

	if (wf32d->detector != IR_DETECTOR) {
	    trlmessage ("");
	    PrSwitch ("shadcorr", wf32d->shadcorr);
	}

	if (extver == 1 && !OmitStep (wf32d->shadcorr)) {

	    PrRefInfo ("shadfile", wf32d->shad.name, wf32d->shad.pedigree,
			wf32d->shad.descrip, "");
	}
}


/* This function verifies whether the overscan region has been
	trimmed from the image during BLEVCORR.  If not, it tells the user,
	and shuts down gracefully.
    *****************
    This test needs to be made robust against sub-arrays!!!!!
    *****************
*/

static int OscnTrimmed (Hdr *phdr, Hdr *hdr) {

	extern int status;

	double ltv1, ltv2;
	int blevcorr;

	int GetKeyDbl (Hdr *, char *, int, double, double *);
	int GetSwitch (Hdr *, char *, int *);

	if (GetSwitch (phdr, "BLEVCORR", &blevcorr))
	    return(status);

	if (GetKeyDbl (hdr, "LTV1", USE_DEFAULT, 0., &ltv1))
	    return (status);
	if (GetKeyDbl (hdr, "LTV2", USE_DEFAULT, 0., &ltv2))
	    return (status);

	/* If there is any overscan region still specified in the
		image header, then we can not process this image yet...

	   WFC3 binned images have ltv values greater than 0 even after
           trimming is performed, so test must be for ltv > 1.
	*/
	if (blevcorr != COMPLETE || (ltv1 > 1. || ltv2 > 1.)) {
	    trlerror ("Overscan region not trimmed.  Please perform BLEVCORR.");
	    status = CAL_STEP_NOT_DONE;
	}

	return (status);
}
