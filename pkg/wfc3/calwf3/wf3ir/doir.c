# include <stdio.h>
# include <string.h>

# include "hstio.h"	/* defines HST I/O functions */
# include "wf3.h"
# include "wf3info.h"
# include "err.h"
# include "trl.h"

extern int status;

/* DOIR: Applies IR calibration steps to an input science data file.
** All calibration steps modify the input image data in-place.
**
** Revision history:
** H.Bushouse	 2-Oct-2000	Initial version copied from CALNICA n_doCalib
**				module and adapted for CALWF3 processing.
** H.Bushouse	24-May-2001	Reordered steps so that doDQIIR is done first;
**				eliminated passing of mask image to doZsig.
** H.Bushouse	16-Nov-2001	Modified photMsg routine to use GRAPHTAB and
**				COMPTAB, instead of PHOTTAB.
** H.Bushouse	10-Apr-2002	Added wf3info to arg list of statcalc in order
**				to be able to skip reference pixels.
** H.Bushouse	03-May-2002	Moved doBlevIR to just before doZsigIR. Also
**				moved copyGroup for zoff image to be just
**				after doBlevIR.
** H.Bushouse	27-Aug-2008	Added call to GetGrp to load LTV offset values.
** H.Bushouse	13-Jan-2010	Change order of processing so that doBlev is
**				after doZsig. Also requires passing zoff image
**				into doBlev so that it can get processed.
** H.Bushouse	14-Jan-2010	Compute zero-read sample time (sampzero) here
**				instead of in zsigcorr routine. (calwf3 v2.0)
** H.Bushouse	27-Apr-2010	Reorder processing steps so that nlincorr is
**				performed before darkcorr. (calwf3 v2.1)
** H.Bushouse	26-Oct-2010	Upgraded crimage header updates to check for
**				flatcorr status when setting BUNIT value. Also
**				modified noisMsg routine to print noiscorr
**				switch value. Both are to support re-entrant
**				processing. (PR 66081)
** H.Bushouse	11-Mar-2011	No longer load dark ref file for zsigcorr.
**				(PR #67728, Trac #681)
** H.Bushouse	 7-Sep-2011	Modified photMsg to work with new imphttab
**				instead of graph and comp tabs.
*/

static void dqiMsg  (WF3Info *);
static void zsigMsg (WF3Info *);
static void zoffMsg (WF3Info *);
static void noisMsg (WF3Info *);
static void darkMsg (WF3Info *);
static void blevMsg (WF3Info *);
static void nlinMsg (WF3Info *);
static void crejMsg (WF3Info *);
static void flatMsg (WF3Info *);
static void photMsg (WF3Info *);
static void unitMsg (WF3Info *);

int DoIR (WF3Info *wf3, MultiNicmosGroup *input, SingleNicmosGroup *crimage) {

/* Arguments:
**	wf3	 i: WF3 info structure
**	input	io: input image data
**	crimage	 o: CR rejected image
*/

	/* Local variables */
	int i;				/* loop index */
	int overscan;			/* indicates whether overscan done */
	char buff[SZ_FITS_REC+1];
	Bool subarray;
	SingleNicmosGroup zoff;		/* original zero-read image */
	static SingleNicmosGroup zsig;	/* zero-read signal image */

	/* Function definitions */
	int getDarkInfo (WF3Info *);
	int copyGroup (SingleNicmosGroup *, SingleNicmosGroup *);
	int GetGrp (WF3Info *, Hdr *);
	int GetCCDTab (WF3Info *, int, int);
	int FindOverscan (WF3Info *, int, int, int *);
	int PutKeyFlt (Hdr *, char *, float, char *);
	int PutKeyStr (Hdr *, char *, char *, char *);
	int GetKeyBool (Hdr *, char *, int, Bool, Bool *);
	void PrRefInfo (char *, char *, char *, char *, char *);

	int doDQIIR  (WF3Info *, MultiNicmosGroup *);
	int doBlevIR (WF3Info *, MultiNicmosGroup *, SingleNicmosGroup *);
	int doZsigIR (WF3Info *, MultiNicmosGroup *, SingleNicmosGroup *);
	int doZoffIR (WF3Info *, MultiNicmosGroup *, SingleNicmosGroup *);
	int doNoisIR (WF3Info *, MultiNicmosGroup *);
	int doDarkIR (WF3Info *, MultiNicmosGroup *);
	int doNlinIR (WF3Info *, MultiNicmosGroup *, SingleNicmosGroup *);
	int doFlatIR (WF3Info *, MultiNicmosGroup *, SingleNicmosGroup *);
	int doUnitIR (WF3Info *, MultiNicmosGroup *);
        int photcalc (WF3Info *, MultiNicmosGroup *);
	int cridcalc (WF3Info *, MultiNicmosGroup *, SingleNicmosGroup *);
	int statcalc (WF3Info *, SingleNicmosGroup *, short);

	int PhotMode (WF3Info *, Hdr *);

	int CCDHistory    (WF3Info *, Hdr *);
	int dqiIRHistory  (WF3Info *, Hdr *);
	int zsigIRHistory (WF3Info *, Hdr *);
	int zoffIRHistory (WF3Info *, Hdr *);
	int noisIRHistory (WF3Info *, Hdr *);
	int darkIRHistory (WF3Info *, Hdr *);
	int blevIRHistory (WF3Info *, Hdr *);
	int nlinIRHistory (WF3Info *, Hdr *);
	int crIRHistory   (WF3Info *, Hdr *);
	int flatIRHistory (WF3Info *, Hdr *);
	int photIRHistory (WF3Info *, Hdr *);
	int unitIRHistory (WF3Info *, Hdr *);

	/* Get header info */
	if (GetGrp (wf3, &input->group[0].sci.hdr))
	    return (status);

	/* Read the SUBARRAY keyword from the input image header */
	if (GetKeyBool (input->group[0].globalhdr, "SUBARRAY", NO_DEFAULT, 0,
			&subarray))
	    return (status);
	if (subarray)
	    wf3->subarray = YES;
	else
	    wf3->subarray = NO;

	/* Load CCD reference table using same function as in WF3CCD. */
	if (GetCCDTab (wf3, input->group[0].sci.data.nx,
			    input->group[0].sci.data.ny))
	    return (status);

	/* Update gain and readnoise keyword values in primary header */
	if (PutKeyFlt(input->group[0].globalhdr,"ATODGNA",wf3->atodgain[0],""))
	    return (status);
	if (PutKeyFlt(input->group[0].globalhdr,"ATODGNB",wf3->atodgain[1],""))
	    return (status);
	if (PutKeyFlt(input->group[0].globalhdr,"ATODGNC",wf3->atodgain[2],""))
	    return (status);
	if (PutKeyFlt(input->group[0].globalhdr,"ATODGND",wf3->atodgain[3],""))
	    return (status);
	if (PutKeyFlt(input->group[0].globalhdr,"READNSEA",wf3->readnoise[0],""))
	    return (status);
	if (PutKeyFlt(input->group[0].globalhdr,"READNSEB",wf3->readnoise[1],""))
	    return (status);
	if (PutKeyFlt(input->group[0].globalhdr,"READNSEC",wf3->readnoise[2],""))
	    return (status);
	if (PutKeyFlt(input->group[0].globalhdr,"READNSED",wf3->readnoise[3],""))
	    return (status);
	trlmessage ("");

	PrRefInfo ("ccdtab", wf3->ccdpar.name, wf3->ccdpar.pedigree,
		   wf3->ccdpar.descrip, wf3->ccdpar.descrip2);

        if (CCDHistory (wf3, input->group[0].globalhdr))
            return (status);

	/* Get overscan region information from OSCNTAB */
	if (FindOverscan (wf3, input->group[0].sci.data.nx,
			  input->group[0].sci.data.ny, &overscan))
	    return (status);

	/* Report readnoise and gain values */
	buff[0] = '\0';
	sprintf (MsgText, "    readnoise =");
	for (i=0; i < NAMPS-1; i++) {
	     if (wf3->readnoise[i] > 0) {
		 sprintf (buff, "%.5g,",wf3->readnoise[i]);
		 strcat (MsgText, buff);
	     }
	}
	if (wf3->readnoise[NAMPS-1] > 0) {
	    sprintf (buff, "%.5g",wf3->readnoise[NAMPS-1]);
	    strcat (MsgText, buff);
	}
	trlmessage (MsgText);

	sprintf (MsgText, "    gain =");
	for (i=0; i < NAMPS-1; i++) {
	     if (wf3->atodgain[i] > 0) {
		 sprintf (buff, "%.5g,",wf3->atodgain[i]);
		 strcat (MsgText, buff);
	     }
	}
	if (wf3->atodgain[NAMPS-1] > 0) {
	    sprintf (buff, "%.5g",wf3->atodgain[NAMPS-1]);
	    strcat (MsgText, buff);
	}
	trlmessage (MsgText);

	/* Load the DARK ref image information */
	if (wf3->darkcorr == PERFORM) {
	    if (getDarkInfo (wf3))
		return (status);
	}

	/* Compute the zero-read sample time, for use later:
	** We compute the effective exposure time of the zeroth read by
	** taking the exposure time of the first read and subtracting the
	** overheads that don't apply to the zeroth read, which are:
	** 20.48 msec TSM expose time
	**    55 usec ADC warm up time
	** (H.Bushouse 01-Dec-2008) */
	wf3->sampzero = wf3->exptime[wf3->ngroups-2] - 0.020535;

	/* Do DQ array initialization */
	dqiMsg (wf3);
	if (doDQIIR (wf3, input))
	    return (status);
	if (dqiIRHistory (wf3, input->group[0].globalhdr))
	    return (status);

	/* Save original zeroth-read image before it gets
	** modified by zsigcorr */
	if (copyGroup (&zoff, &(input->group[wf3->ngroups-1])))
	    return (status);

	/* Do MultiAccum zero-read signal correction */
	zsigMsg (wf3);
	if (doZsigIR (wf3, input, &zsig))
	    return (status);
	if (zsigIRHistory (wf3, input->group[0].globalhdr))
	    return (status);

	/* Do reference pixel subtraction */
	blevMsg (wf3);
	if (doBlevIR (wf3, input, &zoff))
	    return (status);
	if (blevIRHistory (wf3, input->group[0].globalhdr))
	    return (status);

	/* Do MultiAccum zero-read subtraction */
	zoffMsg (wf3);
	if (doZoffIR (wf3, input, &zoff))
	    return (status);
	if (zoffIRHistory (wf3, input->group[0].globalhdr))
	    return (status);
	freeSingleNicmosGroup (&zoff);

	/* Do noise (error) calculation */
	noisMsg (wf3);
	if (doNoisIR (wf3, input))
	    return (status);
	if (noisIRHistory (wf3, input->group[0].globalhdr))
	    return (status);

	/* Do linearity correction */
	nlinMsg (wf3);
	if (doNlinIR (wf3, input, &zsig))
	    return (status);
	if (nlinIRHistory (wf3, input->group[0].globalhdr))
	    return (status);

	/* Do dark subtraction */
	darkMsg (wf3);
	if (doDarkIR (wf3, input))
	    return (status);
	if (darkIRHistory (wf3, input->group[0].globalhdr))
	    return (status);

	/* Free the zsig image data */
	if (wf3->zsigcorr == PERFORM)
	    freeSingleNicmosGroup (&zsig);

        /* Do photometric correction */
	photMsg (wf3);

	/* First, update the PHOTMODE keyword */
	if (PhotMode (wf3, input->group[0].globalhdr))
	    return (status);

	/* Now compute values for the photometry keywords */
        if (photcalc (wf3, input))
            return (status);
	if (photIRHistory (wf3, input->group[0].globalhdr))
	    return (status);

	/* Do units conversion */
	unitMsg (wf3);
	if (doUnitIR (wf3, input))
	    return (status);
	if (unitIRHistory (wf3, input->group[0].globalhdr))
	    return (status);

	/* Do cosmic-ray identification and rejection */
	crejMsg (wf3);
	if (cridcalc (wf3, input, crimage))
	    return (status);
	if (wf3->crcorr == PERFORM) {
	    if (crIRHistory (wf3, crimage->globalhdr))
		return (status);
	}
	if (crIRHistory (wf3, input->group[0].globalhdr))
	    return (status);

	/* Update crimage headers */
	if (wf3->crcorr == PERFORM && wf3->unitcorr == PERFORM) {
	    if (unitIRHistory (wf3, crimage->globalhdr))
		return (status);
	    if (hstio_err())
		return (status = HEADER_PROBLEM);
	    if (wf3->flatcorr == COMPLETE) {
		if (PutKeyStr (&crimage->sci.hdr, "BUNIT", "ELECTRONS/S", ""))
		    return (status);
		if (PutKeyStr (&crimage->err.hdr, "BUNIT", "ELECTRONS/S", ""))
		    return (status);
	    } else {
		if (PutKeyStr (&crimage->sci.hdr, "BUNIT", "COUNTS/S", ""))
		    return (status);
		if (PutKeyStr (&crimage->err.hdr, "BUNIT", "COUNTS/S", ""))
		    return (status);
	    }
	}

	/* Do flat fielding */
	flatMsg (wf3);
	if (doFlatIR (wf3, input, crimage))
	    return (status);
	if (flatIRHistory (wf3, input->group[0].globalhdr))
	    return (status);
	if (wf3->crcorr == PERFORM) {
	    if (flatIRHistory (wf3, crimage->globalhdr))
		return (status);
	}

	/* Compute statistics for calibrated data */
	for (i=0; i < wf3->ngroups; i++) {
	     if (statcalc (wf3, &(input->group[i]), wf3->sdqflags))
		 return (status);
	}
	if (wf3->crcorr == PERFORM) {
	    if (statcalc (wf3, crimage, wf3->sdqflags))
		return (status);
	}

	/* Successful return */
	return (status = 0);
}

static void zsigMsg (WF3Info *wf3) {

	void PrSwitch (char *, int);

	trlmessage ("");
	PrSwitch ("zsigcorr", wf3->zsigcorr);

}

static void zoffMsg (WF3Info *wf3) {

	void PrSwitch (char *, int);

	trlmessage ("");
	PrSwitch ("zoffcorr", wf3->zoffcorr);

}

static void dqiMsg (WF3Info *wf3) {

	int OmitStep (int);
	void PrSwitch (char *, int);
	void PrRefInfo (char *, char *, char *, char *, char *);

	trlmessage ("");
	PrSwitch ("dqicorr", wf3->dqicorr);

	if (!OmitStep (wf3->dqicorr)) {
	    PrRefInfo ("dqitab", wf3->bpix.name, wf3->bpix.pedigree,
			wf3->bpix.descrip, "");
	}
}

static void darkMsg (WF3Info *wf3) {

	int OmitStep (int);
	void PrSwitch (char *, int);
	void PrRefInfo (char *, char *, char *, char *, char *);

	trlmessage ("");
	PrSwitch ("darkcorr", wf3->darkcorr);

	if (!OmitStep (wf3->darkcorr)) {
	    PrRefInfo ("darkfile", wf3->dark.name, wf3->dark.pedigree,
			wf3->dark.descrip, "");
	}
}

static void blevMsg (WF3Info *wf3) {

	void PrSwitch (char *, int);
	void PrRefInfo (char *, char *, char *, char *, char *);
	int OmitStep (int);

	trlmessage ("");
	PrSwitch ("blevcorr", wf3->blevcorr);

	if (!OmitStep (wf3->blevcorr)) {
	    PrRefInfo ("oscntab", wf3->oscn.name, wf3->oscn.pedigree,
		       wf3->oscn.descrip, wf3->oscn.descrip2);
	}
}

static void noisMsg (WF3Info *wf3) {

	void PrSwitch (char *, int);

	trlmessage ("");
	PrSwitch ("noiscorr", wf3->noiscorr);
}

static void nlinMsg (WF3Info *wf3) {

	int OmitStep (int);
	void PrSwitch (char *, int);
	void PrRefInfo (char *, char *, char *, char *, char *);

	trlmessage ("");
	PrSwitch ("nlincorr", wf3->nlincorr);

	if (!OmitStep (wf3->nlincorr)) {
	    PrRefInfo ("nlinfile", wf3->nlin.name, wf3->nlin.pedigree,
			wf3->nlin.descrip, "");
	}
}

static void crejMsg (WF3Info *wf3) {

	void PrSwitch (char *, int);

	trlmessage ("");
	PrSwitch ("crcorr", wf3->crcorr);

}

static void flatMsg (WF3Info *wf3) {

	int GotFileName (char *);
	int OmitStep (int);
	void PrSwitch (char *, int);
	void PrRefInfo (char *, char *, char *, char *, char *);

	trlmessage ("");
	PrSwitch ("flatcorr", wf3->flatcorr);

	if (!OmitStep (wf3->flatcorr)) {

	    if (GotFileName (wf3->pflt.name)) {
		PrRefInfo ("pfltfile", wf3->pflt.name,
			   wf3->pflt.pedigree, wf3->pflt.descrip, "");
	    }

	    if (GotFileName (wf3->dflt.name)) {
		PrRefInfo ("dfltfile", wf3->dflt.name,
			   wf3->dflt.pedigree, wf3->dflt.descrip, "");
	    }

	    if (GotFileName (wf3->lflt.name)) {
		PrRefInfo ("lfltfile", wf3->lflt.name,
			   wf3->lflt.pedigree, wf3->lflt.descrip, "");
	    }
	}
}

static void photMsg (WF3Info *wf3) {

	int OmitStep (int);
	void PrSwitch (char *, int);
	void PrRefInfo (char *, char *, char *, char *, char *);

	trlmessage ("");
	PrSwitch ("photcorr", wf3->photcorr);

	if (!OmitStep (wf3->photcorr)) {

	    PrRefInfo ("imphttab", wf3->phot.name, wf3->phot.pedigree,
			wf3->phot.descrip, wf3->phot.descrip2);
	}
}

static void unitMsg (WF3Info *wf3) {

	void PrSwitch (char *, int);

	trlmessage ("");
	PrSwitch ("unitcorr", wf3->unitcorr);

}

