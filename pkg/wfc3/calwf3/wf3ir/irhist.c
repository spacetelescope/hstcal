# include <stdio.h>
# include <stdlib.h>	/* calloc */
# include <string.h>

# include "hstio.h"
# include "wf3.h"
# include "wf3info.h"
# include "hstcalerr.h"
# include "trl.h"


/* These routines check whether pedigree indicates that the calibration
   file is dummy, and log history info accordingly.  From one to three
   history records will be logged:  the status of the calibration (i.e.
   "complete" or "skipped"), the value of the pedigree keyword, and the
   value of the descrip keyword.  The calibration switch (e.g. FLATCORR)
   in the header will be updated from PERFORM to either COMPLETE or SKIPPED.

   Howard Bushouse, 2001 Apr 13:
	Initial WFC3 version.

   H.Bushouse, 2002 Apr 12:
	Upgraded flatIRHistory to accomodate 3 types of flats.

   H.Bushouse, 2010 Oct 20:
	Upgraded noisIRHistory to first check setting of noiscorr switch
	before adding history keyword, to support re-entrant processing.
	(PR 66081)
   H.Bushouse, 2011 Sep 7:
	Modified photIRHistory to work with new phot table instead of
	graph and comp tabs.
*/

int blevIRHistory (WF3Info *wf3, Hdr *phdr) {

	extern int status;
	int OmitStep (int);
	int PutKeyStr (Hdr *, char *, char *, char *);
	int TabHistory (RefTab *, Hdr *);

	if (OmitStep (wf3->blevcorr))
	    return (status);

	if (PutKeyStr (phdr, "BLEVCORR", "COMPLETE", ""))
	    return (status);
	addHistoryKw (phdr, "BLEVCORR complete.");
	if (hstio_err())
	    return (status = HEADER_PROBLEM);

	addHistoryKw (phdr, "  Overscan region table: ");
	if (TabHistory (&wf3->oscn, phdr))
	    return (status);
	if (hstio_err())
	    return (status = HEADER_PROBLEM);

	return (status);
}

int crIRHistory (WF3Info *wf3, Hdr *phdr) {

	extern int status;
	int OmitStep (int);
	int PutKeyStr (Hdr *, char *, char *, char *);

	if (OmitStep (wf3->crcorr))
	    return (status);

	if (PutKeyStr (phdr, "CRCORR", "COMPLETE", ""))
	    return (status);
	addHistoryKw (phdr, "CRCORR complete.");
	if (hstio_err())
	    return (status = HEADER_PROBLEM);

	return (status);
}

/* Info about CCD parameters table (no specific calibration step). */

int CCDIRHistory (WF3Info *wf3, Hdr *phdr) {

	extern int status;
	int TabHistory (RefTab *, Hdr *);

	addHistoryKw (phdr, "CCD parameters table: ");
	if (TabHistory (&wf3->ccdpar, phdr))
	    return (status);

	return (status);
}

int dqiIRHistory (WF3Info *wf3, Hdr *phdr) {

	extern int status;
	int logit;			/* true if we should log file name */
	int flag;
	int OmitStep (int);
	int TabHistory (RefTab *, Hdr *);
	int UpdateSwitch (char *, int, Hdr *, int *);

	if (OmitStep (wf3->dqicorr))
	    return (status);

	if (wf3->bpix.exists == EXISTS_YES) {
	    if (wf3->bpix.goodPedigree == GOOD_PEDIGREE)
		flag = PERFORM;
	    else if (wf3->bpix.goodPedigree == DUMMY_PEDIGREE)
		flag = DUMMY;
	    else
		flag = IGNORED;
	} else {
	    flag = PERFORM;
	}

	if (UpdateSwitch ("DQICORR", flag, phdr, &logit))
	    return (status);

	if (logit) {
	    if (wf3->bpix.exists == EXISTS_YES) {
		addHistoryKw (phdr, "  DQ array initialized ...");
		if (TabHistory (&(wf3->bpix), phdr))
		    return (status);
	    }
	}

	return (status);
}


int noisIRHistory (WF3Info *wf3, Hdr *phdr) {

	extern int status;
	int OmitStep (int);

	if (OmitStep (wf3->noiscorr))
	    return (status);

	addHistoryKw (phdr, "Uncertainty array initialized.");
	if (hstio_err())
	    return (status = HEADER_PROBLEM);

	return (status);
}

int darkIRHistory (WF3Info *wf3, Hdr *phdr) {

	extern int status;
	int logit;			/* true if we should log file name */
	int OmitStep (int);
	int ImgHistory (RefImage *, Hdr *);
	int UpdateSwitch (char *, int, Hdr *, int *);

	if (OmitStep (wf3->darkcorr))
	    return (status);

	if (UpdateSwitch ("DARKCORR", wf3->darkcorr, phdr, &logit))
	    return (status);

	if (logit) {
	    if (ImgHistory (&(wf3->dark), phdr))
		return (status);
	}

	return (status);
}

int flatIRHistory (WF3Info *wf3, Hdr *phdr) {

	extern int status;

	int logit;			/* true if we should log file name */
	int OmitStep (int);
	int ImgHistory (RefImage *, Hdr *);
	int GotFileName (char *);
	int UpdateSwitch (char *, int, Hdr *, int *);

	if (OmitStep (wf3->flatcorr))
	    return (status);

	if (UpdateSwitch ("FLATCORR", wf3->flatcorr, phdr, &logit))
	    return (status);

	if (logit) {

	    if (GotFileName (wf3->pflt.name)) {
		if (ImgHistory (&(wf3->pflt), phdr))  /* pixel-to-pixel */
		    return (status);
	    }

	    if (GotFileName (wf3->dflt.name)) {
		if (ImgHistory (&(wf3->dflt), phdr))  /* delta flat */
		    return (status);
	    }

	    if (GotFileName (wf3->lflt.name)) {
		if (ImgHistory (&(wf3->lflt), phdr))  /* low-order flat */
		    return (status);
	    }
	}

	return (status);
}

/* Photometry. */

int photIRHistory (WF3Info *wf3, Hdr *phdr) {

	extern int status;
	int logit;			/* true if we should log file name */
	int OmitStep (int);
	int TabHistory (RefTab *, Hdr *);
	int UpdateSwitch (char *, int, Hdr *, int *);

	if (OmitStep (wf3->photcorr))
	    return (status);

	if (UpdateSwitch ("PHOTCORR", wf3->photcorr, phdr, &logit))
	    return (status);

	if (logit) {
	    if (TabHistory (&wf3->phot, phdr))
		return (status);
	}

	return (status);
}

int nlinIRHistory (WF3Info *wf3, Hdr *phdr) {

	extern int status;
	int logit;			/* true if we should log file name */
	int OmitStep (int);
	int ImgHistory (RefImage *, Hdr *);
	int UpdateSwitch (char *, int, Hdr *, int *);

	if (OmitStep (wf3->nlincorr))
	    return (status);

	if (UpdateSwitch ("NLINCORR", wf3->nlincorr, phdr, &logit))
	    return (status);

	if (logit) {
	    if (ImgHistory (&(wf3->nlin), phdr))
		return (status);
	}

	return (status);
}

int unitIRHistory (WF3Info *wf3, Hdr *phdr) {

	extern int status;
	int OmitStep (int);
	int PutKeyStr (Hdr *, char *, char *, char *);

	if (OmitStep (wf3->unitcorr))
	    return (status);

	if (PutKeyStr (phdr, "UNITCORR", "COMPLETE", ""))
	    return (status);
	addHistoryKw (phdr, "UNITCORR complete.");
	if (hstio_err())
	    return (status = HEADER_PROBLEM);

	return (status);
}

int zoffIRHistory (WF3Info *wf3, Hdr *phdr) {

	extern int status;
	int OmitStep (int);
	int PutKeyStr (Hdr *, char *, char *, char *);

	if (OmitStep (wf3->zoffcorr))
	    return (status);

	if (PutKeyStr (phdr, "ZOFFCORR", "COMPLETE", ""))
	    return (status);
	addHistoryKw (phdr, "ZOFFCORR complete.");
	if (hstio_err())
	    return (status = HEADER_PROBLEM);

	return (status);
}

int zsigIRHistory (WF3Info *wf3, Hdr *phdr) {

	extern int status;
	int OmitStep (int);
	int PutKeyStr (Hdr *, char *, char *, char *);

	if (OmitStep (wf3->zsigcorr))
	    return (status);

	if (PutKeyStr (phdr, "ZSIGCORR", "COMPLETE", ""))
	    return (status);
	addHistoryKw (phdr, "ZSIGCORR complete.");
	if (hstio_err())
	    return (status = HEADER_PROBLEM);

	return (status);
}

