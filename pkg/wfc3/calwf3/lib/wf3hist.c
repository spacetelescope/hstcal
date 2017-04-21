# include <stdio.h>
# include <stdlib.h>	/* calloc */
# include <string.h>
# include "hstio.h"
# include "wf3.h"
# include "wf3info.h"
# include "hstcalerr.h"


/* These routines check whether pedigree indicates that the calibration
   file is dummy, and log history info accordingly.  From one to three
   history records will be logged:  the status of the calibration (i.e.
   "complete" or "skipped"), the value of the pedigree keyword, and the
   value of the descrip keyword.  The calibration switch (e.g. FLATCORR)
   in the header will be updated from PERFORM to either COMPLETE or SKIPPED.

   Note that since this info is written to the primary header, these
   routines should only be called for EXTVER = 1.

   Warren Hack, 1998 June 10:
   	Initial ACS version.
	
   Howard Bushouse, 2000 Aug 24:
	Initial WFC3 version.

   H.Bushouse, 2001 May 8:
	Added Post-Flash processing support.

   H.Bushouse, 2001 Nov 16:
	Modified photHistory to write GRAPHTAB and COMPTAB info, instead
	of PHOTTAB.

   H.Bushouse, 2002 June 20:
	Added history reporting on overscan table.

   H.Bushouse, 2011 Sep 7:
	Replaced graph and comp tab with new phot tab in photHistory.
*/

int atodHistory (WF3Info *wf3, Hdr *phdr) {

	extern int status;

	int logit;			/* true if we should log file name */
	int OmitStep (int);
	int TabHistory (RefTab *, Hdr *);
	int UpdateSwitch (char *, int, Hdr *, int *);

	if (OmitStep (wf3->atodcorr))		/* nothing to do */
	    return (status);

	if (UpdateSwitch ("ATODCORR", wf3->atodcorr, phdr, &logit))
	    return (status);

	/* Write history records for the A-to-D table. */
	if (logit) {
	    if (TabHistory (&wf3->atod, phdr))
		return (status);
	}

	return (status);
}

int biasHistory (WF3Info *wf3, Hdr *phdr) {

	extern int status;

	int logit;			/* true if we should log file name */
	int OmitStep (int);
	int ImgHistory (RefImage *, Hdr *);
	int UpdateSwitch (char *, int, Hdr *, int *);

	if (OmitStep (wf3->biascorr))
	    return (status);

	if (UpdateSwitch ("BIASCORR", wf3->biascorr, phdr, &logit))
	    return (status);

	/* Write history records for the bias image. */
	if (logit) {
	    if (ImgHistory (&wf3->bias, phdr))
		return (status);
	}

	return (status);
}

int flashHistory (WF3Info *wf3, Hdr *phdr) {

	extern int status;

	int logit;			/* true if we should log file name */
	int OmitStep (int);
	int ImgHistory (RefImage *, Hdr *);
	int UpdateSwitch (char *, int, Hdr *, int *);

	if (OmitStep (wf3->flashcorr))
	    return (status);

	if (UpdateSwitch ("FLSHCORR", wf3->flashcorr, phdr, &logit))
	    return (status);

	/* Write history records for the flash image. */
	if (logit) {
	    if (ImgHistory (&wf3->flash, phdr))
		return (status);
	}

	return (status);
}

int blevHistory (WF3Info *wf3, Hdr *phdr, int done, int driftcorr) {

	extern int status;
	int OmitStep (int);
	int PutKeyStr (Hdr *, char *, char *, char *);
	int TabHistory (RefTab *, Hdr *);

	if (OmitStep (wf3->blevcorr))
	    return (status);

	if (PutKeyStr (phdr, "BLEVCORR", "COMPLETE", ""))
	    return (status);
	if (done) {
	    addHistoryKw (phdr,
	"BLEVCORR complete; bias level from overscan was subtracted.");
	} else {
	    addHistoryKw (phdr,
	"BLEVCORR complete, but default bias level was subtracted.");
	}
	if (hstio_err())
	    return (status = HEADER_PROBLEM);

	if (driftcorr) {
	    addHistoryKw (phdr,
	"BLEVCORR includes correction for drift along lines.");
	} else {
	    addHistoryKw (phdr,
	"BLEVCORR does not include correction for drift along lines.");
	}

	addHistoryKw (phdr, "  Overscan region table: ");
	if (TabHistory (&wf3->oscn, phdr))
	    return (status);

	if (hstio_err())
	    return (status = HEADER_PROBLEM);

	return (status);
}

/* Info about CCD parameters table (no specific calibration step). */

int CCDHistory (WF3Info *wf3, Hdr *phdr) {

	extern int status;
	int TabHistory (RefTab *, Hdr *);

	addHistoryKw (phdr, "CCD parameters table: ");
	if (TabHistory (&wf3->ccdpar, phdr))
	    return (status);

	return (status);
}


int dqiHistory (WF3Info *wf3, Hdr *phdr) {

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
        if (wf3->detector != IR_DETECTOR) {
	        addHistoryKw (phdr, "  values checked for saturation");
        }
        
	    if (wf3->bpix.exists == EXISTS_YES) {
		    addHistoryKw (phdr, "  DQ array initialized ...");
		    if (TabHistory (&wf3->bpix, phdr))
		        return (status);
	    }
	}

	return (status);
}

/* There are no reference images for assigning initial error values. */

int noiseHistory (Hdr *phdr) {

	extern int status;

	addHistoryKw (phdr, "Uncertainty array initialized.");
	if (hstio_err())
	    return (status = HEADER_PROBLEM);

	return (status);
}

int darkHistory (WF3Info *wf3, Hdr *phdr) {

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
	    if (ImgHistory (&wf3->dark, phdr))
		return (status);
	}

	return (status);
}

int flatHistory (WF3Info *wf3, Hdr *phdr) {

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
		if (ImgHistory (&wf3->pflt, phdr))	/* pixel-to-pixel */
		    return (status);
	    }

	    if (GotFileName (wf3->dflt.name)) {
		if (ImgHistory (&wf3->dflt, phdr))	/* delta flat */
		    return (status);
	    }

	    if (GotFileName (wf3->lflt.name)) {
		if (ImgHistory (&wf3->lflt, phdr))	/* low-order flat */
		    return (status);
	    }
	}

	return (status);
}

int shadHistory (WF3Info *wf3, Hdr *phdr) {

	extern int status;
	int logit;			/* true if we should log file name */
	int OmitStep (int);
	int ImgHistory (RefImage *, Hdr *);
	int UpdateSwitch (char *, int, Hdr *, int *);

	if (OmitStep (wf3->shadcorr))
	    return (status);

	if (UpdateSwitch ("SHADCORR", wf3->shadcorr, phdr, &logit))
	    return (status);

	if (logit) {
	    if (ImgHistory (&wf3->shad, phdr))
		return (status);
	}

	return (status);
}

/* Photometry. */

int photHistory (WF3Info *wf3, Hdr *phdr) {

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

int fluxHistory (WF3Info *wf3, Hdr *phdr) {

	extern int status;
	int logit;			/* true if we should log file name */
	int OmitStep (int);
	int TabHistory (RefTab *, Hdr *);
	int UpdateSwitch (char *, int, Hdr *, int *);

	if (OmitStep (wf3->fluxcorr))
	    return (status);

	if (UpdateSwitch ("FLUXCORR", wf3->fluxcorr, phdr, &logit))
	    return (status);

	return (status);
}

int cteHistory (WF3Info *wf3, Hdr *phdr) {
	extern int status;
	int logit;			/* true if we should log file name */
	int OmitStep (int);
	int TabHistory (RefTab *, Hdr *);
	int UpdateSwitch (char *, int, Hdr *, int *);

	if (OmitStep (wf3->pctecorr))
	    return (status);

	if (UpdateSwitch ("PCTECORR", wf3->pctecorr, phdr, &logit))
	    return (status);

	addHistoryKw (phdr, "CTE parameters table: ");
	if (TabHistory (&wf3->pctetab, phdr))
	    return (status);

	return (status);
}


/* Update the calibration switch in the primary header to COMPLETE
   or SKIPPED.
*/

int UpdateSwitch (char *calSwitch, int flag, Hdr *phdr, int *logit) {

/* arguments:
char *calSwitch   i: name of calibration switch
int flag          i: value of calibration switch
Hdr *phdr         o: primary header
int *logit        o: true if we should log reference file names
*/

	extern int status;

	char *history;
	int PutKeyStr (Hdr *, char *, char *, char *);

	if ((history = (char *) calloc (SZ_LINE+1, sizeof (char))) == NULL)
	    return (status = OUT_OF_MEMORY);

	strcpy (history, calSwitch);

	*logit = 0;
	if (flag == PERFORM) {
	    if (PutKeyStr (phdr, calSwitch, "COMPLETE", ""))
		return (status);
	    strcat (history, " complete ...");
	    addHistoryKw (phdr, history);
	    if (hstio_err())
		return (status = HEADER_PROBLEM);
	    *logit = 1;
	} else if (flag == DUMMY) {
	    if (PutKeyStr (phdr, calSwitch, "SKIPPED", ""))
		return (status);
	    strcat (history, " skipped due to dummy reference file ...");
	    addHistoryKw (phdr,history);
	    if (hstio_err())
		return (status = HEADER_PROBLEM);
	    *logit = 1;
	} else if (flag == IGNORED) {
	    if (PutKeyStr (phdr, calSwitch, "SKIPPED", ""))
		return (status);
	    strcat (history, " not performed ...");  /* for some other reason */
	    addHistoryKw (phdr,history);
	    if (hstio_err())
		return (status = HEADER_PROBLEM);
	    *logit = 1;
	}

	free (history);

	return (status);
}
