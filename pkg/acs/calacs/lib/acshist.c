# include <stdio.h>
# include <stdlib.h>	/* calloc */
# include <string.h>
# include "hstio.h"
# include "acs.h"
# include "acsinfo.h"
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
 Warren Hack, 2000 Sept 12:
 Added Post-FLASH processing support.
 Warren Hack, 2002 March 18:
 Added history reporting on overscan table.
 Pey Lian Lim, 2012 Dec 12:
 Added support for CTE corrected post-flash.
 
 */

int atodHistory (ACSInfo *acs, Hdr *phdr) {
  
	extern int status;
  
	int logit;			/* true if we should log file name */
	int OmitStep (int);
	int TabHistory (RefTab *, Hdr *);
	int UpdateSwitch (char *, int, Hdr *, int *);
  
	if (OmitStep (acs->atodcorr))		/* nothing to do */
    return (status);
  
	if (UpdateSwitch ("ATODCORR", acs->atodcorr, phdr, &logit))
    return (status);
  
	/* Write history records for the A-to-D table. */
	if (logit) {
    if (TabHistory (&acs->atod, phdr))
      return (status);
	}
  
	return (status);
}

int pcteHistory (ACSInfo *acs, Hdr *phdr) {
  
	extern int status;
  
	int logit;			/* true if we should log file name */
	int OmitStep (int);
	int TabHistory (RefTab *, Hdr *);
	int UpdateSwitch (char *, int, Hdr *, int *);
  
	if (OmitStep (acs->pctecorr))		/* nothing to do */
    return (status);
  
	if (UpdateSwitch ("PCTECORR", acs->pctecorr, phdr, &logit))
    return (status);
  
	/* Write history records for the PCTE table. */
	if (logit) {
    if (TabHistory (&acs->pcte, phdr))
      return (status);
	}
  
	return (status);
}

int biasHistory (ACSInfo *acs, Hdr *phdr) {
  
	extern int status;
  
	int logit;			/* true if we should log file name */
	int OmitStep (int);
	int ImgHistory (RefImage *, Hdr *);
	int UpdateSwitch (char *, int, Hdr *, int *);
  
	if (OmitStep (acs->biascorr))
    return (status);
  
	if (UpdateSwitch ("BIASCORR", acs->biascorr, phdr, &logit))
    return (status);
  
	/* Write history records for the bias image. */
	if (logit) {
    if (ImgHistory (&acs->bias, phdr))
      return (status);
	}
  
	return (status);
}

int flashHistory (ACSInfo *acs, Hdr *phdr) {

    extern int status;

    int logit;  /* true if we should log file name */
    int OmitStep (int);
    int ImgHistory (RefImage *, Hdr *);
    int UpdateSwitch (char *, int, Hdr *, int *);
  
    if (OmitStep (acs->flashcorr))
        return (status);

    if (UpdateSwitch ("FLSHCORR", acs->flashcorr, phdr, &logit))
        return (status);

    if (logit) {
        if (acs->pctecorr == PERFORM) {
            if (ImgHistory (&acs->flashcte, phdr))
                return (status);
        } else {
            if (ImgHistory (&acs->flash, phdr))
                return (status);
	}
    }

    return (status);
}

int blevHistory (ACSInfo *acs, Hdr *phdr, int done, int driftcorr) {
  
	extern int status;
	int OmitStep (int);
	int PutKeyStr (Hdr *, char *, char *, char *);
	int TabHistory (RefTab *, Hdr *);
  
	if (OmitStep (acs->blevcorr))
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
	if (TabHistory (&acs->oscn, phdr))
		return (status);
  
	if (hstio_err())
    return (status = HEADER_PROBLEM);
  
	return (status);
}

/* Info about CCD parameters table (no specific calibration step). */

int CCDHistory (ACSInfo *acs, Hdr *phdr) {
  
	extern int status;
	int TabHistory (RefTab *, Hdr *);
  
	addHistoryKw (phdr, "CCD parameters table: ");
	if (TabHistory (&acs->ccdpar, phdr))
    return (status);
  
	return (status);
}

int dqiHistory (ACSInfo *acs, Hdr *phdr) {
  
	extern int status;
	int logit;			/* true if we should log file name */
	int flag;
	int OmitStep (int);
	int TabHistory (RefTab *, Hdr *);
	int UpdateSwitch (char *, int, Hdr *, int *);
  
	if (OmitStep (acs->dqicorr))
    return (status);
  
	if (acs->bpix.exists == EXISTS_YES) {
    if (acs->bpix.goodPedigree == GOOD_PEDIGREE)
      flag = PERFORM;
    else if (acs->bpix.goodPedigree == DUMMY_PEDIGREE)
      flag = DUMMY;
    else
      flag = IGNORED;
	} else {
    flag = PERFORM;
	}
  
	if (UpdateSwitch ("DQICORR", flag, phdr, &logit))
    return (status);
  
	if (logit) {
    if (acs->detector != MAMA_DETECTOR) {
      addHistoryKw (phdr, "  values checked for saturation");
    }
    
    if (acs->bpix.exists == EXISTS_YES) {
      addHistoryKw (phdr, "  DQ array initialized ...");
      if (TabHistory (&acs->bpix, phdr))
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

int darkHistory (ACSInfo *acs, Hdr *phdr) {
  
	extern int status;
	int logit;			/* true if we should log file name */
	int OmitStep (int);
	int ImgHistory (RefImage *, Hdr *);
	int UpdateSwitch (char *, int, Hdr *, int *);
  
	if (OmitStep (acs->darkcorr))
    return (status);
  
	if (UpdateSwitch ("DARKCORR", acs->darkcorr, phdr, &logit))
    return (status);
  
	if (logit && acs->pctecorr == PERFORM) {
    if (ImgHistory (&acs->darkcte, phdr))
      return (status);
  } else if (logit) {
    if (ImgHistory (&acs->dark, phdr))
      return (status);
	}
  
	return (status);
}

int flatHistory (ACSInfo *acs, Hdr *phdr) {
  
	extern int status;
  
	int logit;			/* true if we should log file name */
	int OmitStep (int);
	int ImgHistory (RefImage *, Hdr *);
	int TabHistory (RefTab *, Hdr *);
	int GotFileName (char *);
	int UpdateSwitch (char *, int, Hdr *, int *);
  
	if (OmitStep (acs->flatcorr))
    return (status);
  
	if (UpdateSwitch ("FLATCORR", acs->flatcorr, phdr, &logit))
    return (status);
  
	if (logit) {
    
    if (GotFileName (acs->pflt.name)) {
      if (ImgHistory (&acs->pflt, phdr))	/* pixel-to-pixel */
		    return (status);
    }
    
    if (GotFileName (acs->dflt.name)) {
      if (ImgHistory (&acs->dflt, phdr))	/* delta flat */
		    return (status);
    }
    
    if (GotFileName (acs->lflt.name)) {
      if (ImgHistory (&acs->lflt, phdr))	/* low-order flat */
		    return (status);
    }
    if (GotFileName (acs->cflt.name)) {
      if (ImgHistory (&acs->cflt, phdr))	/* coronographic flat */
		    return (status);
	    if (TabHistory (&acs->spot, phdr)) /* spot reference table */
        return (status);
    }
	}
  
	return (status);
}

int nonlinHistory (ACSInfo *acs, Hdr *phdr) {
  
	extern int status;
	int logglin, loglfl;			/* true if we should log file name */
	int OmitStep (int);
	int TabHistory (RefTab *, Hdr *);
	int UpdateSwitch (char *, int, Hdr *, int *);
  
  /* Flags to record whether GLINCORR and LFLGCORR was performed */
	logglin = 0;
  loglfl = 0;
  
	if (!OmitStep (acs->glincorr))
    if (UpdateSwitch ("GLINCORR", acs->glincorr, phdr, &logglin))
      return (status);
  
	if (!OmitStep (acs->lflgcorr))
    if (UpdateSwitch ("LFLGCORR", acs->lflgcorr, phdr, &loglfl))
      return (status);
  
  /* If either was performed, then record in History comments 
   WJH  27 July 1999
   */
	if (logglin || loglfl)
    if (TabHistory (&acs->mlin, phdr))
      return (status);
  
	return (status);
}

int shadHistory (ACSInfo *acs, Hdr *phdr) {
  
	extern int status;
	int logit;			/* true if we should log file name */
	int OmitStep (int);
	int ImgHistory (RefImage *, Hdr *);
	int UpdateSwitch (char *, int, Hdr *, int *);
  
	if (OmitStep (acs->shadcorr))
    return (status);
  
	if (UpdateSwitch ("SHADCORR", acs->shadcorr, phdr, &logit))
    return (status);
  
	if (logit) {
    if (ImgHistory (&acs->shad, phdr))
      return (status);
	}
  
	return (status);
}

/* Photometry. */

int photHistory (ACSInfo *acs, Hdr *phdr) {
  
	extern int status;
	int logit;			/* true if we should log file name */
	int OmitStep (int);
	int TabHistory (RefTab *, Hdr *);
	int UpdateSwitch (char *, int, Hdr *, int *);
  
	if (OmitStep (acs->photcorr))
    return (status);
  
	if (UpdateSwitch ("PHOTCORR", acs->photcorr, phdr, &logit))
    return (status);
  
	if (logit) {
    if (TabHistory (&acs->phot, phdr))
      return (status);
	}
  
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
  
	if ((history = (char *) calloc (ACS_LINE+1, sizeof (char))) == NULL)
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
