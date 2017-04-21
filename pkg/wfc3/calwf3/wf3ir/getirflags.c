# include <stdio.h>
# include <string.h>		/* for strncmp, strcmp */

# include "hstio.h"
# include "wf3.h"
# include "wf3info.h"
# include "err.h"		/* defines error codes */

static int checkCCD  (Hdr *, WF3Info *, int *);
static int checkDQI  (Hdr *, WF3Info *, int *, int *);
static int checkZoff (Hdr *, WF3Info *, int *, int *);
static int checkDark (Hdr *, WF3Info *, int *, int *);
static int checkBlev (Hdr *, WF3Info *, int *, int *);
static int checkNlin (Hdr *, WF3Info *, int *, int *);
static int checkFlat (Hdr *, WF3Info *, int *, int *);
static int checkPhot (Hdr *, WF3Info *, int *, int *);
static int checkCRRej (Hdr *, WF3Info *, int *, int *);

/* This routine gets the names of reference images and tables from the
   primary header and checks for dummy pedigree.

   Warren Hack, 1998 June 10:
   	Initial ACS version.
   Howard Bushouse, 2000 Aug 29:
	Initial WFC3 UVIS version.
   H.Bushouse, 2001 Apr 19:
	Modified for WFC3 IR use.
   H.Bushouse, 2001 Nov 16:
	Modified checkPhot routine: replaced APERTAB and PHOTTAB with
	GRAPHTAB and COMPTAB (same as CALACS).
   H.Bushouse, 2002 Apr 19:
	Added checkZoff routine.
   H.Bushouse, 2003 Aug 14:
	Corrected logic for reading dark and nlin switches and ref file
	info when zsigcorr is performed.
   H.Bushouse, 2007 Feb 21:
	Added logic to checkDark to turn off zsigcorr if dark=dummy.
   H.Bushouse, 2009 Jan 09:
	Upgraded all checkNNNN routines to verify correct FILETYPE for each
	ref file, as well as matching selection criteria such as DETECTOR
	and FILTER.
   H.Bushouse, 2009 Feb 25:
	Added the checkCRRej routine to check for the existence and correctness
	of the CRREJTAB ref table, for use in CRCORR.
   H.Bushouse, 2010 Oct 21:
	Fixed a reference to dqicorr in checkCRRej routine that should've
	been crcorr.
   H.Bushouse, 2011 Mar 11:
	Removed zsigcorr checks in checkDark, because zsigcorr no longer uses
	the dark. (PR #67728, Trac #681)
   H.Bushouse, 2011 Sep 7:
	Modified checkPhot to work with new imphttab instead of graph and comp
	tabs.
*/

int GetIRFlags (WF3Info *wf3, Hdr *phdr) {

	extern int status;

	int missing = 0;	/* true if any calibration file is missing */
	int nsteps = 0;		/* number of calibration steps to perform */
		
	int GetirSw (WF3Info *, Hdr *);

	/* Get the values for the Calibration Switches from the
	** header for processing. */
	if (GetirSw (wf3, phdr))
	    return (status);

	/* Check each reference file that we need. */
	if (checkCCD (phdr, wf3, &missing))
	    return (status);

	if (checkDQI (phdr, wf3, &missing, &nsteps))
	    return (status);

	if (checkZoff (phdr, wf3, &missing, &nsteps))
	    return (status);

	if (checkDark (phdr, wf3, &missing, &nsteps))
	    return (status);

	if (checkBlev (phdr, wf3, &missing, &nsteps))
	    return (status);

	if (checkNlin (phdr, wf3, &missing, &nsteps))
	    return (status);

	if (checkFlat (phdr, wf3, &missing, &nsteps))
	    return (status);

	if (checkPhot (phdr, wf3, &missing, &nsteps))
	    return (status);

	if (checkCRRej (phdr, wf3, &missing, &nsteps))
	    return (status);

	if (missing) {
	    return (status = CAL_FILE_MISSING);
	} else if (nsteps < 1) {
	    trlwarn ("No calibration switch was set to PERFORM,");
	    trlwarn ("  or all reference files had PEDIGREE = DUMMY.");
	    return (status = NOTHING_TO_DO);
	} else {
	    return (status);
	}
}


static int checkZoff (Hdr *phdr, WF3Info *wf3, int *missing, int *nsteps) {

/* arguments:
Hdr *phdr	i: primary header
WF3Info *wf3	i: switches, file names, etc
int *nsteps	io: incremented if this step can be performed
*/

	extern int status;

	int calswitch;
	int GetSwitch (Hdr *, char *, int *);

	/* Are we supposed to do this step? */
	if (wf3->zoffcorr == PERFORM) {

	    if (GetSwitch (phdr, "ZOFFCORR", &calswitch))
		return (status);

	    if (calswitch == COMPLETE) {
		wf3->zoffcorr = OMIT;
		return (status);
	    }

	    if (wf3->zoffcorr == PERFORM)
		(*nsteps)++;
	}

	return (status);
}


/* We need the CCD parameters table for PHOTCORR.
   This routine checks that the table exists.

   We also need the table for initializing the error array, but we
   don't have a flag for that step.  That's why we need this table
   regardless of which steps are to be performed.
*/

static int checkCCD (Hdr *phdr, WF3Info *wf3, int *missing) {

/* arguments:
Hdr *phdr        i: primary header
WF3Info *wf3     i: switches, file names, etc
int *missing     io: incremented if the table is missing
*/

	extern int status;
	int calswitch;			/* returned by GetTabRef and ignored */
	int GetTabRef (RefFileInfo *, Hdr *, char *, RefTab *, int *);
	void MissingFile (char *, char *, int *);
	int GotFileName (char *);
	void CheckTabType (RefTab *, char *, char *, int *);
	int CheckDetector (char *, int, char *, int *);

	if (GetTabRef (wf3->refnames, phdr, "CCDTAB", &wf3->ccdpar, &calswitch))
	    return (status);

	if (wf3->ccdpar.exists != EXISTS_YES) {

	    MissingFile ("CCDTAB", wf3->ccdpar.name, missing);

	} else if (wf3->ccdpar.goodPedigree != GOOD_PEDIGREE) {

	    (*missing)++;
	    sprintf (MsgText, "CCDTAB `%s' is a dummy table.",wf3->ccdpar.name);
	    trlerror (MsgText);

	} else {

	    /* Is the FILETYPE appropriate for a CCD table? */
	    CheckTabType (&wf3->ccdpar, "CCD PARAMETERS", "CCDTAB", missing);

	    /* Does it have the correct DETECTOR value? */
	    if (CheckDetector(wf3->ccdpar.name, wf3->detector, "DETECTOR",
			      missing))
		return (status);
	}

	/* Load the overscan reference table (OSCNTAB) */
	if (GetTabRef (wf3->refnames, phdr, "OSCNTAB", &wf3->oscn,
		       &wf3->blevcorr))
	    return (status);

	if (wf3->oscn.exists != EXISTS_YES) {
	    if (GotFileName (wf3->oscn.name)) {
		MissingFile ("OSCNTAB", wf3->oscn.name, missing);
	    }

	} else {

	    /* Is the FILETYPE appropriate for an OSCN table? */
	    CheckTabType (&wf3->oscn, "OVERSCAN", "OSCNTAB", missing);

	    /* Does it have the correct DETECTOR value? */
	    if (CheckDetector(wf3->oscn.name, wf3->detector, "DETECTOR",
			      missing))
		return (status);
	}

	return (status);
}

/* If this step is to be performed, check for the existence of the
   dark file.  If it exists, get the pedigree and descrip keyword values.
*/

static int checkDark (Hdr *phdr, WF3Info *wf3, int *missing, int *nsteps) {

/* arguments:
Hdr *phdr        i: primary header
WF3Info *wf3     i: switches, file names, etc
int *missing     io: incremented if the file is missing
int *nsteps      io: incremented if this step can be performed
*/

	extern int status;

	int calswitch;
	int GetSwitch (Hdr *, char *, int *);
	int GetImageRef (RefFileInfo *, Hdr *, char *, RefImage *, int *);
	void MissingFile (char *, char *, int *);
	void CheckImgType (RefImage *, char *, char *, int *);
	int  CheckDetector (char *, int, char *, int *);

	if (wf3->darkcorr == PERFORM) {

	    if (GetSwitch (phdr, "DARKCORR", &calswitch))
		return (status);
	    if (calswitch == COMPLETE) {
		wf3->darkcorr = OMIT;
		return (status);
	    }

	    if (GetImageRef (wf3->refnames, phdr, "DARKFILE", &wf3->dark,
			     &wf3->darkcorr))
		return (status);
	    if (wf3->dark.exists != EXISTS_YES) {
		MissingFile ("DARKFILE", wf3->dark.name, missing);

            } else {

		/* Is the FILETYPE appropriate for a DARK file? */
		CheckImgType (&wf3->dark, "DARK", "DARKFILE", missing);

		/* Does it have the correct DETECTOR value? */
		if (CheckDetector(wf3->dark.name, wf3->detector, "DETECTOR",
				  missing))
		return (status);
	    }

	    if (wf3->darkcorr == PERFORM)
		(*nsteps)++;
	}

	return (status);
}


static int checkBlev (Hdr *phdr, WF3Info *wf3, int *missing, int *nsteps) {

/* arguments:
Hdr *phdr	i: primary header
WF3Info *wf3	i: switches, file names, etc
int *nsteps	io: incremented if this step can be performed
*/

	extern int status;

	int calswitch;
	int GetSwitch (Hdr *, char *, int *);

	/* Are we supposed to do this step? */
	if (wf3->blevcorr == PERFORM) {

	    if (GetSwitch (phdr, "BLEVCORR", &calswitch))
		return (status);

	    if (calswitch == COMPLETE) {
		wf3->blevcorr = OMIT;
		return (status);
	    }

	    if (wf3->blevcorr == PERFORM)
		(*nsteps)++;
	}

	return (status);
}


/* Check whether we should assign initial values to the data quality
   array.  There is a reference table but not an image for this step.

   For the CCD, there are two steps to DQICORR, checking for saturation
   and using the BPIXTAB to initialize the data quality array.  If no
   bad pixel table was specified (name is blank or "N/A"), or if the
   table is dummy, we can still do this calibration step, but it will
   just consist of checking and flagging saturation.  For the IR detector,
   however, if this switch is set to PERFORM, the table must exist.

   Note that, unlike most other switches, dqicorr will not be reset to
   omit if the header value is "COMPLETE", since it would not cause any
   problem to perform this step more than once.  The user might do so
   deliberately in order to accumulate the flags from more than one table.
*/

static int checkDQI (Hdr *phdr, WF3Info *wf3, int *missing, int *nsteps) {

/* arguments:
Hdr *phdr         i: primary header
WF3Info *wf3      i: switches, file names, etc
int *missing     io: incremented if the table is missing
int *nsteps      io: incremented if this step can be performed
*/

	extern int status;

	int GotFileName (char *);
	int GetTabRef    (RefFileInfo *, Hdr *, char *, RefTab *, int *);
	void MissingFile (char *, char *, int *);
	void CheckTabType (RefTab *, char *, char *, int *);
	int  CheckDetector (char *, int, char *, int *);

	if (wf3->dqicorr == PERFORM) {

	    if (GetTabRef (wf3->refnames, phdr, "BPIXTAB", &wf3->bpix,
			   &wf3->dqicorr))
		return (status);

	    if (wf3->bpix.exists != EXISTS_YES) {
		if (GotFileName (wf3->bpix.name)) {
		    MissingFile ("BPIXTAB", wf3->bpix.name, missing);
		}

	    } else {

		/* Is the FILETYPE appropriate for a BPIX table? */
		CheckTabType (&wf3->bpix, "BAD PIXELS", "BPIXTAB", missing);

		/* Does it have the correct DETECTOR value? */
		if (CheckDetector(wf3->bpix.name, wf3->detector, "DETECTOR",
				  missing))
		return (status);
	    }

	    if (wf3->dqicorr == PERFORM)
		(*nsteps)++;
	}

	return (status);
}

/* If this step is to be performed, check for the existence of the
   flat field files.  If they exist, get the pedigree and descrip
   keyword values.  If pedigree is DUMMY, the flag may be reset.
*/

static int checkFlat (Hdr *phdr, WF3Info *wf3, int *missing, int *nsteps) {

/* arguments:
Hdr *phdr        i: primary header
WF3Info *wf3     i: switches, file names, etc
int *missing     io: incremented if a file is missing
int *nsteps      io: incremented if this step can be performed
*/

	extern int status;

	int calswitch;
	int GetSwitch (Hdr *, char *, int *);
	int GotFileName (char *);
	int GetImageRef (RefFileInfo *, Hdr *, char *, RefImage *, int *);
	void MissingFile (char *, char *, int *);
	void CheckImgType (RefImage *, char *, char *, int *);
	int  CheckFilter  (char *, char *, char *, int *);
	int  CheckDetector (char *, int, char *, int *);

	/* Are we supposed to do this step? */
	if (wf3->flatcorr == PERFORM) {
	    if (GetSwitch (phdr, "FLATCORR", &calswitch))
		return (status);
	    if (calswitch == COMPLETE) {
		wf3->flatcorr = OMIT;
		return (status);
	    }

	    /* Initial values; may be reset below. */
	    wf3->pfltcorr = PERFORM;
	    wf3->dfltcorr = PERFORM;
	    wf3->lfltcorr = PERFORM;

	    if (GetImageRef (wf3->refnames, phdr, "PFLTFILE", &wf3->pflt,
			     &wf3->pfltcorr))
		return (status);
	    if (wf3->pflt.exists != EXISTS_YES) {
		if (GotFileName (wf3->pflt.name)) {	/* name specified? */
		    MissingFile ("PFLTFILE", wf3->pflt.name, missing);
		} else {
		    wf3->pfltcorr = OMIT;	/* name was blank or "N/A" */
		}
	    } else {
		/* Is the FILETYPE appropriate for a PFLT file? */
		CheckImgType (&wf3->pflt, "PIXEL-TO-PIXEL FLAT", "PFLTFILE",
			      missing);
		/* Does it have the correct FILTER value? */
		if (CheckFilter(wf3->pflt.name, wf3->filter, "FILTER", missing))
		    return (status);
		/* Does it have the correct DETECTOR value? */
		if (CheckDetector(wf3->pflt.name, wf3->detector, "DETECTOR",
				  missing))
		    return (status);
	    }

	    if (GetImageRef (wf3->refnames, phdr, "DFLTFILE", &wf3->dflt,
			     &wf3->dfltcorr))
		return (status);
	    if (wf3->dflt.exists != EXISTS_YES) {
		if (GotFileName (wf3->dflt.name)) {
		    MissingFile ("DFLTFILE", wf3->dflt.name, missing);
		} else {
		    wf3->dfltcorr = OMIT;
		}
	    } else {
		/* Is the FILETYPE appropriate for a PFLT file? */
		CheckImgType (&wf3->dflt, "DELTA FLAT", "DFLTFILE", missing);
		/* Does it have the correct FILTER value? */
		if (CheckFilter(wf3->dflt.name, wf3->filter, "FILTER",
				missing))
		    return (status);
		/* Does it have the correct DETECTOR value? */
		if (CheckDetector(wf3->dflt.name, wf3->detector, "DETECTOR",
				  missing))
		    return (status);
	    }

	    if (GetImageRef (wf3->refnames, phdr, "LFLTFILE", &wf3->lflt,
			     &wf3->lfltcorr))
		return (status);
	    if (wf3->lflt.exists != EXISTS_YES) {
		if (GotFileName (wf3->lflt.name)) {
		    MissingFile ("LFLTFILE", wf3->lflt.name, missing);
		} else {
		    wf3->lfltcorr = OMIT;
		}
	    } else {
		/* Is the FILETYPE appropriate for a LFLT file? */
		CheckImgType (&wf3->lflt, "LARGE SCALE FLAT", "LFLTFILE",
			      missing);
		/* Does it have the correct FILTER value? */
		if (CheckFilter(wf3->lflt.name, wf3->filter, "FILTER", missing))
		    return (status);
		/* Does it have the correct DETECTOR value? */
		if (CheckDetector(wf3->lflt.name, wf3->detector, "DETECTOR",
				  missing))
		    return (status);
	    }

	    /* If any of the three parts of flat fielding is set to
		PERFORM, then we can do this step.  If not, and if any
		part is DUMMY because of the reference file, reset the
		flat field flag to DUMMY; this will mean that all the
		files that were specified have pedigree=dummy.
	    */
	    if (wf3->pfltcorr == PERFORM || wf3->dfltcorr == PERFORM ||
		wf3->lfltcorr == PERFORM) {
                (*nsteps)++;
	    } else if (wf3->pfltcorr == OMIT && wf3->dfltcorr == OMIT &&
		wf3->lfltcorr == OMIT) {
                (*missing)++;
                trlerror ("PFLTFILE, DFLTFILE, and LFLTFILE are all blank.");
	    } else if (wf3->pfltcorr == DUMMY || wf3->dfltcorr == DUMMY ||
		wf3->lfltcorr == DUMMY) {
                wf3->flatcorr = DUMMY;
	    }
	}

	return (status);
}

/* Check whether we should compute the photometry keyword values.
   There is a reference table but not an image for this step.
   This will only be done for imaging mode, not spectroscopic.
   Also, unlike most other switches, photcorr will not be reset
   to omit if the header value is "COMPLETE", since it would not
   cause any problem to perform this step more than once.
*/

static int checkPhot (Hdr *phdr, WF3Info *wf3, int *missing, int *nsteps) {

/* arguments:
Hdr *phdr        i: primary header
WF3Info *wf3     i: switches, file names, etc
int *missing     io: incremented if the table is missing
int *nsteps      io: incremented if this step can be performed
*/

	extern int status;

	int GetTabRef (RefFileInfo *, Hdr *, char *, RefTab *, int *);
	int GotFileName (char *);
	void MissingFile (char *, char *, int *);

	if (wf3->photcorr == PERFORM) {

	    /* Get the name of the imphttab table. */
	    if (GetTabRef (wf3->refnames, phdr, "IMPHTTAB", &wf3->phot,
			   &wf3->photcorr))
		return (status);
	    if (wf3->phot.exists != EXISTS_YES) {
		MissingFile ("IMPHTTAB", wf3->phot.name, missing);
	    }

	    if (wf3->photcorr == PERFORM)
		(*nsteps)++;
	}

	return (status);
}

/* If this step is to be performed, check for the existence of the
   linearity coefficients file.  If it exists, get the pedigree and descrip
   keyword values.
*/

static int checkNlin (Hdr *phdr, WF3Info *wf3, int *missing, int *nsteps) {

/* arguments:
Hdr *phdr        i: primary header
WF3Info *wf3     i: switches, file names, etc
int *missing     io: incremented if the file is missing
int *nsteps      io: incremented if this step can be performed
*/

	extern int status;

	int calswitch;
	int GetSwitch (Hdr *, char *, int *);
	int GetImageRef (RefFileInfo *, Hdr *, char *, RefImage *, int *);
	void MissingFile (char *, char *, int *);
	void CheckImgType (RefImage *, char *, char *, int *);
	int  CheckDetector (char *, int, char *, int *);

	if (wf3->nlincorr == PERFORM || wf3->zsigcorr == PERFORM) {

	    if (GetSwitch (phdr, "NLINCORR", &calswitch))
		return (status);
	    if (calswitch == COMPLETE) {
		wf3->nlincorr = OMIT;
		return (status);
	    }

	    if (GetImageRef (wf3->refnames, phdr, "NLINFILE", &wf3->nlin,
			     &wf3->nlincorr))
		return (status);
	    if (wf3->nlin.exists != EXISTS_YES) {
		MissingFile ("NLINFILE", wf3->nlin.name, missing);

	    } else {

		/* Is the FILETYPE appropriate for a NLIN file? */
		CheckImgType (&wf3->nlin, "LINEARITY COEFFICIENTS", "NLINFILE",
			      missing);
		/* Does it have the correct DETECTOR value? */
		if (CheckDetector(wf3->nlin.name, wf3->detector, "DETECTOR",
				  missing))
		    return (status);
	    }

	    if (wf3->nlincorr == PERFORM || wf3->zsigcorr == PERFORM)
		(*nsteps)++;
	}

	return (status);
}

/* We need the CRREJ parameters table for CRCORR.
   This routine checks that the table exists and has the proper
   DETECTOR value.
*/

static int checkCRRej (Hdr *phdr, WF3Info *wf3, int *missing, int *nsteps) {

/* arguments:
Hdr *phdr        i: primary header
WF3Info *wf3     i: switches, file names, etc
int *missing     io: incremented if the table is missing
int *nsteps      io: incremented if this step can be performed
*/

	extern int status;
	int calswitch;			/* returned by GetTabRef and ignored */
	int GetTabRef (RefFileInfo *, Hdr *, char *, RefTab *, int *);
	void MissingFile (char *, char *, int *);
	int GotFileName (char *);
	void CheckTabType (RefTab *, char *, char *, int *);
	int CheckDetector (char *, int, char *, int *);

	if (wf3->crcorr == PERFORM) {

	    if (GetTabRef (wf3->refnames, phdr, "CRREJTAB", &wf3->crrej,
			   &calswitch))
		return (status);

	    if (wf3->crrej.exists != EXISTS_YES) {

		MissingFile ("CRREJTAB", wf3->crrej.name, missing);

	    } else if (wf3->crrej.goodPedigree != GOOD_PEDIGREE) {

		(*missing)++;
		sprintf (MsgText,"CRREJTAB '%s' is a dummy table.",
			 wf3->crrej.name);
		trlerror (MsgText);

	    } else {

		/* Is the FILETYPE appropriate for a CRREJ table? */
		CheckTabType (&wf3->crrej, "COSMIC RAY REJECTION", "CRREJTAB",
			      missing);

		/* Does it have the correct DETECTOR value? */
		if (CheckDetector(wf3->crrej.name, wf3->detector, "DETECTOR",
		    missing))
		    return (status);
	    }

	    if (wf3->crcorr == PERFORM)
		(*nsteps)++;
	}

	return (status);
}

