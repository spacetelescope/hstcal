# include <stdio.h>
# include <string.h>		/* for strncmp, strcmp */

# include "hstio.h"

# include "wf3.h"
# include "wf3info.h"
# include "err.h"		/* defines error codes */

static int checkCCD  (Hdr *, WF3Info *, int *);
static int checkDark (Hdr *, WF3Info *, int *, int *);
static int checkDQI  (Hdr *, WF3Info *, int *, int *);
static int checkFlat (Hdr *, WF3Info *, int *, int *);
static int checkPhot (Hdr *, WF3Info *, int *, int *);
static int checkShad (Hdr *, WF3Info *, int *, int *);

/* This routine gets the names of reference images and tables from the
   primary header and checks for dummy pedigree.

   Warren Hack, 1998 June 10:
   	Initial ACS version.
   Howard Bushouse, 2000 Aug 29:
	Initial WFC3 version.
   H. Bushouse, 2001 Nov 16:
	Updated to track CALACS changes - Replaced APERTAB and PHOTTAB
	checks with GRAPHTAB and COMPTAB in checkPhot.
   H. Bushouse, 2009 Jan 09:
	Enhanced all the checkNNNN routines to check for the correct FILETYPE
	for each reference file, as well as matching selection criteria such
	as DETECTOR and FILTER.
   H. Bushouse, 2011 Sep 7:
	Modified checkPhot to work with new imphttab instead of graph and
	comp tabs.
*/

int Get2dFlags (WF3Info *wf32d, Hdr *phdr) {

	extern int status;

	int missing = 0;	/* true if any calibration file is missing */
	int nsteps = 0;		/* number of calibration steps to perform */
		
	/* Check each reference file that we need. */

	if (checkDQI (phdr, wf32d, &missing, &nsteps))
	    return (status);

	if (checkCCD (phdr, wf32d, &missing))
	    return (status);

	if (checkDark (phdr, wf32d, &missing, &nsteps))
	    return (status);

	if (checkFlat (phdr, wf32d, &missing, &nsteps))
	    return (status);

	if (checkShad (phdr, wf32d, &missing, &nsteps))
	    return (status);

	if (checkPhot (phdr, wf32d, &missing, &nsteps))
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


/* If the detector is the CCD, we need the CCD parameters table for
   BIASCORR, DARKCORR, and PHOTCORR.  This routine checks that the table
   exists.

   We also need the table for initializing the error array, but we
   don't have a flag for that step.  That's why we need this table
   regardless of which steps are to be performed.
*/

static int checkCCD (Hdr *phdr, WF3Info *wf32d, int *missing) {

/* arguments:
Hdr *phdr        i: primary header
WF3Info *wf32d   i: switches, file names, etc
int *missing     io: incremented if the table is missing
*/

	extern int status;
	int calswitch;			/* returned by GetTabRef and ignored */
	int GetTabRef (RefFileInfo *, Hdr *, char *, RefTab *, int *);
	void MissingFile (char *, char *, int *);
	void CheckTabType (RefTab *, char *, char *, int *);
	int  CheckDetector (char *, int, char *, int *);

	if (GetTabRef (wf32d->refnames, phdr, "CCDTAB", &wf32d->ccdpar,
		       &calswitch))
	    return (status);

	if (wf32d->ccdpar.exists != EXISTS_YES) {

	    MissingFile ("CCDTAB", wf32d->ccdpar.name, missing);

	} else if (wf32d->ccdpar.goodPedigree != GOOD_PEDIGREE) {

	    (*missing)++;
	    sprintf (MsgText, "CCDTAB `%s' is a dummy table.",
		     wf32d->ccdpar.name);
	    trlerror (MsgText);

	} else {

	    /* Is the FILETYPE appropriate for a CCD table? */
	    CheckTabType (&wf32d->ccdpar, "CCD PARAMETERS", "CCDTAB", missing);

	    /* Does it have the correct DETECTOR value? */
	    if (CheckDetector(wf32d->ccdpar.name, wf32d->detector, "DETECTOR",
		missing))
		return (status);
	}

	return (status);
}

/* If this step is to be performed, check for the existence of the
   dark file.  If it exists, get the pedigree and descrip keyword values.
*/

static int checkDark (Hdr *phdr, WF3Info *wf32d, int *missing, int *nsteps) {

/* arguments:
Hdr *phdr        i: primary header
WF3Info *wf32d   i: switches, file names, etc
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
    char *darktouse;
    char *darktype;
    
	if (wf32d->darkcorr == PERFORM) {
	    if (GetSwitch (phdr, "DARKCORR", &calswitch))
		    return (status);
	    if (calswitch == COMPLETE) {
		    wf32d->darkcorr = OMIT;
		    return (status);
	    } else {
        
            if (GetSwitch (phdr, "PCTECORR", &calswitch))   
                return(status);
            if (calswitch == COMPLETE){
                darktouse="DRKCFILE";
                darktype="CTEDARK";
            } else {
                darktouse="DARKFILE";
                darktype="DARK";
            }
        }

            
	    if (GetImageRef (wf32d->refnames, phdr, darktouse, &wf32d->dark,
			     &wf32d->darkcorr))
		    return (status);

	    if (wf32d->dark.exists != EXISTS_YES) {
		    MissingFile (darktouse, wf32d->dark.name, missing);

	    } else {

		    /* Is the FILETYPE appropriate for a DARK file? */
		    CheckImgType (&wf32d->dark, darktype, darktouse, missing);

		    /* Does it have the correct DETECTOR value? */
		    if (CheckDetector(wf32d->dark.name, wf32d->detector, "DETECTOR",
				  missing))
		        return (status);
	    }

	    if (wf32d->darkcorr == PERFORM)
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
   just consist of checking and flagging saturation.  For IR images,
   however, if this switch is set to PERFORM, the table must exist.

   Note that, unlike most other switches, dqicorr will not be reset to
   omit if the header value is "COMPLETE", since it would not cause any
   problem to perform this step more than once.  The user might do so
   deliberately in order to accumulate the flags from more than one table.
*/

static int checkDQI (Hdr *phdr, WF3Info *wf32d, int *missing, int *nsteps) {

/* arguments:
Hdr *phdr         i: primary header
WF3Info *wf32d    i: switches, file names, etc
int *missing     io: incremented if the table is missing
int *nsteps      io: incremented if this step can be performed
*/

	extern int status;

	int GotFileName (char *);
	int GetTabRef (RefFileInfo *, Hdr *, char *, RefTab *, int *);
	void MissingFile (char *, char *, int *);
	void CheckTabType (RefTab *, char *, char *, int *);
	int  CheckDetector (char *, int, char *, int *);

	if (wf32d->dqicorr == PERFORM) {

	    if (GetTabRef (wf32d->refnames, phdr, "BPIXTAB", &wf32d->bpix,
			   &wf32d->dqicorr))
		return (status);

	    if (wf32d->bpix.exists != EXISTS_YES) {

		if (wf32d->detector == IR_DETECTOR ||
		    GotFileName (wf32d->bpix.name)) {

		    MissingFile ("BPIXTAB", wf32d->bpix.name, missing);
		}

	    } else {

		/* Is the FILETYPE appropriate for a BPIX table? */
		CheckTabType (&wf32d->bpix, "BAD PIXELS", "BPIXTAB", missing);

		/* Does it have the correct DETECTOR value? */
		if (CheckDetector(wf32d->bpix.name, wf32d->detector, "DETECTOR",
				  missing))
		    return (status);
	    }

	    if (wf32d->dqicorr == PERFORM)
		(*nsteps)++;
	}

	return (status);
}

/* If this step is to be performed, check for the existence of the
   flat field files.  If they exist, get the pedigree and descrip
   keyword values.  If pedigree is DUMMY, the flag may be reset.
*/

static int checkFlat (Hdr *phdr, WF3Info *wf32d, int *missing, int *nsteps) {

/* arguments:
Hdr *phdr        i: primary header
WF3Info *wf32d   i: switches, file names, etc
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
	if (wf32d->flatcorr == PERFORM) {

	    if (GetSwitch (phdr, "FLATCORR", &calswitch))
		return (status);
	    if (calswitch == COMPLETE) {
		wf32d->flatcorr = OMIT;
		return (status);
	    }

	    /* Initial values; may be reset below. */
	    wf32d->pfltcorr = PERFORM;
	    wf32d->dfltcorr = PERFORM;
	    wf32d->lfltcorr = PERFORM;

	    if (GetImageRef (wf32d->refnames, phdr, "PFLTFILE", &wf32d->pflt,
			     &wf32d->pfltcorr))
		return (status);
	    if (wf32d->pflt.exists != EXISTS_YES) {
		if (GotFileName (wf32d->pflt.name)) {	/* name specified? */
		    MissingFile ("PFLTFILE", wf32d->pflt.name, missing);
		} else {
		    wf32d->pfltcorr = OMIT;	/* name was blank or "N/A" */
		}
	    } else {
		/* Is the FILETYPE appropriate for a PFLT file? */
		CheckImgType (&wf32d->pflt, "PIXEL-TO-PIXEL FLAT", "PFLTFILE", 
			      missing);
		/* Does it have the correct FILTER value? */
		if (CheckFilter(wf32d->pflt.name, wf32d->filter, "FILTER",
				missing))
		    return (status);
		/* Does it have the correct DETECTOR value? */
		if (CheckDetector(wf32d->pflt.name, wf32d->detector, "DETECTOR",
				  missing))
		    return (status);
	    }

	    if (GetImageRef (wf32d->refnames, phdr, "DFLTFILE", &wf32d->dflt,
			     &wf32d->dfltcorr))
		return (status);
	    if (wf32d->dflt.exists != EXISTS_YES) {
		if (GotFileName (wf32d->dflt.name)) {
		    MissingFile ("DFLTFILE", wf32d->dflt.name, missing);
		} else {
		    wf32d->dfltcorr = OMIT;
		}
	    } else {
		/* Is the FILETYPE appropriate for a DFLT file? */
		CheckImgType (&wf32d->dflt, "DELTA FLAT", "DFLTFILE", 
			      missing);
		/* Does it have the correct FILTER value? */
		if (CheckFilter(wf32d->dflt.name, wf32d->filter, "FILTER",
				missing))
		    return (status);
		/* Does it have the correct DETECTOR value? */
		if (CheckDetector(wf32d->dflt.name, wf32d->detector, "DETECTOR",
				  missing))
		    return (status);
	    }

	    if (GetImageRef (wf32d->refnames, phdr, "LFLTFILE", &wf32d->lflt,
			     &wf32d->lfltcorr))
		return (status);
	    if (wf32d->lflt.exists != EXISTS_YES) {
		if (GotFileName (wf32d->lflt.name)) {
		    MissingFile ("LFLTFILE", wf32d->lflt.name, missing);
		} else {
		    wf32d->lfltcorr = OMIT;
		}
	    } else {
		/* Is the FILETYPE appropriate for a LFLT file? */
		CheckImgType (&wf32d->lflt, "LARGE SCALE FLAT", "LFLTFILE", 
			      missing);
		/* Does it have the correct FILTER value? */
		if (CheckFilter(wf32d->lflt.name, wf32d->filter, "FILTER",
				missing))
		    return (status);
		/* Does it have the correct DETECTOR value? */
		if (CheckDetector(wf32d->lflt.name, wf32d->detector, "DETECTOR",
				  missing))
		    return (status);
	    }

	    /* If any of the three parts of flat fielding is set to
		PERFORM, then we can do this step.  If not, and if any
		part is DUMMY because of the reference file, reset the
		flat field flag to DUMMY; this will mean that all the
		files that were specified have pedigree=dummy.
	    */
	    if (wf32d->pfltcorr == PERFORM || wf32d->dfltcorr == PERFORM ||
		wf32d->lfltcorr == PERFORM) {
                (*nsteps)++;
	    } else if (wf32d->pfltcorr == OMIT && wf32d->dfltcorr == OMIT &&
		wf32d->lfltcorr == OMIT) {
                (*missing)++;
                trlerror ("PFLTFILE, DFLTFILE, and LFLTFILE are all blank.");
	    } else if (wf32d->pfltcorr == DUMMY || wf32d->dfltcorr == DUMMY ||
		wf32d->lfltcorr == DUMMY) {
                wf32d->flatcorr = DUMMY;
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

static int checkPhot (Hdr *phdr, WF3Info *wf32d, int *missing, int *nsteps) {

/* arguments:
Hdr *phdr        i: primary header
WF3Info *wf32d   i: switches, file names, etc
int *missing     io: incremented if the table is missing
int *nsteps      io: incremented if this step can be performed
*/

	extern int status;

	int GetTabRef (RefFileInfo *, Hdr *, char *, RefTab *, int *);
	int GotFileName (char *);
	void MissingFile (char *, char *, int *);

	if (wf32d->photcorr == PERFORM) {

	    /* Get the name of the imphttab table. */
	    if (GetTabRef (wf32d->refnames, phdr, "IMPHTTAB", &wf32d->phot,
			   &wf32d->photcorr))
		return (status);
	    if (wf32d->phot.exists != EXISTS_YES) {
		MissingFile ("IMPHTTAB", wf32d->phot.name, missing);
	    }

	    if (wf32d->photcorr == PERFORM)
		(*nsteps)++;
	}

	return (status);
}

/* If this step is to be performed, check for the existence of the
   shutter shading file.  If it exists, get the pedigree and descrip
   keyword values.
*/

static int checkShad (Hdr *phdr, WF3Info *wf32d, int *missing, int *nsteps) {

/* arguments:
Hdr *phdr        i: primary header
WF3Info *wf32d   i: switches, file names, etc
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

	if (wf32d->shadcorr == PERFORM) {

	    if (GetSwitch (phdr, "SHADCORR", &calswitch))
		return (status);
	    if (calswitch == COMPLETE) {
		wf32d->shadcorr = OMIT;
		return (status);
	    }

	    if (GetImageRef (wf32d->refnames, phdr, "SHADFILE", &wf32d->shad,
			     &wf32d->shadcorr))
		return (status);
	    if (wf32d->shad.exists != EXISTS_YES) {
		MissingFile ("SHADFILE", wf32d->shad.name, missing);

	    } else {

		/* Is the FILETYPE appropriate for a SHAD file? */
		CheckImgType (&wf32d->shad, "SHUTTER SHADING", "SHADFILE",
			      missing);

		/* Does it have the correct DETECTOR value? */
		if (CheckDetector(wf32d->shad.name, wf32d->detector, "DETECTOR",
				  missing))
		    return (status);
	    }

	    if (wf32d->shadcorr == PERFORM)
		(*nsteps)++;
	}

	return (status);
}
