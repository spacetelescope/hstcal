# include <stdio.h>
# include <string.h>		/* for strncmp, strcmp */

#include "hstcal.h"
# include "hstio.h"
# include "wf3.h"
# include "msg.h"
# include "wf3info.h"
# include "hstcalerr.h"		/* defines error codes */

static int checkAtoD (Hdr *, WF3Info *, int *, int *);
static int checkBias (Hdr *, WF3Info *, int *, int *);
static int checkBlev (Hdr *, WF3Info *, int *, int *);
static int checkCCD  (Hdr *, WF3Info *, int *);
static int checkDQI  (Hdr *, WF3Info *, int *, int *);
static int checkFlash(Hdr *, WF3Info *, int *, int *);
int GetImageRef (RefFileInfo *, Hdr *, char *, RefImage *, int *);

/* This routine gets the names of reference images and tables from the
   primary header and checks for dummy pedigree.

   Warren Hack, 1998 June 8:
    Original ACS version based on Phil Hodge's CALSTIS routine...

   Howard Bushouse, 2000 Aug 29:
    Original WFC3 version based on Warren Hack's CALACS routine.

   H.Bushouse, 2001 Nov 16:
    Updates to track CALACS changes - finished revisions for supporting
    post-flash processing.

   H.Bushouse, 2009 Jan 08:
    Enhanced all the checkNNNN routines to check for the correct FILETYPE
    for each reference file, as well as verifying correct selection criteria
    such as DETECTOR, FILTER, and CCDGAIN.

   M. De La Pena, 2022 February
    Added new SATUFILE: Full-well saturation image.

   M. De La Pena, 2023 May
    Only try to access the SATUFILE keyword if it is actually available in
    the header.  If the keyword is missing or does not contain a filename,
    the algorithm will indicate the original method of flagging saturated
    pixels by using a single value threshold should be used.

   M. De La Pena, 2023 October
    In addition to the current check of the BIASFILE not existing, added
    handling for the situation when the BIASFILE has a bad pedigree. In this
    instance BIASCORR cannot be performed.  As a consequence, full-well
    saturation must fall back to be applied as a scalar threshold.
*/

int GetFlags (WF3Info *wf3, Hdr *phdr) {

	extern int status;

	int missing = 0;	/* true if any calibration file is missing */
	int nsteps = 0;		/* number of calibration steps to perform */

	int GetccdSw (WF3Info *, Hdr *);

	/* Get the values for the Calibration Switches from the
	**	header for processing.  */
	if (GetccdSw (wf3, phdr) )
	    return(status);

	/* Check each reference file that we need. */
	if (checkDQI (phdr, wf3, &missing, &nsteps))
	    return (status);

	if (checkAtoD (phdr, wf3, &missing, &nsteps))
	    return (status);

	if (checkBlev (phdr, wf3, &missing, &nsteps))	/* no reference file */
	    return (status);

	if (checkCCD (phdr, wf3, &missing))
	    return (status);

	if (checkBias (phdr, wf3, &missing, &nsteps))
	    return (status);

	if (checkFlash (phdr, wf3, &missing, &nsteps))
	    return (status);

	if (missing) {
	    return (status = CAL_FILE_MISSING);
	} else if (nsteps < 1) {
	    trlwarn("No calibration switch was set to PERFORM, ");
	    trlwarn("            or all reference files had PEDIGREE = DUMMY.");
	    return (status = NOTHING_TO_DO);
	} else {
	    return (status);
	}
}

/* If this step is to be performed, check for the existence of the
   atod table.
*/

static int checkAtoD (Hdr *phdr, WF3Info *wf3, int *missing, int *nsteps) {

/* arguments:
Hdr *phdr         i: primary header
WF3Info *wf3      i: switches, file names, etc
int *missing     io: incremented if the file is missing or wrong type
int *nsteps      io: incremented if this step can be performed
*/

	extern int status;

	int calswitch;
	int GetSwitch (Hdr *, char *, int *);
	int GetTabRef (RefFileInfo *, Hdr *, char *, RefTab *, int *);
	void MissingFile (char *, char *, int *);
	void CheckTabType (RefTab *, char *, char *, int *);

	/* Are we supposed to do this step? */
	if (wf3->atodcorr == PERFORM) {

	    if (GetSwitch (phdr, "ATODCORR", &calswitch))
		return (status);
	    if (calswitch == COMPLETE) {
		wf3->atodcorr = OMIT;
		return (status);
	    }

	    /* Get the table name, check that the file exists, and get
	    ** pedigree and descrip.  */
	    if (GetTabRef (wf3->refnames, phdr, "ATODTAB", &wf3->atod,
			   &wf3->atodcorr))
		return (status);
	    if (wf3->atod.exists != EXISTS_YES) {
		MissingFile ("ATODTAB", wf3->atod.name, missing);

	    } else {

		/* Is the FILETYPE appropriate for an AtoD table? */
		CheckTabType (&wf3->atod, "ANALOG-TO-DIGITAL", "ATODTAB",
			      missing);
	    }

	    if (wf3->atodcorr == PERFORM)
		(*nsteps)++;
	}

	return (status);
}

/* If this step is to be performed, check for the existence of the
   bias file.  If it exists, get the pedigree and descrip keyword values.
*/

static int checkBias (Hdr *phdr, WF3Info *wf3, int *missing, int *nsteps) {

/* arguments:
Hdr *phdr         i: primary header
WF3Info *wf3      i: switches, file names, etc
int *missing     io: incremented if the file is missing
int *nsteps      io: incremented if this step can be performed
*/

	extern int status;

    int saveBiasCorr = GOOD_PEDIGREE;

	int calswitch;
	int GetSwitch (Hdr *, char *, int *);
	void MissingFile (char *, char *, int *);
	void CheckImgType (RefImage *, char *, char *, int *);
	int CheckGain (char *, float, char *, int *);

	/* Are we supposed to do this step? */
	if (wf3->biascorr == PERFORM) {

	    if (GetSwitch (phdr, "BIASCORR", &calswitch))
		return (status);
	    if (calswitch == COMPLETE) {
		wf3->biascorr = OMIT;
		return (status);
	    }

        if (GetImageRef (wf3->refnames, phdr, "BIASFILE", &wf3->bias, &wf3->biascorr))
            return (status);

        /*
            If the BIASFILE has a DUMMY pedigree, the GetImageRef command does not error.
            However, the wf3->biacorr will be set to DUMMY, and this will cause the bias
            correction to be skipped.  Due to the new implementation which uses a
            full-well saturation *image*, both BLEVCORR *and* BIASCORR must be performed.
        */
        if (wf3->biascorr != PERFORM) {
            wf3->scalar_satflag = True;
            trlwarn("There is an issue with the BIASFILE, so BIASCORR will not be performed.\n");
            trlmessage("A single threshold value will be used for full-well saturation flagging.");
            return (status);
        }

        if (wf3->bias.exists != EXISTS_YES) {
            MissingFile ("BIASFILE", wf3->bias.name, missing);
        } else {

            /* Is the FILETYPE appropriate for a BIAS file? */
            CheckImgType (&wf3->bias, "BIAS", "BIASFILE", missing);

            /* Does it have the correct GAIN value? */
            if (CheckGain(wf3->bias.name, wf3->ccdgain, "CCDGAIN", missing))
                return (status);
        }

	    if (wf3->biascorr == PERFORM)
		(*nsteps)++;

        /* Save the value for recovery */
        saveBiasCorr = wf3->biascorr;

        /*
          Also check for the new full-well saturation image which is
          applied after BLEVCORR and BIASCORR are done. Since the reference
          file is not associated with its own "calibration step keyword"
          (e.g., SATUCORR), just using the BIASCORR key as a standin here -
          make sure the BIASCORR retains its value as set in the above code.

          This is a kludge.
       */
        if (GetImageRef (wf3->refnames, phdr,
                         "SATUFILE", &wf3->satmap, &wf3->biascorr))
        {
            wf3->scalar_satflag = True;
	        trlerror("SATUFILE not found or cannot be opened.");
	        trlmessage("A single threshold value will be used for full-well saturation flagging.");
            return (status);
        }

        /* Recover the biascorr setting */
        wf3->biascorr = saveBiasCorr;

        /* Accommodate a missing SATUFILE keyword or associated value */
        if (wf3->satmap.exists != EXISTS_YES) {
            wf3->scalar_satflag = True;
            MissingFile ("SATUFILE", wf3->satmap.name, missing);
            *missing = 0;
	        trlmessage("A single threshold value will be used for full-well saturation flagging.");
        }
    /*
      At the least BIASCORR is not set to PERFORM, so issue a message and set the scalar_satflag so
      the 2D saturation image will not be used.
    */
	} else {
        wf3->scalar_satflag = True;
	    trlmessage("BIASCORR is *NOT* set to PERFORM.");
	    trlmessage("A single threshold value will be used for full-well saturation flagging.\n");
    }

	return (status);
}


static int checkBlev (Hdr *phdr, WF3Info *wf3, int *missing, int *nsteps) {

/* arguments:
Hdr *phdr         i: primary header
WF3Info *wf3      i: switches, file names, etc
int *nsteps      io: incremented if this step can be performed
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

/* If this step is to be performed, check for the existence of the
   post-flash file.  If it exists, get the pedigree and descrip keyword values.
*/

static int checkFlash (Hdr *phdr, WF3Info *wf3, int *missing, int *nsteps) {

/* arguments:
Hdr *phdr        i: primary header
WF3Info *wf3     i: switches, file names, etc
int *missing     io: incremented if the file is missing
int *nsteps      io: incremented if this step can be performed
*/

	extern int status;

	int calswitch;
	int GetSwitch (Hdr *, char *, int *);
	void MissingFile (char *, char *, int *);
	void CheckImgType (RefImage *, char *, char *, int *);

	/* Are we supposed to do this step? */
	if (wf3->flashcorr == PERFORM) {

	    if (GetSwitch (phdr, "FLSHCORR", &calswitch))
		return (status);
	    if (calswitch == COMPLETE) {
		wf3->flashcorr = OMIT;
		return (status);
	    }

	    if (GetImageRef (wf3->refnames, phdr, "FLSHFILE", &wf3->flash,
			     &wf3->flashcorr))
		return (status);
	    if (wf3->flash.exists != EXISTS_YES) {
		MissingFile ("FLSHFILE", wf3->flash.name, missing);

	    } else {

		/* Is the FILETYPE appropriate for a FLASH file? */
		CheckImgType (&wf3->flash, "POST FLASH", "FLSHFILE", missing);
	    }

	    if (wf3->flashcorr == PERFORM)
		(*nsteps)++;
	}

	return (status);
}

/* If the detector is the CCD, we need the CCD parameters table for
   BIASCORR, DARKCORR, and PHOTCORR.  This routine checks that the table
   exists.

   We also need the table for initializing the error array, but we
   don't have a flag for that step.  That's why we need this table
   regardless of which steps are to be performed.
*/

static int checkCCD (Hdr *phdr, WF3Info *wf3, int *missing) {

/* arguments:
Hdr *phdr         i: primary header
WF3Info *wf3      i: switches, file names, etc
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
	    trlerror("CCDTAB `%s' is a dummy table.",wf3->ccdpar.name);
	} else {

	    /* Is the FILETYPE appropriate for a CCD table? */
	    CheckTabType (&wf3->ccdpar, "CCD PARAMETERS", "CCDTAB", missing);

	    /* Does it have the correct DETECTOR value? */
	    if (CheckDetector(wf3->ccdpar.name, wf3->detector, "DETECTOR",
			      missing))
		return (status);
	}

	if (GetTabRef (wf3->refnames, phdr,
			"OSCNTAB", &wf3->oscn, &wf3->blevcorr))
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



/* Check whether we should assign initial values to the data quality
   array.  There is a reference table but not an image for this step.

   For the CCD, there are two steps to DQICORR, checking for saturation
   and using the BPIXTAB to initialize the data quality array.  If no
   bad pixel table was specified (name is blank or "N/A"), or if the
   table is dummy, we can still do this calibration step, but it will
   just consist of checking and flagging saturation.

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
	int GetTabRef (RefFileInfo *, Hdr *, char *, RefTab *, int *);
	void MissingFile (char *, char *, int *);
	void CheckTabType (RefTab *, char *, char *, int *);
	int  CheckDetector (char *, int, char *, int *);

	if (wf3->dqicorr == PERFORM) {

	    if (GetImageRef (wf3->refnames, phdr, "SNKCFILE", &wf3->sink,
			     &wf3->dqicorr))
		    return (status);

	    if (wf3->sink.exists != EXISTS_YES)
		    MissingFile ("SNKCFILE", wf3->sink.name, missing);


	    if (GetTabRef (wf3->refnames, phdr,
				"BPIXTAB", &wf3->bpix, &wf3->dqicorr))
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

