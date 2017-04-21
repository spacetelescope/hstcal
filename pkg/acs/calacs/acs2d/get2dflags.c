# include <stdio.h>
# include <string.h>		/* for strncmp, strcmp */

# include "hstio.h"
# include "acs.h"
# include "acsinfo.h"
# include "hstcalerr.h"		/* defines error codes */

static int checkCCD (Hdr *, ACSInfo *, int *);
static int checkDark (Hdr *, ACSInfo *, int *, int *);
static int checkDarkCTE (Hdr *, ACSInfo *, int *, int *);
static int checkFlash (Hdr *, ACSInfo *, int *, int *);
static int checkFlashCTE (Hdr *, ACSInfo *, int *, int *);
static int checkDQI (Hdr *, ACSInfo *, int *, int *);
static int checkFlat (Hdr *, ACSInfo *, int *, int *);
static int checkNonLin (Hdr *, ACSInfo *, int *, int *);
static int checkPhot (Hdr *, ACSInfo *, int *, int *);
static int checkShad (Hdr *, ACSInfo *, int *, int *);

/* This routine gets the names of reference images and tables from the
 primary header and checks for dummy pedigree.

 Warren Hack, 1998 June 10:
 Initial ACS version.
 Warren Hack, 2001 Oct 17:
 Replaced APERTAB and PHOTTAB checks with GRAPHTAB and COMPTAB
 in checkPhot.
 Pey Lian Lim, 2012 Dec 12:
 Moved FLASHCORR from ACSCCD. Removed POSTFLSH support per PR 44872.
 */

int Get2dFlags (ACSInfo *acs2d, Hdr *phdr) {

    extern int status;

    int missing = 0;	/* true if any calibration file is missing */
    int nsteps = 0;		/* number of calibration steps to perform */

    /* Check each reference file that we need. */

    if (checkDQI (phdr, acs2d, &missing, &nsteps))
        return (status);

    if (acs2d->detector != MAMA_DETECTOR) {
        if (checkCCD (phdr, acs2d, &missing))
            return (status);
    }

    if (checkNonLin (phdr, acs2d, &missing, &nsteps))
        return (status);

    if (checkDark (phdr, acs2d, &missing, &nsteps))
        return (status);

    if (checkDarkCTE (phdr, acs2d, &missing, &nsteps))
        return (status);

    if (checkFlash (phdr, acs2d, &missing, &nsteps))
        return (status);

    if (checkFlashCTE (phdr, acs2d, &missing, &nsteps))
        return (status);

    if (checkFlat (phdr, acs2d, &missing, &nsteps))
        return (status);

    if (checkShad (phdr, acs2d, &missing, &nsteps))
        return (status);

    if (checkPhot (phdr, acs2d, &missing, &nsteps))
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

static int checkCCD (Hdr *phdr, ACSInfo *acs2d, int *missing) {

  /* arguments:
   Hdr *phdr        i: primary header
   ACSInfo *acs2d   i: switches, file names, etc
   int *missing     io: incremented if the table is missing
   */

    extern int status;
    int calswitch;			/* returned by GetTabRef and ignored */
    int GetTabRef (RefFileInfo *, Hdr *, char *, RefTab *, int *);
    void MissingFile (char *, char *, int *);

    if (acs2d->detector == MAMA_DETECTOR)
        return (status);

    if (GetTabRef (acs2d->refnames, phdr,
                   "CCDTAB", &acs2d->ccdpar, &calswitch))
        return (status);

    if (acs2d->ccdpar.exists != EXISTS_YES) {

        MissingFile ("CCDTAB", acs2d->ccdpar.name, missing);

    } else if (acs2d->ccdpar.goodPedigree != GOOD_PEDIGREE) {

        (*missing)++;
        sprintf (MsgText, "CCDTAB `%s' is a dummy table.", acs2d->ccdpar.name);
        trlerror (MsgText);
    }

    return (status);
}

/* If this step is to be performed, check for the existence of the
 dark file.  If it exists, get the pedigree and descrip keyword values.
 */

static int checkDark (Hdr *phdr, ACSInfo *acs2d, int *missing, int *nsteps) {

    /* arguments:
       Hdr *phdr        i: primary header
       ACSInfo *acs2d   i: switches, file names, etc
       int *missing     io: incremented if the file is missing
       int *nsteps      io: incremented if this step can be performed
     */

    extern int status;

    int calswitch;
    int GetSwitch (Hdr *, char *, int *);
    int GetImageRef (RefFileInfo *, Hdr *, char *, RefImage *, int *);
    void MissingFile (char *, char *, int *);

    if (acs2d->darkcorr == PERFORM && acs2d->pctecorr != PERFORM) {

        if (GetSwitch (phdr, "DARKCORR", &calswitch))
            return (status);

        if (calswitch == COMPLETE) {
            acs2d->darkcorr = OMIT;
            return (status);
        }

    if (GetImageRef (acs2d->refnames, phdr,
                     "DARKFILE", &acs2d->dark, &acs2d->darkcorr))
        return (status);

    if (acs2d->dark.exists != EXISTS_YES)
        MissingFile ("DARKFILE", acs2d->dark.name, missing);

    if (acs2d->darkcorr == PERFORM)
       (*nsteps)++;
    }

    return (status);
}

/* If this step is to be performed, check for the existence of the
 dark file.  If it exists, get the pedigree and descrip keyword values.
 */

static int checkDarkCTE (Hdr *phdr, ACSInfo *acs2d, int *missing, int *nsteps) {

    /* arguments:
       Hdr *phdr        i: primary header
       ACSInfo *acs2d   i: switches, file names, etc
       int *missing     io: incremented if the file is missing
       int *nsteps      io: incremented if this step can be performed
     */

    extern int status;

    int calswitch;
    int GetSwitch (Hdr *, char *, int *);
    int GetImageRef (RefFileInfo *, Hdr *, char *, RefImage *, int *);
    void MissingFile (char *, char *, int *);

    if (acs2d->darkcorr == PERFORM && acs2d->pctecorr == PERFORM) {

        if (GetSwitch (phdr, "DARKCORR", &calswitch))
            return (status);

        if (calswitch == COMPLETE) {
            acs2d->darkcorr = OMIT;
            return (status);
        }

        if (GetImageRef (acs2d->refnames, phdr,
                         "DRKCFILE", &acs2d->darkcte, &acs2d->darkcorr))
            return (status);

        if (acs2d->darkcte.exists != EXISTS_YES)
            MissingFile ("DRKCFILE", acs2d->darkcte.name, missing);

        if (acs2d->darkcorr == PERFORM)
            (*nsteps)++;
    }

    return (status);
}

/* If this step is to be performed, check for the existence of the
 post-flash file.  If it exists, get the pedigree and descrip keyword values.
 */

static int checkFlash (Hdr *phdr, ACSInfo *acs2d, int *missing, int *nsteps) {

    /* arguments:
       Hdr *phdr        i: primary header
       ACSInfo *acs2d   i: switches, file names, etc
       int *missing     io: incremented if the file is missing
       int *nsteps      io: incremented if this step can be performed
     */

    extern int status;

    int calswitch;
    int GetSwitch (Hdr *, char *, int *);
    int GetImageRef (RefFileInfo *, Hdr *, char *, RefImage *, int *);
    void MissingFile (char *, char *, int *);

    SingleGroup y;
    double ltm[2], ltv[2];
    const int extver=1;
    int GetLT (Hdr *, double *, double *);

    if (acs2d->flashcorr == PERFORM && acs2d->pctecorr != PERFORM) {

        if (GetSwitch (phdr, "FLSHCORR", &calswitch))
            return (status);

        if (calswitch == COMPLETE) {
            acs2d->flashcorr = OMIT;
            return (status);
        }

        if (GetImageRef (acs2d->refnames, phdr,
                         "FLSHFILE", &acs2d->flash, &acs2d->flashcorr))
            return (status);

        if (acs2d->flash.exists != EXISTS_YES)
            MissingFile ("FLSHFILE", acs2d->flash.name, missing);

        /* Error message if LTV not 0 */
        initSingleGroup (&y);
	getSingleGroup (acs2d->flash.name, extver, &y);
        if (hstio_err())
            return (status = OPEN_FAILED);
        if (GetLT (&y.sci.hdr, ltm, ltv)) {
            freeSingleGroup(&y);
            return (status);
	}
        if ((ltv[0] != 0) || (ltv[1] != 0)) {
            sprintf(MsgText, "FLSHFILE `%s' has untrimmed overscans.",
                    acs2d->flash.name);
            trlerror (MsgText);
            freeSingleGroup(&y);
            return (status = SIZE_MISMATCH);
	}
        freeSingleGroup(&y);

        if (acs2d->flashcorr == PERFORM)
            (*nsteps)++;
    }

    return (status);
}

/* If this step is to be performed, check for the existence of the
 post-flash file.  If it exists, get the pedigree and descrip keyword values.
 */

static int checkFlashCTE (Hdr *phdr, ACSInfo *acs2d, int *missing, int *nsteps) {

    /* This is basically the same as checkFlash().
       It was written to support CTE-corrected post-flash reference file.
       But then ACS Team decided to use non-corrected post-flash instead.
     */

    /* arguments:
       Hdr *phdr        i: primary header
       ACSInfo *acs2d   i: switches, file names, etc
       int *missing     io: incremented if the file is missing
       int *nsteps      io: incremented if this step can be performed
     */

    extern int status;

    int calswitch;
    int GetSwitch (Hdr *, char *, int *);
    int GetImageRef (RefFileInfo *, Hdr *, char *, RefImage *, int *);
    void MissingFile (char *, char *, int *);

    SingleGroup y;
    double ltm[2], ltv[2];
    const int extver=1;
    int GetLT (Hdr *, double *, double *);

    if (acs2d->flashcorr == PERFORM && acs2d->pctecorr == PERFORM) {

        if (GetSwitch (phdr, "FLSHCORR", &calswitch))
            return (status);

        if (calswitch == COMPLETE) {
            acs2d->flashcorr = OMIT;
            return (status);
        }

        if (GetImageRef (acs2d->refnames, phdr,
                         "FLSHFILE", &acs2d->flashcte, &acs2d->flashcorr))
            return (status);

        if (acs2d->flashcte.exists != EXISTS_YES)
            MissingFile ("FLSHFILE", acs2d->flashcte.name, missing);

        /* Error message if LTV not 0 */
        initSingleGroup (&y);
	getSingleGroup (acs2d->flashcte.name, extver, &y);
        if (hstio_err())
            return (status = OPEN_FAILED);
        if (GetLT (&y.sci.hdr, ltm, ltv)) {
            freeSingleGroup(&y);
            return (status);
	}
        if ((ltv[0] != 0) || (ltv[1] != 0)) {
            sprintf(MsgText, "FLSHFILE `%s' has untrimmed overscans.",
                    acs2d->flashcte.name);
            trlerror (MsgText);
            freeSingleGroup(&y);
            return (status = SIZE_MISMATCH);
	}
        freeSingleGroup(&y);

        if (acs2d->flashcorr == PERFORM)
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
 just consist of checking and flagging saturation.  For the MAMAs,
 however, if this switch is set to PERFORM, the table must exist.

 Note that, unlike most other switches, dqicorr will not be reset to
 omit if the header value is "COMPLETE", since it would not cause any
 problem to perform this step more than once.  The user might do so
 deliberately in order to accumulate the flags from more than one table.
 */

static int checkDQI (Hdr *phdr, ACSInfo *acs2d, int *missing, int *nsteps) {

  /* arguments:
   Hdr *phdr        i: primary header
   ACSInfo *acs2d   i: switches, file names, etc
   int *missing     io: incremented if the table is missing
   int *nsteps      io: incremented if this step can be performed
   */

    extern int status;

    int GotFileName (char *);
    int GetTabRef (RefFileInfo *, Hdr *, char *, RefTab *, int *);
    void MissingFile (char *, char *, int *);

    if (acs2d->dqicorr == PERFORM) {

        if (GetTabRef (acs2d->refnames, phdr,
                       "BPIXTAB", &acs2d->bpix, &acs2d->dqicorr))
            return (status);

        if (acs2d->bpix.exists != EXISTS_YES) {

            if (acs2d->detector == MAMA_DETECTOR ||
                GotFileName (acs2d->bpix.name)) {

                MissingFile ("BPIXTAB", acs2d->bpix.name, missing);
            }
        }

        if (acs2d->dqicorr == PERFORM)
            (*nsteps)++;
    }

    return (status);
}

/* If this step is to be performed, check for the existence of the
 flat field files.  If they exist, get the pedigree and descrip
 keyword values.  If pedigree is DUMMY, the flag may be reset.
 */

static int checkFlat (Hdr *phdr, ACSInfo *acs2d, int *missing, int *nsteps) {

  /* arguments:
   Hdr *phdr        i: primary header
   ACSInfo *acs2d   i: switches, file names, etc
   int *missing     io: incremented if a file is missing
   int *nsteps      io: incremented if this step can be performed
   */

    extern int status;

    int calswitch;
    int GetSwitch (Hdr *, char *, int *);
    int GotFileName (char *);
    int GetImageRef (RefFileInfo *, Hdr *, char *, RefImage *, int *);
    int GetTabRef (RefFileInfo *, Hdr *, char *, RefTab *, int *);
    void MissingFile (char *, char *, int *);

    /* Are we supposed to do this step? */
    if (acs2d->flatcorr == PERFORM) {

        if (GetSwitch (phdr, "FLATCORR", &calswitch))
            return (status);
        if (calswitch == COMPLETE) {
            acs2d->flatcorr = OMIT;
            return (status);
        }

        /* Initial values; may be reset below. */
        acs2d->pfltcorr = PERFORM;
        acs2d->dfltcorr = PERFORM;
        acs2d->lfltcorr = PERFORM;
        acs2d->cfltcorr = PERFORM;

        if (GetImageRef (acs2d->refnames, phdr,
                         "PFLTFILE", &acs2d->pflt, &acs2d->pfltcorr))
            return (status);
        if (acs2d->pflt.exists != EXISTS_YES) {
            if (GotFileName (acs2d->pflt.name)) {	/* name specified? */
                MissingFile ("PFLTFILE", acs2d->pflt.name, missing);
            } else {
                acs2d->pfltcorr = OMIT;	/* name was blank or "N/A" */
            }
        }

        if (GetImageRef (acs2d->refnames, phdr,
                         "DFLTFILE", &acs2d->dflt, &acs2d->dfltcorr))
            return (status);
        if (acs2d->dflt.exists != EXISTS_YES) {
            if (GotFileName (acs2d->dflt.name)) {
                MissingFile ("DFLTFILE", acs2d->dflt.name, missing);
            } else {
                acs2d->dfltcorr = OMIT;
            }
        }

        if (GetImageRef (acs2d->refnames, phdr,
                         "LFLTFILE", &acs2d->lflt, &acs2d->lfltcorr))
            return (status);
        if (acs2d->lflt.exists != EXISTS_YES) {
            if (GotFileName (acs2d->lflt.name)) {
                MissingFile ("LFLTFILE", acs2d->lflt.name, missing);
            } else {
                acs2d->lfltcorr = OMIT;
            }
        }

        if (GetImageRef (acs2d->refnames, phdr,
                         "CFLTFILE", &acs2d->cflt, &acs2d->cfltcorr))
            return (status);
        if (acs2d->cflt.exists != EXISTS_YES) {
            if (GotFileName (acs2d->cflt.name)) {
                MissingFile ("CFLTFILE", acs2d->cflt.name, missing);
            } else {
                acs2d->cfltcorr = OMIT;
            }
        }

        if (GetTabRef (acs2d->refnames, phdr,
                       "SPOTTAB", &acs2d->spot, &acs2d->cfltcorr))
            return (status);

        if (acs2d->spot.exists != EXISTS_YES) {

            if (GotFileName (acs2d->spot.name)) {
                MissingFile ("SPOTTAB", acs2d->spot.name, missing);
            } else {
                acs2d->cfltcorr = OMIT;
	    }
        }

        /* If any of the three parts of flat fielding is set to
           PERFORM, then we can do this step.  If not, and if any
           part is DUMMY because of the reference file, reset the
           flat field flag to DUMMY; this will mean that all the
           files that were specified have pedigree=dummy.
        */
        if (acs2d->pfltcorr == PERFORM || acs2d->dfltcorr == PERFORM ||
            acs2d->lfltcorr == PERFORM) {
            (*nsteps)++;
        } else if (acs2d->pfltcorr == OMIT && acs2d->dfltcorr == OMIT &&
                   acs2d->lfltcorr == OMIT) {
            (*missing)++;
            trlerror ("PFLTFILE, DFLTFILE, and LFLTFILE are all blank.");
        } else if (acs2d->pfltcorr == DUMMY || acs2d->dfltcorr == DUMMY ||
                   acs2d->lfltcorr == DUMMY) {
            acs2d->flatcorr = DUMMY;
        }
    }

    return (status);
}

/* Check whether we should check for and possibly correct nonlinearity.
 There is a reference table but not an image for this step.  The table
 is opened, and the values are read into the ACSInfo structure.
 */

static int checkNonLin (Hdr *phdr, ACSInfo *acs2d, int *missing, int *nsteps) {

  /* arguments:
   Hdr *phdr        i: primary header
   ACSInfo *acs2d   i: switches, file names, etc
   int *missing     io: incremented if the table is missing
   int *nsteps      io: incremented if this step can be performed
   */

    extern int status;

    int calswitch;
    int GetSwitch (Hdr *, char *, int *);
    int GetTabRef (RefFileInfo *, Hdr *, char *, RefTab *, int *);
    void MissingFile (char *, char *, int *);

    if (acs2d->glincorr == PERFORM || acs2d->lflgcorr == PERFORM) {

        if (GetSwitch (phdr, "GLINCORR", &calswitch))
            return (status);
        if (calswitch == COMPLETE)
            acs2d->glincorr = OMIT;

        if (GetSwitch (phdr, "LFLGCORR", &calswitch))
            return (status);
        if (calswitch == COMPLETE)
            acs2d->lflgcorr = OMIT;

        if (acs2d->glincorr != PERFORM && acs2d->lflgcorr != PERFORM)
            return (status);

        /* Get the name of the MAMA linearity information table,
           and check whether it exists.  calswitch is ignored;
           glincorr and lflgcorr may be reset to DUMMY below.
        */
        if (GetTabRef (acs2d->refnames, phdr,
                       "MLINTAB", &acs2d->mlin, &calswitch))
            return (status);
        if (acs2d->mlin.exists != EXISTS_YES)
            MissingFile ("MLINTAB", acs2d->mlin.name, missing);

        if (acs2d->mlin.goodPedigree != GOOD_PEDIGREE) {
            if (acs2d->glincorr == PERFORM)
                acs2d->glincorr = DUMMY;
            if (acs2d->lflgcorr == PERFORM)
                acs2d->lflgcorr = DUMMY;
        }

        if (acs2d->glincorr == PERFORM || acs2d->lflgcorr == PERFORM)
            (*nsteps)++;
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

static int checkPhot (Hdr *phdr, ACSInfo *acs2d, int *missing, int *nsteps) {

  /* arguments:
   Hdr *phdr        i: primary header
   ACSInfo *acs2d   i: switches, file names, etc
   int *missing     io: incremented if the table is missing
   int *nsteps      io: incremented if this step can be performed
   */

    extern int status;

    int GetTabRef (RefFileInfo *, Hdr *, char *, RefTab *, int *);
    int GotFileName (char *);
    void MissingFile (char *, char *, int *);

    if (acs2d->photcorr == PERFORM) {

        /* Get the name of the graphtab table. */
        if (GetTabRef (acs2d->refnames, phdr,
                       "IMPHTTAB", &acs2d->phot, &acs2d->photcorr))
            return (status);
        if (acs2d->phot.exists != EXISTS_YES) {
            MissingFile ("IMPHTTAB", acs2d->phot.name, missing);
        }

        if (acs2d->photcorr == PERFORM)
            (*nsteps)++;
    }

    return (status);
}

/* If this step is to be performed, check for the existence of the
 shutter shading file.  If it exists, get the pedigree and descrip
 keyword values.
 */

static int checkShad (Hdr *phdr, ACSInfo *acs2d, int *missing, int *nsteps) {

  /* arguments:
   Hdr *phdr        i: primary header
   ACSInfo *acs2d   i: switches, file names, etc
   int *missing     io: incremented if the file is missing
   int *nsteps      io: incremented if this step can be performed
   */

    extern int status;

    int calswitch;
    int GetSwitch (Hdr *, char *, int *);
    int GetImageRef (RefFileInfo *, Hdr *, char *, RefImage *, int *);
    void MissingFile (char *, char *, int *);

    if (acs2d->shadcorr == PERFORM) {

        if (GetSwitch (phdr, "SHADCORR", &calswitch))
            return (status);
        if (calswitch == COMPLETE) {
            acs2d->shadcorr = OMIT;
            return (status);
        }

        if (GetImageRef (acs2d->refnames, phdr,
                         "SHADFILE", &acs2d->shad, &acs2d->shadcorr))
            return (status);
        if (acs2d->shad.exists != EXISTS_YES)
            MissingFile ("SHADFILE", acs2d->shad.name, missing);

        if (acs2d->shadcorr == PERFORM)
            (*nsteps)++;
    }

    return (status);
}
