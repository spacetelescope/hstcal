# include <stdio.h>
# include <string.h>		/* for strncmp, strcmp, strlen */

# include "c_iraf.h"
# include "hstio.h"
# include "stis.h"
# include "calstis1.h"
# include "hstcalerr.h"		/* defines error codes */
# include "stisdef.h"

static int checkAtoD (Hdr *, StisInfo1 *, int *, int *);
static int checkBias (Hdr *, StisInfo1 *, int *, int *);
static int checkBlev (Hdr *, StisInfo1 *, int *);
static int checkCCD (Hdr *, StisInfo1 *, int *);
static int checkDark (Hdr *, StisInfo1 *, int *, int *);
static int checkDQI (Hdr *, StisInfo1 *, int *, int *);
static int checkFlat (Hdr *, StisInfo1 *, int *, int *);
static int checkNonLin (Hdr *, StisInfo1 *, int *, int *);
static int checkPhot (Hdr *, StisInfo1 *, int *, int *);
static int checkShad (Hdr *, StisInfo1 *, int *, int *);
static int GetImageRef (RefFileInfo *, Hdr *, char *, RefImage *, int *);
static int GetTabRef (RefFileInfo *, Hdr *, char *, RefTab *, int *);
static void MissingFile (char *, char *, int *);

/* This routine gets the names of reference images and tables from the
   primary header and checks for dummy pedigree.

   Phil Hodge, 1997 Dec 10:
	Get apertab in checkPhot;
	add sts->refnames to calling sequence of GetRefName.

   Phil Hodge, 1998 Jan 14:
	APERTAB for PHOTCORR is optional.  Add GetImageRef, GetTabRef, and
	MissingFile, based on code in cs7/getflags7.c.

   Phil Hodge, 1998 May 19:
	Don't get OBSTYPE from the header, because it's now in sts.

   Ivo Busko, 2002 Mar 29:
        NUV dark time-dependency.

   Paul Barrett, 2003 Sep 18:
        Read EPCTAB keyword from the header.

   Paul Barrett, 2003 Sep 25:
        Read TDSTAB keyword from the header.

   Paul Barrett, 2004 Jun 25:
        Removed reading of EPCTAB keyword from header, moved to GetGrpInfo1

   Phil Hodge, 2004 Aug 2:
	In checkDark, there was an "} else" that apparently should have
	been "} else {" with a closing brace later, based on the indentation.
	This actually did the right thing, but fix it anyway.

  Phil Hodge, 2011 May 9:
	Modify checkPhot to replace keyword name PHOTTAB with IMPHTTAB;
	don't get APERTAB or TDSTAB; don't set sts->filtcorr or sts->tdscorr.

  Phil Hodge, 2011 Nov 17:
	Modify checkPhot to include the TDSTAB.

   Phil Hodge, 2012 Oct 15:
	In checkDark, don't allow the TDCTAB to be optional.
*/

int GetFlags1 (StisInfo1 *sts, Hdr *phdr) {

	int status;

	int missing = 0;	/* true if any calibration file is missing */
	int nsteps = 0;		/* number of calibration steps to perform */

	/* Check each reference file that we need. */

	if ((status = checkDQI (phdr, sts, &missing, &nsteps)))
	    return (status);

	if ((status = checkAtoD (phdr, sts, &missing, &nsteps)))
	    return (status);

	if ((status = checkBlev (phdr, sts, &nsteps))) /* no reference file */
	    return (status);

	if (sts->detector == CCD_DETECTOR) {
	    if ((status = checkCCD (phdr, sts, &missing)))
		return (status);
	}

	if (sts->lorscorr == PERFORM)
	    nsteps++;

	if ((status = checkNonLin (phdr, sts, &missing, &nsteps)))
	    return (status);

	if ((status = checkBias (phdr, sts, &missing, &nsteps)))
	    return (status);

	if ((status = checkDark (phdr, sts, &missing, &nsteps)))
	    return (status);

	if ((status = checkFlat (phdr, sts, &missing, &nsteps)))
	    return (status);

	if ((status = checkShad (phdr, sts, &missing, &nsteps)))
	    return (status);

	if ((status = checkPhot (phdr, sts, &missing, &nsteps)))
	    return (status);

	if (missing) {
	    return (CAL_FILE_MISSING);
	} else if (nsteps < 1) {
	    printf ("Warning  No calibration switch was set to PERFORM, \\\n");
	    printf ("Warning  or all reference files had PEDIGREE = DUMMY.\n");
	    return (NOTHING_TO_DO);
	} else {
	    return (0);
	}
}

/* If this step is to be performed, check for the existence of the
   atod table.
*/

static int checkAtoD (Hdr *phdr, StisInfo1 *sts, int *missing, int *nsteps) {

/* arguments:
Hdr *phdr        i: primary header
StisInfo1 *sts   i: switches, file names, etc
int *missing     io: incremented if the file is missing
int *nsteps      io: incremented if this step can be performed
*/

	int status;

	int calswitch;

	/* Are we supposed to do this step? */
	if (sts->atodcorr == PERFORM) {

	    if ((status = GetSwitch (phdr, "ATODCORR", &calswitch)))
		return (status);
	    if (calswitch == COMPLETE) {
		sts->atodcorr = OMIT;
		return (0);
	    }

	    /* Get the table name, check that the file exists, and get
		pedigree and descrip.
	    */
	    if ((status = GetTabRef (sts->refnames, phdr,
                                     "ATODTAB", &sts->atod, &sts->atodcorr)))
		return (status);
	    if (sts->atod.exists != EXISTS_YES)
		MissingFile ("ATODTAB", sts->atod.name, missing);

	    if (sts->atodcorr == PERFORM)
		(*nsteps)++;
	}

	return (0);
}

/* If this step is to be performed, check for the existence of the
   bias file.  If it exists, get the pedigree and descrip keyword values.
*/

static int checkBias (Hdr *phdr, StisInfo1 *sts, int *missing, int *nsteps) {

/* arguments:
Hdr *phdr        i: primary header
StisInfo1 *sts   i: switches, file names, etc
int *missing     io: incremented if the file is missing
int *nsteps      io: incremented if this step can be performed
*/

	int status;

	int calswitch;

	/* Are we supposed to do this step? */
	if (sts->biascorr == PERFORM) {

	    if ((status = GetSwitch (phdr, "BIASCORR", &calswitch)))
		return (status);
	    if (calswitch == COMPLETE) {
		sts->biascorr = OMIT;
		return (0);
	    }

	    if ((status = GetImageRef (sts->refnames, phdr,
                                       "BIASFILE", &sts->bias, &sts->biascorr)))
		return (status);
	    if (sts->bias.exists != EXISTS_YES)
		MissingFile ("BIASFILE", sts->bias.name, missing);

	    if (sts->biascorr == PERFORM)
		(*nsteps)++;
	}

	return (0);
}

/* There's no reference file for this step. */

static int checkBlev (Hdr *phdr, StisInfo1 *sts, int *nsteps) {

/* arguments:
Hdr *phdr        i: primary header
StisInfo1 *sts   i: switches, file names, etc
int *nsteps      io: incremented if this step can be performed
*/

	int status;

	int calswitch;

	/* Are we supposed to do this step? */
	if (sts->blevcorr == PERFORM) {

	    if ((status = GetSwitch (phdr, "BLEVCORR", &calswitch)))
		return (status);
	    if (calswitch == COMPLETE)
		sts->blevcorr = OMIT;

	    if (sts->blevcorr == PERFORM)
		(*nsteps)++;
	}

	return (0);
}

/* If the detector is the CCD, we need the CCD parameters table for
   BIASCORR, DARKCORR, and PHOTCORR.  This routine checks that the table
   exists.

   We also need the table for initializing the error array, but we
   don't have a flag for that step.  That's why we need this table
   regardless of which steps are to be performed.
*/

static int checkCCD (Hdr *phdr, StisInfo1 *sts, int *missing) {

/* arguments:
Hdr *phdr        i: primary header
StisInfo1 *sts   i: switches, file names, etc
int *missing     io: incremented if the table is missing
*/

	int status;
	int calswitch;			/* returned by GetTabRef and ignored */

	if (sts->detector != CCD_DETECTOR)
	    return (0);

	if ((status = GetTabRef (sts->refnames, phdr, "CCDTAB",
                                 &sts->ccdpar, &calswitch)))
	    return (status);

	if (sts->ccdpar.exists != EXISTS_YES) {

	    MissingFile ("CCDTAB", sts->ccdpar.name, missing);

	} else if (sts->ccdpar.goodPedigree != GOOD_PEDIGREE) {

	    (*missing)++;
	    printf (
		"ERROR    CCDTAB `%s' is a dummy table.\n", sts->ccdpar.name);
	}
	return (0);
}

/* If this step is to be performed, check for the existence of the
   dark file.  If it exists, get the pedigree and descrip keyword values.
   This function also checks for the existence of the tdc table in the case
   of a NUV MAMA detector.
*/

static int checkDark (Hdr *phdr, StisInfo1 *sts, int *missing, int *nsteps) {

/* arguments:
Hdr *phdr        i: primary header
StisInfo1 *sts   i: switches, file names, etc
int *missing     io: incremented if the file is missing
int *nsteps      io: incremented if this step can be performed
*/

	int status, calswitch;

	if (sts->darkcorr == PERFORM) {

	    if ((status = GetSwitch (phdr, "DARKCORR", &calswitch)))
		return (status);
	    if (calswitch == COMPLETE) {
		sts->darkcorr = OMIT;
		return (status = 0);
	    }

	    if ((status = GetImageRef (sts->refnames, phdr,
                                       "DARKFILE", &sts->dark, &sts->darkcorr)))
		return (status);
	    if (sts->dark.exists != EXISTS_YES)
		MissingFile ("DARKFILE", sts->dark.name, missing);

            if (sts->detector == NUV_MAMA_DETECTOR) {

		if ((status = GetTabRef (sts->refnames, phdr,
                                         "TDCTAB", &sts->tdctab, &calswitch)))
		    return (status);

		if (sts->tdctab.exists != EXISTS_YES) {

		    MissingFile ("TDCTAB", sts->tdctab.name, missing);

		} else if (sts->tdctab.goodPedigree != GOOD_PEDIGREE) {

                    (*missing)++;
                    printf ("ERROR    TDCTAB `%s' is a dummy table.\n",
                         sts->tdctab.name);
		}
	    }

	    if (sts->darkcorr == PERFORM)
		(*nsteps)++;
	}

	return (0);
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

static int checkDQI (Hdr *phdr, StisInfo1 *sts, int *missing, int *nsteps) {

/* arguments:
Hdr *phdr        i: primary header
StisInfo1 *sts   i: switches, file names, etc
int *missing     io: incremented if the table is missing
int *nsteps      io: incremented if this step can be performed
*/

	int status;

	if (sts->dqicorr == PERFORM) {

	    if ((status = GetTabRef (sts->refnames, phdr,
                                     "BPIXTAB", &sts->bpix, &sts->dqicorr)))
		return (status);

	    if (sts->bpix.exists != EXISTS_YES) {

		if (sts->detector != CCD_DETECTOR ||
			GotFileName (sts->bpix.name)) {

		    MissingFile ("BPIXTAB", sts->bpix.name, missing);
		}
	    }

	    if (sts->dqicorr == PERFORM)
		(*nsteps)++;
	}

	return (0);
}

/* If this step is to be performed, check for the existence of the
   flat field files.  If they exist, get the pedigree and descrip
   keyword values.  If pedigree is DUMMY, the flag may be reset.
*/

static int checkFlat (Hdr *phdr, StisInfo1 *sts, int *missing, int *nsteps) {

/* arguments:
Hdr *phdr        i: primary header
StisInfo1 *sts   i: switches, file names, etc
int *missing     io: incremented if a file is missing
int *nsteps      io: incremented if this step can be performed
*/

	int status;

	int calswitch;

	/* Are we supposed to do this step? */
	if (sts->flatcorr == PERFORM) {

	    if ((status = GetSwitch (phdr, "FLATCORR", &calswitch)))
		return (status);
	    if (calswitch == COMPLETE) {
		sts->flatcorr = OMIT;
		return (status = 0);
	    }

	    /* Initial values; may be reset below. */
	    sts->pfltcorr = PERFORM;
	    sts->dfltcorr = PERFORM;
	    sts->lfltcorr = PERFORM;

	    if ((status = GetImageRef (sts->refnames, phdr,
                                       "PFLTFILE", &sts->pflt, &sts->pfltcorr)))
		return (status);
	    if (sts->pflt.exists != EXISTS_YES) {
		if (GotFileName (sts->pflt.name)) {	/* name specified? */
		    MissingFile ("PFLTFILE", sts->pflt.name, missing);
		} else {
		    sts->pfltcorr = OMIT;	/* name was blank or "N/A" */
		}
	    }

	    if ((status = GetImageRef (sts->refnames, phdr,
                                       "DFLTFILE", &sts->dflt, &sts->dfltcorr)))
		return (status);
	    if (sts->dflt.exists != EXISTS_YES) {
		if (GotFileName (sts->dflt.name)) {
		    MissingFile ("DFLTFILE", sts->dflt.name, missing);
		} else {
		    sts->dfltcorr = OMIT;
		}
	    }

	    if ((status = GetImageRef (sts->refnames, phdr,
                                       "LFLTFILE", &sts->lflt, &sts->lfltcorr)))
		return (status);
	    if (sts->lflt.exists != EXISTS_YES) {
		if (GotFileName (sts->lflt.name)) {
		    MissingFile ("LFLTFILE", sts->lflt.name, missing);
		} else {
		    sts->lfltcorr = OMIT;
		}
	    }

	    /* If any of the three parts of flat fielding is set to
		PERFORM, then we can do this step.  If not, and if any
		part is DUMMY because of the reference file, reset the
		flat field flag to DUMMY; this will mean that all the
		files that were specified have pedigree=dummy.
	    */
	    if (sts->pfltcorr == PERFORM || sts->dfltcorr == PERFORM ||
		sts->lfltcorr == PERFORM) {
		(*nsteps)++;
	    } else if (sts->pfltcorr == OMIT && sts->dfltcorr == OMIT &&
		sts->lfltcorr == OMIT) {
		(*missing)++;
		printf (
		"ERROR    PFLTFILE, DFLTFILE, and LFLTFILE are all blank.\n");
	    } else if (sts->pfltcorr == DUMMY || sts->dfltcorr == DUMMY ||
		sts->lfltcorr == DUMMY) {
		sts->flatcorr = DUMMY;
	    }
	}

	return (0);
}

/* Check whether we should check for and possibly correct nonlinearity.
   There is a reference table but not an image for this step.  The table
   is opened, and the values are read into the StisInfo1 structure.
*/

static int checkNonLin (Hdr *phdr, StisInfo1 *sts, int *missing, int *nsteps) {

/* arguments:
Hdr *phdr        i: primary header
StisInfo1 *sts   i: switches, file names, etc
int *missing     io: incremented if the table is missing
int *nsteps      io: incremented if this step can be performed
*/

	int status;

	int calswitch;

	if (sts->glincorr == PERFORM || sts->lflgcorr == PERFORM) {

	    if ((status = GetSwitch (phdr, "GLINCORR", &calswitch)))
		return (status);
	    if (calswitch == COMPLETE)
		sts->glincorr = OMIT;

	    if ((status = GetSwitch (phdr, "LFLGCORR", &calswitch)))
		return (status);
	    if (calswitch == COMPLETE)
		sts->lflgcorr = OMIT;

	    if (sts->glincorr != PERFORM && sts->lflgcorr != PERFORM)
		return (0);

	    /* Get the name of the MAMA linearity information table,
		and check whether it exists.  calswitch is ignored;
		glincorr and lflgcorr may be reset to DUMMY below.
	    */
	    if ((status = GetTabRef (sts->refnames, phdr,
                                     "MLINTAB", &sts->mlin, &calswitch)))
		return (status);
	    if (sts->mlin.exists != EXISTS_YES)
		MissingFile ("MLINTAB", sts->mlin.name, missing);

	    if (sts->mlin.goodPedigree != GOOD_PEDIGREE) {
		if (sts->glincorr == PERFORM)
		    sts->glincorr = DUMMY;
		if (sts->lflgcorr == PERFORM)
		    sts->lflgcorr = DUMMY;
	    }

	    if (sts->glincorr == PERFORM || sts->lflgcorr == PERFORM)
		(*nsteps)++;
	}

	return (0);
}

/* Check whether we should compute the photometry keyword values.
   There is a reference table but not an image for this step.
   This will only be done for imaging mode, not spectroscopic.
   Also, unlike most other switches, photcorr will not be reset
   to omit if the header value is "COMPLETE", since it would not
   cause any problem to perform this step more than once.
*/

static int checkPhot (Hdr *phdr, StisInfo1 *sts, int *missing, int *nsteps) {

/* arguments:
Hdr *phdr        i: primary header
StisInfo1 *sts   i: switches, file names, etc
int *missing     io: incremented if the table is missing
int *nsteps      io: incremented if this step can be performed
*/

	int status;

	if (sts->photcorr == PERFORM) {

	    if (strcmp (sts->obstype, "IMAGING") != 0) {
		sts->photcorr = OMIT;
		return (0);
	    }

	    /* Get the name of the _pht throughput table. */
	    if ((status = GetTabRef (sts->refnames, phdr,
                                     "IMPHTTAB", &sts->phot, &sts->photcorr)))
		return (status);
	    if (sts->phot.exists != EXISTS_YES) {
		MissingFile ("IMPHTTAB", sts->phot.name, missing);
		sts->photcorr = OMIT;
	    }

	    /* Get the name of the _tds sensitivity table */
	    if ((status = GetTabRef(sts->refnames, phdr, "TDSTAB",
                                    &sts->tdstab, &sts->tdscorr)))
		return (status);
	    if (sts->tdstab.exists != EXISTS_YES) {
		if (GotFileName (sts->tdstab.name))
		    MissingFile("TDSTAB", sts->tdstab.name, missing);
		sts->tdscorr = OMIT;
	    }

	    if (sts->photcorr == PERFORM)
		(*nsteps)++;
	}

	return (0);
}

/* If this step is to be performed, check for the existence of the
   shutter shading file.  If it exists, get the pedigree and descrip
   keyword values.
*/

static int checkShad (Hdr *phdr, StisInfo1 *sts, int *missing, int *nsteps) {

/* arguments:
Hdr *phdr        i: primary header
StisInfo1 *sts   i: switches, file names, etc
int *missing     io: incremented if the file is missing
int *nsteps      io: incremented if this step can be performed
*/

	int status;

	int calswitch;

	if (sts->shadcorr == PERFORM) {

	    if ((status = GetSwitch (phdr, "SHADCORR", &calswitch)))
		return (status);
	    if (calswitch == COMPLETE) {
		sts->shadcorr = OMIT;
		return (0);
	    }

	    if ((status = GetImageRef (sts->refnames, phdr,
                                       "SHADFILE", &sts->shad, &sts->shadcorr)))
		return (status);
	    if (sts->shad.exists != EXISTS_YES)
		MissingFile ("SHADFILE", sts->shad.name, missing);

	    if (sts->shadcorr == PERFORM)
		(*nsteps)++;
	}

	return (0);
}

/* This routine gets the name of a reference image, checks whether it
   exists, and gets pedigree and descrip if they are present.  The image
   name will be null if the keyword is not present in the header; this is
   not an error.
*/

static int GetImageRef (RefFileInfo *refnames, Hdr *phdr,
		char *keyword, RefImage *image, int *calswitch) {

	int status;

	/* Get the reference image name. */
	if ((status = GetRefName (refnames, phdr, keyword, image->name)))
	    return (status);

	/* ImgPedigree opens the image to verify that it exists, and if so,
	   gets pedigree & descrip.
	*/
	if ((status = ImgPedigree (image)))
	    return (status);
	if (image->exists == EXISTS_YES) {
	    if (image->goodPedigree != GOOD_PEDIGREE)
		*calswitch = DUMMY;
	}

	return (0);
}

/* This routine gets the name of a reference table, checks whether it
   exists, and gets pedigree and descrip if they are present.  The table
   name will be null if the keyword is not present in the header; this is
   not an error.
*/

static int GetTabRef (RefFileInfo *refnames, Hdr *phdr,
		char *keyword, RefTab *table, int *calswitch) {

	int status;

	/* Get the reference table name. */
	if ((status = GetRefName (refnames, phdr, keyword, table->name)))
	    return (status);

	/* TabPedigree opens the table to verify that it exists, and if so,
	   gets pedigree & descrip.
	*/
	if ((status = TabPedigree (table)))
	    return (status);

	if (table->exists == EXISTS_YES) {
	    if (table->goodPedigree != GOOD_PEDIGREE)
		*calswitch = DUMMY;
	}

	return (0);
}

static void MissingFile (char *keyword, char *filename, int *missing) {

	printf ("ERROR    %s `%s' not found or can't open.\n",
				keyword, filename);
	(*missing)++;
}
