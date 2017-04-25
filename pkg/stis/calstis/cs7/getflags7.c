# include <stdio.h>
# include <string.h>

# include "c_iraf.h"
# include "hstio.h"

# include "stis.h"
# include "calstis7.h"
# include "hstcalerr.h"	/* defines error codes */
# include "stisdef.h"

static int CheckFlux (Hdr *, StisInfo7 *, int *, int *);
static int CheckSmGeo (Hdr *, StisInfo7 *, int *, int *);
static int CheckWave (Hdr *, StisInfo7 *);
static int CheckWX2D (Hdr *, StisInfo7 *);
static int CheckX2D (Hdr *, StisInfo7 *, int *, int *);
static int GetCheckRef (RefFileInfo *, Hdr *, char *, RefTab *, int *);
static void MissingFile (char *, char *, int *);

/* This routine gets the names of reference images and tables from the
   primary header.  A routine will be called to check the existence of
   the reference files and to get pedigree and descrip.

   Note that the calibration switches themselves are not gotten, except
   for WAVECORR, and that's just for information.

   Phil Hodge, 1997 Dec 11:
	Get pctab in CheckFlux; add sts->refnames to calling sequence
	of GetCheckRef, and change that to call GetRefName.

   Phil Hodge, 1998 Jan 12:
	Modify GetCheckRef and add MissingFile; the latter now prints
	the error message and increments missing.  MissingFile is not
	called for pctab if a file name was not specified.

   Phil Hodge, 2000 June 30:
	In CheckX2D, move the code for apdestab to within the section
	for spectroscopic mode; that is, don't look for apdestab for
	imaging data.

   Phil Hodge, 2000 Aug 9:
	Remove the check on the MAMA offset table (mofftab) in CheckX2D.

   Ivo Busko, 2002 Jan 10:
	Include check on tdstab.

   Phil Hodge, 2006 Sept 21:
	Include a test on WX2DCORR.
*/

int GetFlags7 (StisInfo7 *sts, Hdr *phdr) {

	int status;

	int missing = 0;	/* true if any calibration file is missing */
	int nsteps = 0;		/* number of calibration steps to perform */

	/* Check each reference file. */

	if ((status = CheckX2D (phdr, sts, &missing, &nsteps)))
	    return (status);

	if ((status = CheckSmGeo (phdr, sts, &missing, &nsteps)))
	    return (status);

	if (sts->obstype == SPECTROSCOPIC_TYPE) {

	    if ((status = CheckWave (phdr, sts))) /* just check flag */
		return (status);

	    if ((status = CheckWX2D (phdr, sts))) /* wx2dcorr done? */
		return (status);

	    if (sts->heliocorr == PERFORM)
		nsteps++;

	    if ((status = CheckFlux (phdr, sts, &missing, &nsteps)))
		return (status);

	} else {

	    sts->heliocorr = OMIT;
	    sts->fluxcorr = OMIT;
	    sts->wx2dcorr = OMIT;
	}

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

/* Check whether the wavecal has been used to update the coordinates. */

static int CheckWave (Hdr *phdr, StisInfo7 *sts) {

/* arguments:
Hdr *phdr       i: primary header
StisInfo7 *sts  i: switches, file names, etc
*/

	int status;

	if ((status = GetSwitch (phdr, "WAVECORR", &sts->wavecorr)))
	    return (status);

	return (0);
}

/* Check whether the input file is the output of wx2d (i.e. whether the
   image has already been corrected for the trace).
*/

static int CheckWX2D (Hdr *phdr, StisInfo7 *sts) {

/* arguments:
Hdr *phdr       i: primary header
StisInfo7 *sts  i: switches, file names, etc
*/

	int status;

	if ((status = GetSwitch (phdr, "WX2DCORR", &sts->wx2dcorr)))
	    return (status);

	return (0);
}

/* Get the names of the reference files needed. */

static int CheckX2D (Hdr *phdr, StisInfo7 *sts, int *missing, int *nsteps) {

/* arguments:
Hdr *phdr       i: primary header
StisInfo7 *sts  i: switches, file names, etc
int *missing    io: incremented if the file is missing
int *nsteps     io: incremented if this step can be performed
*/

	int status;

	/* Of course we're supposed to do X2DCORR. */

	if (sts->obstype == SPECTROSCOPIC_TYPE) {

	    /* distntab = SDCTAB */
	    if ((status = GetCheckRef (sts->refnames, phdr,
                    "SDCTAB", &sts->distntab, &sts->x2dcorr)))
		return (status);
	    if (sts->distntab.exists != EXISTS_YES)
		MissingFile ("SDCTAB", sts->distntab.name, missing);

	    /* Aperture description table. */
	    if ((status = GetCheckRef (sts->refnames, phdr,
                    "APDESTAB", &sts->apdestab, &sts->x2dcorr)))
		return (status);
	    if (sts->apdestab.exists != EXISTS_YES)
		MissingFile ("APDESTAB", sts->apdestab.name, missing);

	    /* Dispersion coefficients. */
	    if ((status = GetCheckRef (sts->refnames, phdr,
                    "DISPTAB", &sts->disptab, &sts->x2dcorr)))
		return (status);
	    if (sts->disptab.exists != EXISTS_YES)
		MissingFile ("DISPTAB", sts->disptab.name, missing);

	    /* Incidence-angle correction table. */
	    if ((status = GetCheckRef (sts->refnames, phdr,
                    "INANGTAB", &sts->inangtab, &sts->x2dcorr)))
		return (status);
	    if (sts->inangtab.exists != EXISTS_YES)
		MissingFile ("INANGTAB", sts->inangtab.name, missing);

	    /* Spectrum trace table. */
	    if ((status = GetCheckRef (sts->refnames, phdr,
                    "SPTRCTAB", &sts->sptrctab, &sts->x2dcorr)))
		return (status);
	    if (sts->sptrctab.exists != EXISTS_YES)
		MissingFile ("SPTRCTAB", sts->sptrctab.name, missing);

	} else {			/* imaging mode */

	    /* distntab = IDCTAB */
	    if ((status = GetCheckRef (sts->refnames, phdr,
                    "IDCTAB", &sts->distntab, &sts->x2dcorr)))
                return (status);
	    if (sts->distntab.exists != EXISTS_YES)
		MissingFile ("IDCTAB", sts->distntab.name, missing);
	}

	/* Has the switch not been reset to DUMMY? */
	if (sts->x2dcorr == PERFORM) {
	    (*nsteps)++;
	} else if (sts->x2dcorr == DUMMY) {
	    printf ("Warning  Dummy reference file encountered.\n");
	    return (NOTHING_TO_DO);
	} else {
	    printf ("Warning  X2DCORR not set to PERFORM.\n");
	    return (NOTHING_TO_DO);
	}

	return (0);
}

/* If this step is to be performed, check for the existence of the
   small-scale distortion file sdstfile, and get pedigree and descrip.
*/

static int CheckSmGeo (Hdr *phdr, StisInfo7 *sts, int *missing, int *nsteps) {

/* arguments:
Hdr *phdr       i: primary header
StisInfo7 *sts  i: switches, file names, etc
int *missing    io: incremented if the file is missing
int *nsteps     io: incremented if this step can be performed
*/

	int status;

	/* Are we supposed to do this step? */
	if (sts->sgeocorr == PERFORM) {

	    /* Small-scale distortion file. */
	    if ((status = GetRefName (sts->refnames, phdr, "SDSTFILE",
                                      sts->sdstfile.name)))
		return (status);

	    /* Open the image to verify that it exists, and if it does,
		get pedigree & descrip.
	    */
	    if ((status = ImgPedigree (&sts->sdstfile)))
		return (status);
	    if (sts->sdstfile.exists != EXISTS_YES) {
		(*missing)++;
		printf ("ERROR    SDSTFILE `%s' not found.\n",
			sts->sdstfile.name);
	    }
	    if (sts->sdstfile.goodPedigree != GOOD_PEDIGREE)
		sts->sgeocorr = DUMMY;

	    if (sts->sgeocorr == PERFORM)
		(*nsteps)++;
	}

	return (0);
}

/* Check whether we should convert to absolute flux units.  If so,
   check for the existence of the absolute sensitivity table and
   aperture throughput table.
   Also check for the photometric correction table (pctab), but
   if that table is missing, the pctcorr switch will be reset to
   indicate that that particular correction should not be included.
*/

static int CheckFlux (Hdr *phdr, StisInfo7 *sts, int *missing, int *nsteps) {

/* arguments:
Hdr *phdr       i: primary header
StisInfo7 *sts  i: switches, file names, etc
int *missing    io: incremented if a table is missing
int *nsteps     io: incremented if this step can be performed
*/

	int status;

	if (sts->fluxcorr == PERFORM) {

	    /* Throughput table. */
	    if ((status = GetCheckRef (sts->refnames, phdr,
                    "PHOTTAB", &sts->phottab, &sts->fluxcorr)))
		return (status);
	    if (sts->phottab.exists != EXISTS_YES)
		MissingFile ("PHOTTAB", sts->phottab.name, missing);

	    /* Relative aperture throughput table. */
	    if ((status = GetCheckRef (sts->refnames, phdr,
                    "APERTAB", &sts->apertab, &sts->fluxcorr)))
		return (status);
	    if (sts->apertab.exists != EXISTS_YES)
		MissingFile ("APERTAB", sts->apertab.name, missing);

	    /* Table for corrections for slit height.  Note that we pass
		pctcorr instead of fluxcorr.
	    */
	    if ((status = GetCheckRef (sts->refnames, phdr,
                    "PCTAB", &sts->pctab, &sts->pctcorr)))
		return (status);
	    if (sts->pctab.exists != EXISTS_YES) {
		/* Only print an error message if a name was specified. */
		if (GotFileName (sts->pctab.name))
		    MissingFile ("PCTAB", sts->pctab.name, missing);
		sts->pctcorr = OMIT;
	    }

	    /* Time-dependent sensitivity table. Handled like pct. */
	    sts->tdscorr = PERFORM;
	    if ((status = GetCheckRef (sts->refnames, phdr,
                    "TDSTAB", &sts->tdstab, &sts->tdscorr)))
		return (status);
	    if (sts->tdstab.exists != EXISTS_YES) {
		/* Only print an error message if a name was specified. */
		if (GotFileName (sts->tdstab.name))
		    MissingFile ("TDSTAB", sts->tdstab.name, missing);
		sts->tdscorr = OMIT;
            }

            if (sts->fluxcorr == PERFORM)
                (*nsteps)++;
	}

	return (0);
}

static int GetCheckRef (RefFileInfo *refnames, Hdr *phdr,
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
