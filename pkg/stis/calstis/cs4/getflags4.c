/* This file contains:
	GetFlags4
	CheckWave
	CheckX2D
	GetCheckRef4
*/

# include <stdio.h>

# include "c_iraf.h"
# include "hstio.h"

# include "stis.h"
# include "calstis4.h"
# include "err.h"		/* defines error codes */
# include "stisdef.h"

static int CheckWave (Hdr *, StisInfo4 *, int *, int *);
static int CheckX2D (Hdr *, StisInfo4 *);
static int GetCheckRef4 (RefFileInfo *, Hdr *, char *, RefTab *, int *, int *);

/* This routine gets the names of reference tables from the primary
   header, checks that they exist, and gets pedigree and descrip.  If
   pedigree indicates that a reference table is DUMMY, the wavecorr
   switch will be reset to DUMMY.

   The X2DCORR switch is checked.  It must be COMPLETE for first-order
   data, and it must not be COMPLETE for echelle or prism data.

   Phil Hodge, 1997 Dec 11:
	Call GetRefName instead of GetKeyS in CheckWave.

   Phil Hodge, 1998 Dec 11:
	Get WCPTAB name.

   Phil Hodge, 2000 Jan 5:
	Get table names for global echelle offset.
	Include GetCheckRef4.

   Phil Hodge, 2001 Mar 7:
	In CheckX2D, check for prism.  In CheckWave, get disptab, etc.,
	for first order as well as for echelle.  Get sdctab for prism.
*/

int GetFlags4 (StisInfo4 *sts, Hdr *phdr) {

	int status;

	int missing = 0;	/* true if any calibration file is missing */
	int nsteps = 0;		/* number of calibration steps to perform */

	sts->wavecorr = PERFORM;	/* initial value; may be reset below */

	/* Check the X2DCORR switch. */
	if ((status = CheckX2D (phdr, sts)))
	    return (status);

	if (sts->wavecorr != PERFORM)
	    return (NOTHING_TO_DO);

	/* Check the reference files. */
	if ((status = CheckWave (phdr, sts, &missing, &nsteps)))
	    return (status);

	if (missing) {
	    return (CAL_FILE_MISSING);
	} else if (nsteps < 1) {
	    return (NOTHING_TO_DO);
	} else {
	    return (0);
	}
}

/* This routine gets the names of the reference files needed, depending
   on whether we have first-order or echelle data.  The tables are opened
   to verify that they exist and that they are not dummy.
*/

static int CheckWave (Hdr *phdr, StisInfo4 *sts, int *missing, int *nsteps) {

/* arguments:
Hdr *phdr        i: primary header
StisInfo4 *sts   i: switches, file names, etc
int *missing    io: incremented if the file is missing
int *nsteps     io: incremented if this step can be performed
*/

	int status;

	/* Parameters that control wavecal processing.
	   We're not using GetCheckRef4 for wcptab because it's OK for
	   the keyword to not exist or for the name to be left blank.
	*/
	if ((status = GetRefName (sts->refnames, phdr,
                                  "WCPTAB", sts->wcptab.name)))
	    return (status);

	/* Open the table to verify that it exists, and if it does,
	    get pedigree & descrip.
	*/
	if ((status = TabPedigree (&sts->wcptab)))
	    return (status);
	if (sts->wcptab.exists != EXISTS_YES) {
	    if (GotFileName (sts->wcptab.name)) {
		(*missing)++;
		printf ("ERROR    WCPTAB `%s' not found.\n", sts->wcptab.name);
	    }
	} else if (sts->wcptab.goodPedigree != GOOD_PEDIGREE) {
	    printf ("Warning  WCPTAB has PEDIGREE = DUMMY; \\\n");
	    printf ("Warning  default parameters will be used.\n");
	    sts->wcptab.exists = EXISTS_NO;
	}

	/* Spectrum of calibration lamp. */
	if ((status = GetCheckRef4 (sts->refnames, phdr,
                "LAMPTAB", &sts->lamptab, &sts->wavecorr, missing)))
	    return (status);

	/* Aperture description table. */
	if ((status = GetCheckRef4 (sts->refnames, phdr,
                "APDESTAB", &sts->apdestab, &sts->wavecorr, missing)))
	    return (status);

	/* We need additional reference tables for echelle or prism data. */
	if (sts->disp_type == ECHELLE_DISP || sts->disp_type == PRISM_DISP) {

	    /* Dispersion coefficients. */
	    if ((status = GetCheckRef4 (sts->refnames, phdr,
                    "DISPTAB", &sts->disptab, &sts->wavecorr, missing)))
		return (status);

	    /* Incidence-angle correction table. */
	    if ((status = GetCheckRef4 (sts->refnames, phdr,
                    "INANGTAB", &sts->inangtab, &sts->wavecorr, missing)))
		return (status);

	    /* 1-D spectrum trace table. */
	    if ((status = GetCheckRef4 (sts->refnames, phdr,
                    "SPTRCTAB", &sts->sptrctab, &sts->wavecorr, missing)))
		return (status);
	}

	if (sts->disp_type == PRISM_DISP) {
	    if ((status = GetCheckRef4 (sts->refnames, phdr,
                    "SDCTAB", &sts->sdctab, &sts->wavecorr, missing)))
		return (status);
	}

	/* Has the switch not been reset to DUMMY? */
	if (sts->wavecorr == PERFORM)
	    (*nsteps)++;

	return (0);
}

/* Check whether the input image has been 2-D rectified.  If not,
   reset the wavecorr switch.  This is to protect against the user
   unintentionally running cs4 on a flt or crj file.
*/

static int CheckX2D (Hdr *phdr, StisInfo4 *sts) {

/* arguments:
Hdr *phdr       i: primary header
StisInfo4 *sts  i: switches, file names, etc
*/

	int status;

	int x2dcorr;

	if ((status = GetSwitch (phdr, "X2DCORR", &x2dcorr)))
	    return (status);

	if (sts->disp_type == ECHELLE_DISP) {
	    if (x2dcorr == COMPLETE) {
		printf (
	"ERROR    Input file has already been 2-D rectified; \\\n");
		printf (
	"ERROR    for echelle data the input should be the _flt file.\n");
		sts->wavecorr = OMIT;
	    }
	} else if (sts->disp_type == PRISM_DISP) {
	    if (x2dcorr == COMPLETE) {
		printf (
	"ERROR    Input file has already been 2-D rectified; \\\n");
		printf (
	"ERROR    for prism data the input should be the _flt file.\n");
		sts->wavecorr = OMIT;
	    }
	} else {
	    if (x2dcorr != COMPLETE) {
		printf (
	"ERROR    Input file has not been 2-D rectified; \\\n");
		printf (
	"ERROR    for first-order data you must first run x2d.\n");
		sts->wavecorr = OMIT;
	    }
	}

	return (0);
}

/* This routine checks for the existence of a reference table.  If it
   does exist, pedigree and descrip are gotten.  This version differs
   from some others in that we only call this one if the reference table
   must exist; if it doesn't exist or is dummy, a message will be printed
   and *missing will be incremented.
*/

static int GetCheckRef4 (RefFileInfo *refnames, Hdr *phdr,
	char *keyword, RefTab *table, int *calswitch, int *missing) {

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

	    if (table->goodPedigree != GOOD_PEDIGREE) {
		printf ("ERROR    %s has PEDIGREE = DUMMY.\n", table->name);
		*calswitch = DUMMY;
		(*missing)++;
	    }

	} else {

	    printf ("ERROR    %s `%s' not found or can't open.\n",
				keyword, table->name);
	    (*missing)++;
	}

	return (0);
}
