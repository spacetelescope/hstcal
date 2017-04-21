# include <stdio.h>
# include <stdlib.h>
# include <string.h>

# include "c_iraf.h"
# include "hstio.h"
# include "xtables.h"

# include "stis.h"
# include "calstis7.h"
# include "err.h"

typedef struct {
	IRAFPointer tp;			/* pointer to table descriptor */
	IRAFPointer cp_aperture;	/* column descriptors */
	IRAFPointer cp_offset[2];
	IRAFPointer cp_pedigree;
	int nrows;			/* number of rows in table */
} TblInfo;

typedef struct {
	char aperture[STIS_CBUF+1];
	double offset[2];		/* offset from ref_aper */
} TblRow;

static int OpenApTab (char *, TblInfo *);
static int ReadApTab (TblInfo *, int, TblRow *);
static int CheckPedigree (TblInfo *, int, int *);
static int CloseApTab (TblInfo *);

/* This routine gets the offset from the aperture used to measure the
   dispersion relation to the aperture used for the current observation.

   NOTE that this ref_aper will in general be different from the position
   reference aperture.

   The aperture description table should contain the following:
	columns:
		APERTURE:  aperture name (string)
		OFFSET1, OFFSET2:  offset from nominal position (float)

   The table is read to find the row for which the value of APERTURE
   is the same as in the input image header.  For that row, OFFSET1
   is read; note that this value is still in arcseconds, not pixels.

   Phil Hodge, 2000 Jan 13:
	Add one to aperture buffer size.
*/

int GetApOffset (StisInfo7 *sts, ApInfo *slit, char *ref_aper, double *delta) {

/* arguments:
StisInfo7 *sts   i: calibration switches and info
ApInfo *slit     i: description of slit (for ap_offset)
char *ref_aper   i: name of reference aperture for dispersion relation coeff
double *delta    o: offset in dispersion direction, arcseconds
*/

	int status;

	TblInfo tabinfo;	/* pointer to table descriptor, etc */
	TblRow tabrow;		/* values read from a table row */

	int row;		/* number of rows, and loop index */
	int pedigree;		/* for selected row */
	int foundit = 0;	/* true if aperture found in table */

	/* Open the aperture description table. */
	if ((status = OpenApTab (sts->apdestab.name, &tabinfo)))
	    return (status);

	/* Check each row for a match with ref_aper. */

	for (row = 1;  row <= tabinfo.nrows;  row++) {

	    if ((status = ReadApTab (&tabinfo, row, &tabrow)))
		return (status);

	    if (SameString (tabrow.aperture, ref_aper)) {

		foundit = 1;
		if ((status = CheckPedigree (&tabinfo, row, &pedigree)))
		    return (status);
		if (pedigree == DUMMY_PEDIGREE) {
		    printf ("Warning  DUMMY pedigree in row %d of %s.\n",
			row, sts->apdestab.name);
		    *delta = 0.;
		    sts->x2dcorr_o = DUMMY;
		    CloseApTab (&tabinfo);
		    return (0);
		} else {
		    if (sts->dispaxis == 1)
			*delta = slit->ap_offset[0] - tabrow.offset[0];
		    else if (sts->dispaxis == 2)
			*delta = slit->ap_offset[1] - tabrow.offset[1];
		    else
			*delta = 0.;
		}

		break;
	    }
	}

	if ((status = CloseApTab (&tabinfo)))
	    return (status);

	if (!foundit) {
	    printf ("Warning  APERTURE %s not found in APDESTAB %s.\n",
		ref_aper, sts->apdestab.name);
	    sts->x2dcorr_o = OMIT;
	}

	return (0);
}

/* This routine opens the aperture description table, finds the columns
   that we need, and gets the total number of rows in the table.
*/

static int OpenApTab (char *tname, TblInfo *tabinfo) {

	tabinfo->tp = c_tbtopn (tname, IRAF_READ_ONLY, 0);
	if (c_iraferr()) {
	    printf ("ERROR    APDESTAB `%s' not found.\n", tname);
	    return (OPEN_FAILED);
	}

	tabinfo->nrows = c_tbpsta (tabinfo->tp, TBL_NROWS);

	/* Find the columns. */
	c_tbcfnd1 (tabinfo->tp, "APERTURE", &tabinfo->cp_aperture);
	c_tbcfnd1 (tabinfo->tp, "OFFSET1", &tabinfo->cp_offset[0]);
	c_tbcfnd1 (tabinfo->tp, "OFFSET2", &tabinfo->cp_offset[1]);
	if (tabinfo->cp_aperture == 0 ||
	    tabinfo->cp_offset[0] == 0 || tabinfo->cp_offset[1] == 0) {
	    printf ("ERROR    Column not found in APDESTAB.\n");
	    c_tbtclo (tabinfo->tp);
	    return (COLUMN_NOT_FOUND);
	}

	/* Pedigree is not a required column. */
	c_tbcfnd1 (tabinfo->tp, "PEDIGREE", &tabinfo->cp_pedigree);

	return (0);
}

/* This routine reads the column (aperture name) used to select the
   correct row.
*/

static int ReadApTab (TblInfo *tabinfo, int row, TblRow *tabrow) {

	c_tbegtt (tabinfo->tp, tabinfo->cp_aperture, row,
			tabrow->aperture, STIS_CBUF);
	c_tbegtd (tabinfo->tp, tabinfo->cp_offset[0], row, &tabrow->offset[0]);
	c_tbegtd (tabinfo->tp, tabinfo->cp_offset[1], row, &tabrow->offset[1]);
	if (c_iraferr())
	    return (TABLE_ERROR);

	return (0);
}

/* This routine reads the pedigree from the current row, if the
   column exists, and checks for DUMMY.
*/

static int CheckPedigree (TblInfo *tabinfo, int row, int *pedigree) {

	char *str_pedigree;

	if (tabinfo->cp_pedigree > 0) {

	    if ((str_pedigree = calloc (STIS_LINE+1, sizeof(char))) == NULL) {
		printf ("ERROR    Out of memory.\n");
		return (OUT_OF_MEMORY);
	    }
	    c_tbegtt (tabinfo->tp, tabinfo->cp_pedigree, row,
			str_pedigree, STIS_LINE);
	    if (c_iraferr())
		return (TABLE_ERROR);

	    if (strncmp (str_pedigree, "DUMMY", 5) == 0)
		*pedigree = DUMMY_PEDIGREE;
	    else
		*pedigree = GOOD_PEDIGREE;
	    free (str_pedigree);

	} else {

	    *pedigree = GOOD_PEDIGREE;
	}

	return (0);
}

/* This routine closes the apdestab table. */

static int CloseApTab (TblInfo *tabinfo) {

	c_tbtclo (tabinfo->tp);
	if (c_iraferr())
	    return (TABLE_ERROR);

	return (0);
}
