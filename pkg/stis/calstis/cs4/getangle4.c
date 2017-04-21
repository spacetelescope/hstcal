# include <stdio.h>
# include <stdlib.h>
# include <string.h>

# include "c_iraf.h"
# include "hstio.h"
# include "xtables.h"

# include "stis.h"
# include "calstis4.h"
# include "err.h"
# include "stisdef.h"

typedef struct {
	IRAFPointer tp;			/* pointer to table descriptor */
	IRAFPointer cp_aperture;	/* column descriptors */
	IRAFPointer cp_offset;
	IRAFPointer cp_pedigree;
	IRAFPointer cp_descrip;
	int nrows;			/* number of rows in table */
} TblInfo;

typedef struct {
	char aperture[STIS_CBUF+1];
	double offset;			/* offset from ref_aper */
} TblRow;

static int OpenApTab (char *, TblInfo *);
static int ReadApTab (TblInfo *, int, TblRow *);
static int CloseApTab (TblInfo *);

/* This routine gets the offset from the aperture used to measure the
   dispersion relation to the aperture used for the current observation.

   NOTE that this ref_aper will in general be different from the position
   reference aperture.

   The aperture description table should contain the following:
	columns:
		APERTURE:  aperture name (string)
		OFFSET1:   offset from nominal position (float)

   The table is read to find the row for which the value of APERTURE
   is the same as in the input image header.  For that row, OFFSET1
   is read; note that this value is still in arcseconds, not pixels.
   Since only OFFSET1 is gotten, this cannot be used for data with
   dispaxis = 2.

   Phil Hodge, 2000 Jan 14:
	Created from cs7/getapoffset.c.
*/

int GetAngle4 (StisInfo4 *sts, char *ref_aper, double *angle) {

/* arguments:
StisInfo4 *sts   i: calibration switches and info
char *ref_aper   i: name of reference aperture for dispersion relation coeff
double *angle    o: incidence angle, in arcseconds
*/

	int status;

	TblInfo tabinfo;	/* pointer to table descriptor, etc */
	TblRow tabrow;		/* values read from a table row */

	double ap_offset;	/* offset of aperture used for observation */
	double ref_offset;	/* offset of reference aperture */
	int row;		/* number of rows, and loop index */
	int ap_found;		/* true if aperture found in table */
	int ref_found;		/* true if reference aperture found in table */

	/* Open the aperture description table. */
	if ((status = OpenApTab (sts->apdestab.name, &tabinfo)))
	    return (status);

	/* Check each row for a match with aperture or ref_aper. */

	ap_found = 0;
	ref_found = 0;
	for (row = 1;  row <= tabinfo.nrows;  row++) {

	    if ((status = ReadApTab (&tabinfo, row, &tabrow)))
		return (status);

	    if (SameString (tabrow.aperture, sts->aperture)) {
		ap_found = 1;
		ap_offset = tabrow.offset;
		/* Get pedigree & descrip from this row. */
		if ((status = RowPedigree (&sts->apdestab, row,
                        tabinfo.tp, tabinfo.cp_pedigree, tabinfo.cp_descrip)))
		    return (status);
		if (sts->apdestab.goodPedigree == DUMMY_PEDIGREE)
		    printf ("Warning  APDESTAB has PEDIGREE = DUMMY.\n");
	    }

	    if (SameString (tabrow.aperture, ref_aper)) {
		ref_found = 1;
		ref_offset = tabrow.offset;
	    }
	    if (ap_found && ref_found)
		break;
	}

	if ((status = CloseApTab (&tabinfo)))
	    return (status);

	if (!ap_found || !ref_found) {
	    if (!ap_found)
		printf ("APERTURE %s not found in APDESTAB.\n", sts->aperture);
	    if (!ref_found)
		printf ("REF_APER %s not found in APDESTAB.\n", ref_aper);
	    printf ("  APDESTAB = %s\n", sts->apdestab.name);
	    return (TABLE_ERROR);
	}

	*angle = ap_offset - ref_offset;

	return (0);
}

/* This routine opens the aperture description table, finds the columns
   that we need, and gets the total number of rows in the table.
*/

static int OpenApTab (char *tname, TblInfo *tabinfo) {

	tabinfo->tp = c_tbtopn (tname, IRAF_READ_ONLY, 0);
	if (c_iraferr()) {
	    printf ("APDESTAB `%s' not found.\n", tname);
	    return (OPEN_FAILED);
	}

	tabinfo->nrows = c_tbpsta (tabinfo->tp, TBL_NROWS);

	/* Find the columns. */
	c_tbcfnd1 (tabinfo->tp, "APERTURE", &tabinfo->cp_aperture);
	c_tbcfnd1 (tabinfo->tp, "OFFSET1", &tabinfo->cp_offset);
	if (tabinfo->cp_aperture == 0 || tabinfo->cp_offset == 0) {
	    printf ("Column not found in APDESTAB.\n");
	    c_tbtclo (tabinfo->tp);
	    return (COLUMN_NOT_FOUND);
	}

	/* Pedigree and descrip are optional columns. */
	c_tbcfnd1 (tabinfo->tp, "PEDIGREE", &tabinfo->cp_pedigree);
	c_tbcfnd1 (tabinfo->tp, "DESCRIP", &tabinfo->cp_descrip);

	return (0);
}

/* This routine reads the column (aperture name) used to select the
   correct row.
*/

static int ReadApTab (TblInfo *tabinfo, int row, TblRow *tabrow) {

	c_tbegtt (tabinfo->tp, tabinfo->cp_aperture, row,
			tabrow->aperture, STIS_CBUF);
	c_tbegtd (tabinfo->tp, tabinfo->cp_offset, row, &tabrow->offset);
	if (c_iraferr())
	    return (TABLE_ERROR);

	return (0);
}

/* This routine closes the apdestab table. */

static int CloseApTab (TblInfo *tabinfo) {

	c_tbtclo (tabinfo->tp);
	if (c_iraferr())
	    return (TABLE_ERROR);

	return (0);
}
