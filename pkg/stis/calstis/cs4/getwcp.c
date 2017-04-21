# include <stdio.h>

# include "c_iraf.h"
# include "hstio.h"
# include "xtables.h"

# include "stis.h"
# include "calstis4.h"
# include "hstcalerr.h"
# include "stisdef.h"

typedef struct {
	IRAFPointer tp;			/* pointer to table descriptor */
				/* column descriptors */
	IRAFPointer cp_detector;	/* detector name */
	IRAFPointer cp_opt_elem;	/* name of grating */
	IRAFPointer cp_wl_trim1;
	IRAFPointer cp_wl_trim2;
	IRAFPointer cp_sp_trim1;
	IRAFPointer cp_sp_trim2;
	IRAFPointer cp_wl_range;
	IRAFPointer cp_sp_range;
	IRAFPointer cp_nsigma_cr;
	IRAFPointer cp_nsigma_illum;
	IRAFPointer cp_mad_reject;
	IRAFPointer cp_min_mad;
	IRAFPointer cp_pedigree;
	IRAFPointer cp_descrip;
	int nrows;			/* number of rows in table */
} TblInfo;

typedef struct {
	char detector[STIS_CBUF+1];	/* detector name */
	char opt_elem[STIS_CBUF+1];	/* grating name */
} TblRow;

static int OpenWCPTab (char *, TblInfo *);
static int ReadWCPTab (TblInfo *, int, TblRow *);
static int ReadWCPArray (TblInfo *, int, StisInfo4 *);
static int CloseWCPTab (TblInfo *);

/* This routine gets the wavecal parameters from the WCPTAB.

   The wavecal parameters table should contain the following:
	header parameters:
		none needed
	columns:
		DETECTOR:  detector name (string)
		OPT_ELEM:  grating name (string)
		WL_TRIM1:  pixels to trim from first axis (int)
		WL_TRIM2:  pixels to trim from second axis (int)
		SP_TRIM1:  pixels to trim from first axis (int)
		SP_TRIM2:  pixels to trim from second axis (int)
		WL_RANGE:  search range for wavelength (int)
		SP_RANGE:  search range in spatial direction (int)
		NSIG_CR:   number of sigma for rejecting cosmic rays (double)
		NSIG_ILL:  number of sigma in illuminated regions (double)
		MAD_REJ::  number of MAD for rejecting outliers (double)
		MIN_MAD:   minimum value of MAD (double)

   The table is read to find the row for which the detector and grating
   are the same as that used for the observation.

   Phil Hodge, 1998 Dec 11:
	Function created.

   Phil Hodge, 2000 Jan 13:
	Add one to detector and opt_elem buffer sizes.
*/

int GetWCP (StisInfo4 *sts) {

/* argument:
StisInfo4 *sts    i: calibration switches and info
*/

	int status;

	TblInfo tabinfo;	/* pointer to table descriptor, etc */
	TblRow tabrow;		/* values read from a table row */

	int row;		/* loop index */
	int foundit;		/* true if detector name found in table */

	/* Open the WCPTAB table. */
	if ((status = OpenWCPTab (sts->wcptab.name, &tabinfo)))
	    return (status);

	foundit = 0;
	for (row = 1;  row <= tabinfo.nrows;  row++) {

	    if ((status = ReadWCPTab (&tabinfo, row, &tabrow)))
		return (status);

	    if (SameString (tabrow.detector, sts->det) &&
		SameString (tabrow.opt_elem, sts->opt_elem)) {

		foundit = 1;

		/* Get pedigree & descrip from the row. */
		if ((status = RowPedigree (&sts->wcptab, row,
                        tabinfo.tp, tabinfo.cp_pedigree, tabinfo.cp_descrip)))
		    return (status);
		if (sts->wcptab.goodPedigree == DUMMY_PEDIGREE) {
		    printf ("Warning  WCPTAB has PEDIGREE = DUMMY; \\\n");
		    printf ("Warning  default parameters will be used.\n");
		    break;
		}

		/* Read parameters into sts. */
		if ((status = ReadWCPArray (&tabinfo, row, sts)))
		    return (status);

		break;
	    }
	}

	if ((status = CloseWCPTab (&tabinfo)))
	    return (status);

	if (!foundit) {
	    printf ("ERROR    DETECTOR %s, OPT_ELEM %s not found in %s.\n",
		sts->det, sts->opt_elem, sts->wcptab.name);
	    return (GENERIC_ERROR_CODE);
	}

	return (0);
}

/* This routine opens the wavecal parameters table, finds the columns
   that we need, and gets the total number of rows in the table.
*/

static int OpenWCPTab (char *tname, TblInfo *tabinfo) {

	tabinfo->tp = c_tbtopn (tname, IRAF_READ_ONLY, 0);
	if (c_iraferr()) {
	    printf ("ERROR    WCPTAB `%s' not found.\n", tname);
	    return (OPEN_FAILED);
	}

	tabinfo->nrows = c_tbpsta (tabinfo->tp, TBL_NROWS);

	/* Find the columns. */
	c_tbcfnd1 (tabinfo->tp, "DETECTOR", &tabinfo->cp_detector);
	c_tbcfnd1 (tabinfo->tp, "OPT_ELEM", &tabinfo->cp_opt_elem);
	c_tbcfnd1 (tabinfo->tp, "WL_TRIM1", &tabinfo->cp_wl_trim1);
	c_tbcfnd1 (tabinfo->tp, "WL_TRIM2", &tabinfo->cp_wl_trim2);
	c_tbcfnd1 (tabinfo->tp, "SP_TRIM1", &tabinfo->cp_sp_trim1);
	c_tbcfnd1 (tabinfo->tp, "SP_TRIM2", &tabinfo->cp_sp_trim2);
	c_tbcfnd1 (tabinfo->tp, "WL_RANGE", &tabinfo->cp_wl_range);
	c_tbcfnd1 (tabinfo->tp, "SP_RANGE", &tabinfo->cp_sp_range);
	c_tbcfnd1 (tabinfo->tp, "NSIG_CR", &tabinfo->cp_nsigma_cr);
	c_tbcfnd1 (tabinfo->tp, "NSIG_ILL", &tabinfo->cp_nsigma_illum);
	c_tbcfnd1 (tabinfo->tp, "MAD_REJ", &tabinfo->cp_mad_reject);
	c_tbcfnd1 (tabinfo->tp, "MIN_MAD", &tabinfo->cp_min_mad);

	if (tabinfo->cp_detector == 0 ||
	    tabinfo->cp_opt_elem == 0 ||
	    tabinfo->cp_wl_trim1 == 0 ||
	    tabinfo->cp_wl_trim2 == 0 ||
	    tabinfo->cp_sp_trim1 == 0 ||
	    tabinfo->cp_sp_trim2 == 0 ||
	    tabinfo->cp_wl_range == 0 ||
	    tabinfo->cp_sp_range == 0 ||
	    tabinfo->cp_nsigma_cr == 0 ||
	    tabinfo->cp_nsigma_illum == 0 ||
	    tabinfo->cp_mad_reject == 0 ||
	    tabinfo->cp_min_mad == 0) {
	    printf ("ERROR    Column not found in WCPTAB.\n");
	    c_tbtclo (tabinfo->tp);
	    return (COLUMN_NOT_FOUND);
	}

	/* Pedigree and descrip are optional columns. */
	c_tbcfnd1 (tabinfo->tp, "PEDIGREE", &tabinfo->cp_pedigree);
	c_tbcfnd1 (tabinfo->tp, "DESCRIP", &tabinfo->cp_descrip);

	return (0);
}

/* This routine reads the DETECTOR and OPT_ELEM, which are used to select
   the correct row.
*/

static int ReadWCPTab (TblInfo *tabinfo, int row, TblRow *tabrow) {

	c_tbegtt (tabinfo->tp, tabinfo->cp_detector, row,
			tabrow->detector, STIS_CBUF);
	if (c_iraferr())
	    return (TABLE_ERROR);

	c_tbegtt (tabinfo->tp, tabinfo->cp_opt_elem, row,
			tabrow->opt_elem, STIS_CBUF);
	if (c_iraferr())
	    return (TABLE_ERROR);

	return (0);
}

/* Read the info into the sts structure.
   This routine doesn't read any arrays; the name is used for consistency
   with other functions that read the contents of a selected row.
*/

static int ReadWCPArray (TblInfo *tabinfo, int row, StisInfo4 *sts) {

	c_tbegti (tabinfo->tp, tabinfo->cp_wl_trim1, row, &sts->wl_trim1);
	if (c_iraferr())
	    return (TABLE_ERROR);

	c_tbegti (tabinfo->tp, tabinfo->cp_wl_trim2, row, &sts->wl_trim2);
	if (c_iraferr())
	    return (TABLE_ERROR);

	c_tbegti (tabinfo->tp, tabinfo->cp_sp_trim1, row, &sts->sp_trim1);
	if (c_iraferr())
	    return (TABLE_ERROR);

	c_tbegti (tabinfo->tp, tabinfo->cp_sp_trim2, row, &sts->sp_trim2);
	if (c_iraferr())
	    return (TABLE_ERROR);

	c_tbegti (tabinfo->tp, tabinfo->cp_wl_range, row, &sts->wl_range);
	if (c_iraferr())
	    return (TABLE_ERROR);

	c_tbegti (tabinfo->tp, tabinfo->cp_sp_range, row, &sts->sp_range);
	if (c_iraferr())
	    return (TABLE_ERROR);

	c_tbegtd (tabinfo->tp, tabinfo->cp_nsigma_cr, row, &sts->nsigma_cr);
	if (c_iraferr())
	    return (TABLE_ERROR);

	c_tbegtd (tabinfo->tp, tabinfo->cp_nsigma_illum, row,
		&sts->nsigma_illum);
	if (c_iraferr())
	    return (TABLE_ERROR);

	c_tbegtd (tabinfo->tp, tabinfo->cp_mad_reject, row, &sts->mad_reject);
	if (c_iraferr())
	    return (TABLE_ERROR);

	c_tbegtd (tabinfo->tp, tabinfo->cp_min_mad, row, &sts->min_mad);
	if (c_iraferr())
	    return (TABLE_ERROR);

	return (0);
}

/* This routine closes the wcptab table. */

static int CloseWCPTab (TblInfo *tabinfo) {

	c_tbtclo (tabinfo->tp);
	if (c_iraferr())
	    return (TABLE_ERROR);

	return (0);
}
