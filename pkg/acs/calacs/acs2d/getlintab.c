# include <stdio.h>

# include "hstio.h"
# include "xtables.h"
# include "acs.h"
# include "acsinfo.h"
# include "hstcalerr.h"

typedef struct {
	IRAFPointer tp;			/* pointer to table descriptor */
	IRAFPointer cp_det;		/* column descriptors */
	IRAFPointer cp_global;
	IRAFPointer cp_local;
	IRAFPointer cp_tau;
	IRAFPointer cp_expand;
	IRAFPointer cp_pedigree;
	IRAFPointer cp_descrip;
	int nrows;			/* number of rows in table */
} TblInfo;

typedef struct {
	char detector[ACS_CBUF];	/* detector name */
	double global_limit;
	double local_limit;
	double tau;
	float expand;
} TblRow;

static int OpenLinTab (char *, TblInfo *);
static int ReadLinTab (TblInfo *, int, TblRow *);
static int CloseLinTab (TblInfo *);

/* This routine gets information from the MAMA linearity table.

   The MAMA linearity correction table should contain the following:
   floating point header parameter:
		none needed
   and the following double precision columns:
		DETECTOR:  "NUV_MAMA" or "FUV_MAMA"
		GLOBAL_LIMIT:  count rate resulting in 10% global nonlinearity
		LOCAL_LIMIT:  count rate resulting in 10% local nonlinearity
		TAU:  time constant in global nonlinearity expression
		EXPAND:  radius in high-res pixels
   The table is read to find the row that matches the current detector,
   and that row is read to obtain the maximum allowed count rates
   (counts/second) for the global and local linearity limits and the time
   constant TAU, which is used for computing the global nonlinearity.
   When flagging local nonlinearity, pixels within a radius of EXPAND
   of any nonlinear pixel will also be flagged.
*/

int GetLinTab (ACSInfo *acs2d) {

/* arguments:
ACS2dInfo *acs2d     io: calibration switches, etc
*/

	extern int status;

	TblInfo tabinfo;	/* pointer to table descriptor, etc */
	TblRow tabrow;		/* values read from a table row */

	int row;			/* loop index */
	int foundit;			/* detector name found in table? */
	int RowPedigree (RefTab *, int, IRAFPointer, IRAFPointer, IRAFPointer);
	int SameString (char *, char *);

	/* Open the MAMA linearity table and find columns. */
	if (OpenLinTab (acs2d->mlin.name, &tabinfo))
	    return (status);

	foundit = 0;

	/* Check each row for a match with detector, and get the info
	   from the matching row.
	*/

	for (row = 1;  row <= tabinfo.nrows;  row++) {

	    /* Read the current row into tabrow. */
	    if (ReadLinTab (&tabinfo, row, &tabrow))
		return (status);

	    if (SameString (tabrow.detector, acs2d->det)) {

		foundit = 1;

		/* Get pedigree & descrip from the row. */
		if (RowPedigree (&acs2d->mlin, row,
			tabinfo.tp, tabinfo.cp_pedigree, tabinfo.cp_descrip))
		    return (status);
		if (acs2d->mlin.goodPedigree == DUMMY_PEDIGREE) {
		    if (acs2d->glincorr == PERFORM)
			acs2d->glincorr = DUMMY;
		    if (acs2d->lflgcorr == PERFORM)
			acs2d->lflgcorr = DUMMY;
		    CloseLinTab (&tabinfo);
		    return (status);
		}

		acs2d->global_limit = tabrow.global_limit;
		acs2d->local_limit = tabrow.local_limit;
		acs2d->tau = tabrow.tau;
		acs2d->expand = tabrow.expand;
		break;
	    }
	}

	if (!foundit) {
	    sprintf (MsgText, "Detector `%s' not found in MLINTAB `%s'.",
		acs2d->det, acs2d->mlin.name);
	    trlerror (MsgText);
		CloseLinTab (&tabinfo);
	    return (status = TABLE_ERROR);
	}

	if (CloseLinTab (&tabinfo))		/* close the table */
	    return (status);

	return (status);
}

/* This routine opens the MAMA linearity table, finds the columns that we
   need, gets one header keyword, and gets the total number of rows in the
   table.  The columns are DETECTOR, GLOBAL_LIMIT, LOCAL_LIMIT, and TAU.
*/

static int OpenLinTab (char *tname, TblInfo *tabinfo) {

	extern int status;

	tabinfo->tp = c_tbtopn (tname, IRAF_READ_ONLY, 0);
	if (c_iraferr()) {
	    sprintf (MsgText, "MLINTAB `%s' not found.", tname);
	    trlerror (MsgText);
		return (status = OPEN_FAILED);
	}

	tabinfo->nrows = c_tbpsta (tabinfo->tp, TBL_NROWS);

	/* Find the columns. */
	c_tbcfnd1 (tabinfo->tp, "DETECTOR", &tabinfo->cp_det);
	c_tbcfnd1 (tabinfo->tp, "GLOBAL_LIMIT", &tabinfo->cp_global);
	c_tbcfnd1 (tabinfo->tp, "LOCAL_LIMIT", &tabinfo->cp_local);
	c_tbcfnd1 (tabinfo->tp, "TAU", &tabinfo->cp_tau);
	c_tbcfnd1 (tabinfo->tp, "EXPAND", &tabinfo->cp_expand);
	if (tabinfo->cp_det == 0 ||
	    tabinfo->cp_global == 0 ||
	    tabinfo->cp_local == 0 ||
	    tabinfo->cp_tau == 0 ||
	    tabinfo->cp_expand == 0) {
	    trlerror ("Column not found in MLINTAB.");
	    c_tbtclo (tabinfo->tp);
	    return (status = COLUMN_NOT_FOUND);
	}

	/* Pedigree and descrip are optional columns. */
	c_tbcfnd1 (tabinfo->tp, "PEDIGREE", &tabinfo->cp_pedigree);
	c_tbcfnd1 (tabinfo->tp, "DESCRIP", &tabinfo->cp_descrip);

	return (status);
}

/* This routine reads the relevant data from one row.  The detector name,
   global and local linearity limits, and tau are gotten.
*/

static int ReadLinTab (TblInfo *tabinfo, int row, TblRow *tabrow) {

	extern int status;

	c_tbegtt (tabinfo->tp, tabinfo->cp_det, row,
			tabrow->detector, ACS_CBUF-1);
	if (c_iraferr())
	    return (status = TABLE_ERROR);

	c_tbegtd (tabinfo->tp, tabinfo->cp_global, row, &tabrow->global_limit);
	if (c_iraferr())
	    return (status = TABLE_ERROR);

	c_tbegtd (tabinfo->tp, tabinfo->cp_local, row, &tabrow->local_limit);
	if (c_iraferr())
	    return (status = TABLE_ERROR);

	c_tbegtd (tabinfo->tp, tabinfo->cp_tau, row, &tabrow->tau);
	if (c_iraferr())
	    return (status = TABLE_ERROR);

	c_tbegtr (tabinfo->tp, tabinfo->cp_expand, row, &tabrow->expand);
	if (c_iraferr())
	    return (status = TABLE_ERROR);

	return (status);
}

/* This routine closes the mlintab table. */

static int CloseLinTab (TblInfo *tabinfo) {

	extern int status;

	c_tbtclo (tabinfo->tp);
	if (c_iraferr())
	    return (status = TABLE_ERROR);

	return (status);
}
