# include <stdio.h>

# include "c_iraf.h"
# include "hstio.h"
# include "xtables.h"
# include "stis.h"
# include "calstis1.h"
# include "err.h"
# include "stisdef.h"

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
	char detector[STIS_CBUF+1];	/* detector name */
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

   Phil Hodge, 2000 Jan 13:
	Add one to detector buffer size.
*/

int GetLinTab (StisInfo1 *sts) {

/* arguments:
StisInfo1 *sts     io: calibration switches, etc
*/

	int status;

	TblInfo tabinfo;	/* pointer to table descriptor, etc */
	TblRow tabrow;		/* values read from a table row */

	int row;			/* loop index */
	int foundit;			/* detector name found in table? */

	/* Open the MAMA linearity table and find columns. */
	if ((status = OpenLinTab (sts->mlin.name, &tabinfo)))
	    return (status);

	foundit = 0;

	/* Check each row for a match with detector, and get the info
	   from the matching row.
	*/

	for (row = 1;  row <= tabinfo.nrows;  row++) {

	    /* Read the current row into tabrow. */
	    if ((status = ReadLinTab (&tabinfo, row, &tabrow)))
		return (status);

	    if (SameString (tabrow.detector, sts->det)) {

		foundit = 1;

		/* Get pedigree & descrip from the row. */
		if ((status = RowPedigree (&sts->mlin, row,
                        tabinfo.tp, tabinfo.cp_pedigree, tabinfo.cp_descrip)))
		    return (status);
		if (sts->mlin.goodPedigree == DUMMY_PEDIGREE) {
		    if (sts->glincorr == PERFORM)
			sts->glincorr = DUMMY;
		    if (sts->lflgcorr == PERFORM)
			sts->lflgcorr = DUMMY;
		    CloseLinTab (&tabinfo);
		    return (0);
		}

		sts->global_limit = tabrow.global_limit;
		sts->local_limit = tabrow.local_limit;
		sts->tau = tabrow.tau;
		sts->expand = tabrow.expand;
		break;
	    }
	}

	if (!foundit) {
	    printf ("ERROR    Detector `%s' not found in MLINTAB `%s'.\n",
		sts->det, sts->mlin.name);
	    CloseLinTab (&tabinfo);
	    return (TABLE_ERROR);
	}

	if ((status = CloseLinTab (&tabinfo)))  /* close the table */
	    return (status);

	return (0);
}

/* This routine opens the MAMA linearity table, finds the columns that we
   need, gets one header keyword, and gets the total number of rows in the
   table.  The columns are DETECTOR, GLOBAL_LIMIT, LOCAL_LIMIT, and TAU.
*/

static int OpenLinTab (char *tname, TblInfo *tabinfo) {

	tabinfo->tp = c_tbtopn (tname, IRAF_READ_ONLY, 0);
	if (c_iraferr()) {
	    printf ("ERROR    MLINTAB `%s' not found.\n", tname);
	    return (OPEN_FAILED);
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
	    printf ("ERROR    Column not found in MLINTAB.\n");
	    c_tbtclo (tabinfo->tp);
	    return (COLUMN_NOT_FOUND);
	}

	/* Pedigree and descrip are optional columns. */
	c_tbcfnd1 (tabinfo->tp, "PEDIGREE", &tabinfo->cp_pedigree);
	c_tbcfnd1 (tabinfo->tp, "DESCRIP", &tabinfo->cp_descrip);

	return (0);
}

/* This routine reads the relevant data from one row.  The detector name,
   global and local linearity limits, and tau are gotten.
*/

static int ReadLinTab (TblInfo *tabinfo, int row, TblRow *tabrow) {

	c_tbegtt (tabinfo->tp, tabinfo->cp_det, row,
			tabrow->detector, STIS_CBUF);
	if (c_iraferr())
	    return (TABLE_ERROR);

	c_tbegtd (tabinfo->tp, tabinfo->cp_global, row, &tabrow->global_limit);
	if (c_iraferr())
	    return (TABLE_ERROR);

	c_tbegtd (tabinfo->tp, tabinfo->cp_local, row, &tabrow->local_limit);
	if (c_iraferr())
	    return (TABLE_ERROR);

	c_tbegtd (tabinfo->tp, tabinfo->cp_tau, row, &tabrow->tau);
	if (c_iraferr())
	    return (TABLE_ERROR);

	c_tbegtr (tabinfo->tp, tabinfo->cp_expand, row, &tabrow->expand);
	if (c_iraferr())
	    return (TABLE_ERROR);

	return (0);
}

/* This routine closes the mlintab table. */

static int CloseLinTab (TblInfo *tabinfo) {

	c_tbtclo (tabinfo->tp);
	if (c_iraferr())
	    return (TABLE_ERROR);

	return (0);
}
