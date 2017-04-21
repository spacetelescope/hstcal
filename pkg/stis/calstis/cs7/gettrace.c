# include <stdio.h>
# include <stdlib.h>
# include <string.h>

# include "c_iraf.h"
# include "hstio.h"
# include "xtables.h"

# include "stis.h"
# include "calstis7.h"
# include "err.h"
# include "stisdef.h"

typedef struct {
	IRAFPointer tp;			/* pointer to table descriptor */
	IRAFPointer cp_opt_elem;
	IRAFPointer cp_sporder;
	IRAFPointer cp_cenwave;

	IRAFPointer cp_a1center;
	IRAFPointer cp_a2center;
	IRAFPointer cp_nelem;
	IRAFPointer cp_a2displ;

	IRAFPointer cp_pedigree;
	IRAFPointer cp_descrip;
	IRAFPointer cp_mjd;
	IRAFPointer cp_degperyr;
	int nrows;			/* number of rows in table */
} TblInfo;

typedef struct {
	char opt_elem[STIS_CBUF+1];	/* optical element name */
	int sporder;			/* spectral order */
	int cenwave;			/* central wavelength */
} TblRow;

static int OpenTraceTab (char *, TblInfo *);
static int ReadTraceTab (TblInfo *, int, TblRow *);
static int ReadTraceArray (TblInfo *, int, StisInfo7*,  SpTrace **);
static int CloseTraceTab (TblInfo *);

/* This routine reads the coordinate information from the spectrum trace
   table SPTRCTAB.  This is only used for obstype=SPECTROSCOPIC.

   The spectrum trace table should contain the following:
	header parameters:
		none needed
	columns:
		OPT_ELEM:  grating (or mirror) name (string)
		CENWAVE:  central wavelength (int)
		SPORDER:  order number (int)
		A1CENTER:  X location on detector corresponding to CRPIX1
		A2CENTER:  Y location on detector corresponding to CRPIX2
		NELEM:  size (max 1024) of A2DISPL array (int)
		A2DISPL:  spectrum trace (array of double)
		MJD:  MJD (double)
                DEGPERYR: rate of trace rotation (double)


   The table is read to find all rows for which the values of OPT_ELEM,
   SPORDER, and CENWAVE are the same as in the input image header.  There
   can be several such rows, each with a different value of A2CENTER.
   All these rows are read into memory, pointed to by SpTrace.
   It is an error to have duplicate values of A2CENTER in the SPTRCTAB table.
   The SPTRCTAB table need not be sorted.

   If no matching row is found in the table, sts->x2dcorr_o will be set
   to OMIT, and a warning will be printed.  If a matching row is found,
   but pedigree for that row begins with "DUMMY", sts->x2dcorr_o will be
   set to DUMMY, and a warning will be printed.  Note that status = 0 for
   both these cases, even though the data will not have been gotten.

   Note:
	Memory is allocated for the SpTrace list; it should be freed by
	calling FreeTrace.

   Phil Hodge, 2000 Jan 13:
	In ReadTraceArray, newrec is now saved by NewTrace, and freed only
	by FreeTrace.
	Add one to opt_elem buffer size.

   Phil Hodge, 2006 Sept 21:
	In ReadTraceArray, if sts->wx2dcorr is COMPLETE, set the trace
	(a2displ) to zero and set the rotation (sts->trace_rotation) to zero.
*/

int GetTrace (StisInfo7 *sts, int sporder, SpTrace **trace) {

/* arguments:
StisInfo7 *sts    i: calibration switches and info
int sporder       i: spectral order number
SpTrace **trace   o: size and coordinate info for output
*/

	int status;

	TblInfo tabinfo;	/* pointer to table descriptor, etc */
	TblRow tabrow;		/* values read from a table row */

	int row;		/* loop index */
	int CheckTrace (SpTrace **);
	void FreeTrace (SpTrace **);

	/* Open the spectrum trace table. */
	if ((status = OpenTraceTab (sts->sptrctab.name, &tabinfo)))
	    return (status);

	for (row = 1;  row <= tabinfo.nrows;  row++) {

	    if ((status = ReadTraceTab (&tabinfo, row, &tabrow)))
		return (status);

	    /* Check for a match with opt_elem, cenwave, and sporder. */

	    if (SameString (tabrow.opt_elem, sts->opt_elem) &&
		SameInt (tabrow.cenwave, sts->cenwave) &&
		SameInt (tabrow.sporder, sporder)) {

		/* Get pedigree & descrip from the row. */
		if ((status = RowPedigree (&sts->sptrctab, row,
                        tabinfo.tp, tabinfo.cp_pedigree, tabinfo.cp_descrip)))
		    return (status);
		if (sts->sptrctab.goodPedigree == DUMMY_PEDIGREE) {
		    printf ("Warning  DUMMY pedigree in row %d of %s.\n",
			row, sts->sptrctab.name);
		    sts->x2dcorr_o = DUMMY;
		    CloseTraceTab (&tabinfo);
		    return (0);
		}

		/* Read data from this row. */
		if ((status = ReadTraceArray (&tabinfo, row, sts, trace)))
		    return (status);
	    }
	}



	/* Get for duplicate a2center or non-duplicate a1center. */
	if ((status = CheckTrace (trace))) {
	    FreeTrace (trace);
	    if (status < 0) {
		printf ("Warning  Matching row not found in SPTRCTAB %s; \\\n",
				sts->sptrctab.name);
		printf ("Warning  OPT_ELEM %s, SPORDER %d, CENWAVE %d.\n",
			sts->opt_elem, sporder, sts->cenwave);
		sts->x2dcorr_o = OMIT;
	    } else {
		return (status);
	    }
	}

	if ((status = CloseTraceTab (&tabinfo)))
	    return (status);

	return (0);
}

/* This routine opens the spectrum trace table, finds the columns
   that we need, and gets the total number of rows in the table.
*/

static int OpenTraceTab (char *tname, TblInfo *tabinfo) {

	tabinfo->tp = c_tbtopn (tname, IRAF_READ_ONLY, 0);
	if (c_iraferr()) {
	    printf ("ERROR    SPTRCTAB `%s' not found.\n", tname);
	    return (OPEN_FAILED);
	}

	tabinfo->nrows = c_tbpsta (tabinfo->tp, TBL_NROWS);

	/* Find the columns. */

	c_tbcfnd1 (tabinfo->tp, "OPT_ELEM", &tabinfo->cp_opt_elem);
	c_tbcfnd1 (tabinfo->tp, "CENWAVE", &tabinfo->cp_cenwave);
	c_tbcfnd1 (tabinfo->tp, "SPORDER", &tabinfo->cp_sporder);

	c_tbcfnd1 (tabinfo->tp, "A1CENTER", &tabinfo->cp_a1center);
	c_tbcfnd1 (tabinfo->tp, "A2CENTER", &tabinfo->cp_a2center);
	c_tbcfnd1 (tabinfo->tp, "NELEM", &tabinfo->cp_nelem);
	c_tbcfnd1 (tabinfo->tp, "A2DISPL", &tabinfo->cp_a2displ);

	if (tabinfo->cp_opt_elem == 0 || tabinfo->cp_cenwave == 0 ||
	    tabinfo->cp_sporder == 0  ||
	    tabinfo->cp_a1center == 0  || tabinfo->cp_a2center == 0 ||
	    tabinfo->cp_nelem == 0   || tabinfo->cp_a2displ == 0) {

	    c_tbtclo (tabinfo->tp);
	    printf ("ERROR    Column not found in SPTRCTAB.\n");
	    return (COLUMN_NOT_FOUND);
	}

	/* Pedigree and descrip are optional columns. */
	c_tbcfnd1 (tabinfo->tp, "PEDIGREE", &tabinfo->cp_pedigree);
	c_tbcfnd1 (tabinfo->tp, "DESCRIP", &tabinfo->cp_descrip);

	/* mjd an degperyr are optional columns for trace rotation. */
        c_tbcfnd1 (tabinfo->tp, "MJD", &tabinfo->cp_mjd);
	c_tbcfnd1 (tabinfo->tp, "DEGPERYR", &tabinfo->cp_degperyr);

	return (0);
}

/* This routine reads the columns (OPT_ELEM, SPORDER, and CENWAVE) used to
   select the correct row.
*/

static int ReadTraceTab (TblInfo *tabinfo, int row, TblRow *tabrow) {

	c_tbegtt (tabinfo->tp, tabinfo->cp_opt_elem, row,
			tabrow->opt_elem, STIS_CBUF);
	if (c_iraferr())
	    return (TABLE_ERROR);

	c_tbegti (tabinfo->tp, tabinfo->cp_sporder, row, &tabrow->sporder);
	if (c_iraferr())
	    return (TABLE_ERROR);

	c_tbegti (tabinfo->tp, tabinfo->cp_cenwave, row, &tabrow->cenwave);
	if (c_iraferr())
	    return (TABLE_ERROR);

	return (0);
}

/* This routine reads the data from one row into the SpTrace structure.
   Several scalar column values and one array are gotten.

   Note that newrec is allocated but not freed here; NewTrace adds
   this to the list of traces.
*/

static int ReadTraceArray (TblInfo *tabinfo, int row, StisInfo7 *sts,
		SpTrace **trace) {

	int status;

	int nelem;		/* number of elements read from table */
	double mjd;
	double degperyr;

	SpTrace *newrec;	/* new trace record */
	int NewTrace (SpTrace **, SpTrace *);

	if ((newrec = malloc (sizeof (SpTrace))) == NULL) {
	    printf ("ERROR    Can't allocate memory in GetTrace.\n");
	    return (OUT_OF_MEMORY);
	}
	newrec->next = NULL;

	/* Get spectrum trace and other info. */
	c_tbegtd (tabinfo->tp, tabinfo->cp_a1center, row, &newrec->a1center);
	c_tbegtd (tabinfo->tp, tabinfo->cp_a2center, row, &newrec->a2center);

	c_tbegti (tabinfo->tp, tabinfo->cp_nelem, row, &newrec->nelem);
	if (newrec->nelem > MAX_SP_TRACE) {
	    printf ("ERROR    Spectrum trace in SPTRCTAB is too large.\n");
	    return (TABLE_ERROR);
	}
	if (sts->wx2dcorr == COMPLETE) {
	    /* set the trace to zero (appropriate for the output of wx2d) */
	    memset (newrec->a2displ, 0, (newrec->nelem)*sizeof(double));
	    nelem = newrec->nelem;
	} else {
	    /* read the trace from the table */
	    nelem = c_tbagtd (tabinfo->tp, tabinfo->cp_a2displ, row,
				newrec->a2displ, 1, newrec->nelem);
	    if (c_iraferr())
		return (TABLE_ERROR);
	}

	if (tabinfo->cp_mjd == 0 || sts->wx2dcorr == COMPLETE) {
	    sts->trace_rotation = 0.;
	} else {
	    c_tbegtd (tabinfo->tp, tabinfo->cp_mjd, row, &mjd);
	    c_tbegtd (tabinfo->tp, tabinfo->cp_degperyr, row, &degperyr);
	    sts->trace_rotation = rotatetrace(sts->expstart, mjd, degperyr,
			newrec->a2displ, nelem);
	}

	/* Convert a1center and a2center to zero-indexed. */
	newrec->a1center--;
	newrec->a2center--;

	if (nelem < newrec->nelem) {
	    c_tbtclo (tabinfo->tp);
	    printf ("ERROR    Not all elements were read from SPTRCTAB.\n");
	    return (TABLE_ERROR);
	}

	/* Insert newrec into the SpTrace list. */
	if ((status = NewTrace (trace, newrec)))
	    return (status);

	return (0);
}

/* This routine closes the SPTRCTAB table. */

static int CloseTraceTab (TblInfo *tabinfo) {

	c_tbtclo (tabinfo->tp);
	if (c_iraferr())
	    return (TABLE_ERROR);

	return (0);
}
