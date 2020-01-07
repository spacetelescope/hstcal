# include <stdio.h>
# include <stdlib.h>

# include "xtables.h"
# include "hstio.h"

# include "stis.h"
# include "stisdef.h"
# include "calstis6.h"
# include "hstcalerr.h"

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
	char opt_elem[STIS_CBUF];	/* optical element name */
	int sporder;			/* spectral order */
	int cenwave;			/* central wavelength */
} TblRow;

static int OpenTraceTab (char *, TblInfo *);
static int ReadTraceTab (TblInfo *, int, TblRow *);
static int ReadTraceArray (TblInfo *, int, StisInfo6*, SpTrace **);
static int CloseTraceTab (TblInfo *);


/* This routine reads the coordinate information from the spectrum trace
   table SPTRCTAB.

   The spectrum trace table should contain the following:
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

   Note:
	Memory is allocated for the SpTrace list; it should be freed by
	calling FreeTrace.




   Revision history:
   ----------------
   20 Feb 97  -  Borrowed from calstis7 (I.Busko)
   20 Feb 97  -  Empty trace table is not an error (IB)
   24 Feb 97  -  Rename routine to avoid conflict with cs7 (IB)
   10 Apr 97  -  Changes after code review (IB):
                 - explicit cast for malloc-returned pointer.
                 - rename "new" variable to avoid conflict in C++
   02 May 97  -  Set x1d_o flag for rows with DUMMY pedigree (IB)
   08 May 97  -  Conform to new _trl standard (IB)
   24 Jun 97  -  Modified logic that skips entries in trace table (IB)
   10 Apr 06  - Implemented trace rotation (ND)
*/

int GetTrace6 (StisInfo6 *sts, int sporder, SpTrace **trace) {

/* arguments:
StisInfo6 *sts    i: calibration switches and info
int sporder       i: spectral order number
SpTrace **trace   o: size and coordinate info for output
*/

	int status;

	TblInfo tabinfo;	/* pointer to table descriptor, etc */
	TblRow tabrow;		/* values read from a table row */

	int i;
	int found_trace;
	int row;		/* loop index */
	int CheckTrace6 (SpTrace **);
	void FreeTrace6 (SpTrace **);

	/* Open the spectrum trace table. */
	if ((status = OpenTraceTab (sts->sptrctab.name, &tabinfo)))
	    return (status);

	found_trace = 0;
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
		    sts->x1d_o = DUMMY;
		    CloseTraceTab (&tabinfo);
		    return (0);
		}

		/* Read data from this row. */
		if ((status = ReadTraceArray (&tabinfo, row, sts, trace)))
		    return (status);
	        found_trace = 1;
	    }
	}


	/* If nothing found, return silently. */
	if (!found_trace) {
	    if ((status = CloseTraceTab (&tabinfo)))
	        return (status);
	    return (NO_TRACE);
	}

	/* Get for duplicate a2center or non-duplicate a1center.
           An empty trace table is not an error.
        */
	i = CheckTrace6 (trace);
	if (i == NO_TRACE) {
	    if ((status = CloseTraceTab (&tabinfo)))
	        return (status);
	    return (NO_TRACE);
	} else if (i == ERROR_TRACE) {
	    FreeTrace6 (trace);
	    return (ERROR_TRACE);
	}

	if ((status = CloseTraceTab (&tabinfo)))
	    return (status);

	return (0);
}



/* This routine opens the spectrum trace table, finds the columns
   that we need, and gets the total number of rows in the table.
*/

static int OpenTraceTab (char *tname, TblInfo *tabinfo) {

        clear_cvoserr();
	tabinfo->tp = c_tbtopn (tname, IRAF_READ_ONLY, 0);
	if (c_iraferr()) {
	    printf ("ERROR    SPTRCTAB `%s' not found\n", tname);
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
	    printf ("ERROR    Column not found in SPTRCTAB\n");
	    return (COLUMN_NOT_FOUND);
	}

	/* Pedigree and descrip are optional columns. */
	c_tbcfnd1 (tabinfo->tp, "PEDIGREE", &tabinfo->cp_pedigree);
	c_tbcfnd1 (tabinfo->tp, "DESCRIP", &tabinfo->cp_descrip);

	/* mjd and degperyr are optional columns used for trace rotation. */
	c_tbcfnd1 (tabinfo->tp, "MJD", &tabinfo->cp_mjd);
	c_tbcfnd1 (tabinfo->tp, "DEGPERYR", &tabinfo->cp_degperyr);

	return (0);
}



/* This routine reads the columns (OPT_ELEM, SPORDER, and CENWAVE) used to
   select the correct row.
*/

static int ReadTraceTab (TblInfo *tabinfo, int row, TblRow *tabrow) {

	c_tbegtt (tabinfo->tp, tabinfo->cp_opt_elem, row,
			tabrow->opt_elem, STIS_CBUF-1);
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
*/

static int ReadTraceArray (TblInfo *tabinfo, int row, StisInfo6* sts, SpTrace **trace) {

	int status;

	int nelem;		/* number of elements read from table */
	double mjd;             /* MJD */
	double degperyr;        /* rate of trace rotation */

	SpTrace *newd;
	int NewTrace6 (SpTrace **, SpTrace *);

	if ((newd = (SpTrace *) malloc (sizeof (SpTrace))) == NULL) {
	    printf ("ERROR    Can't allocate memory.\n");
	    return (OUT_OF_MEMORY);
	}
	newd->next = NULL;

	/* Get spectrum trace and other info. */
	c_tbegtd (tabinfo->tp, tabinfo->cp_a1center, row, &newd->a1center);
	c_tbegtd (tabinfo->tp, tabinfo->cp_a2center, row, &newd->a2center);

	c_tbegti (tabinfo->tp, tabinfo->cp_nelem, row, &newd->nelem);
	if (newd->nelem > MAX_SP_TRACE) {
	    printf ("ERROR    Spectrum trace in SPTRCTAB is too large.\n");
	    return (TABLE_ERROR);
	}
	nelem = c_tbagtd (tabinfo->tp, tabinfo->cp_a2displ, row,
			newd->a2displ, 1, newd->nelem);
	if (c_iraferr())
	    return (TABLE_ERROR);

	if (tabinfo->cp_mjd != 0) {
	    c_tbegtd (tabinfo->tp, tabinfo->cp_mjd, row, &mjd);
	    c_tbegtd (tabinfo->tp, tabinfo->cp_degperyr, row, &degperyr);
	    sts->trace_rotation = rotatetrace(sts->expstart, mjd, degperyr, newd->a2displ, nelem);
	}

	/* Convert a1center and a2center to zero-indexed. */
	newd->a1center--;
	newd->a2center--;

	if (nelem < newd->nelem) {
	    c_tbtclo (tabinfo->tp);
	    printf ("ERROR    Not all elements were read from SPTRCTAB\n");
	    return (TABLE_ERROR);
	}

	/* Insert newd into the SpTrace list. */
	if ((status = NewTrace6 (trace, newd)))
	    return (status);

	free (newd);

	return (0);
}



/* This routine closes the SPTRCTAB table. */

static int CloseTraceTab (TblInfo *tabinfo) {

	c_tbtclo (tabinfo->tp);
	if (c_iraferr())
	    return (TABLE_ERROR);

	return (0);
}
