# include <stdio.h>
# include <stdlib.h>

# include "c_iraf.h"
# include "hstio.h"
# include "xtables.h"

# include "stis.h"
# include "err.h"
# include "stisdef.h"

typedef struct {
	IRAFPointer tp;			/* pointer to table descriptor */
	IRAFPointer cp_opt_elem;
	IRAFPointer cp_cenwave;

	IRAFPointer cp_mref;
	IRAFPointer cp_yref;
	IRAFPointer cp_a4corr;

	IRAFPointer cp_pedigree;
	IRAFPointer cp_descrip;
	int nrows;			/* number of rows in table */
} TblInfo;

typedef struct {
	char opt_elem[STIS_CBUF+1];	/* optical element name */
	int cenwave;			/* central wavelength */
} TblRow;

static int OpenMOCTab (char *, TblInfo *);
static int ReadMOCTab (TblInfo *, int, TblRow *);
static int ReadMOCArray (TblInfo *, int, int *, double *, double *);
static int CloseMOCTab (TblInfo *);

/* This routine reads the MAMA offset parameters from the dispersion
   coefficients table DISPTAB.  This is only used for echelle data.

   The dispersion coefficients table should contain the following:
	header parameters:
		none needed
	columns:
		OPT_ELEM:  grating (or mirror) name (string)
		CENWAVE:  central wavelength (int)
		MREF:  spectral order number of reference order (int)
		YREF:  Y location of spectral order MREF (double)
		A4CORR:  correction factor (double)

   The table is read to find the row for which the values of OPT_ELEM
   and CENWAVE are the same as in the input image header.
   This row is read into memory.

   If no matching row is found in the table, a4corr will be set
   to zero, and a warning will be printed.  If a matching row is found,
   but pedigree for that row begins with "DUMMY", a4corr will be
   set to zero, and a warning will be printed.  Note that status = 0 for
   both these cases, even though the data will not have been gotten.

   Phil Hodge, 2001 May 4:
	Function created, based on the old getmoc.
*/

int GetMOC (RefTab *disptab, char *opt_elem, int cenwave,
		int *mref, double *yref, double *a4corr) {

/* arguments:
RefTab *disptab    i: disptab
char *opt_elem     i: grating name (selection criterion)
int cenwave        i: central wavelength (selection criterion)
int *mref          o: spectral order number of reference order
double *yref       o: Y location of spectral order mref
double *a4corr     o: correction factor read from table
*/

	int status;

	TblInfo tabinfo;	/* pointer to table descriptor, etc */
	TblRow tabrow;		/* values read from a table row */

	int row;		/* loop index */
	int foundit = 0;

	/* initial values, in case of a non-fatal error */
	*mref = 1;
	*yref = 512.;
	*a4corr = 0.;

	/* Open the dispersion coefficients table. */
	if ((status = OpenMOCTab (disptab->name, &tabinfo))) {
	    if (status < 0)
		status = 0;
	    return (status);
	}

	for (row = 1;  row <= tabinfo.nrows;  row++) {

	    if ((status = ReadMOCTab (&tabinfo, row, &tabrow)) != 0) {
		if (status < 0)
		    status = 0;
		return (status);
	    }

	    /* Check for a match with opt_elem, and cenwave. */

	    if (SameString (tabrow.opt_elem, opt_elem) &&
		SameInt (tabrow.cenwave, cenwave)) {

		foundit = 1;

		/* Get pedigree & descrip from the row. */
		if ((status = RowPedigree (disptab, row,
                        tabinfo.tp, tabinfo.cp_pedigree, tabinfo.cp_descrip)))
		    return (status);
		if (disptab->goodPedigree == DUMMY_PEDIGREE) {
		    printf ("Warning  DUMMY pedigree in row %d of %s. \\\n",
			row, disptab->name);
		    printf ("Warning  MAMA offset coefficient set to zero.\n");
		    CloseMOCTab (&tabinfo);
		    return (0);
		}

		/* Read data from this row. */
		if ((status = ReadMOCArray (&tabinfo, row, mref, yref, a4corr)))
		    return (status);
	    }
	}

	if (!foundit) {
	    printf ("Warning  Matching row not found in %s \\\n",
			disptab->name);
	    printf ("Warning  OPT_ELEM %s, CENWAVE %d \\\n",
			opt_elem, cenwave);
	    printf ("Warning  MAMA offset coefficient set to zero.\n");
	}

	if ((status = CloseMOCTab (&tabinfo)))
	    return (status);

	return (0);
}

/* This routine opens the dispersion coefficients table, finds the
   columns that we need, and gets the total number of rows in the table.
*/

static int OpenMOCTab (char *tname, TblInfo *tabinfo) {

	if (!GotFileName (tname)) {
	    printf ("Warning  DISPTAB = `%s' \\\n", tname);
	    printf ("Warning  MAMA offset coefficient set to zero.\n");
	    return (-1);
	}

	tabinfo->tp = c_tbtopn (tname, IRAF_READ_ONLY, 0);
	if (c_iraferr()) {
	    printf ("ERROR:  Can't open `%s'\n", tname);
	    return (OPEN_FAILED);
	}

	tabinfo->nrows = c_tbpsta (tabinfo->tp, TBL_NROWS);

	/* Find the columns. */

	c_tbcfnd1 (tabinfo->tp, "OPT_ELEM", &tabinfo->cp_opt_elem);
	c_tbcfnd1 (tabinfo->tp, "CENWAVE", &tabinfo->cp_cenwave);

	c_tbcfnd1 (tabinfo->tp, "MREF", &tabinfo->cp_mref);
	c_tbcfnd1 (tabinfo->tp, "YREF", &tabinfo->cp_yref);
	c_tbcfnd1 (tabinfo->tp, "A4CORR", &tabinfo->cp_a4corr);

	if (tabinfo->cp_opt_elem == 0 || tabinfo->cp_cenwave == 0) {

	    c_tbtclo (tabinfo->tp);
	    printf ("ERROR:  Column not found in %s\n", tname);
	    return (COLUMN_NOT_FOUND);
	}

	/* The original disptab did not include these columns. */
	if (tabinfo->cp_mref == 0 || tabinfo->cp_yref == 0 ||
	    tabinfo->cp_a4corr == 0) {
	    c_tbtclo (tabinfo->tp);
	    printf ("Warning  DISPTAB appears to be the old format; \\\n");
	    printf ("Warning  MAMA offset coefficient set to zero.\n");
	    return (-1);	/* not a fatal error */
	}

	/* Pedigree and descrip are optional columns. */
	c_tbcfnd1 (tabinfo->tp, "PEDIGREE", &tabinfo->cp_pedigree);
	c_tbcfnd1 (tabinfo->tp, "DESCRIP", &tabinfo->cp_descrip);

	return (0);
}

/* This routine reads the columns (OPT_ELEM and CENWAVE) used to
   select the correct row.
*/

static int ReadMOCTab (TblInfo *tabinfo, int row, TblRow *tabrow) {

	c_tbegtt (tabinfo->tp, tabinfo->cp_opt_elem, row,
			tabrow->opt_elem, STIS_CBUF);
	if (c_iraferr())
	    return (TABLE_ERROR);

	c_tbegti (tabinfo->tp, tabinfo->cp_cenwave, row, &tabrow->cenwave);
	if (c_iraferr())
	    return (TABLE_ERROR);

	return (0);
}

/* This routine reads mref, yref, and a4corr from the current row. */

static int ReadMOCArray (TblInfo *tabinfo, int row,
		int *mref, double *yref, double *a4corr) {

	c_tbegti (tabinfo->tp, tabinfo->cp_mref, row, mref);
	if (c_iraferr())
	    return (TABLE_ERROR);

	c_tbegtd (tabinfo->tp, tabinfo->cp_yref, row, yref);
	if (c_iraferr())
	    return (TABLE_ERROR);

	/* convert to zero-indexing */
	yref--;

	c_tbegtd (tabinfo->tp, tabinfo->cp_a4corr, row, a4corr);
	if (c_iraferr())
	    return (TABLE_ERROR);

	return (0);
}

/* This routine closes the DISPTAB table. */

static int CloseMOCTab (TblInfo *tabinfo) {

	c_tbtclo (tabinfo->tp);
	if (c_iraferr())
	    return (TABLE_ERROR);

	return (0);
}
