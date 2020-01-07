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
	IRAFPointer cp_cenwave;

	IRAFPointer cp_a2center;
	IRAFPointer cp_ncoeff;
	IRAFPointer cp_coeff;
	IRAFPointer cp_ref_aper;

	IRAFPointer cp_pedigree;
	IRAFPointer cp_descrip;
	int nrows;			/* number of rows in table */
} TblInfo;

typedef struct {
	char opt_elem[STIS_CBUF];	/* optical element name */
	int cenwave;			/* central wavelength */
} TblRow;

static int OpenDSPTab (char *, TblInfo *);
static int ReadDSPTab (TblInfo *, int, TblRow *);
static int ReadDSPArray (TblInfo *, int, DispRelation **);
static int CloseDSPTab (TblInfo *);

/* This routine reads the coordinate information from the dispersion
   coefficients table DISPTAB.  This is only used for obstype=SPECTROSCOPIC.

   The dispersion coefficients table should contain the following:
	header parameters:
		none needed
	columns:
		OPT_ELEM:  grating (or mirror) name (string)
		CENWAVE:  central wavelength (int)
		REF_APER:  name of the reference aperture (string)
		A2CENTER:  Y location on detector corresponding to CRPIX2
		NCOEFF:  size of COEFF array (int)
		COEFF:  array (max 10) of dispersion coefficients (double)

   The table is read to find all rows for which the values of OPT_ELEM and
   CENWAVE are the same as in the input image header.  There can be several
   such rows, each with a different value of A2CENTER. All these rows are
   read into memory, pointed to by disp. It is an error to have duplicate
   values of A2CENTER in the DISPTAB table. The DISPTAB table need not be
   sorted.

   Note:
	Memory is allocated for the disp list; it should be freed by
	calling FreeDisp6.



   Revision history:
   ----------------
   20 Feb 97  -  Borrowed from calstis7 (I.Busko)
   24 Feb 97  -  Rename routine to avoid conflict with cs7 (IB)
   10 Apr 97  -  Changes after code review (IB):
                 - replaced literal by ROW_NOT_FOUND constant
                 - explicit cast for malloc-returned pointer.
                 - rename "new" variable to avoid conflict in C++
   02 May 97  -  Set x1d_o flag for rows with DUMMY pedigree (IB)
   08 May 97  -  Conform to new _trl standard (IB)
   09 May 01  -  Do not use SPORDER to select in table (IB)

*/

int GetDisp6 (StisInfo6 *sts, DispRelation **disp) {

/* arguments:
StisInfo6 *sts       i: calibration switches and info
DispRelation **disp  o: size and coordinate info for output
*/

	int status;

	TblInfo tabinfo;	/* pointer to table descriptor, etc */
	TblRow tabrow;		/* values read from a table row */

	int row;		/* loop index */
	int CheckDisp6 (DispRelation **);
	void FreeDisp6 (DispRelation **);

	/* Open the dispersion coefficients table. */
	if ((status = OpenDSPTab (sts->disptab.name, &tabinfo)))
	    return (status);

	for (row = 1;  row <= tabinfo.nrows;  row++) {

	    if ((status = ReadDSPTab (&tabinfo, row, &tabrow)))
		return (status);

	    /* Check for a match with opt_elem and cenwave. */

	    if (SameString (tabrow.opt_elem, sts->opt_elem) &&
		SameInt (tabrow.cenwave, sts->cenwave)) {

		/* Get pedigree & descrip from the row. */
		if ((status = RowPedigree (&sts->disptab, row,
                        tabinfo.tp, tabinfo.cp_pedigree, tabinfo.cp_descrip)))
		    return (status);
		if (sts->disptab.goodPedigree == DUMMY_PEDIGREE) {
		    sts->x1d_o = DUMMY;
		    CloseDSPTab (&tabinfo);
		    return (0);
		}

		/* Read data from this row. */
		if ((status = ReadDSPArray (&tabinfo, row, disp)))
		    return (status);
	    }
	}

	/* Get for duplicate a2center or non-duplicate ref_aper. */
	if ((status = CheckDisp6 (disp))) {
	    FreeDisp6 (disp);
	    if (status < 0) {
		printf ("ERROR    Matching row not found in DISPTAB %s\n",
				sts->disptab.name);
		printf ("ERROR    OPT_ELEM %s, CENWAVE %d\n",
		    sts->opt_elem, sts->cenwave);
		return (ROW_NOT_FOUND);
	    } else {
		return (status);
	    }
	}

	if ((status = CloseDSPTab (&tabinfo)))
	    return (status);

	return (0);
}

/* This routine opens the dispersion coefficients table, finds the columns
   that we need, and gets the total number of rows in the table.
*/

static int OpenDSPTab (char *tname, TblInfo *tabinfo) {

	tabinfo->tp = c_tbtopn (tname, IRAF_READ_ONLY, 0);
	if (c_iraferr()) {
	    printf ("ERROR    DISPTAB `%s' not found\n", tname);
	    return (OPEN_FAILED);
	}

	tabinfo->nrows = c_tbpsta (tabinfo->tp, TBL_NROWS);

	/* Find the columns. */

	c_tbcfnd1 (tabinfo->tp, "OPT_ELEM", &tabinfo->cp_opt_elem);
	c_tbcfnd1 (tabinfo->tp, "CENWAVE", &tabinfo->cp_cenwave);

	c_tbcfnd1 (tabinfo->tp, "A2CENTER", &tabinfo->cp_a2center);
	c_tbcfnd1 (tabinfo->tp, "NCOEFF", &tabinfo->cp_ncoeff);
	c_tbcfnd1 (tabinfo->tp, "COEFF", &tabinfo->cp_coeff);
	c_tbcfnd1 (tabinfo->tp, "REF_APER", &tabinfo->cp_ref_aper);

	if (tabinfo->cp_opt_elem == 0 || tabinfo->cp_cenwave == 0 ||
	    tabinfo->cp_a2center == 0 ||
	    tabinfo->cp_ncoeff == 0   || tabinfo->cp_coeff == 0 ||
	    tabinfo->cp_ref_aper == 0) {

	    c_tbtclo (tabinfo->tp);
	    printf ("ERROR    Column not found in DISPTAB\n");
	    return (COLUMN_NOT_FOUND);
	}

	/* Pedigree and descrip are optional columns. */
	c_tbcfnd1 (tabinfo->tp, "PEDIGREE", &tabinfo->cp_pedigree);
	c_tbcfnd1 (tabinfo->tp, "DESCRIP", &tabinfo->cp_descrip);

	return (0);
}

/* This routine reads the columns (OPT_ELEM and CENWAVE) used to
   select the correct row.
*/

static int ReadDSPTab (TblInfo *tabinfo, int row, TblRow *tabrow) {

	c_tbegtt (tabinfo->tp, tabinfo->cp_opt_elem, row,
			tabrow->opt_elem, STIS_CBUF-1);
	if (c_iraferr())
	    return (TABLE_ERROR);

	c_tbegti (tabinfo->tp, tabinfo->cp_cenwave, row, &tabrow->cenwave);
	if (c_iraferr())
	    return (TABLE_ERROR);

	return (0);
}

/* This routine reads the data from one row into the disp structure.
   Several scalar column values and one array are gotten.
*/

static int ReadDSPArray (TblInfo *tabinfo, int row, DispRelation **disp) {

	int status;
	int ncoeff;		/* number of coefficients read from table */
	DispRelation *newd;
	int NewDisp6 (DispRelation **, DispRelation *);

	if ((newd = (DispRelation *) malloc (sizeof (DispRelation))) == NULL) {
	    printf ("ERROR    Can't allocate memory.\n");
	    return (OUT_OF_MEMORY);
	}
	newd->next = NULL;

	/* Get dispersion coefficients and other info. */
	c_tbegtd (tabinfo->tp, tabinfo->cp_a2center, row, &newd->a2center);
	c_tbegtt (tabinfo->tp, tabinfo->cp_ref_aper, row, newd->ref_aper,
		STIS_CBUF-1);
	c_tbegti (tabinfo->tp, tabinfo->cp_ncoeff, row, &newd->ncoeff);
	if (newd->ncoeff > MAX_DISP_COEFF) {
	    printf (
                "ERROR    Too many dispersion coefficients %d in DISPTAB.\n",
		newd->ncoeff);
	    return (TABLE_ERROR);
	}
	ncoeff = c_tbagtd (tabinfo->tp, tabinfo->cp_coeff, row,
			newd->coeff, 1, newd->ncoeff);
	if (c_iraferr())
	    return (TABLE_ERROR);

	/* Convert a2center to zero-indexed. */
	newd->a2center--;

	if (ncoeff < newd->ncoeff) {
	    c_tbtclo (tabinfo->tp);
	    printf ("ERROR    Not all coefficients were read from DISPTAB\n");
	    return (TABLE_ERROR);
	}

	/* Insert newd into the disp list. */
	if ((status = NewDisp6 (disp, newd)))
	    return (status);

	free (newd);

	return (0);
}

/* This routine closes the DISPTAB table. */

static int CloseDSPTab (TblInfo *tabinfo) {

	c_tbtclo (tabinfo->tp);
	if (c_iraferr())
	    return (TABLE_ERROR);

	return (0);
}
