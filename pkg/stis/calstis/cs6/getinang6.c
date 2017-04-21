# include <stdio.h>
# include <stdlib.h>

# include "xtables.h"
# include "hstio.h"

# include "stis.h"
# include "stisdef.h"
# include "calstis6.h"
# include "err.h"

typedef struct {
	IRAFPointer tp;			/* pointer to table descriptor */
	IRAFPointer cp_opt_elem;
	IRAFPointer cp_sporder;
	IRAFPointer cp_cenwave;

	IRAFPointer cp_ncoeff1;
	IRAFPointer cp_coeff1;
	IRAFPointer cp_ncoeff2;
	IRAFPointer cp_coeff2;

	IRAFPointer cp_pedigree;
	IRAFPointer cp_descrip;
	int nrows;			/* number of rows in table */
} TblInfo;

typedef struct {
	char opt_elem[STIS_CBUF];	/* optical element name */
	int sporder;			/* spectral order */
	int cenwave;			/* central wavelength */
} TblRow;

static int OpenIACTab (char *, TblInfo *);
static int ReadIACTab (TblInfo *, int, TblRow *);
static int ReadIACArray (TblInfo *, int, InangInfo *);
static int CloseIACTab (TblInfo *);

/* This routine reads the coordinate information from the incidence-
   angle coefficients table INANGTAB.

   The incidence-angle coefficients table should contain the following:
	header parameters:
		none needed
	columns:
		OPT_ELEM:  grating (or mirror) name (string)
		CENWAVE:  central wavelength (int)
		SPORDER:  order number (int)
		NCOEFF1:  size of COEFF1 arrays (int)
		NCOEFF2:  size of COEFF2 arrays (int)
		COEFF1:  array of incidence-angle coeff (double)
		COEFF2:  array of incidence-angle coeff (double)

   The table is read to find the row for which the values of OPT_ELEM,
   SPORDER, and CENWAVE are the same as in the input image header.
   This row is read into memory.

   Notes:
   Memory is allocated for the arrays of coefficients.  This memory
   should be freed by calling FreeInang.




   Revision history:
   ----------------
   20 Feb 97  -  Borrowed from calstis7 (I.Busko)
   24 Feb 97  -  Rename routine to avoid conflict with cs7 (IB)
   10 Apr 97  -  Changes after code review (IB):
                 - replaced literal by ROW_NOT_FOUND constant
                 - explicit cast for calloc and malloc-returned pointer.
   02 May 97  -  Set x1d_o flag for rows with DUMMY pedigree (IB)
   08 May 97  -  Conform to new _trl standard (IB)

*/

int GetInang6 (StisInfo6 *sts, RefTab *table, int sporder, InangInfo *iac) {

/* arguments:
StisInfo6 *sts     i: calibration switches and info
RefTab *table      i: inangtab or mofftab
int sporder        i: spectral order number
InangInfo *iac     o: incidence-angle info
*/

	int status;

	TblInfo tabinfo;	/* pointer to table descriptor, etc */
	TblRow tabrow;		/* values read from a table row */

	int row;		/* loop index */
	int foundit = 0;

	/* Open the incidence-angle coefficients table. */
	if ((status = OpenIACTab (table->name, &tabinfo)))
	    return (status);

	for (row = 1;  row <= tabinfo.nrows;  row++) {

	    if ((status = ReadIACTab (&tabinfo, row, &tabrow)))
		return (status);

	    /* Check for a match with opt_elem, cenwave, and sporder. */

	    if (SameString (tabrow.opt_elem, sts->opt_elem) &&
		SameInt (tabrow.cenwave, sts->cenwave) &&
		SameInt (tabrow.sporder, sporder)) {

		foundit = 1;

		/* Get pedigree & descrip from the row. */
		if ((status = RowPedigree (table, row,
                        tabinfo.tp, tabinfo.cp_pedigree, tabinfo.cp_descrip)))
		    return (status);
		if (table->goodPedigree == DUMMY_PEDIGREE) {
		    sts->x1d_o = DUMMY;
		    CloseIACTab (&tabinfo);
		    return (0);
		}

		/* Read data from this row. */
		if ((status = ReadIACArray (&tabinfo, row, iac)))
		    return (status);
	    }
	}

	if (!foundit) {
	    printf ("ERROR    Matching row not found in %s\n", table->name);
	    printf ("ERROR    OPT_ELEM %s, SPORDER %d, CENWAVE %d\n",
		    sts->opt_elem, sporder, sts->cenwave);
	    return (ROW_NOT_FOUND);
	}

	if ((status = CloseIACTab (&tabinfo)))
	    return (status);

	return (0);
}

/* This routine opens the incidence-angle coefficients table, finds the
   columns that we need, and gets the total number of rows in the table.
*/

static int OpenIACTab (char *tname, TblInfo *tabinfo) {

	tabinfo->tp = c_tbtopn (tname, IRAF_READ_ONLY, 0);
	if (c_iraferr()) {
	    printf ("ERROR    Can't open `%s'\n", tname);
	    return (OPEN_FAILED);
	}

	tabinfo->nrows = c_tbpsta (tabinfo->tp, TBL_NROWS);

	/* Find the columns. */

	c_tbcfnd1 (tabinfo->tp, "OPT_ELEM", &tabinfo->cp_opt_elem);
	c_tbcfnd1 (tabinfo->tp, "CENWAVE", &tabinfo->cp_cenwave);
	c_tbcfnd1 (tabinfo->tp, "SPORDER", &tabinfo->cp_sporder);

	c_tbcfnd1 (tabinfo->tp, "NCOEFF1", &tabinfo->cp_ncoeff1);
	c_tbcfnd1 (tabinfo->tp, "COEFF1", &tabinfo->cp_coeff1);
	c_tbcfnd1 (tabinfo->tp, "NCOEFF2", &tabinfo->cp_ncoeff2);
	c_tbcfnd1 (tabinfo->tp, "COEFF2", &tabinfo->cp_coeff2);

	if (tabinfo->cp_opt_elem == 0 || tabinfo->cp_cenwave == 0 ||
	    tabinfo->cp_sporder == 0  ||
	    tabinfo->cp_ncoeff1 == 0   || tabinfo->cp_coeff1 == 0 ||
	    tabinfo->cp_ncoeff2 == 0   || tabinfo->cp_coeff2 == 0) {

	    c_tbtclo (tabinfo->tp);
	    printf ("ERROR    Column not found in %s\n", tname);
	    return (COLUMN_NOT_FOUND);
	}

	/* Pedigree and descrip are optional columns. */
	c_tbcfnd1 (tabinfo->tp, "PEDIGREE", &tabinfo->cp_pedigree);
	c_tbcfnd1 (tabinfo->tp, "DESCRIP", &tabinfo->cp_descrip);

	return (0);
}

/* This routine reads the columns (OPT_ELEM, SPORDER, and CENWAVE) used to
   select the correct row.
*/

static int ReadIACTab (TblInfo *tabinfo, int row, TblRow *tabrow) {

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

/* This routine reads the data from one row into the iac structure.
   ncoeff1 and ncoeff2 are gotten, and if they are greater than zero,
   memory is allocated and the arrays of coefficients are gotten.
*/

static int ReadIACArray (TblInfo *tabinfo, int row, InangInfo *iac) {

	int ncoeff1, ncoeff2;	/* number of coefficients read from table */
	char *tname;		/* for possible error message */

	c_tbegti (tabinfo->tp, tabinfo->cp_ncoeff1, row, &ncoeff1);
	c_tbegti (tabinfo->tp, tabinfo->cp_ncoeff2, row, &ncoeff2);
	iac->ncoeff1 = ncoeff1;
	iac->ncoeff2 = ncoeff2;

	if (ncoeff1 > 0) {
	    if ((iac->coeff1 = (double *) malloc (iac->ncoeff1 *
                                                  sizeof(double))) == NULL)
		return (OUT_OF_MEMORY);
	    /* replace ncoeff1 with actual number of elements read */
	    ncoeff1 = c_tbagtd (tabinfo->tp, tabinfo->cp_coeff1, row,
			iac->coeff1, 1, iac->ncoeff1);
	    if (c_iraferr())
		return (TABLE_ERROR);
	}
	if (ncoeff2 > 0) {
	    if ((iac->coeff2 = (double *) malloc (iac->ncoeff2 *
                                                  sizeof(double))) == NULL)
		return (OUT_OF_MEMORY);
	    ncoeff2 = c_tbagtd (tabinfo->tp, tabinfo->cp_coeff2, row,
			iac->coeff2, 1, iac->ncoeff2);
	    if (c_iraferr())
		return (TABLE_ERROR);
	}
	iac->allocated = 1;

	if (ncoeff1 < iac->ncoeff1 || ncoeff2 < iac->ncoeff2) {
	    if ((tname = (char *) calloc (STIS_LINE+1, sizeof (char))) == NULL)
		return (OUT_OF_MEMORY);
	    c_tbtnam (tabinfo->tp, tname, STIS_LINE);
	    c_tbtclo (tabinfo->tp);
	    printf (
	"ERROR    Not all coefficients were read from %s\n", tname);
	    free (tname);
	    return (TABLE_ERROR);
	}

	return (0);
}

/* This routine closes the INANGTAB table. */

static int CloseIACTab (TblInfo *tabinfo) {

	c_tbtclo (tabinfo->tp);
	if (c_iraferr())
	    return (TABLE_ERROR);

	return (0);
}
