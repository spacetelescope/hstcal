# include <stdio.h>
# include <stdlib.h>

# include "c_iraf.h"
# include "hstio.h"
# include "xtables.h"

# include "stis.h"
# include "calstis4.h"
# include "hstcalerr.h"
# include "stisdef.h"

typedef struct {
	IRAFPointer tp;			/* pointer to table descriptor */
	IRAFPointer cp_opt_elem;
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
	char opt_elem[STIS_CBUF+1];	/* optical element name */
	int cenwave;			/* central wavelength */
} TblRow;

static int OpenIACTab (char *, TblInfo *);
static int ReadIACTab (TblInfo *, int, TblRow *);
static int ReadIACArray (TblInfo *, int, InangInfo *);
static int CloseIACTab (TblInfo *);
static void AdjustDisp (DispRelation *, double, InangInfo *);
static void FreeInang (InangInfo *);

/* This routine reads the coefficients from the incidence-angle table INANGTAB
   and modifies the dispersion relation coefficients accordingly.

   The incidence-angle coefficients table should contain the following:
	header parameters:
		none needed
	columns:
		OPT_ELEM:  grating (or mirror) name (string)
		CENWAVE:  central wavelength (int)
		NCOEFF1:  size of COEFF1 array (int)
		COEFF1:  array of incidence-angle coeff (double)

   The table is read to find the row for which the values of OPT_ELEM
   and CENWAVE are the same as in the input image header.  Note that
   we assume that the row does not depend on SPORDER.  The matching row
   is read into memory, and those coefficients are applied, using the
   incidence angle, to correct the dispersion coefficients.

   Phil Hodge, 2000 Jan 13:
	Created from cs7/getinang.c.
*/

int GetInang4 (StisInfo4 *sts, DispRelation *disp, double angle) {

/* arguments:
StisInfo4 *sts       i: calibration switches and info
DispRelation *disp  io: dispersion relation; coeffients will be modified
double angle         i: incidence angle
*/

	int status;

	TblInfo tabinfo;	/* pointer to table descriptor, etc */
	TblRow tabrow;		/* values read from a table row */
	InangInfo iac;		/* incidence-angle info */

	int row;		/* loop index */
	int foundit = 0;

	/* Open the incidence-angle coefficients table. */
	if ((status = OpenIACTab (sts->inangtab.name, &tabinfo)))
	    return (status);

	iac.allocated = 0;
	for (row = 1;  row <= tabinfo.nrows;  row++) {

	    if ((status = ReadIACTab (&tabinfo, row, &tabrow)))
		return (status);

	    /* Check for a match with opt_elem, cenwave. */

	    if (SameString (tabrow.opt_elem, sts->opt_elem) &&
		SameInt (tabrow.cenwave, sts->cenwave)) {

		foundit = 1;

		/* Get pedigree & descrip from the row. */
		if ((status = RowPedigree (&sts->inangtab, row,
                        tabinfo.tp, tabinfo.cp_pedigree, tabinfo.cp_descrip)))
		    return (status);
		if (sts->inangtab.goodPedigree == DUMMY_PEDIGREE) {
		    printf ("Warning  DUMMY pedigree in row %d of %s.\n",
			row, sts->inangtab.name);
		}

		/* Read data from this row. */
		if ((status = ReadIACArray (&tabinfo, row, &iac)))
		    return (status);

		/* Modify the dispersion coefficients. */
		AdjustDisp (disp, angle, &iac);

		break;
	    }
	}

	if (!foundit) {
	    printf ("Matching row not found in INANGTAB %s;\n",
			sts->inangtab.name);
	    printf ("  OPT_ELEM %s, CENWAVE %d.\n",
			sts->opt_elem, sts->cenwave);
	    return (TABLE_ERROR);
	}

	if ((status = CloseIACTab (&tabinfo)))
	    return (status);

	FreeInang (&iac);

	return (0);
}

/* This routine opens the incidence-angle coefficients table, finds the
   columns that we need, and gets the total number of rows in the table.
*/

static int OpenIACTab (char *tname, TblInfo *tabinfo) {

	tabinfo->tp = c_tbtopn (tname, IRAF_READ_ONLY, 0);
	if (c_iraferr()) {
	    printf ("Can't open `%s'.\n", tname);
	    return (OPEN_FAILED);
	}

	tabinfo->nrows = c_tbpsta (tabinfo->tp, TBL_NROWS);

	/* Find the columns. */

	c_tbcfnd1 (tabinfo->tp, "OPT_ELEM", &tabinfo->cp_opt_elem);
	c_tbcfnd1 (tabinfo->tp, "CENWAVE", &tabinfo->cp_cenwave);

	c_tbcfnd1 (tabinfo->tp, "NCOEFF1", &tabinfo->cp_ncoeff1);
	c_tbcfnd1 (tabinfo->tp, "COEFF1", &tabinfo->cp_coeff1);
	c_tbcfnd1 (tabinfo->tp, "NCOEFF2", &tabinfo->cp_ncoeff2);
	c_tbcfnd1 (tabinfo->tp, "COEFF2", &tabinfo->cp_coeff2);

	if (tabinfo->cp_opt_elem == 0 || tabinfo->cp_cenwave == 0 ||
	    tabinfo->cp_ncoeff1 == 0   || tabinfo->cp_coeff1 == 0 ||
	    tabinfo->cp_ncoeff2 == 0   || tabinfo->cp_coeff2 == 0) {

	    c_tbtclo (tabinfo->tp);
	    printf ("Column not found in %s.\n", tname);
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

static int ReadIACTab (TblInfo *tabinfo, int row, TblRow *tabrow) {

	c_tbegtt (tabinfo->tp, tabinfo->cp_opt_elem, row,
			tabrow->opt_elem, STIS_CBUF);
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
	    if ((iac->coeff1 = malloc (iac->ncoeff1 * sizeof(double))) == NULL)
		return (OUT_OF_MEMORY);
	    /* replace ncoeff1 with actual number of elements read */
	    ncoeff1 = c_tbagtd (tabinfo->tp, tabinfo->cp_coeff1, row,
			iac->coeff1, 1, iac->ncoeff1);
	    if (c_iraferr())
		return (TABLE_ERROR);
	}
	if (ncoeff2 > 0) {
	    if ((iac->coeff2 = malloc (iac->ncoeff2 * sizeof(double))) == NULL)
		return (OUT_OF_MEMORY);
	    ncoeff2 = c_tbagtd (tabinfo->tp, tabinfo->cp_coeff2, row,
			iac->coeff2, 1, iac->ncoeff2);
	    if (c_iraferr())
		return (TABLE_ERROR);
	}
	iac->allocated = 1;

	if (ncoeff1 < iac->ncoeff1 || ncoeff2 < iac->ncoeff2) {
	    if ((tname = calloc (STIS_FNAME+1, sizeof (char))) == NULL)
		return (OUT_OF_MEMORY);
	    c_tbtnam (tabinfo->tp, tname, STIS_FNAME);
	    c_tbtclo (tabinfo->tp);
	    printf (
	"Not all coefficients were read from %s.\n", tname);
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

/* This routine applies the incidence-angle correction to the dispersion
   coefficients.
*/

static void AdjustDisp (DispRelation *disp, double angle, InangInfo *iac) {

/* arguments:
DispRelation *disp   io: dispersion relation interpolated
double angle          i: offset of slit from reference slit
InangInfo *iac        i: incidence-angle correction coefficients
*/

	int ncoeff;
	int i;

	if (disp->ncoeff < iac->ncoeff1) {
	    ncoeff = disp->ncoeff;
	    printf (
	"Warning  %d dispersion coefficients, but %d incidence-angle coeff.\n",
			disp->ncoeff, iac->ncoeff1);
	} else {
	    ncoeff = iac->ncoeff1;
	}

	/* first coefficients */
	for (i = 0;  i < ncoeff;  i++)
	    disp->coeff[i] += iac->coeff1[i] * angle;

	/* second coefficients */
	if (iac->ncoeff2 > 0)
	    disp->coeff[0] += iac->coeff2[0] * angle;

	if (iac->ncoeff2 > 1)
	    disp->coeff[0] += iac->coeff2[1] * angle * angle;

	if (iac->ncoeff2 > 2) {
	    printf (
	"Warning  %d incidence-angle second coefficents, limit is 2;\n",
			iac->ncoeff2);
	    printf (
	"  the remaining coefficents will not be applied.\n");
	}
}

/* This routine frees memory for the incidence-angle coefficients. */

static void FreeInang (InangInfo *iac) {

	if (iac->allocated) {
	    if (iac->ncoeff1 > 0)
		free (iac->coeff1);
	    if (iac->ncoeff2 > 0)
		free (iac->coeff2);
	    iac->allocated = 0;
	}
}
