# include <stdio.h>
# include <stdlib.h>
# include <math.h>		/* for fabs */

# include "c_iraf.h"
# include "hstio.h"
# include "xtables.h"

# include "stis.h"
# include "calstis4.h"
# include "err.h"
# include "stisdef.h"

# define MIDDLE_LINE (512.)	/* middle line of detector, in Y direction */

typedef struct {
	IRAFPointer tp;			/* pointer to table descriptor */
	IRAFPointer cp_opt_elem;
	IRAFPointer cp_cenwave;
	IRAFPointer cp_a2center;

	IRAFPointer cp_ncoeff;
	IRAFPointer cp_coeff;
	IRAFPointer cp_ref_aper;

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
	double a2center;		/* y location */
} TblRow;

static int OpenDSPTab (char *, TblInfo *);
static int ReadDSPTab (TblInfo *, int, TblRow *);
static int ReadDSPArray (TblInfo *, int, DispRelation *, char *);
static int CloseDSPTab (TblInfo *);

/* This routine reads the coordinate information from the dispersion
   coefficients table DISPTAB.  This is only used for obstype=SPECTROSCOPIC,
   and currently it is used only for echelle and prism data.

   The dispersion coefficients table should contain the following:
	header parameters:
		none needed
	columns:
		OPT_ELEM:  grating name (string)
		CENWAVE:  central wavelength (int)
		REF_APER:  name of the reference aperture (string)
		A2CENTER:  Y location on detector
		NCOEFF:  size of COEFF array (int)
		COEFF:  array (max 10) of dispersion coefficients (double)
		MREF:  spectral order number of reference order (int)
		YREF:  Y location of spectral order MREF (double)
		A4CORR:  correction factor (double)

   The table is read to find the row for which the values of OPT_ELEM
   and CENWAVE are the same as in the input image header.  If more than
   one row matches, the one with A2CENTER closest to 512 is used.  The
   spectral order number is not expected to be needed, since it is a
   part of the dispersion relation.

   It is a fatal error if no matching row is found in the table.

   Phil Hodge, 2000 Jan 13:
	File created, based on cs7/getdisp.c.

   Phil Hodge, 2001 Feb 23:
	Select the row with a2center closest to 512.

   Phil Hodge, 2011 Jan 5:
	Copy code from ../lib/getmoc.c to here, to populate mref, yref and
	a4corr in the DispRelation struct.
*/

int GetDisp4 (StisInfo4 *sts, DispRelation *disp, char *ref_aper) {

/* arguments:
StisInfo4 *sts       i: calibration switches and info
DispRelation *disp   o: dispersion relation
char *ref_aper       o: name of aperture used to measure dispersion relation
*/

	int status;

	TblInfo tabinfo;	/* pointer to table descriptor, etc */
	TblRow tabrow;		/* values read from a table row */

	int row;		/* loop index */
	int best_row;		/* read info from this row */
	double a2diff;		/* difference between A2CENTER and middle */
	double min_a2diff;	/* min value of a2diff */

	/* Open the dispersion coefficients table. */
	if ((status = OpenDSPTab (sts->disptab.name, &tabinfo)))
	    return (status);

	best_row = -1;			/* initial values */
	min_a2diff = 2. * MIDDLE_LINE;

	for (row = 1;  row <= tabinfo.nrows;  row++) {

	    if ((status = ReadDSPTab (&tabinfo, row, &tabrow)))
		return (status);

	    /* Check for a match with opt_elem and cenwave. */

	    if (SameString (tabrow.opt_elem, sts->opt_elem) &&
		SameInt (tabrow.cenwave, sts->cenwave)) {

		a2diff = fabs (tabrow.a2center - MIDDLE_LINE);
		if (a2diff < min_a2diff) {
		    min_a2diff = a2diff;
		    best_row = row;
		}
	    }
	}

	if (best_row < 1) {
	    printf ("Matching row not found in DISPTAB %s:\n",
				sts->disptab.name);
	    printf ("    OPT_ELEM %s, CENWAVE %d\n",
			sts->opt_elem, sts->cenwave);
	    return (TABLE_ERROR);
	}

	/* Read data from the most appropriate row. */
	if ((status = ReadDSPArray (&tabinfo, best_row, disp, ref_aper)))
	    return (status);

	/* Get pedigree & descrip from this row. */
	if ((status = RowPedigree (&sts->disptab, best_row,
                tabinfo.tp, tabinfo.cp_pedigree, tabinfo.cp_descrip)))
	    return (status);
	if (sts->disptab.goodPedigree == DUMMY_PEDIGREE) {
	    printf ("Warning  DUMMY pedigree in row %d of %s.\n",
		best_row, sts->disptab.name);
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
	    printf ("DISPTAB `%s' not found.\n", tname);
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

	/* for the a4corr for echelle data */
	c_tbcfnd1 (tabinfo->tp, "MREF", &tabinfo->cp_mref);
	c_tbcfnd1 (tabinfo->tp, "YREF", &tabinfo->cp_yref);
	c_tbcfnd1 (tabinfo->tp, "A4CORR", &tabinfo->cp_a4corr);

	if (tabinfo->cp_opt_elem == 0 || tabinfo->cp_cenwave == 0 ||
	    tabinfo->cp_a2center == 0 ||
	    tabinfo->cp_ncoeff == 0   || tabinfo->cp_coeff == 0 ||
	    tabinfo->cp_ref_aper == 0) {

	    c_tbtclo (tabinfo->tp);
	    printf ("Column not found in DISPTAB.\n");
	    return (COLUMN_NOT_FOUND);
	}

	/* Pedigree and descrip are optional columns. */
	c_tbcfnd1 (tabinfo->tp, "PEDIGREE", &tabinfo->cp_pedigree);
	c_tbcfnd1 (tabinfo->tp, "DESCRIP", &tabinfo->cp_descrip);

	if (tabinfo->cp_mref == 0 || tabinfo->cp_yref == 0 ||
	    tabinfo->cp_a4corr == 0) {
	    printf ("Warning  A4CORR not found.\n");
	    tabinfo->cp_a4corr = 0;	/* so we can test on just one value */
	}

	return (0);
}

/* This routine reads the columns (OPT_ELEM, CENWAVE, A2CENTER) used to
   select the correct row.
*/

static int ReadDSPTab (TblInfo *tabinfo, int row, TblRow *tabrow) {

	c_tbegtt (tabinfo->tp, tabinfo->cp_opt_elem, row,
			tabrow->opt_elem, STIS_CBUF);
	if (c_iraferr())
	    return (TABLE_ERROR);

	c_tbegti (tabinfo->tp, tabinfo->cp_cenwave, row, &tabrow->cenwave);
	if (c_iraferr())
	    return (TABLE_ERROR);

	c_tbegtd (tabinfo->tp, tabinfo->cp_a2center, row, &tabrow->a2center);
	if (c_iraferr())
	    return (TABLE_ERROR);

	return (0);
}

/* This routine reads the data from one row into the disp structure.
   The reference aperture name is also read.
*/

static int ReadDSPArray (TblInfo *tabinfo, int row,
		DispRelation *disp, char *ref_aper) {

	int ncoeff;		/* number of coefficients read from table */
	int i;

	for (i = 0;  i < MAX_DISP_COEFF;  i++)
	    disp->coeff[i] = 0.;

	/* Get name of aperture used to measure dispersion relation. */
	c_tbegtt (tabinfo->tp, tabinfo->cp_ref_aper, row, ref_aper, STIS_CBUF);

	/* Number of coefficients in dispersion relation. */
	c_tbegti (tabinfo->tp, tabinfo->cp_ncoeff, row, &disp->ncoeff);
	if (disp->ncoeff > MAX_DISP_COEFF) {
	    printf ("Too many dispersion coefficients %d in DISPTAB.\n",
		disp->ncoeff);
	    return (TABLE_ERROR);
	}

	/* Get the coefficients themselves. */
	ncoeff = c_tbagtd (tabinfo->tp, tabinfo->cp_coeff, row,
			disp->coeff, 1, disp->ncoeff);
	if (c_iraferr())
	    return (TABLE_ERROR);

	if (ncoeff < disp->ncoeff) {
	    c_tbtclo (tabinfo->tp);
	    printf ("Not all coefficients were read from DISPTAB.\n");
	    return (TABLE_ERROR);
	}

	/* Get the a4corr info. */
	if (tabinfo->cp_a4corr == 0) {
	    disp->mref = 0;
	    disp->yref = 0.;
	    disp->a4corr = 0.;
	} else {
	    double yref;

	    c_tbegti (tabinfo->tp, tabinfo->cp_mref, row, &disp->mref);
	    if (c_iraferr())
		return (TABLE_ERROR);

	    c_tbegtd (tabinfo->tp, tabinfo->cp_yref, row, &yref);
	    if (c_iraferr())
		return (TABLE_ERROR);
	    /* convert to zero-indexing */
	    disp->yref = yref - 1.;

	    c_tbegtd (tabinfo->tp, tabinfo->cp_a4corr, row, &disp->a4corr);
	    if (c_iraferr())
		return (TABLE_ERROR);
	}

	return (0);
}

/* This routine closes the DISPTAB table. */

static int CloseDSPTab (TblInfo *tabinfo) {

	c_tbtclo (tabinfo->tp);
	if (c_iraferr())
	    return (TABLE_ERROR);

	return (0);
}
