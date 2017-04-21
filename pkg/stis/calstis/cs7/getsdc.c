# include <stdio.h>
# include <stdlib.h>	/* malloc */
# include <string.h>	/* strcpy */

# include "c_iraf.h"
# include "hstio.h"
# include "xtables.h"

# include "stis.h"
# include "calstis7.h"
# include "err.h"
# include "stisdef.h"

typedef struct {
	IRAFPointer tp;			/* pointer to table descriptor */
	IRAFPointer cp_aperture;
	IRAFPointer cp_opt_elem;
	IRAFPointer cp_cenwave;

	IRAFPointer cp_sporder;
	IRAFPointer cp_a2center;
	IRAFPointer cp_npix[2];
	IRAFPointer cp_crpix[2];
	IRAFPointer cp_crval[2];
	IRAFPointer cp_cdelt[2];

	IRAFPointer cp_pedigree;
	IRAFPointer cp_descrip;
	int nrows;			/* number of rows in table */
} TblInfo;

typedef struct {
	char aperture[STIS_CBUF+1];	/* slit name */
	char opt_elem[STIS_CBUF+1];	/* optical element name */
	int cenwave;			/* central wavelength */
} TblRow;

static int OpenSDistTab (char *, TblInfo *);
static int ReadSDistTab (TblInfo *, int, TblRow *);
static int ReadSDistArray (TblInfo *, int, CoordInfo **);
static int CloseSDistTab (TblInfo *);

/* This routine reads the coordinate information from the spectroscopic
   distortion table SDCTAB.  This is only used for obstype=SPECTROSCOPIC.

   The distortion information table should contain the following:
	header parameters:
		none needed
	columns:
		APERTURE:  name of the slit (string)
		OPT_ELEM:  grating (or mirror) name (string)
		CENWAVE:  central wavelength (int)
		SPORDER:  order number (int)
		A2CENTER:  Y location on detector corresponding to CRPIX2
		NPIX1, NPIX2:  size for output image (int)
		CRPIX1, CRPIX2:  reference pixel (double)
		CRVAL1, CRVAL2:  coordinates at reference pixel (double)
		CDELT1, CDELT2:  increment per pixel (double)

   The table is read to find all rows for which the values of OPT_ELEM
   and CENWAVE are the same as in the input image header.  If an APERTURE
   column exists, it will also be used to select the row.  There should be
   one or more matching rows, corresponding to different values of spectral
   order SPORDER.  All these rows are read into memory, pointed to by
   coords.  The minimum and maximum spectral order numbers are returned.
   It is an error to have duplicate values of SPORDER in the SDCTAB table,
   and all values between minorder and maxorder must be present.  The
   SDCTAB table need not be sorted.

   Since APERTURE was added as a column and a selection criterion, there
   can be two matching rows with the same spectral order number if the
   aperture used was one of the CCD pseudo-apertures (sporder=1); the
   two matching rows would have APERTURE equal to the actual name and
   APERTURE = "ANY".  In this case we must stop looking for matching rows
   after finding the pseudo-aperture name, not only because that is the
   row that contains the correct information for the aperture in use, but
   also to prevent an error due to duplicate SPORDER values.  For this
   reason, the rows for APERTURE = "ANY" must follow the rows for specific
   apertures.  There is an explicit test for this in the loop over rows:
   if a matching row has been found and the spectral order number is one,
   we break out of the loop.

   Note:
	Memory is allocated for the coords list; it should be freed
	by calling FreeCoord.
	coords should have been initialized to NULL.

   Phil Hodge, 2000 Jan 13:
	Add one to opt_elem buffer size.

   Phil Hodge, 2000 Sept 5:
	Check for the existence of an APERTURE column, and if it exists,
	use it as a selection criterion.  For first-order data, break out
	of the loop after finding one matching row; we do this because
	if the aperture is one of the CCD pseudo-apertures (...E1), then
	that aperture name and ANY would both match, and we must take only
	the row for the pseudo-aperture.
*/

int GetSDC (StisInfo7 *sts, CoordInfo **coords, int *minorder, int *maxorder) {

/* arguments:
StisInfo7 *sts             i: calibration switches and info
CoordInfo **coords         o: size and coordinate info for output
int *minorder, *maxorder   o: minimum and maximum values of SPORDER
*/

	int status;

	TblInfo tabinfo;	/* pointer to table descriptor, etc */
	TblRow tabrow;		/* values read from a table row */

	int row;		/* loop index */
	int RangeCoord (CoordInfo **, int *, int *);

	/* Open the spectroscopic distortion information table. */
	if ((status = OpenSDistTab (sts->distntab.name, &tabinfo)))
	    return (status);

	for (row = 1;  row <= tabinfo.nrows;  row++) {

	    if ((status = ReadSDistTab (&tabinfo, row, &tabrow)))
		return (status);

	    /* Check for a match with aperture, opt_elem and cenwave. */

	    if (SameString (tabrow.aperture, sts->aperture) &&
		SameString (tabrow.opt_elem, sts->opt_elem) &&
		SameInt (tabrow.cenwave, sts->cenwave)) {

		/* Get pedigree & descrip from the row. */
		if ((status = RowPedigree (&sts->distntab, row,
                        tabinfo.tp, tabinfo.cp_pedigree, tabinfo.cp_descrip)))
		    return (status);
		if (sts->distntab.goodPedigree == DUMMY_PEDIGREE) {
		    printf ("Warning  DUMMY pedigree in row %d of %s.\n",
			row, sts->distntab.name);
		    sts->x2dcorr_o = DUMMY;
		    CloseSDistTab (&tabinfo);
		    return (0);
		}

		/* Read data from this row. */
		if ((status = ReadSDistArray (&tabinfo, row, coords)))
		    return (status);

		/* It's only for echelle data that we need multiple rows.
		   This test was added because there may be two matching
		   rows if one of the CCD pseudo-apertures was used, the
		   actual aperture name and aperture=ANY, and in this case
		   we must stop looking for matching rows after finding
		   the pseudo-aperture name.
		*/
		if ((*coords)->sporder == 1)
		    break;
	    }
	}

	/* Get the range of order numbers. */
	if ((status = RangeCoord (coords, minorder, maxorder))) {
	    if (status < 0) {
		printf ("Warning  Matching row not found in SDCTAB %s; \\\n",
				sts->distntab.name);
		printf ("Warning  OPT_ELEM %s, CENWAVE %d.\n",
		    sts->opt_elem, sts->cenwave);
		sts->x2dcorr_o = OMIT;
	    } else {
		return (status);
	    }
	}

	if ((status = CloseSDistTab (&tabinfo)))
	    return (status);

	return (0);
}

/* This routine opens the 2-D distortion table, finds the columns
   that we need, and gets the total number of rows in the table.
*/

static int OpenSDistTab (char *tname, TblInfo *tabinfo) {

	tabinfo->tp = c_tbtopn (tname, IRAF_READ_ONLY, 0);
	if (c_iraferr()) {
	    printf ("ERROR    SDCTAB `%s' not found.\n", tname);
	    return (OPEN_FAILED);
	}

	tabinfo->nrows = c_tbpsta (tabinfo->tp, TBL_NROWS);

	/* Find the columns. */

	c_tbcfnd1 (tabinfo->tp, "APERTURE", &tabinfo->cp_aperture);
	c_tbcfnd1 (tabinfo->tp, "OPT_ELEM", &tabinfo->cp_opt_elem);
	c_tbcfnd1 (tabinfo->tp, "CENWAVE", &tabinfo->cp_cenwave);

	c_tbcfnd1 (tabinfo->tp, "SPORDER", &tabinfo->cp_sporder);
	c_tbcfnd1 (tabinfo->tp, "A2CENTER", &tabinfo->cp_a2center);
	c_tbcfnd1 (tabinfo->tp, "NPIX1", &tabinfo->cp_npix[0]);
	c_tbcfnd1 (tabinfo->tp, "NPIX2", &tabinfo->cp_npix[1]);
	c_tbcfnd1 (tabinfo->tp, "CRPIX1", &tabinfo->cp_crpix[0]);
	c_tbcfnd1 (tabinfo->tp, "CRPIX2", &tabinfo->cp_crpix[1]);
	c_tbcfnd1 (tabinfo->tp, "CRVAL1", &tabinfo->cp_crval[0]);
	c_tbcfnd1 (tabinfo->tp, "CRVAL2", &tabinfo->cp_crval[1]);
	c_tbcfnd1 (tabinfo->tp, "CDELT1", &tabinfo->cp_cdelt[0]);
	c_tbcfnd1 (tabinfo->tp, "CDELT2", &tabinfo->cp_cdelt[1]);

	/* aperture is optional */

	if (tabinfo->cp_opt_elem == 0 || tabinfo->cp_cenwave == 0 ||
	    tabinfo->cp_sporder == 0  || tabinfo->cp_a2center == 0 ||
	    tabinfo->cp_npix[0] == 0  || tabinfo->cp_npix[1] == 0 ||
	    tabinfo->cp_crpix[0] == 0 || tabinfo->cp_crpix[1] == 0 ||
	    tabinfo->cp_crval[0] == 0 || tabinfo->cp_crval[1] == 0 ||
	    tabinfo->cp_cdelt[0] == 0 || tabinfo->cp_cdelt[1] == 0) {

	    c_tbtclo (tabinfo->tp);
	    printf ("ERROR    Column not found in SDCTAB.\n");
	    return (COLUMN_NOT_FOUND);
	}

	/* Pedigree and descrip are optional columns. */
	c_tbcfnd1 (tabinfo->tp, "PEDIGREE", &tabinfo->cp_pedigree);
	c_tbcfnd1 (tabinfo->tp, "DESCRIP", &tabinfo->cp_descrip);

	return (0);
}

/* This routine reads the columns (OPT_ELEM and CENWAVE) used to select
   the correct row.
*/

static int ReadSDistTab (TblInfo *tabinfo, int row, TblRow *tabrow) {

	/* The APERTURE column is optional. */
	if (tabinfo->cp_aperture == 0) {
	    strcpy (tabrow->aperture, "ANY");
	} else {
	    c_tbegtt (tabinfo->tp, tabinfo->cp_aperture, row,
			tabrow->aperture, STIS_CBUF);
	    if (c_iraferr())
		return (TABLE_ERROR);
	}

	c_tbegtt (tabinfo->tp, tabinfo->cp_opt_elem, row,
			tabrow->opt_elem, STIS_CBUF);
	if (c_iraferr())
	    return (TABLE_ERROR);

	c_tbegti (tabinfo->tp, tabinfo->cp_cenwave, row, &tabrow->cenwave);
	if (c_iraferr())
	    return (TABLE_ERROR);

	return (0);
}

/* This routine reads the data from one row into the coords structure.
   Several scalar column values are gotten.  We're not actually reading
   arrays; the name ReadSDistArray was chosen for consistency with other
   similar routines.
*/

static int ReadSDistArray (TblInfo *tabinfo, int row, CoordInfo **coords) {

	int status;

	CoordInfo *newrec;
	int NewCoord (CoordInfo **, CoordInfo *);

	if ((newrec = malloc (sizeof (CoordInfo))) == NULL) {
	    printf ("ERROR    (GetSDC) can't allocate memory.\n");
	    return (OUT_OF_MEMORY);
	}
	newrec->next = NULL;

	/* Get size info and coordinate parameters. */
	c_tbegti (tabinfo->tp, tabinfo->cp_sporder, row, &newrec->sporder);
	c_tbegtd (tabinfo->tp, tabinfo->cp_a2center, row, &newrec->a2center);
	c_tbegti (tabinfo->tp, tabinfo->cp_npix[0], row, &newrec->npix[0]);
	c_tbegti (tabinfo->tp, tabinfo->cp_npix[1], row, &newrec->npix[1]);
	c_tbegtd (tabinfo->tp, tabinfo->cp_crpix[0], row, &newrec->crpix[0]);
	c_tbegtd (tabinfo->tp, tabinfo->cp_crpix[1], row, &newrec->crpix[1]);
	c_tbegtd (tabinfo->tp, tabinfo->cp_crval[0], row, &newrec->crval[0]);
	c_tbegtd (tabinfo->tp, tabinfo->cp_crval[1], row, &newrec->crval[1]);
	c_tbegtd (tabinfo->tp, tabinfo->cp_cdelt[0], row, &newrec->cdelt[0]);
	c_tbegtd (tabinfo->tp, tabinfo->cp_cdelt[1], row, &newrec->cdelt[1]);
	if (c_iraferr())
	    return (TABLE_ERROR);

	/* Convert crpix and a2center to zero-indexed. */
	newrec->crpix[0]--;
	newrec->crpix[1]--;
	newrec->a2center--;

	/* Insert in the coords list. */
	if ((status = NewCoord (coords, newrec)))
	    return (status);

	free (newrec);

	return (0);
}

/* This routine closes the SDCTAB table. */

static int CloseSDistTab (TblInfo *tabinfo) {

	c_tbtclo (tabinfo->tp);
	if (c_iraferr())
	    return (TABLE_ERROR);

	return (0);
}
