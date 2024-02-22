# include <stdio.h>
# include <stdlib.h>	/* malloc */

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

	IRAFPointer cp_sporder;
	IRAFPointer cp_a2center;
	IRAFPointer cp_cdelt2;
	IRAFPointer cp_npix;

	IRAFPointer cp_pedigree;
	IRAFPointer cp_descrip;
	int nrows;			/* number of rows in table */
} TblInfo;

typedef struct {
	char opt_elem[STIS_CBUF];	/* optical element name */
	int cenwave;			/* central wavelength */
} TblRow;

static int OpenSDistTab (char *, TblInfo *);
static int ReadSDistTab (TblInfo *, int, TblRow *);
static int ReadSDistArray (TblInfo *, int, CoordInfo **);
static int CloseSDistTab (TblInfo *);

/* This routine reads the coordinate information from the spectroscopic
   distortion table SDCTAB.

   The distortion information table should contain the following:
	header parameters:
		none needed
	columns:
		OPT_ELEM:  grating (or mirror) name (string)
		CENWAVE:   central wavelength (int)
		SPORDER:   order number (int)
		A2CENTER:  Y location on detector corresponding to CRPIX2
		CDELT2:    increment per pixel (double)
		NPIX2:     size of rectified image (int)

   The table is read to find all rows for which the values of OPT_ELEM
   and CENWAVE are the same as in the input image header.  There should be
   one or more such rows, corresponding to different values of spectral
   order SPORDER.  All these rows are read into memory, pointed to by
   coords.  The minimum and maximum spectral order numbers are returned.
   It is an error to have duplicate values of SPORDER in the SDCTAB table,
   and all values between minorder and maxorder must be present.  The
   SDCTAB table need not be sorted.

   Calstis6 basically needs from the SDCTAB table two quantities: the plate
   scale value CDELT2 in the spatial (slit) direction, and the size of the
   rectified image in the A2 direction created by calstis7. The plate scale
   is used for correcting offsets due to POSTARG2. The CDELT1 plate scale
   is not needed because offsets in the dispersion direction must be input
   directly in arcsec, the unit used by the incidence angle correction
   formula. The rectified image size is used by the scattered light
   correction algorithm (IB).

   Note:
	Memory is allocated for the coords list; it should be freed
	by calling FreeCoord6.
	coords should have been initialized to NULL.




   Revision history:
   ----------------
   14 Nov 97  -  Borrowed from calstis7 (I.Busko)
   14 Nov 97  -  Removed input from some columns used in calstis7 (IB)
   14 Nov 97  -  Rename routine to avoid conflict with cs7 (IB)
   07 Jan 00  -  Read NPIX2 (IB)
   05 Sep 00  -  For first-order data, break out of the loop after finding
                 any matching row (PH)
*/

int GetSDC6 (StisInfo6 *sts, CoordInfo **coords, int *minorder, int *maxorder)
{

/* arguments:
StisInfo6 *sts             i: calibration switches and info
CoordInfo **coords         o: size and coordinate info for output
int *minorder, *maxorder   o: minimum and maximum values of SPORDER
*/

	int status;

	TblInfo tabinfo;	/* pointer to table descriptor, etc */
	TblRow tabrow;		/* values read from a table row */

	int row;		/* loop index */
	int RangeCoord6 (CoordInfo **, int *, int *);

	/* Open the spectroscopic distortion information table. */
	if ((status = OpenSDistTab (sts->distntab.name, &tabinfo)))
	    return (status);

	for (row = 1;  row <= tabinfo.nrows;  row++) {

	    if ((status = ReadSDistTab (&tabinfo, row, &tabrow)))
		return (status);

	    /* Check for a match with opt_elem and cenwave. */
	    if (SameString (tabrow.opt_elem, sts->opt_elem) &&
		SameInt (tabrow.cenwave, sts->cenwave)) {

		/* Get pedigree & descrip from the row. */
		if ((status = RowPedigree (&sts->distntab, row,
                       tabinfo.tp, tabinfo.cp_pedigree, tabinfo.cp_descrip)))
		    return (status);
		if (sts->distntab.goodPedigree == DUMMY_PEDIGREE) {
		    printf ("Warning  DUMMY pedigree in row %d of %s.\n",
			row, sts->distntab.name);
		    sts->x1d_o = DUMMY;
		    CloseSDistTab (&tabinfo);
		    return (0);
		}

		/* Read data from this row. */
		if ((status = ReadSDistArray (&tabinfo, row, coords)))
		    return (status);

		/* This was added because the SDCTAB may contain multiple
		   entries for the CCD pseudo-apertures.  In this case, all
		   matching rows contain the same NPIX2 value, so we can
		   take the value from the first matching row.
		*/
		if ((*coords)->sporder == 1)
		    break;
	    }
	}

	/* Get the range of order numbers. */
	if ((status = RangeCoord6 (coords, minorder, maxorder))) {
	    if (status < 0) {
		printf ("Warning  Matching row not found in SDCTAB %s; \\\n",
				sts->distntab.name);
		printf ("Warning  OPT_ELEM %s, CENWAVE %d.\n",
		    sts->opt_elem, sts->cenwave);
		sts->x1d_o = OMIT;
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

	c_tbcfnd1 (tabinfo->tp, "OPT_ELEM", &tabinfo->cp_opt_elem);
	c_tbcfnd1 (tabinfo->tp, "CENWAVE", &tabinfo->cp_cenwave);

	c_tbcfnd1 (tabinfo->tp, "SPORDER", &tabinfo->cp_sporder);
	c_tbcfnd1 (tabinfo->tp, "A2CENTER", &tabinfo->cp_a2center);
	c_tbcfnd1 (tabinfo->tp, "CDELT2", &tabinfo->cp_cdelt2);
	c_tbcfnd1 (tabinfo->tp, "NPIX2", &tabinfo->cp_npix);

	if (tabinfo->cp_opt_elem == 0 || tabinfo->cp_cenwave == 0 ||
	    tabinfo->cp_sporder == 0  || tabinfo->cp_a2center == 0 ||
	    tabinfo->cp_cdelt2 == 0   || tabinfo->cp_npix == 0) {

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

	c_tbegtt (tabinfo->tp, tabinfo->cp_opt_elem, row,
			tabrow->opt_elem, STIS_CBUF-1);
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
	int NewCoord6 (CoordInfo **, CoordInfo *);

	if ((newrec = malloc (sizeof (CoordInfo))) == NULL) {
	    printf ("ERROR    (GetSDC) can't allocate memory.\n");
	    return (OUT_OF_MEMORY);
	}
	newrec->next = NULL;

	/* Get size info and coordinate parameters. */
	c_tbegti (tabinfo->tp, tabinfo->cp_sporder, row, &newrec->sporder);
	c_tbegtd (tabinfo->tp, tabinfo->cp_a2center, row, &newrec->a2center);
	c_tbegtd (tabinfo->tp, tabinfo->cp_cdelt2, row, &newrec->cdelt2);
	c_tbegti (tabinfo->tp, tabinfo->cp_npix, row, &newrec->npix);
	if (c_iraferr())
	    return (TABLE_ERROR);

	/* Insert in the coords list. */
	if ((status = NewCoord6 (coords, newrec)))
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

