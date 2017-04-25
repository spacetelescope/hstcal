# include <stdio.h>
# include <stdlib.h>		/* calloc */
# include <string.h>
# include <math.h>		/* fabs, sqrt */

# include "c_iraf.h"
# include "hstio.h"
# include "xtables.h"

# include "stis.h"
# include "calstis7.h"
# include "hstcalerr.h"
# include "stisdef.h"

typedef struct {
	IRAFPointer tp;			/* pointer to table descriptor */
	IRAFPointer cp_direction;	/* column descriptors */
	IRAFPointer cp_filter;

	IRAFPointer cp_scale;
	IRAFPointer cp_xref;
	IRAFPointer cp_yref;
	IRAFPointer cp_cxsize;
	IRAFPointer cp_cysize;
	IRAFPointer cp_cxref;
	IRAFPointer cp_cyref;

	IRAFPointer cp_pedigree;
	IRAFPointer cp_descrip;
	int nrows;			/* number of rows in table */
} TblInfo;

typedef struct {
	char direction[STIS_CBUF+1];	/* forward or inverse mapping */
	char filter[STIS_CBUF+1];	/* name of filter, or CLEAR */
} TblRow;

static int OpenDistTab (char *, TblInfo *);
static int ReadDistTab (TblInfo *, int, TblRow *);
static int ReadCoefficients (TblInfo *, int, CoordInfo **, DistInfo *);
static int CloseDistTab (TblInfo *);

/* This routine reads the distortion information from the imaging
   distortion table IDCTAB.  This is only used for obstype=IMAGING.

   The distortion information table should contain the following:
	header parameters:
		NORDER:  polynomial order (int)
	columns:
		DIRECTION:  direction of mapping (string)
		FILTER:  name of filter (string) (not required)
		XREF, YREF:  zero point in input image (double)
		CXREF, CYREF:  zero point in output image (double)
		CXSIZE, CYSIZE:  size for output image (int)
			(saved in both dist and coords, as npix[0] & npix[1])
		SCALE:  arcseconds / pixel (double)
		CX00, CY00, etc:  arrays of coefficients (double)

   The table is read to find the first row for which the value of
   DIRECTION is "INVERSE".  If an FILTER column exists, that column
   will also be used to select the row.

   Note that for all columns that include pixel in the units, the size is
   that of a reference pixel, not necessarily an image pixel.

   The output image size is read from CXSIZE & CYSIZE, and this is the
   only (useful) information saved in the CoordInfo structure.  It's
   copied there because Do2Dx uses it when creating the output image.

   Note:
	Memory is allocated for the coords list; it should be freed
	by calling FreeCoord.
	coords should have been initialized to NULL.
	Memory is allocated for the distortion coefficients, dist->xcoeff
	and dist->ycoeff.  The memory should be freed by calling freeDist.

   Phil Hodge, 2000 Jan 13:
	Add one to aperture and opt_elem buffer sizes.

   Phil Hodge, 2001 Jan 3:
	Remove references to aperture; select row on DIRECTION = "INVERSE"
	and OPT_ELEM (except that the latter is optional).
	Get columns CX00, CX10, ... CY00, CY10, ... .

   Phil Hodge, 2001 Mar 9:
	The columns for the linear coefficients are now required to exist;
	these are CX10, CX11, CY10, CY11.  All other coefficients can take
	the default value of zero.

   Phil Hodge, 2001 Aug 8:
	Select the row based on FILTER instead of OPT_ELEM.
*/

int GetIDC (StisInfo7 *sts, CoordInfo **coords, DistInfo *dist) {

/* arguments:
StisInfo7 *sts      i: calibration switches and info
CoordInfo **coords  o: coordinate parameters for output image
DistInfo *dist      o: distortion coefficients
*/

	int status;

	TblInfo tabinfo;	/* pointer to table descriptor, etc */
	TblRow tabrow;		/* values read from a table row */

	int row;		/* loop index */
	int foundit = 0;	/* true if parameters found in table */

	/* Open the distortion information table. */
	if ((status = OpenDistTab (sts->distntab.name, &tabinfo)))
	    return (status);

	/* Find the row for the inverse mapping. */

	for (row = 1;  row <= tabinfo.nrows;  row++) {

	    if ((status = ReadDistTab (&tabinfo, row, &tabrow)))
		return (status);

	    if (SameString (tabrow.direction, "INVERSE") &&
		SameString (tabrow.filter, sts->filter)) {

		/* This is the row; read the info. */
		foundit = 1;

		/* Get pedigree & descrip from the row. */
		if ((status = RowPedigree (&sts->distntab, row,
                        tabinfo.tp, tabinfo.cp_pedigree, tabinfo.cp_descrip)))
		    return (status);
		if (sts->distntab.goodPedigree == DUMMY_PEDIGREE) {
		    printf ("Warning  DUMMY pedigree in row %d of %s.\n",
			row, sts->distntab.name);
		    sts->x2dcorr_o = DUMMY;
		    CloseDistTab (&tabinfo);
		    return (0);
		}

		if ((status = ReadCoefficients (&tabinfo, row, coords, dist)))
		    return (status);

		break;
	    }
	}

	if (!foundit) {
	    printf ("Warning  Matching row not found in IDCTAB %s \\\n",
				sts->distntab.name);

	    printf ("Warning  for DIRECTION = 'INVERSE')\n");
	    sts->x2dcorr_o = OMIT;
	}

	if ((status = CloseDistTab (&tabinfo)))
	    return (status);

	return (0);
}

/* This routine opens the 2-D distortion table, finds the columns
   that we need, and gets the total number of rows in the table.
*/

static int OpenDistTab (char *tname, TblInfo *tabinfo) {

	IRAFPointer cp[4];	/* for checking whether columns exist */

	tabinfo->tp = c_tbtopn (tname, IRAF_READ_ONLY, 0);
	if (c_iraferr()) {
	    printf ("ERROR    IDCTAB `%s' not found.\n", tname);
	    return (OPEN_FAILED);
	}

	tabinfo->nrows = c_tbpsta (tabinfo->tp, TBL_NROWS);

	/* Find the columns. */

	c_tbcfnd1 (tabinfo->tp, "DIRECTION", &tabinfo->cp_direction);
	c_tbcfnd1 (tabinfo->tp, "FILTER", &tabinfo->cp_filter);

	c_tbcfnd1 (tabinfo->tp, "SCALE", &tabinfo->cp_scale);
	c_tbcfnd1 (tabinfo->tp, "XREF", &tabinfo->cp_xref);
	c_tbcfnd1 (tabinfo->tp, "YREF", &tabinfo->cp_yref);
	c_tbcfnd1 (tabinfo->tp, "CXSIZE", &tabinfo->cp_cxsize);
	c_tbcfnd1 (tabinfo->tp, "CYSIZE", &tabinfo->cp_cysize);
	c_tbcfnd1 (tabinfo->tp, "CXREF", &tabinfo->cp_cxref);
	c_tbcfnd1 (tabinfo->tp, "CYREF", &tabinfo->cp_cyref);

	if (tabinfo->cp_direction == 0 ||
	    tabinfo->cp_scale == 0 ||
	    tabinfo->cp_xref == 0 || tabinfo->cp_yref == 0 ||
	    tabinfo->cp_cxsize == 0 || tabinfo->cp_cysize == 0 ||
	    tabinfo->cp_cxref == 0 || tabinfo->cp_cyref == 0) {

	    c_tbtclo (tabinfo->tp);
	    printf ("ERROR    Column not found in IDCTAB.\n");
	    return (COLUMN_NOT_FOUND);
	}

	/* Look for these columns just to make sure they're present. */
	c_tbcfnd1 (tabinfo->tp, "CX10", &cp[0]);
	c_tbcfnd1 (tabinfo->tp, "CX11", &cp[1]);
	c_tbcfnd1 (tabinfo->tp, "CY10", &cp[2]);
	c_tbcfnd1 (tabinfo->tp, "CY11", &cp[3]);
	if (cp[0] == 0 || cp[1] == 0 || cp[2] == 0 || cp[3] == 0) {
	    c_tbtclo (tabinfo->tp);
	    printf (
"ERROR    Columns for linear transformation coefficients are required; \\\n");
	    printf ("ERROR    these columns are CX10, CX11, CY10, CY11.\n");
	    return (COLUMN_NOT_FOUND);
	}

	/* Pedigree and descrip are optional columns. */
	c_tbcfnd1 (tabinfo->tp, "PEDIGREE", &tabinfo->cp_pedigree);
	c_tbcfnd1 (tabinfo->tp, "DESCRIP", &tabinfo->cp_descrip);

	return (0);
}

/* This routine reads the columns (DIRECTION and optionally FILTER)
   used to select the correct row.
*/

static int ReadDistTab (TblInfo *tabinfo, int row, TblRow *tabrow) {

	c_tbegtt (tabinfo->tp, tabinfo->cp_direction, row,
			tabrow->direction, STIS_CBUF);
	if (c_iraferr())
	    return (TABLE_ERROR);

	if (tabinfo->cp_filter == 0) {

	    /* default value, if column was not present in IDC table */
	    strcpy (tabrow->filter, "ANY");

	} else {

	    c_tbegtt (tabinfo->tp, tabinfo->cp_filter, row,
			tabrow->filter, STIS_CBUF);
	    if (c_iraferr())
		return (TABLE_ERROR);
	}

	return (0);
}

/* This routine reads the data from one row into the coords and dist
   structures.  The size and coordinate parameters are gotten into coords.
   The number of elements in the arrays of distortion coefficients is
   gotten into dist, memory is allocated, and the array values are read.
   The order of the polynomial coefficients is gotten from the NORDER
   header keyword.
*/

static int ReadCoefficients (TblInfo *tabinfo, int row,
		CoordInfo **coords, DistInfo *dist) {

	char colname[SZ_COLNAME];	/* a column name */
	IRAFPointer cp;		/* for getting coefficients */
	int i, j, k;
	int kx, ky;		/* indexes for coefficients of x and y */
	int npix[2];		/* the size of the output image */
	int status;

	CoordInfo *newrec;
	int NewCoord (CoordInfo **, CoordInfo *);

	if ((newrec = malloc (sizeof (CoordInfo))) == NULL) {
	    printf ("ERROR    Can't allocate memory.\n");
	    return (OUT_OF_MEMORY);
	}
	newrec->sporder = 1;	/* use nominal values for imaging mode */
	newrec->a2center = 1.;
	newrec->next = NULL;

	/* Get size info and coordinate parameters. */
	c_tbegti (tabinfo->tp, tabinfo->cp_cxsize, row, &npix[0]);
	c_tbegti (tabinfo->tp, tabinfo->cp_cysize, row, &npix[1]);
	if (c_iraferr())
	    return (TABLE_ERROR);

	newrec->npix[0] = npix[0];
	newrec->npix[1] = npix[1];

	/* dummy values (won't be used) */
	newrec->crpix[0] = 0.;
	newrec->crpix[1] = 0.;
	newrec->cdelt[0] = 1.;
	newrec->cdelt[1] = 1.;
	newrec->crval[0] = 0.;
	newrec->crval[1] = 0.;

	/* Insert in the coords list. */
	if ((status = NewCoord (coords, newrec)))
	    return (status);

	free (newrec);

	/* Allocate memory for the coefficients. */
	dist->xcoeff = calloc (MAX_NCOEFF, sizeof(double));
	dist->ycoeff = calloc (MAX_NCOEFF, sizeof(double));
	if (dist->xcoeff == NULL || dist->ycoeff == NULL) {
	    c_tbtclo (tabinfo->tp);
	    return (OUT_OF_MEMORY);
	}
	dist->allocated = 1;			/* set flag */

	/* Initialize.  It is not required that all columns be present;
	   missing columns result in default values for coefficients.
	*/

	for (k = 0;  k < MAX_NCOEFF;  k++) {
	    dist->xcoeff[k] = 0.;
	    dist->ycoeff[k] = 0.;
	}
	kx = WHICH_COEFF (1, 1);	/* CX11, coefficient of x */
	ky = WHICH_COEFF (1, 0);	/* CY10, coefficient of y */
	dist->xcoeff[kx] = 1.;		/* note:  just default values */
	dist->ycoeff[ky] = 1.;

	/* Copy image size. */
	dist->npix[0] = npix[0];
	dist->npix[1] = npix[1];

	/* The polynomial order is a header keyword. */
	dist->norder = c_tbhgti (tabinfo->tp, "NORDER");
	if (c_iraferr()) {
	    printf ("Warning  Can't read NORDER from IDCTAB table header.\n");
	    dist->norder = MAX_ORDER;
	    clear_cvoserr();
	}
	if (dist->norder > MAX_ORDER) {
	    printf (
	"ERROR    Polynomial order = %d in IDCTAB exceeds maximum of %d\n",
		dist->norder, MAX_ORDER);
	    return (GENERIC_ERROR_CODE);
	}

	/* Get the reference points for the input and output images (the
	   zero points for the mapping), and convert to zero indexing.
	*/

	/* input, distorted image */
	c_tbegtd (tabinfo->tp, tabinfo->cp_xref, row, &dist->xref);
	if (c_iraferr())
	    return (TABLE_ERROR);
	dist->xref--;
	c_tbegtd (tabinfo->tp, tabinfo->cp_yref, row, &dist->yref);
	if (c_iraferr())
	    return (TABLE_ERROR);
	dist->yref--;

	/* output, corrected image */
	c_tbegtd (tabinfo->tp, tabinfo->cp_cxref, row, &dist->cxref);
	if (c_iraferr())
	    return (TABLE_ERROR);
	dist->cxref--;
	c_tbegtd (tabinfo->tp, tabinfo->cp_cyref, row, &dist->cyref);
	if (c_iraferr())
	    return (TABLE_ERROR);
	dist->cyref--;

	/* Get the scale, because the CXij and CYij coefficients map from
	   arcseconds to pixels, and we need pixels to pixels.
	*/
	c_tbegtd (tabinfo->tp, tabinfo->cp_scale, row, &dist->scale);
	if (c_iraferr())
	    return (TABLE_ERROR);

	/* Get the coefficients. */

	for (i = 0;  i <= dist->norder;  i++) {

	    for (j = 0;  j <= i;  j++) {

		k = WHICH_COEFF (i, j);

		/* CX coefficient for x coordinate. */

		sprintf (colname, "CX%d%d", i, j);

		c_tbcfnd1 (tabinfo->tp, colname, &cp);
		if (cp != 0) {
		    c_tbegtd (tabinfo->tp, cp, row, &dist->xcoeff[k]);
		    if (c_iraferr())
			return (TABLE_ERROR);
		}

		/* CY coefficient for y coordinate. */

		sprintf (colname, "CY%d%d", i, j);

		c_tbcfnd1 (tabinfo->tp, colname, &cp);
		if (cp != 0) {
		    c_tbegtd (tabinfo->tp, cp, row, &dist->ycoeff[k]);
		    if (c_iraferr())
			return (TABLE_ERROR);
		}
	    }
	}

	return (0);
}

/* This routine closes the IDCTAB table. */

static int CloseDistTab (TblInfo *tabinfo) {

	c_tbtclo (tabinfo->tp);
	if (c_iraferr())
	    return (TABLE_ERROR);

	return (0);
}
