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
	IRAFPointer cp_aperture;
	IRAFPointer cp_opt_elem;
	IRAFPointer cp_cenwave;

	IRAFPointer cp_sporder;
	IRAFPointer cp_extrsize;
	IRAFPointer cp_bksize[2];
	IRAFPointer cp_bkoffset[2];
	IRAFPointer cp_ncoeffsl;
	IRAFPointer cp_sltcoeff;
	IRAFPointer cp_ncoeffbk;
	IRAFPointer cp_bktcoeff;
	IRAFPointer cp_backord;
	IRAFPointer cp_xtracalg;
	IRAFPointer cp_maxsearch;

	IRAFPointer cp_pedigree;
	IRAFPointer cp_descrip;
	int nrows;			/* number of rows in table */
} TblInfo;

typedef struct {
	char aperture[STIS_CBUF];	/* aperture name */
	char opt_elem[STIS_CBUF];	/* optical element name */
	int cenwave;			/* central wavelength */
} TblRow;


static int OpenXtractTab (char *, TblInfo *);
static int ReadXtractTab (TblInfo *, int, TblRow *);
static int ReadXtractArray (TblInfo *, int, XtractInfo **);
static int CloseXtractTab (TblInfo *);


/* This routine reads the coordinate information from the 1-D extraction
   table XTRACTAB.

   The extraction information table should contain the following:
	header parameters:
		none needed
	columns:
		APERTURE: aperture name (string)
		OPT_ELEM: grating (or mirror) name (string)
		CENWAVE: central wavelength (int)
		SPORDER: order number (int)
		EXTRSIZE: size of spectrum extraction box (float)
		BK1SIZE, BK2SIZE: size of background extraction boxes (float)
		BK1OFFST, BK2OFFST: offset of background extraction boxes
				      from spectrum extraction box (float)
	        SLTCOEFF: coefficient array that describes spectrum box tilt
	        NCOEFFSL: size of above array
	        BKTCOEFF: coefficient array that describes backgr. box tilt
	        NCOEFFBK: size of above array
		BACKORD: order of polynomial to fit to background (short)
		XTRACALG: extraction algorithm (string)
		MAXSRCH:  maximum search range for cross correlation (short)

   The table is read to find all rows for which the values of APERTURE,
   OPT_ELEM and CENWAVE are the same as in the input image header. There
   should be one or more such rows, corresponding to different values of
   spectral order SPORDER.  All these rows are read into memory, pointed
   to by extract.  The minimum and maximum spectral order numbers are returned.
   It is an error to have duplicate values of SPORDER in the XTRACTAB table,
   and all values between minorder and maxorder must be present. The XTRACTAB
   table need not be sorted.

   Note:
	Memory is allocated for the extract list; it should be freed
	by calling FreeXtract.
	extract should have been initialized to NULL.




   Revision history:
   ----------------
   20 Feb 97  -  Implemented using getsdc in calstis7 as model (I.Busko)
   10 Apr 97  -  Changes after code review (IB):
                 - replaced literal by ROW_NOT_FOUND constant
                 - explicit cast for malloc-returned pointer.
                 - rename "new" variable to avoid conflict in C++
   16 Apr 97  -  Replaced scalar bktilt by a polynomial description. Also
                 implemented tilted spectrum extraction box (what for ?) (IB)
   22 Apr 97  -  Changed column name from MAXSEARCH to MAXSRCH (IB)
   02 May 97  -  Set x1d_o flag for rows with DUMMY pedigree (IB)
   08 May 97  -  Conform to new _trl standard (IB)
   01 Dec 00  -  For first-order data, break out of the loop after finding
                 any matching row (PH, IB)
*/

int GetXtract (StisInfo6 *sts, XtractInfo **extract, int *minorder,
               int *maxorder) {

/* arguments:
StisInfo6 *sts             i: calibration switches and info
XtractInfo **extract       o: size and coordinate info for output
int *minorder, *maxorder   o: minimum and maximum values of SPORDER
*/

	int status;

	TblInfo tabinfo;	/* pointer to table descriptor, etc */
	TblRow tabrow;		/* values read from a table row */

	int row;		/* loop index */
	int RangeXtract (XtractInfo **, int *, int *);

	/* Open the 1-D extraction table. */
	if ((status = OpenXtractTab (sts->xtrctab.name, &tabinfo)))
	    return (status);

	for (row = 1;  row <= tabinfo.nrows;  row++) {

	    if ((status = ReadXtractTab (&tabinfo, row, &tabrow)))
		return (status);

	    /* Check for a match with aperture, opt_elem and cenwave. */

	    if (SameString (tabrow.aperture, sts->aperture) &&
	        SameString (tabrow.opt_elem, sts->opt_elem) &&
		SameInt (tabrow.cenwave, sts->cenwave)) {

		/* Get pedigree & descrip from the row. */
		if ((status = RowPedigree (&sts->xtrctab, row,
                        tabinfo.tp, tabinfo.cp_pedigree, tabinfo.cp_descrip)))
		    return (status);
		if (sts->xtrctab.goodPedigree == DUMMY_PEDIGREE) {
		    sts->x1d_o = DUMMY;
		    CloseXtractTab (&tabinfo);
		    return (0);
		}

		/* Read data from this row. */
		if ((status = ReadXtractArray (&tabinfo, row, extract)))
		    return (status);

		/* This was added because the XTRACTAB may contain multiple
		   entries for the CCD pseudo-apertures.  In this case, all
		   matching rows contain the same NPIX2 value, so we can
		   take the value from the first matching row.
		*/
		if ((*extract)->sporder == 1)
		    break;
	    }
	}

	/* Get the range of order numbers. */
	if ((status = RangeXtract (extract, minorder, maxorder))) {
	    if (status < 0) {
		printf ("ERROR    Matching row not found in XTRACTAB %s\n",
				sts->xtrctab.name);
		printf ("ERROR    APERTURE %s, OPT_ELEM %s, CENWAVE %d\n",
		    sts->aperture, sts->opt_elem, sts->cenwave);
		return (ROW_NOT_FOUND);
	    } else {
		return (status);
	    }
	}

	if ((status = CloseXtractTab (&tabinfo)))
	    return (status);

	return (0);
}



/* This routine opens the 1-D extraction table, finds the columns
   that we need, and gets the total number of rows in the table.
*/

static int OpenXtractTab (char *tname, TblInfo *tabinfo) {

	tabinfo->tp = c_tbtopn (tname, IRAF_READ_ONLY, 0);
	if (c_iraferr()) {
	    printf ("ERROR    XTRACTAB `%s' not found\n", tname);
	    return (OPEN_FAILED);
	}

	tabinfo->nrows = c_tbpsta (tabinfo->tp, TBL_NROWS);

	/* Find the columns. */

	c_tbcfnd1 (tabinfo->tp, "APERTURE", &tabinfo->cp_aperture);
	c_tbcfnd1 (tabinfo->tp, "OPT_ELEM", &tabinfo->cp_opt_elem);
	c_tbcfnd1 (tabinfo->tp, "CENWAVE",  &tabinfo->cp_cenwave);
	c_tbcfnd1 (tabinfo->tp, "SPORDER",  &tabinfo->cp_sporder);

	c_tbcfnd1 (tabinfo->tp, "EXTRSIZE", &tabinfo->cp_extrsize);
	c_tbcfnd1 (tabinfo->tp, "BK1SIZE",  &tabinfo->cp_bksize[0]);
	c_tbcfnd1 (tabinfo->tp, "BK2SIZE",  &tabinfo->cp_bksize[1]);
	c_tbcfnd1 (tabinfo->tp, "BK1OFFST", &tabinfo->cp_bkoffset[0]);
	c_tbcfnd1 (tabinfo->tp, "BK2OFFST", &tabinfo->cp_bkoffset[1]);
	c_tbcfnd1 (tabinfo->tp, "BKTCOEFF", &tabinfo->cp_bktcoeff);
	c_tbcfnd1 (tabinfo->tp, "NCOEFFBK", &tabinfo->cp_ncoeffbk);
	c_tbcfnd1 (tabinfo->tp, "SLTCOEFF", &tabinfo->cp_sltcoeff);
	c_tbcfnd1 (tabinfo->tp, "NCOEFFSL", &tabinfo->cp_ncoeffsl);
	c_tbcfnd1 (tabinfo->tp, "BACKORD",  &tabinfo->cp_backord);
	c_tbcfnd1 (tabinfo->tp, "XTRACALG", &tabinfo->cp_xtracalg);
	c_tbcfnd1 (tabinfo->tp, "MAXSRCH",  &tabinfo->cp_maxsearch);

	if (tabinfo->cp_opt_elem    == 0 ||
            tabinfo->cp_cenwave     == 0 ||
	    tabinfo->cp_sporder     == 0 ||
            tabinfo->cp_aperture    == 0 ||
	    tabinfo->cp_extrsize    == 0 ||
            tabinfo->cp_bktcoeff    == 0 ||
            tabinfo->cp_ncoeffbk    == 0 ||
            tabinfo->cp_sltcoeff    == 0 ||
            tabinfo->cp_ncoeffsl    == 0 ||
	    tabinfo->cp_bksize[0]   == 0 ||
            tabinfo->cp_bksize[1]   == 0 ||
	    tabinfo->cp_bksize[0]   == 0 ||
            tabinfo->cp_bksize[1]   == 0 ||
	    tabinfo->cp_bkoffset[0] == 0 ||
            tabinfo->cp_bkoffset[1] == 0 ||
	    tabinfo->cp_backord     == 0 ||
            tabinfo->cp_xtracalg    == 0 ||
	    tabinfo->cp_maxsearch   == 0) {
	    c_tbtclo (tabinfo->tp);
	    printf ("ERROR    Column not found in XTRACTAB\n");
	    return (COLUMN_NOT_FOUND);
	}

	/* Pedigree and descrip are optional columns. */
	c_tbcfnd1 (tabinfo->tp, "PEDIGREE", &tabinfo->cp_pedigree);
	c_tbcfnd1 (tabinfo->tp, "DESCRIP",  &tabinfo->cp_descrip);

	return (0);
}



/* This routine reads the columns (APERTURE, OPT_ELEM and CENWAVE) used
   to select the correct row.
*/

static int ReadXtractTab (TblInfo *tabinfo, int row, TblRow *tabrow) {

	c_tbegtt (tabinfo->tp, tabinfo->cp_aperture, row, tabrow->aperture,
                  STIS_CBUF-1);
	if (c_iraferr())
	    return (TABLE_ERROR);

	c_tbegtt (tabinfo->tp, tabinfo->cp_opt_elem, row, tabrow->opt_elem,
                  STIS_CBUF-1);
	if (c_iraferr())
	    return (TABLE_ERROR);

	c_tbegti (tabinfo->tp, tabinfo->cp_cenwave, row, &tabrow->cenwave);
	if (c_iraferr())
	    return (TABLE_ERROR);

	return (0);
}



/* This routine reads the data from one row into the extract structure.
   Several scalar columns and two arrays are gotten.
*/

static int ReadXtractArray (TblInfo *tabinfo, int row, XtractInfo **extract) {

	int status;

	int ncoeff;		/* number of coefficients read from table */
	XtractInfo *newd;
	int NewXtract (XtractInfo **, XtractInfo *);

	if ((newd = (XtractInfo *) malloc (sizeof (XtractInfo))) == NULL) {
	    printf ("ERROR    Can't allocate memory.\n");
	    return (OUT_OF_MEMORY);
	}
	newd->next = NULL;

	/* Get extraction parameters. */
	c_tbegts (tabinfo->tp, tabinfo->cp_sporder, row, &newd->sporder);
	c_tbegtr (tabinfo->tp, tabinfo->cp_extrsize, row, &newd->extrsize);
	c_tbegtr (tabinfo->tp, tabinfo->cp_bksize[0], row, &newd->bksize[0]);
	c_tbegtr (tabinfo->tp, tabinfo->cp_bksize[1], row, &newd->bksize[1]);
	c_tbegtr (tabinfo->tp, tabinfo->cp_bkoffset[0], row,
                  &newd->bkoffset[0]);
	c_tbegtr (tabinfo->tp, tabinfo->cp_bkoffset[1], row,
                  &newd->bkoffset[1]);
	c_tbegti (tabinfo->tp, tabinfo->cp_ncoeffsl, row, &newd->ncoeffsl);

	if (newd->ncoeffsl > MAX_SLIT_COEFF) {
	    printf (
               "ERROR    Too many slit tilt coefficients %d in XTRACTAB.\n",
		newd->ncoeffsl);
	    return (TABLE_ERROR);
	}
	ncoeff = c_tbagtd (tabinfo->tp, tabinfo->cp_sltcoeff, row,
                           newd->sltcoeff, 1, newd->ncoeffsl);
	if (c_iraferr())
	    return (TABLE_ERROR);

	c_tbegti (tabinfo->tp, tabinfo->cp_ncoeffbk, row, &newd->ncoeffbk);
	if (newd->ncoeffbk > MAX_BACK_COEFF) {
	    printf (
            "ERROR    Too many background tilt coefficients %d in XTRACTAB.\n",
		newd->ncoeffbk);
	    return (TABLE_ERROR);
	}
	ncoeff = c_tbagtd (tabinfo->tp, tabinfo->cp_bktcoeff, row,
                           newd->bktcoeff, 1, newd->ncoeffbk);
	if (c_iraferr())
	    return (TABLE_ERROR);

	c_tbegts (tabinfo->tp, tabinfo->cp_backord, row, &newd->backord);
	c_tbegtt (tabinfo->tp, tabinfo->cp_xtracalg, row, newd->xtracalg,
                  STIS_CBUF);
	c_tbegts (tabinfo->tp, tabinfo->cp_maxsearch, row, &newd->maxsearch);
	if (c_iraferr())
	    return (TABLE_ERROR);

	/* Insert in the extract list. */
	if ((status = NewXtract (extract, newd)))
	    return (status);

	free (newd);

	return (0);
}



/* This routine closes the XTRACTAB table. */

static int CloseXtractTab (TblInfo *tabinfo) {

	c_tbtclo (tabinfo->tp);
	if (c_iraferr())
	    return (TABLE_ERROR);

	return (0);
}
