# include <stdio.h>
# include <stdlib.h>	/* calloc */

# include "c_iraf.h"
# include "hstio.h"
# include "xtables.h"

# include "stis.h"
# include "calstis7.h"
# include "err.h"
# include "stisdef.h"

typedef struct {
	IRAFPointer tp;			/* pointer to table descriptor */
	IRAFPointer cp_aperture;	/* column descriptors */
	IRAFPointer cp_nelem;
	IRAFPointer cp_wl;
	IRAFPointer cp_thr;
	IRAFPointer cp_pedigree;
	IRAFPointer cp_descrip;
	int nrows;			/* number of rows in table */
} TblInfo;

typedef struct {
	char aperture[STIS_CBUF+1];
} TblRow;

static int OpenThruTab (char *, TblInfo *);
static int ReadThruTab (TblInfo *, int, TblRow *);
static int ReadThruArray (TblInfo *, int, ApInfo *);
static int CloseThruTab (TblInfo *);

/* This routine gets aperture throughput from APERTAB (_apt) and saves
   the info in the aperture information structure.

   The aperture info table should contain the following:
	header parameters:
		none needed
	columns:
		APERTURE:  aperture name (string)
		NELEM:  actual number of elements in arrays (int)
		WAVELENGTH:  wavelengths in Angstroms (double)
		THROUGHPUT:  the fraction of light passed by the slit (float)
			(corresponding to WAVELENGTH)

   Rows are selected on APERTURE.  If a matching row is found, the sizes
   of the wavelength and throughput arrays are gotten from column NELEM,
   memory is allocated, and WAVELENGTH and THROUGHPUT are read.

   When done, the memory should be freed by calling freeThroughput.

   Phil Hodge, 2000 Jan 13:
	Add one to aperture buffer size.
*/

int GetApThr (StisInfo7 *sts, ApInfo *slit) {

/* arguments:
StisInfo7 *sts   i: calibration switches and info
ApInfo *slit     o: description of slit
*/

	int status;

	TblInfo tabinfo;	/* pointer to table descriptor, etc */
	TblRow tabrow;		/* values read from a table row */

	int row;		/* loop index */
	int foundit = 0;	/* true if parameters found in table */

	/* Open the aperture throughput table. */
	if ((status = OpenThruTab (sts->apertab.name, &tabinfo)))
	    return (status);

	for (row = 1;  row <= tabinfo.nrows;  row++) {

	    if ((status = ReadThruTab (&tabinfo, row, &tabrow)))
		return (status);

	    if (SameString (tabrow.aperture, sts->aperture)) {

		foundit = 1;

		/* Get pedigree & descrip from the row. */
		if ((status = RowPedigree (&sts->apertab, row,
                        tabinfo.tp, tabinfo.cp_pedigree, tabinfo.cp_descrip)))
		    return (status);
                if (sts->apertab.goodPedigree == DUMMY_PEDIGREE) {
                    printf ("Warning  DUMMY pedigree in row %d of %s.\n",
                            row, sts->apertab.name);
		    sts->x2dcorr_o = DUMMY;
		    CloseThruTab (&tabinfo);
		    return (0);
		}

		/* Read wavelengths and throughputs into slit structure. */
                if ((status = ReadThruArray (&tabinfo, row, slit)))
                    return (status);
	    }
	}

	if ((status = CloseThruTab (&tabinfo)))
	    return (status);

	if (!foundit) {
	    printf ("Warning  Matching row not found in APERTAB %s; \\\n",
			sts->apertab.name);
	    printf ("Warning  APERTURE %s.\n", sts->aperture);
	    sts->x2dcorr_o = OMIT;
	}

	return (0);
}

/* This routine opens the aperture throughput table, finds the columns
   that we need, and gets the total number of rows in the table.
*/

static int OpenThruTab (char *tname, TblInfo *tabinfo) {

	tabinfo->tp = c_tbtopn (tname, IRAF_READ_ONLY, 0);
	if (c_iraferr()) {
	    printf ("ERROR    APERTAB `%s' not found.\n", tname);
	    return (OPEN_FAILED);
	}

	tabinfo->nrows = c_tbpsta (tabinfo->tp, TBL_NROWS);

	/* Find the columns. */
	c_tbcfnd1 (tabinfo->tp, "APERTURE", &tabinfo->cp_aperture);
	c_tbcfnd1 (tabinfo->tp, "NELEM", &tabinfo->cp_nelem);
	c_tbcfnd1 (tabinfo->tp, "WAVELENGTH", &tabinfo->cp_wl);
	c_tbcfnd1 (tabinfo->tp, "THROUGHPUT", &tabinfo->cp_thr);
	if (tabinfo->cp_aperture == 0 ||
	    tabinfo->cp_nelem == 0 ||
	    tabinfo->cp_wl == 0 ||
	    tabinfo->cp_thr == 0) {
	    printf ("ERROR    Column not found in APERTAB.\n");
	    c_tbtclo (tabinfo->tp);
	    return (COLUMN_NOT_FOUND);
	}

	/* Pedigree and descrip are optional columns. */
	c_tbcfnd1 (tabinfo->tp, "PEDIGREE", &tabinfo->cp_pedigree);
	c_tbcfnd1 (tabinfo->tp, "DESCRIP", &tabinfo->cp_descrip);

	return (0);
}

/* This routine reads the column (aperture name) used to select the
   correct row.
*/

static int ReadThruTab (TblInfo *tabinfo, int row, TblRow *tabrow) {

	c_tbegtt (tabinfo->tp, tabinfo->cp_aperture, row,
			tabrow->aperture, STIS_CBUF);
	if (c_iraferr())
	    return (TABLE_ERROR);

	return (0);
}

/* This routine reads the array data from one row.  The number of elements
   in the arrays is gotten, the arrays are allocated, and the wavelengths
   and throughputs are read into the arrays.
*/

static int ReadThruArray (TblInfo *tabinfo, int row, ApInfo *slit) {

	int nwl, nthru;		/* number of elements actually read */

	/* Find out how many elements there are in the throughput arrays,
	   and allocate space for the arrays to be read from the table.
	*/
	c_tbegti (tabinfo->tp, tabinfo->cp_nelem, row, &slit->nelem);
	if (c_iraferr())
	    return (TABLE_ERROR);

	/* Allocate memory. */
	slit->wl = calloc (slit->nelem, sizeof(double));
	slit->thr = calloc (slit->nelem, sizeof(double));
	if (slit->wl == NULL || slit->thr == NULL) {
	    CloseThruTab (tabinfo);
	    return (OUT_OF_MEMORY);
	}

	nwl = c_tbagtd (tabinfo->tp, tabinfo->cp_wl, row,
			slit->wl, 1, slit->nelem);
	if (c_iraferr())
	    return (TABLE_ERROR);

	nthru = c_tbagtd (tabinfo->tp, tabinfo->cp_thr, row,
			slit->thr, 1, slit->nelem);
	if (c_iraferr())
	    return (TABLE_ERROR);

	if (nwl < slit->nelem || nthru < slit->nelem) {
	    c_tbtclo (tabinfo->tp);
	    free (slit->wl);
	    free (slit->thr);
	    printf ("ERROR    not all coefficients were read from APERTAB.\n");
	    return (TABLE_ERROR);
	}

	slit->allocated = 1;			/* set flag */

	return (0);
}

/* This routine closes the APERTAB table. */

static int CloseThruTab (TblInfo *tabinfo) {

	c_tbtclo (tabinfo->tp);
	if (c_iraferr())
	    return (TABLE_ERROR);

	return (0);
}
