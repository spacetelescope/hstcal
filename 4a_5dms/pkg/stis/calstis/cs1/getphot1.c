# include <stdio.h>
# include <stdlib.h>	/* malloc */

# include "c_iraf.h"
# include "hstio.h"
# include "xtables.h"

# include "stis.h"
# include "calstis1.h"
# include "stiserr.h"
# include "stisdef.h"

typedef struct {
	IRAFPointer tp;			/* pointer to table descriptor */
					/* column descriptors */
	IRAFPointer cp_opt_elem;	/* opt_elem (mirror or grating) */
	IRAFPointer cp_nelem;
	IRAFPointer cp_wl;
	IRAFPointer cp_thru;
	IRAFPointer cp_pedigree;
	IRAFPointer cp_descrip;
	int nrows;			/* number of rows in table */
} TblInfo;

typedef struct {
	char opt_elem[STIS_CBUF+1];
} TblRow;

static int OpenPhotTab (char *, TblInfo *);
static int ReadPhotTab (TblInfo *, int, TblRow *);
static int ReadPhotArray (TblInfo *, int, PhotInfo *);
static int ClosePhotTab (TblInfo *);

/* This routine gets the absolute flux conversion from PHOTTAB and saves
   the info in the photometry information structure.

   The absolute-flux table should contain the following:
	header parameters:
		none used
	columns:
		OPT_ELEM:  mirror name (string)
		NELEM:  actual number of elements in arrays (int)
		WAVELENGTH:  array of wavelengths (read as float)
		THROUGHPUT:  array of throughputs (float)

   The table is read to find the row for which the value of OPT_ELEM
   is the same as in the input image header.  For that row, the number
   of elements NELEM is read, and the arrays of wavelength and throughput
   are read.  The synphot routine phopar is then called to determine
   the inverse sensitivity, reference magnitude (actually a constant),
   pivot wavelength, and RMS bandwidth.  These are written to keywords
   in the primary header.

   Phil Hodge, 1997 Nov 13:
	Extracted from dophot.c.

   Phil Hodge, 1998 Feb 6:
	Remove the include for "hstio.h".

   Phil Hodge, 1998 Oct 5:
	Change status value 1010 to GENERIC_ERROR_CODE.

   Phil Hodge, 2000 Jan 13:
	Add one to opt_elem buffer size.
*/

int GetPhot1 (StisInfo1 *sts, PhotInfo *phot) {

/* arguments:
StisInfo1 *sts   i: calibration switches, etc
PhotInfo *phot   o: QE throughput values are allocated and assigned
*/

	int status;

	TblInfo tabinfo;	/* pointer to table descriptor, etc */
	TblRow tabrow;		/* values read from a table row */

	int row;		/* loop index */
	int foundit;		/* true if parameters found in table */

	/* Open the photometry table. */
	if (status = OpenPhotTab (sts->phot.name, &tabinfo))
	    return (status);

	/* Check each row for a match with keyword values, then read
	   the arrays of wavelength and throughput if there's a match.
	*/

	foundit = 0;
	for (row = 1;  row <= tabinfo.nrows;  row++) {

	    if (status = ReadPhotTab (&tabinfo, row, &tabrow))
		return (status);

	    if (SameString (tabrow.opt_elem, sts->opt_elem)) {

		foundit = 1;

		/* Get pedigree & descrip from the row. */
		if (status = RowPedigree (&sts->phot, row,
			tabinfo.tp, tabinfo.cp_pedigree, tabinfo.cp_descrip))
		    return (status);
		if (sts->phot.goodPedigree == DUMMY_PEDIGREE) {
		    sts->photcorr = DUMMY;
		    ClosePhotTab (&tabinfo);
		    return (0);
		}

		/* Read wavelengths and throughputs into phot. */
		if (status = ReadPhotArray (&tabinfo, row, phot))
		    return (status);

		break;
	    }
	}

	if (status = ClosePhotTab (&tabinfo))
	    return (status);

	if (!foundit) {
	    printf ("ERROR    OPT_ELEM %s not found in PHOTTAB %s\n",
			sts->opt_elem, sts->phot.name);
	    return (GENERIC_ERROR_CODE);
	}

	return (0);
}

/* This routine opens the throughput table, finds the columns that we
   need, and gets the total number of rows in the table.
*/

static int OpenPhotTab (char *tname, TblInfo *tabinfo) {

	tabinfo->tp = c_tbtopn (tname, IRAF_READ_ONLY, 0);
	if (c_iraferr()) {
	    printf ("ERROR    PHOTTAB `%s' not found.\n", tname);
	    return (OPEN_FAILED);
	}

	tabinfo->nrows = c_tbpsta (tabinfo->tp, TBL_NROWS);

	/* Find the columns. */
	c_tbcfnd1 (tabinfo->tp, "OPT_ELEM", &tabinfo->cp_opt_elem);
	c_tbcfnd1 (tabinfo->tp, "NELEM", &tabinfo->cp_nelem);
	c_tbcfnd1 (tabinfo->tp, "WAVELENGTH", &tabinfo->cp_wl);
	c_tbcfnd1 (tabinfo->tp, "THROUGHPUT", &tabinfo->cp_thru);
	if (tabinfo->cp_opt_elem == 0 ||
	    tabinfo->cp_nelem == 0 ||
	    tabinfo->cp_wl == 0 ||
	    tabinfo->cp_thru == 0) {
	    printf ("ERROR    Column not found in PHOTTAB.\n");
	    c_tbtclo (tabinfo->tp);
	    return (COLUMN_NOT_FOUND);
	}

	/* Pedigree and descrip are optional columns. */
	c_tbcfnd1 (tabinfo->tp, "PEDIGREE", &tabinfo->cp_pedigree);
	c_tbcfnd1 (tabinfo->tp, "DESCRIP", &tabinfo->cp_descrip);

	return (0);
}

/* This routine reads the mirror name, which is used to select the
   correct row.
*/

static int ReadPhotTab (TblInfo *tabinfo, int row, TblRow *tabrow) {

	c_tbegtt (tabinfo->tp, tabinfo->cp_opt_elem, row,
			tabrow->opt_elem, STIS_CBUF);
	if (c_iraferr())
	    return (TABLE_ERROR);

	return (0);
}

/* This routine reads the array data from one row.  The number of elements
   in the arrays is gotten, the array is allocated, and the wavelengths
   and throughputs are read into the arrays.
*/

static int ReadPhotArray (TblInfo *tabinfo, int row, PhotInfo *phot) {

	int nwl, nthru;		/* number of elements actually read */

	/* Find out how many elements there are in the throughput arrays,
	   and allocate space for the arrays to be read from the table.
	*/
	c_tbegti (tabinfo->tp, tabinfo->cp_nelem, row, &phot->p_nelem);
	if (c_iraferr())
	    return (TABLE_ERROR);

	phot->p_wl = calloc (phot->p_nelem, sizeof(float));
	phot->p_thru = calloc (phot->p_nelem, sizeof(float));
	if (phot->p_wl == NULL || phot->p_thru == NULL)
	    return (OUT_OF_MEMORY);

	nwl = c_tbagtr (tabinfo->tp, tabinfo->cp_wl, row,
			phot->p_wl, 1, phot->p_nelem);
	if (c_iraferr())
	    return (TABLE_ERROR);

	nthru = c_tbagtr (tabinfo->tp, tabinfo->cp_thru, row,
			phot->p_thru, 1, phot->p_nelem);
	if (c_iraferr())
	    return (TABLE_ERROR);

	if (nwl < phot->p_nelem || nthru < phot->p_nelem) {
	    c_tbtclo (tabinfo->tp);
	    printf ("ERROR    Not all coefficients were read from PHOTTAB.\n");
	    return (TABLE_ERROR);
	}

	return (0);
}

/* This routine closes the phottab table. */

static int ClosePhotTab (TblInfo *tabinfo) {

	c_tbtclo (tabinfo->tp);
	if (c_iraferr())
	    return (TABLE_ERROR);

	return (0);
}
