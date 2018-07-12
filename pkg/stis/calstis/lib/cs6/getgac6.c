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
					/* column descriptors */
	IRAFPointer cp_opt_elem;	/* mirror or grating */
	IRAFPointer cp_cenwave;
	IRAFPointer cp_aperture;
	IRAFPointer cp_nelem;
	IRAFPointer cp_wl;
	IRAFPointer cp_thr;
	IRAFPointer cp_pedigree;
	IRAFPointer cp_descrip;
	int nrows;			/* number of rows in table */
} TblInfo;

typedef struct {
	char opt_elem[STIS_CBUF];
	char aperture[STIS_CBUF];
	int cenwave;
} TblRow;


static int OpenGacTab (StisInfo6 *, TblInfo *);
static int ReadGacTab (TblInfo *, int, TblRow *);
static int ReadGacData (TblInfo *, int, ApInfo *);
static int CloseGacTab (TblInfo *);


/* This routine gets the absolute flux conversion from GACTAB and saves
   the info in the aperture information structure.

   The grating-aperture correction table should contain the following:
	header parameters:
		none used
	columns:
		OPT_ELEM:  grating name (string)
		CENWAVE:  central wavelength (int)
		APERTURE:  aperture name (string)
		NELEM:  number of elements in arrays (int)
		WAVELENGTH:  array of wavelengths (double)
		THROUGHPUT:  array of correction factors (double)

   The table is read to find the row for which the values of OPT_ELEM,
   CENWAVE and APERTURE are the same as in the input image header.
   For that row, the number of elements NELEM is read, memory is
   allocated in the slit structure for the arrays of wavelength and
   throughput, and those arrays are read.

   When done, the memory should be freed by calling FreeThroughput6.

   Revision history:
   ----------------
   08 Apr 2005  -  Initial version (PEH)
*/

int GetGAC6 (StisInfo6 *sts, ApInfo *slit) {

/* arguments:
StisInfo *sts   i: calibration switches and info
ApInfo *slit    o: includes grating-aperture correction factors
*/

	int status;

	TblInfo tabinfo;	/* pointer to table descriptor, etc */
	TblRow tabrow;		/* values read from a table row */

	int row;		/* loop index */
	int foundit = 0;	/* true if parameters found in table */

	/* Open the grating-aperture correction table. */
	if ((status = OpenGacTab (sts, &tabinfo)))
	    return (status);

	/* Check each row for a match with keyword values, then read
	   the arrays of wavelength and throughput if there's a match.
	*/

	for (row = 1;  row <= tabinfo.nrows;  row++) {

	    if ((status = ReadGacTab (&tabinfo, row, &tabrow)))
		return (status);

	    if (SameString (tabrow.opt_elem, sts->opt_elem) &&
		SameInt (tabrow.cenwave, sts->cenwave) &&
                SameString (tabrow.aperture, sts->aperture)) {

		foundit = 1;

		/* Get pedigree & descrip from the row. */
		if ((status = RowPedigree (&sts->gactab, row,
                        tabinfo.tp, tabinfo.cp_pedigree, tabinfo.cp_descrip))) {
		    return (status);
	        }

		if (sts->gactab.goodPedigree == DUMMY_PEDIGREE) {
		    sts->gaccorr = OMIT;
		    CloseGacTab (&tabinfo);
		    return (0);
		}

		/* Read wavelengths and throughputs into slit. */
		if ((status = ReadGacData (&tabinfo, row, slit)))
		    return (status);

		break;
	    }
	}

	if ((status = CloseGacTab (&tabinfo)))
	    return (status);

	if (!foundit) {
	    /* silently turn off the grating-aperture correction */
	    sts->gaccorr = OMIT;
	    slit->gac_nelem = 0;
	}

	return (0);
}



/* This routine opens the grating-aperture correction table, finds the
   columns that we need, and gets the total number of rows in the table.
*/

static int OpenGacTab (StisInfo6 *sts, TblInfo *tabinfo) {

	tabinfo->tp = c_tbtopn (sts->gactab.name, IRAF_READ_ONLY, 0);

	if (c_iraferr()) {
	    printf ("ERROR    GACTAB `%s' not found\n", sts->gactab.name);
	    return (OPEN_FAILED);
	}

	tabinfo->nrows = c_tbpsta (tabinfo->tp, TBL_NROWS);

	/* Find the mandatory columns. */
	c_tbcfnd1 (tabinfo->tp, "OPT_ELEM",   &tabinfo->cp_opt_elem);
	c_tbcfnd1 (tabinfo->tp, "CENWAVE", &tabinfo->cp_cenwave);
	c_tbcfnd1 (tabinfo->tp, "APERTURE", &tabinfo->cp_aperture);
	c_tbcfnd1 (tabinfo->tp, "NELEM",      &tabinfo->cp_nelem);
	c_tbcfnd1 (tabinfo->tp, "WAVELENGTH", &tabinfo->cp_wl);
	c_tbcfnd1 (tabinfo->tp, "THROUGHPUT", &tabinfo->cp_thr);
	if (tabinfo->cp_opt_elem == 0 ||
	    tabinfo->cp_cenwave  == 0 ||
	    tabinfo->cp_aperture == 0 ||
	    tabinfo->cp_nelem    == 0 ||
	    tabinfo->cp_wl       == 0 ||
	    tabinfo->cp_thr      == 0) {
	    printf ("ERROR    Column not found in GACTAB\n");
	    c_tbtclo (tabinfo->tp);
	    return (COLUMN_NOT_FOUND);
	}

	/* Pedigree and descrip are optional columns. */
	c_tbcfnd1 (tabinfo->tp, "PEDIGREE", &tabinfo->cp_pedigree);
	c_tbcfnd1 (tabinfo->tp, "DESCRIP", &tabinfo->cp_descrip);

	return (0);
}



/* This routine reads the columns used to select the correct row.
   The grating name, central wavelength and aperture are gotten.
*/

static int ReadGacTab (TblInfo *tabinfo, int row, TblRow *tabrow) {

	c_tbegtt (tabinfo->tp, tabinfo->cp_opt_elem, row,
			tabrow->opt_elem, STIS_CBUF-1);
	if (c_iraferr())
	    return (TABLE_ERROR);

	c_tbegti (tabinfo->tp, tabinfo->cp_cenwave, row, &tabrow->cenwave);
	if (c_iraferr())
	    return (TABLE_ERROR);

	c_tbegtt (tabinfo->tp, tabinfo->cp_aperture, row,
			tabrow->aperture, STIS_CBUF-1);
	if (c_iraferr())
	    return (TABLE_ERROR);

	return (0);
}



/* This routine reads all relevant data from one row.  The number of elements
   in the arrays is gotten, the arrays are allocated, and the wavelengths
   and throughput factors are read into the arrays.
*/

static int ReadGacData (TblInfo *tabinfo, int row, ApInfo *slit) {

	int nwl, nthr;		/* number of elements actually read */

	/* Find out how many elements there are in the throughput arrays,
	   and allocate space for the arrays to be read from the table.
	*/
	c_tbegti (tabinfo->tp, tabinfo->cp_nelem, row, &slit->gac_nelem);
	if (c_iraferr())
	    return (TABLE_ERROR);

	slit->gac_wl = (double *) malloc (slit->gac_nelem * sizeof(double));
	slit->gac_thr = (double *) malloc (slit->gac_nelem * sizeof(double));
	if (slit->gac_wl == NULL || slit->gac_thr == NULL)
	    return (OUT_OF_MEMORY);

	nwl = c_tbagtd (tabinfo->tp, tabinfo->cp_wl, row,
	                slit->gac_wl, 1, slit->gac_nelem);
	if (c_iraferr())
	    return (TABLE_ERROR);

	nthr = c_tbagtd (tabinfo->tp, tabinfo->cp_thr, row,
	                  slit->gac_thr, 1, slit->gac_nelem);
	if (c_iraferr())
	    return (TABLE_ERROR);

	if (nwl < slit->gac_nelem || nthr < slit->gac_nelem) {
	    c_tbtclo (tabinfo->tp);
	    free (slit->gac_wl);
	    free (slit->gac_thr);
	    printf ("ERROR    Not all coefficients were read from GACTAB\n");
	    return (TABLE_ERROR);
	}

	/* set flag to indicate that memory has been allocated */
	slit->gac_allocated = 1;

	return (0);
}



/* This routine closes the gactab table. */

static int CloseGacTab (TblInfo *tabinfo) {

	c_tbtclo (tabinfo->tp);
	if (c_iraferr())
	    return (TABLE_ERROR);

	return (0);
}
