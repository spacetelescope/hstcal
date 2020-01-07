# include <stdio.h>
# include <stdlib.h>	/* malloc */
# include <string.h>	/* strcpy */

# include "c_iraf.h"
# include "hstio.h"
# include "xtables.h"

# include "stis.h"
# include "calstis7.h"
# include "hstcalerr.h"
# include "stisdef.h"
# include "stispht.h"

typedef struct {
	IRAFPointer tp;			/* pointer to table descriptor */
					/* column descriptors */
	IRAFPointer cp_opt_elem;	/* mirror or grating */
	IRAFPointer cp_cenwave;
	IRAFPointer cp_sporder;
	IRAFPointer cp_nelem;
	IRAFPointer cp_wl;
	IRAFPointer cp_thru;
	IRAFPointer cp_mref;
	IRAFPointer cp_wref;
	IRAFPointer cp_yref;
	IRAFPointer cp_mjd;
	IRAFPointer cp_mx;
	IRAFPointer cp_my;
	IRAFPointer cp_mt;
	IRAFPointer cp_m0;
	IRAFPointer cp_pedigree;
	IRAFPointer cp_descrip;
	int nrows;			/* number of rows in table */
} TblInfo;

typedef struct {
	char opt_elem[STIS_CBUF+1];
	int cenwave;
	int sporder;
} TblRow;

static int OpenPhotTab (StisInfo7 *, TblInfo *, PhotInfo *, int *);
static int ReadPhotTab (TblInfo *, int, TblRow *);
static int ReadPhotData (TblInfo *, int, PhotInfo *);
static int ClosePhotTab (TblInfo *);

/* This routine gets the absolute flux conversion from PHOTTAB and saves
   the info in the photometry information structure.

   The absolute-flux table should contain the following:
	header parameters:
		none used
	columns:
		OPT_ELEM:  grating name (string)
		CENWAVE:  central wavelength (int)
		SPORDER:  spectral order (int)
		NELEM:  actual number of elements in arrays (int)
		WAVELENGTH:  array of wavelengths (double)
		THROUGHPUT:  array of throughputs (float)
		REFORD: reference order # of sens. function (int)
		REFWAV: central wavelength of mref order (double)
		REFY:  A2 position of mref order (double)
		REFMJD: observation date of blaze function (double)
		BSHIFT_VS_X: x shift scale factor (double
		BSHIFT_VS_Y: y shift scale factor (double)
		BSHIFT_VS_T: time shift scale factor (double)
		BSHIFT_OFFSET: zero point offset (double)

   The table is read to find the row for which the value of OPT_ELEM
   is the same as in the input image header.  For spectroscopic data,
   CENWAVE is also compared with the header, and SPORDER must be the
   current spectral order being processed.  For that row, the number of
   elements NELEM is read, memory is allocated in the phot structure for
   the arrays of wavelength and throughput, and those arrays are read.
   If existent, values associated with the blaze shift correction are
   then read. When done, the memory should be freed by calling freePhot.

   Phil Hodge, 2000 Jan 13:
	Add one to opt_elem buffer size.

   Ivo Busko, 2002 Feb 04
       Blaze shift correction.

   Paul Barrett, 2004 Jul 27
       Remove blaze shift warning for first order mode.

   Phil Hodge, 2005 Feb 15:
	Initialize phot->blazecorr to OMIT.

   Phil Hodge, 2006 June 27:
	Also get BSHIFT_OFFSET and include in phot as phot->m0.
*/

int GetAbsPhot (StisInfo7 *sts, int sporder, PhotInfo *phot, int print,
                int *warn) {

/* arguments:
StisInfo *sts   i: calibration switches and info
int sporder     i: current order number
PhotInfo *phot  o: factors to convert to absolute flux
int print	i: if set to zero, turn off error message printing
int *warn	io: if set to zero, turn off blaze shift warning
*/

	int status;

	TblInfo tabinfo;	/* pointer to table descriptor, etc */
	TblRow tabrow;		/* values read from a table row */

	int row;		/* loop index */
	int foundit = 0;	/* true if parameters found in table */

	/* Assign the units for the calibrated data. */
	strcpy (phot->bunit, "erg /s /cm**2 /angstrom /arcsec**2");

	/* Open the photometry table. */
	if ((status = OpenPhotTab (sts, &tabinfo, phot, warn)))
	    return (status);

	/* Check each row for a match with keyword values, then read
	   the arrays of wavelength and throughput if there's a match.
	*/

	for (row = 1;  row <= tabinfo.nrows;  row++) {

	    if ((status = ReadPhotTab (&tabinfo, row, &tabrow)))
		return (status);

	    if (sts->obstype == SPECTROSCOPIC_TYPE) {
		if (!SameInt (tabrow.cenwave, sts->cenwave) ||
		    !SameInt (tabrow.sporder, sporder))
		    continue;
	    }

	    if (SameString (tabrow.opt_elem, sts->opt_elem)) {

		foundit = 1;

		/* Get pedigree & descrip from the row. */
		if ((status = RowPedigree (&sts->phottab, row,
                        tabinfo.tp, tabinfo.cp_pedigree, tabinfo.cp_descrip)))
		    return (status);
		if (sts->phottab.goodPedigree == DUMMY_PEDIGREE) {
		    printf ("Warning  DUMMY pedigree in row %d of %s.\n",
			row, sts->phottab.name);
		    sts->x2dcorr_o = DUMMY;
		    ClosePhotTab (&tabinfo);
		    return (0);
		}

		/* Read wavelengths and throughputs into phot. */
		if ((status = ReadPhotData (&tabinfo, row, phot)))
		    return (status);

		break;
	    }
	}

	if ((status = ClosePhotTab (&tabinfo)))
	    return (status);

	if (!foundit) {
	    if (print) {
	        printf ("Warning  Matching row not found in PHOTTAB %s; \\\n",
			 sts->phottab.name);
	        if (sts->obstype == SPECTROSCOPIC_TYPE) {
		    printf ("Warning  OPT_ELEM %s, CENWAVE %d, SPORDER %d.\n",
			    sts->opt_elem, sts->cenwave, sporder);
	        } else {
		    printf ("Warning  OPT_ELEM %s\n", sts->opt_elem);
	        }
	        sts->x2dcorr_o = OMIT;
	    } else
	        return (ROW_NOT_FOUND);
	}

	return (0);
}

/* This routine opens the throughput table, finds the columns that we
   need, and gets the total number of rows in the table.

   It is not an error if the columns associated with the blaze function
   are not found.
*/

static int OpenPhotTab (StisInfo7 *sts, TblInfo *tabinfo, PhotInfo *phot,
			int *warn) {

	tabinfo->tp = c_tbtopn (sts->phottab.name, IRAF_READ_ONLY, 0);
	if (c_iraferr()) {
	    printf ("ERROR    PHOTTAB `%s' not found.\n", sts->phottab.name);
	    return (OPEN_FAILED);
	}

	tabinfo->nrows = c_tbpsta (tabinfo->tp, TBL_NROWS);

	/* Find the mandatory columns. */
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
	if (sts->obstype == SPECTROSCOPIC_TYPE) {
	    c_tbcfnd1 (tabinfo->tp, "CENWAVE", &tabinfo->cp_cenwave);
	    c_tbcfnd1 (tabinfo->tp, "SPORDER", &tabinfo->cp_sporder);
	    if (tabinfo->cp_cenwave == 0 || tabinfo->cp_sporder == 0) {
		printf (
	"ERROR    Column (CENWAVE or SPORDER) not found in PHOTTAB.\n");
		c_tbtclo (tabinfo->tp);
		return (COLUMN_NOT_FOUND);
	    }
	} else {
	    /* Flag the fact that we didn't try to find these columns. */
	    tabinfo->cp_cenwave = 0;
	    tabinfo->cp_sporder = 0;
	}

	/* Now look for the blaze shift columns. It's not an error if
	   they are not found.
	*/
	phot->blazecorr = OMIT;
	if (!(sts->first_order)) {
	    phot->blazecorr = PERFORM;
	    c_tbcfnd1 (tabinfo->tp, "REFORD",      &tabinfo->cp_mref);
	    c_tbcfnd1 (tabinfo->tp, "REFWAV",      &tabinfo->cp_wref);
	    c_tbcfnd1 (tabinfo->tp, "REFY",        &tabinfo->cp_yref);
	    c_tbcfnd1 (tabinfo->tp, "REFMJD",      &tabinfo->cp_mjd);
	    c_tbcfnd1 (tabinfo->tp, "BSHIFT_VS_X", &tabinfo->cp_mx);
	    c_tbcfnd1 (tabinfo->tp, "BSHIFT_VS_Y", &tabinfo->cp_my);
	    c_tbcfnd1 (tabinfo->tp, "BSHIFT_VS_T", &tabinfo->cp_mt);
	    /* bshift_offset is an optional column */
	    c_tbcfnd1 (tabinfo->tp, "BSHIFT_OFFSET", &tabinfo->cp_m0);
	    if (tabinfo->cp_mref == 0 ||
		tabinfo->cp_wref == 0 ||
		tabinfo->cp_yref == 0 ||
		tabinfo->cp_mjd  == 0 ||
		tabinfo->cp_mx   == 0 ||
		tabinfo->cp_my   == 0 ||
		tabinfo->cp_mt   == 0) {

		phot->blazecorr = OMIT;

		if (*warn)
		    printf ("Warning  PHOTTAB does not "
			    "contain blaze shift information.\n");
		*warn = 0;
	    }
        }

	/* Pedigree and descrip are optional columns. */
	c_tbcfnd1 (tabinfo->tp, "PEDIGREE", &tabinfo->cp_pedigree);
	c_tbcfnd1 (tabinfo->tp, "DESCRIP", &tabinfo->cp_descrip);

	return (0);
}

/* This routine reads the columns used to select the correct row.
   The grating name, central wavelength, and spectral order number
   are gotten.
*/

static int ReadPhotTab (TblInfo *tabinfo, int row, TblRow *tabrow) {

	c_tbegtt (tabinfo->tp, tabinfo->cp_opt_elem, row,
			tabrow->opt_elem, STIS_CBUF);
	if (c_iraferr())
	    return (TABLE_ERROR);

	if (tabinfo->cp_cenwave != 0) {
	    c_tbegti (tabinfo->tp, tabinfo->cp_cenwave, row, &tabrow->cenwave);
	    if (c_iraferr())
		return (TABLE_ERROR);
	    c_tbegti (tabinfo->tp, tabinfo->cp_sporder, row, &tabrow->sporder);
	    if (c_iraferr())
		return (TABLE_ERROR);
	} else {
	    tabrow->cenwave = 0;
	    tabrow->sporder = 0;
	}

	return (0);
}

/* This routine reads all relevant data from one row. The number of elements
   in the arrays is gotten, the arrays are allocated, and the wavelengths
   and throughputs are read into the arrays. Then, if blazecorr is enabled,
   the associated parameters are gotten as well
*/

static int ReadPhotData (TblInfo *tabinfo, int row, PhotInfo *phot) {

	int nwl, nthru;		/* number of elements actually read */

	/* Find out how many elements there are in the throughput arrays,
	   and allocate space for the arrays to be read from the table.
	*/
	c_tbegti (tabinfo->tp, tabinfo->cp_nelem, row, &phot->nelem);
	if (c_iraferr())
	    return (TABLE_ERROR);

	phot->wl = (double *) malloc (phot->nelem * sizeof(double));
	phot->thru = (double *) malloc (phot->nelem * sizeof(double));
	if (phot->wl == NULL || phot->thru == NULL)
	    return (OUT_OF_MEMORY);

	nwl = c_tbagtd (tabinfo->tp, tabinfo->cp_wl, row,
			phot->wl, 1, phot->nelem);
	if (c_iraferr())
	    return (TABLE_ERROR);

	nthru = c_tbagtd (tabinfo->tp, tabinfo->cp_thru, row,
			phot->thru, 1, phot->nelem);
	if (c_iraferr())
	    return (TABLE_ERROR);

	if (nwl < phot->nelem || nthru < phot->nelem) {
	    c_tbtclo (tabinfo->tp);
	    free (phot->wl);
	    free (phot->thru);
	    printf ("ERROR    Not all elements were read from PHOTTAB.\n");
	    return (TABLE_ERROR);
	}

	phot->allocated = 1;			/* set flag */

	if (phot->blazecorr == PERFORM) {

	    c_tbegti (tabinfo->tp, tabinfo->cp_mref, row, &phot->mref);
	    c_tbegtd (tabinfo->tp, tabinfo->cp_wref, row, &phot->wref);
	    c_tbegtd (tabinfo->tp, tabinfo->cp_yref, row, &phot->yref);
	    c_tbegtd (tabinfo->tp, tabinfo->cp_mjd,  row, &phot->mjd);
	    c_tbegtd (tabinfo->tp, tabinfo->cp_mx,   row, &phot->mx);
	    c_tbegtd (tabinfo->tp, tabinfo->cp_my,   row, &phot->my);
	    c_tbegtd (tabinfo->tp, tabinfo->cp_mt,   row, &phot->mt);
	    if (tabinfo->cp_m0 == 0)
		phot->m0 = 0.;
	    else
		c_tbegtd (tabinfo->tp, tabinfo->cp_m0, row, &phot->m0);

	    if (c_iraferr()) {
	        free (phot->wl);
	        free (phot->thru);
	        free (phot->error);
	        return (TABLE_ERROR);
	    }

	    /* Reference data is 1-indexed ! */
	    phot->yref -= 1.0;
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
