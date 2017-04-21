# include <stdio.h>
# include <stdlib.h>	/* malloc */
# include <string.h>	/* strcpy */

# include "c_iraf.h"
# include "hstio.h"
# include "xtables.h"

# include "stis.h"
# include "calstis4.h"
# include "hstcalerr.h"
# include "stiswild.h"	/* for STRING_WILDCARD */
# include "stisdef.h"

typedef struct {
	IRAFPointer tp;			/* pointer to table descriptor */
				/* column descriptors */
	IRAFPointer cp_sclamp;		/* name of lamp */
	IRAFPointer cp_lampset;		/* current through lamp */
	IRAFPointer cp_opt_elem;	/* grating or prism name */
	IRAFPointer cp_nelem;
	IRAFPointer cp_wl;
	IRAFPointer cp_flux;
	IRAFPointer cp_pedigree;
	IRAFPointer cp_descrip;
	int nrows;			/* number of rows in table */
} TblInfo;

typedef struct {
	char sclamp[STIS_CBUF+1];	/* name of lamp */
	char lampset[STIS_CBUF+1];	/* current through lamp */
	char opt_elem[STIS_CBUF+1];	/* grating or prism name */
} TblRow;

static int OpenLampTab (char *, TblInfo *);
static int ReadLampTab (TblInfo *, int, TblRow *);
static int ReadLampArray (TblInfo *, int, LampInfo *);
static int CloseLampTab (TblInfo *);

/* This routine gets the calibration lamp spectrum from LAMPTAB.

   The spectral calibration lamp table should contain the following:
	header parameters:
		none needed
	columns:
		SCLAMP:  name of lamp (string)
		LAMPSET:  current through lamp (string)
		OPT_ELEM:  grating or prism name (string)
		NELEM:  actual number of elements in arrays (int)
		WAVELENGTH:  array of wavelengths in Angstroms (double)
		FLUX:  array of intensities (float)

   The table is read to find the row for which the value of LAMPSET is
   the same as in the input image header.  If the input table contains
   the SCLAMP and OPT_ELEM columns (they may be absent, for backward
   compatibility), those columns will also be used to select the row.
   For the selected row, the number of elements NELEM is read, memory is
   allocated, and the arrays of wavelength and flux are read into the
   allocated memory.  These arrays in the table are in columns that
   contain fixed-length arrays, with NELEM values of data followed
   by fill values to the fixed length.

   After reading the values from the table, the wavelengths will be modified
   to give the wavelengths at the midpoints between values in the flux
   array.  flux[i] is the flux for wavelengths from wl[i] to wl[i+1].
   Thus, the final wavelength array will have one more value than NELEM.

   When done, the memory should be freed by calling FreeLampSpec.

   Phil Hodge, 1997 Dec 4:
	Include SCLAMP as a selection criterion.

   Phil Hodge, 1998 Sept 15:
	Initialize foundit; include lamp name in error message if no
	matching row was found.

   Phil Hodge, 1998 Oct 5:
	Change status value 1010 to GENERIC_ERROR_CODE.

   Phil Hodge, 2000 Jan 13:
	Change sclampCurrent to lampset.
	Add one to sclamp and lampset buffer sizes.

   Phil Hodge, 2001 Mar 7:
	Also select on opt_elem, if the column exists.
*/

int GetLamp (StisInfo4 *sts, LampInfo *lamp) {

/* arguments:
StisInfo4 *sts    i: calibration switches and info
LampInfo *lamp    o: spectrum of calibration lamp
*/

	int status;

	TblInfo tabinfo;	/* pointer to table descriptor, etc */
	TblRow tabrow;		/* values read from a table row */

	int row;		/* loop index */
	int foundit;		/* true if lamp name found in table */

	/* Open the LAMPTAB table. */
	if ((status = OpenLampTab (sts->lamptab.name, &tabinfo)))
	    return (status);

	foundit = 0;
	for (row = 1;  row <= tabinfo.nrows;  row++) {

	    if ((status = ReadLampTab (&tabinfo, row, &tabrow)))
		return (status);

	    if (SameString (tabrow.lampset, sts->lampset) &&
		SameString (tabrow.sclamp, sts->sclamp) &&
		SameString (tabrow.opt_elem, sts->opt_elem)) {

		foundit = 1;

		/* Get pedigree & descrip from the row. */
		if ((status = RowPedigree (&sts->lamptab, row,
                        tabinfo.tp, tabinfo.cp_pedigree, tabinfo.cp_descrip)))
		    return (status);
		if (sts->lamptab.goodPedigree == DUMMY_PEDIGREE) {
		    sts->wavecorr = DUMMY;
		    CloseLampTab (&tabinfo);
		    printf ("Warning  LAMPTAB has PEDIGREE = DUMMY.\n");
		    return (NOTHING_TO_DO);
		}

		/* Read wavelengths and throughputs into lamp structure. */
		if ((status = ReadLampArray (&tabinfo, row, lamp)))
		    return (status);

		break;
	    }
	}

	if ((status = CloseLampTab (&tabinfo)))
	    return (status);

	if (!foundit) {
	    printf ("ERROR    LAMP %s, LAMPSET %s not found in %s.\n",
		sts->sclamp, sts->lampset, sts->lamptab.name);
	    return (GENERIC_ERROR_CODE);
	}

	return (0);
}

/* This routine opens the template lamp spectrum table, finds the columns
   that we need, and gets the total number of rows in the table.
*/

static int OpenLampTab (char *tname, TblInfo *tabinfo) {

	tabinfo->tp = c_tbtopn (tname, IRAF_READ_ONLY, 0);
	if (c_iraferr()) {
	    printf ("ERROR    LAMPTAB `%s' not found.\n", tname);
	    return (OPEN_FAILED);
	}

	tabinfo->nrows = c_tbpsta (tabinfo->tp, TBL_NROWS);

	/* Find the columns. */
	c_tbcfnd1 (tabinfo->tp, "LAMPSET", &tabinfo->cp_lampset);
	c_tbcfnd1 (tabinfo->tp, "NELEM", &tabinfo->cp_nelem);
	c_tbcfnd1 (tabinfo->tp, "WAVELENGTH", &tabinfo->cp_wl);
	c_tbcfnd1 (tabinfo->tp, "FLUX", &tabinfo->cp_flux);

	if (tabinfo->cp_lampset == 0 ||
	    tabinfo->cp_nelem == 0 ||
	    tabinfo->cp_wl == 0 || tabinfo->cp_flux == 0) {
	    printf ("ERROR    Column not found in LAMPTAB.\n");
	    c_tbtclo (tabinfo->tp);
	    return (COLUMN_NOT_FOUND);
	}

	/* Pedigree and descrip are optional columns. */
	c_tbcfnd1 (tabinfo->tp, "PEDIGREE", &tabinfo->cp_pedigree);
	c_tbcfnd1 (tabinfo->tp, "DESCRIP", &tabinfo->cp_descrip);

	/* SCLAMP and OPT_ELEM are optional, for backward compatibility. */
	c_tbcfnd1 (tabinfo->tp, "SCLAMP", &tabinfo->cp_sclamp);
	c_tbcfnd1 (tabinfo->tp, "OPT_ELEM", &tabinfo->cp_opt_elem);

	return (0);
}

/* This routine reads the columns used to select the correct row.
   These columns are LAMPSET and optionally SCLAMP and OPT_ELEM.
   If SCLAMP and/or OPT_ELEM are not found in the input table (they
   were not present in earlier versions), "ANY" will be copied to
   tabrow->sclamp and/or tabrow->opt_elem to effectively disable
   selection on those columns.
*/

static int ReadLampTab (TblInfo *tabinfo, int row, TblRow *tabrow) {

	c_tbegtt (tabinfo->tp, tabinfo->cp_lampset, row,
			tabrow->lampset, STIS_CBUF);
	if (c_iraferr())
	    return (TABLE_ERROR);

	if (tabinfo->cp_sclamp == 0) {
	    strcpy (tabrow->sclamp, STRING_WILDCARD);
	} else {
	    c_tbegtt (tabinfo->tp, tabinfo->cp_sclamp, row,
			tabrow->sclamp, STIS_CBUF);
	    if (c_iraferr())
		return (TABLE_ERROR);
	}

	if (tabinfo->cp_opt_elem == 0) {
	    strcpy (tabrow->opt_elem, STRING_WILDCARD);
	} else {
	    c_tbegtt (tabinfo->tp, tabinfo->cp_opt_elem, row,
			tabrow->opt_elem, STIS_CBUF);
	    if (c_iraferr())
		return (TABLE_ERROR);
	}

	return (0);
}

/* Read the info into memory. */

static int ReadLampArray (TblInfo *tabinfo, int row, LampInfo *lamp) {

	double wsave, wlast, wnext;
	int nwl, nflux;		/* actual number of elements read */
	int j;

	c_tbegti (tabinfo->tp, tabinfo->cp_nelem, row, &lamp->nelem);
	if (c_iraferr())
	    return (TABLE_ERROR);

	lamp->wl = malloc ((1 + lamp->nelem) * sizeof(double));
	lamp->flux = malloc (lamp->nelem * sizeof(double));
	if (lamp->wl == NULL || lamp->flux == NULL)
	    return (OUT_OF_MEMORY);
	lamp->allocated = 1;			/* set flag */

	nwl = c_tbagtd (tabinfo->tp, tabinfo->cp_wl, row, lamp->wl,
			1, lamp->nelem);
	nflux = c_tbagtd (tabinfo->tp, tabinfo->cp_flux, row, lamp->flux,
			1, lamp->nelem);
	if (c_iraferr())
	    return (TABLE_ERROR);

	if (nwl < lamp->nelem || nflux < lamp->nelem) {
	    c_tbtclo (tabinfo->tp);
	    printf ("ERROR    Not all elements were read from LAMPTAB.\n");
	    return (TABLE_ERROR);
	}

	/* Change the wavelength values to be the wavelengths at the
	   midpoints between elements.
	*/
	wsave = lamp->wl[0] - (lamp->wl[1] - lamp->wl[0]) / 2.;
	wlast = lamp->wl[nwl-1] + (lamp->wl[nwl-1] - lamp->wl[nwl-2]) / 2.;
	for (j = 1;  j < nwl;  j++) {
	    wnext = (lamp->wl[j-1] + lamp->wl[j]) / 2.;
	    lamp->wl[j-1] = wsave;	/* left edge of pixel j-1 */
	    wsave = wnext;		/* left edge of pixel j */
	}
	lamp->wl[nwl-1] = wnext;	/* left edge of last pixel */
	lamp->wl[nwl] = wlast;		/* right edge of last pixel */

	return (0);
}

/* This routine closes the lamptab table. */

static int CloseLampTab (TblInfo *tabinfo) {

	c_tbtclo (tabinfo->tp);
	if (c_iraferr())
	    return (TABLE_ERROR);

	return (0);
}
