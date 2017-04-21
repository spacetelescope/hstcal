# include <stdio.h>
# include <stdlib.h>	/* malloc, abs */
# include <string.h>	/* strcpy */

# include "xtables.h"
# include "hstio.h"

# include "stis.h"
# include "stisdef.h"
# include "calstis6.h"
# include "stispht.h"
# include "err.h"

typedef struct {
	IRAFPointer tp;			/* pointer to table descriptor */
					/* column descriptors */
	IRAFPointer cp_aperture;
	IRAFPointer cp_cenwave;
	IRAFPointer cp_extrheight;
	IRAFPointer cp_nelem;
	IRAFPointer cp_wl;
	IRAFPointer cp_pc;
	IRAFPointer cp_pedigree;
	IRAFPointer cp_descrip;
	int maxhght;			/* extrheight to use for "infinity" */
	int nrows;			/* number of rows in table */
} TblInfo;

typedef struct {
	char aperture[STIS_CBUF];
	int cenwave;
	int extrheight;
} TblRow;

static int OpenPCTab (char *, TblInfo *);
static int ReadPCTab (TblInfo *, int, TblRow *);
static int ReadPCArray (TblInfo *, int, PhotInfo *);
static int ClosePCTab (TblInfo *);
static int PCDummy (PhotInfo *);

/* This routine gets the PCT correction to the absolute sensitivity,
   the factor (as a function of wavelength) to correct to an "infinite"
   extraction aperture from the default extraction aperture.
   the info in the photometry information structure.

   The PCT table should contain the following:
	header parameters:
		MAXHGHT:  value of EXTRHEIGHT for "infinite" height
	columns:
		APERTURE:  aperture name (string)
		CENWAVE:  central wavelength (int)
		EXTRHEIGHT:  height of spectrum extraction box (int)
		NELEM:  actual number of elements in arrays (int)
		WAVELENGTH:  array of wavelengths (double)
		THROUGHPUT:  array of factors (float)


   The table is read to find the row for which both CENWAVE and APERTURE
   are the same as in the input image header. The third row matching
   criterion is based on the 'extrsize' argument. If set to zero, then the
   row is selected that has EXTRHEIGHT equal to the table header keyword
   MAXHGHT. If set to anything larger than zero, then the row is selected
   that has EXTRHEIGHT closest to an entry in the table. In this case a
   warning message is issued.

   For the selected row, the number of elements NELEM is read, and the
   arrays of wavelength and correction factor are read.

   These are only coarsely sampled values; spline interpolation is then
   used to obtain a correction factor at each of the wavelengths in the
   phot->wl array.  Note that this means that GetAbsPhot must have been
   called before this function is called.

   It is not a fatal error for the PCTAB table to not exist, as long as
   this is flagged by the name of the PCTAB name being null.  In this case,
   a dummy phot->pcorr array will be allocated, and the values will be set
   to one.

   Memory allocated by this routine will be freed along with the
   photometry info when FreePhot is called.

   Phil Hodge, 1997 Nov 13:
	Function created.




   Revision history:
   ----------------
   26 Jan 98  -  Borrowed from calstis7 (I.Busko)
   02 Feb 98  -  Read both the maximum and the current height data (IB)

*/

int GetPCT6 (StisInfo6 *sts, PhotInfo *phot, double extrsize, int warn) {

/* arguments:
StisInfo6 *sts  i: calibration switches and info
PhotInfo *phot  o: phot->pcorr is the array of factors to correct
                   absolute flux for finite extraction box size
double extrsize i: the extraction height size.
int warn	i: warning message control
*/

	int status;

	TblInfo tabinfo;	/* pointer to table descriptor, etc */
	TblRow tabrow;		/* values read from a table row */

	int row;		/* loop index */
	int foundit = -1;	/* status of found entry */
	int mheight;		/* the height to look for */
	int mindist;		/* minimum distance in set of matching rows */

	/* Check that the photometry info has in fact been gotten;
	   it's phot->wl that we need.
	*/
	if (!phot->allocated)
	    return (INTERNAL_ERROR);

	if (sts->pctcorr != PERFORM) {
	    /* Allocate memory for phot->pcorr and set the values to 1. */
	    if ((status = PCDummy (phot)))
		return (status);
	    else {
	        if (warn == WARN) {
printf ("Warning  No PCTAB information in input file.\n");
printf ("Warning  No extraction box photometric correction is used to\n");
printf ("Warning  correct the flux-calibrated result.\n");
                }
		return (0);
	    }
	}

	/* Open the PCTAB table. */
	if ((status = OpenPCTab (sts->pctab.name, &tabinfo)))
	    return (status);

	/* Select proper height to look for. */
	mheight = ((int)extrsize > 0) ? (int)extrsize : tabinfo.maxhght;

	/* Check each row for a match with keyword values, then read
	   the arrays of wavelength and throughput if there's a match.
	*/
	mindist = 2000;
	for (row = 1;  row <= tabinfo.nrows;  row++) {

	    if ((status = ReadPCTab (&tabinfo, row, &tabrow)))
		return (status);

	    if (SameString (tabrow.aperture, sts->aperture) &&
		SameInt (tabrow.cenwave, sts->cenwave) &&
		SameInt (tabrow.extrheight, mheight)) {

		foundit = 0;  /* this is an exact match */

		/* Get pedigree & descrip from the row. */
		if ((status = RowPedigree (&sts->pctab, row,
                        tabinfo.tp, tabinfo.cp_pedigree, tabinfo.cp_descrip)))
		    return (status);
		if (sts->pctab.goodPedigree == DUMMY_PEDIGREE) {
		    printf ("Warning  DUMMY pedigree in row %d of %s.\n",
		            row, sts->pctab.name);
		    sts->pctcorr = DUMMY;
		    if ((status = PCDummy (phot)))
			return (status);
		    ClosePCTab (&tabinfo);
		    return (0);
		}

		/* Read wavelengths and throughputs into phot. */
		if ((status = ReadPCArray (&tabinfo, row, phot)))
		    return (status);

		break;

	    /* This handles the case of an approximate match. */

	    } else if (SameString (tabrow.aperture, sts->aperture) &&
	               SameInt (tabrow.cenwave, sts->cenwave)) {
	        if (mindist > abs (tabrow.extrheight - mheight)) {
	            mindist = abs (tabrow.extrheight - mheight);
		    if ((status = ReadPCArray (&tabinfo, row, phot)))
		        return (status);
	            foundit = tabrow.extrheight;
	        }
	    }
	}

	if ((status = ClosePCTab (&tabinfo)))
	    return (status);

	if (foundit == -1) {
	    printf ("Warning  No appropriate row found in PCTAB %s; \\\n",
			sts->pctab.name);
	    printf ("Warning  APERTURE %s, CENWAVE %d, EXTRHEIGHT %d.\n",
		sts->aperture, sts->cenwave, tabinfo.maxhght);
	    sts->pctcorr = OMIT;
	    if ((status = PCDummy (phot)))
		return (status);
	} else if (foundit > 0) {
	    printf
            ("Warning  No exact extraction box height match in PCTAB.\n");
	    printf
            ("Warning  Using entry for height %d\n", foundit);
	}

	return (0);
}

/* This routine opens the throughput table, finds the columns that we
   need, and gets the total number of rows in the table.
*/

static int OpenPCTab (char *tname, TblInfo *tabinfo) {

	tabinfo->tp = c_tbtopn (tname, IRAF_READ_ONLY, 0);
	if (c_iraferr()) {
	    printf ("ERROR    PCTAB `%s' not found.\n", tname);
	    return (OPEN_FAILED);
	}

	tabinfo->nrows = c_tbpsta (tabinfo->tp, TBL_NROWS);

	tabinfo->maxhght = c_tbhgti (tabinfo->tp, "MAXHGHT");

	/* Find the columns. */
	c_tbcfnd1 (tabinfo->tp, "APERTURE", &tabinfo->cp_aperture);
	c_tbcfnd1 (tabinfo->tp, "CENWAVE", &tabinfo->cp_cenwave);
	c_tbcfnd1 (tabinfo->tp, "EXTRHEIGHT", &tabinfo->cp_extrheight);
	c_tbcfnd1 (tabinfo->tp, "NELEM", &tabinfo->cp_nelem);
	c_tbcfnd1 (tabinfo->tp, "WAVELENGTH", &tabinfo->cp_wl);
	c_tbcfnd1 (tabinfo->tp, "THROUGHPUT", &tabinfo->cp_pc);
	if (tabinfo->cp_aperture == 0 ||
	    tabinfo->cp_cenwave == 0 ||
	    tabinfo->cp_extrheight == 0 ||
	    tabinfo->cp_nelem == 0 ||
	    tabinfo->cp_wl == 0 ||
	    tabinfo->cp_pc == 0) {
	    printf ("ERROR    Column not found in PCTAB.\n");
	    c_tbtclo (tabinfo->tp);
	    return (COLUMN_NOT_FOUND);
	}

	/* Pedigree and descrip are optional columns. */
	c_tbcfnd1 (tabinfo->tp, "PEDIGREE", &tabinfo->cp_pedigree);
	c_tbcfnd1 (tabinfo->tp, "DESCRIP", &tabinfo->cp_descrip);

	return (0);
}

/* This routine reads the columns used to select the correct row.
   The grating name, aperture, central wavelength, and extraction box
   height are gotten.
*/

static int ReadPCTab (TblInfo *tabinfo, int row, TblRow *tabrow) {

	c_tbegtt (tabinfo->tp, tabinfo->cp_aperture, row,
			tabrow->aperture, STIS_CBUF-1);
	if (c_iraferr())
	    return (TABLE_ERROR);

	c_tbegti (tabinfo->tp, tabinfo->cp_cenwave, row, &tabrow->cenwave);
	if (c_iraferr())
	    return (TABLE_ERROR);
	c_tbegti (tabinfo->tp, tabinfo->cp_extrheight, row,
			&tabrow->extrheight);
	if (c_iraferr())
	    return (TABLE_ERROR);

	return (0);
}

/* This routine reads the array data from one row.  The number of elements
   in the arrays is gotten, the arrays are allocated, and the wavelengths
   and throughputs are read into the arrays.
*/

static int ReadPCArray (TblInfo *tabinfo, int row, PhotInfo *phot) {

	int status;
	double *wl, *pc;	/* for values read from table */
	int nelem;		/* number of elements we expect to read */
	int nret_wl, nret_pc;	/* number of elements actually read */

	/* Find out how many elements there are in the table row. */
	c_tbegti (tabinfo->tp, tabinfo->cp_nelem, row, &nelem);
	if (c_iraferr())
	    return (TABLE_ERROR);

	/* Allocate space for the arrays to be read from the table.
	   Note that nelem will likely be much smaller than phot->nelem.
	*/
	pc = malloc (nelem * sizeof(double));
	wl = malloc (nelem * sizeof(double));
	if (pc == NULL || wl == NULL)
	    return (OUT_OF_MEMORY);

	/* Read the wavelength and correction factors from the table. */

	nret_wl = c_tbagtd (tabinfo->tp, tabinfo->cp_wl, row,
			wl, 1, nelem);
	if (c_iraferr()) {
	    free (wl);
	    free (pc);
	    return (TABLE_ERROR);
	}

	nret_pc = c_tbagtd (tabinfo->tp, tabinfo->cp_pc, row,
			pc, 1, nelem);
	if (c_iraferr()) {
	    free (wl);
	    free (pc);
	    return (TABLE_ERROR);
	}

	if (nret_wl < nelem || nret_pc < nelem) {
	    c_tbtclo (tabinfo->tp);
	    printf ("ERROR    Not all coefficients were read from PCTAB.\n");
	    free (wl);
	    free (pc);
	    return (TABLE_ERROR);
	}

	/* Allocate space for the interpolated correction factors.
           If the array was allocated in a previous execution of this
           function, de-allocate it first. This can happen when multiple
           table rows satisfy the approximate distance criterion.
        */
	if (phot->pcorr != NULL)
	    free (phot->pcorr);
	phot->pcorr = malloc (phot->nelem * sizeof(double));
	if (phot->pcorr == NULL) {
	    free (wl);
	    free (pc);
	    return (OUT_OF_MEMORY);
	}

	/* Interpolate. */
	if ((status = splint_nr (wl, pc, nelem,
                                 phot->wl, phot->pcorr, phot->nelem))) {
	    free (wl);
	    free (pc);
	    return (status);
	}

	free (wl);
	free (pc);

	return (0);
}

/* This routine closes the PCTAB table. */

static int ClosePCTab (TblInfo *tabinfo) {

	c_tbtclo (tabinfo->tp);
	if (c_iraferr())
	    return (TABLE_ERROR);

	return (0);
}

/* This routine allocates the pcorr array and fills it with ones. */

static int PCDummy (PhotInfo *phot) {

	int i;

	phot->pcorr = malloc (phot->nelem * sizeof(double));
	if (phot->pcorr == NULL)
	    return (OUT_OF_MEMORY);

	for (i = 0;  i < phot->nelem;  phot->pcorr[i++] = 1.);

	return (0);
}
