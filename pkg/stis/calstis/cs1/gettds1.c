# include <stdio.h>
# include <string.h>	/* strcpy */
# include <stdlib.h>	/* calloc */

# include "c_iraf.h"
# include "xtables.h"
# include "hstio.h"

# include "stis.h"
# include "stisdef.h"
# include "err.h"
# include "stistds.h"


typedef struct {
	IRAFPointer tp;			/* pointer to table descriptor */
	IRAFPointer cp_wl;		/* column descriptors */
	IRAFPointer cp_nwl;
	IRAFPointer cp_optelem;
	IRAFPointer cp_pedigree;
	IRAFPointer cp_descrip;
	IRAFPointer cp_reftemp;
	IRAFPointer cp_tempsens;
	int nrows;			/* number of rows in table */
} TblInfo;

typedef struct {
	char opt_elem[STIS_CBUF];
} TblRow;


static int OpenTdsTab1 (char *, TblInfo *, TdsInfo *);
static int ReadTdsTab1 (TblInfo *, int, TblRow *);
static int ReadTdsArray1 (TblInfo *, int, TdsInfo *);
static int CloseTdsTab1 (TblInfo *);

/* This routine gets time-dependent sensitivity info from TDSTAB (_tds)
   and saves it in the time-dependent sensitivity information structure.

   The time-dependent sensitivity table should contain the following:
        header parameters:
          none needed
        columns:
          OPT_ELEM:   Grating name.
          WAVELENGTH: Array of wavelengths (Angstroms) for correction factors.
          NWL:        Number of wavelengths in wavelength array.
          REFTEMP:    Reference temperature.
          TEMPSENS:   Array (NWL) of temperature sensitivity factors.

   Rows are selected on OPT_ELEM. If a matching row is found, the array sizes
   are gotten, memory is allocated, and the arrays are read from the table.

   When done, the memory should be freed by calling FreeTds1.

   Revision history:
   ----------------
   2011 Nov 17  -  Copied from ../lib/gettds.c, modified to compute just
		   the temperature correction, not time-dependent corr.
   2011 Dec 12  -  In GetTds1 and ReadTdsArray1, set status (but ignore it)
		   from the value returned by CloseTdsTab1.
*/

int GetTds1 (char *tabname, char *opt_elem, TdsInfo *tds) {

/* arguments:
char *tabname    i: table name
char *opt_elem   i: grating name
TdsInfo *tds     o: time-dependent sensitivity info
*/

	int status;

	TblInfo tabinfo;	/* pointer to table descriptor, etc */
	TblRow tabrow;		/* values read from a table row */

	int row;		/* loop index */
	int foundit = 0;	/* true if parameters found in table */

	RefTab tempreftab;	/* holding place for table name */
	void FreeTds1 (TdsInfo *);

	/* Open the time-dependent sensitivity table. */
	if ((status = OpenTdsTab1 (tabname, &tabinfo, tds)))
	    return (status);

	for (row = 1;  row <= tabinfo.nrows;  row++) {

	    if ((status = ReadTdsTab1 (&tabinfo, row, &tabrow)))

		return (status);

	    if (SameString (tabrow.opt_elem, opt_elem)) {

		foundit = 1;

		/* Get pedigree & descrip from the row. */
	        InitRefTab (&tempreftab);
	        strcpy (tempreftab.name, tabname);
		if ((status = RowPedigree (&tempreftab, row,
                        tabinfo.tp, tabinfo.cp_pedigree, tabinfo.cp_descrip)))
		    return (status);
		if (tempreftab.goodPedigree == DUMMY_PEDIGREE) {
		    status = CloseTdsTab1 (&tabinfo);
		    return (DUMMY);
		}

		/* Read temperature-dependent sensitivity info into the
		   tds structure.
		*/
		if ((status = ReadTdsArray1 (&tabinfo, row, tds)))
		    return (status);
	    }
	}

	if ((status = CloseTdsTab1 (&tabinfo)))
	    return (status);

	if (!foundit) {
	    printf ("Warning  Matching row not found in TDSTAB %s\n", tabname);
	    printf ("Warning  OPT_ELEM %s\n", opt_elem);
	    printf ("Warning  Skipping TDS correction.\n");
	    return (ROW_NOT_FOUND);
	}

	return (0);
}

/* This routine opens the TDS table, finds the columns for temperature
   dependence, and gets the total number of rows in the table.
*/

static int OpenTdsTab1 (char *tname, TblInfo *tabinfo, TdsInfo *tds) {

	tabinfo->tp = c_tbtopn (tname, IRAF_READ_ONLY, 0);
	if (c_iraferr()) {
	    printf ("ERROR    TDSTAB `%s' not found\n", tname);
	    return (OPEN_FAILED);
	}

	tabinfo->nrows = c_tbpsta (tabinfo->tp, TBL_NROWS);

	c_tbcfnd1 (tabinfo->tp, "OPT_ELEM",   &tabinfo->cp_optelem);
	c_tbcfnd1 (tabinfo->tp, "WAVELENGTH", &tabinfo->cp_wl);
	c_tbcfnd1 (tabinfo->tp, "NWL",        &tabinfo->cp_nwl);
	if (tabinfo->cp_optelem == 0 ||
	    tabinfo->cp_wl      == 0 ||
	    tabinfo->cp_nwl     == 0) {
	    printf ("ERROR    Column not found in TDSTAB\n");
	    c_tbtclo (tabinfo->tp);
	    return (COLUMN_NOT_FOUND);
	}

	/* Find the columns for correcting for temperature dependence. */
	c_tbcfnd1 (tabinfo->tp, "REFTEMP",  &tabinfo->cp_reftemp);
	c_tbcfnd1 (tabinfo->tp, "TEMPSENS", &tabinfo->cp_tempsens);
	if (tabinfo->cp_reftemp == 0 ||tabinfo->cp_tempsens == 0) {
	    printf ("Warning  Column REFTEMP or TEMPSENS not found in %s;\n",
			tname);
	    printf (
	"Warning  no temperature correction applied to sensitivity\n");
	}

	/* Pedigree and descrip are optional columns. */
	c_tbcfnd1 (tabinfo->tp, "PEDIGREE", &tabinfo->cp_pedigree);
	c_tbcfnd1 (tabinfo->tp, "DESCRIP", &tabinfo->cp_descrip);

	return (0);
}

/* This routine reads the column (mirror or grating name) used to select the
   correct row.
*/

static int ReadTdsTab1 (TblInfo *tabinfo, int row, TblRow *tabrow) {

	c_tbegtt (tabinfo->tp, tabinfo->cp_optelem, row,
                  tabrow->opt_elem, STIS_CBUF-1);
	if (c_iraferr())
	    return (TABLE_ERROR);

	return (0);
}



/* This routine reads array data from one row. */

static int ReadTdsArray1 (TblInfo *tabinfo, int row, TdsInfo *tds) {

	int nwl, ntemp, i;
        /*int nt, ns, ini, ndim, dim[2];*/
	int status = 0;

	/* Find out how many elements there are in the arrays. */

	c_tbegti (tabinfo->tp, tabinfo->cp_nwl, row, &tds->nwl);
	if (c_iraferr())
	    return (TABLE_ERROR);

	/* Allocate memory. */

	tds->wl        = (double *) calloc (tds->nwl, sizeof(double));
	tds->temp_sens = (double *) calloc (tds->nwl, sizeof(double));
	if (tds->temp_sens == NULL || tds->wl == NULL) {
	    status = CloseTdsTab1 (tabinfo);
	    return (OUT_OF_MEMORY);
	}
	tds->format = COS_TDS_FORMAT;		/* this is not relevant */
	tds->time = NULL;
	tds->slope = NULL;
	tds->intercept = NULL;
	tds->nt = 0;
	tds->ref_time = 0.;

	tds->allocated = 1;

	nwl = c_tbagtd (tabinfo->tp, tabinfo->cp_wl, row, tds->wl, 1, tds->nwl);
	if (c_iraferr())
	    return (TABLE_ERROR);

	/* Check if all elements that were expected were actually read. */

	if (nwl < tds->nwl) {
	    c_tbtclo (tabinfo->tp);
	    printf ("ERROR    Not all values were read from TDSTAB\n");
	    return (TABLE_ERROR);
	}

	/* Get the reference temperature and temperature sensitivity. */
	if (tabinfo->cp_reftemp == 0 || tabinfo->cp_tempsens == 0) {
	    /* set the temperature sensitivity to zero */
	    tds->ref_temp = -1.;
	    for (i = 0; i < tds->nwl; i++)
		tds->temp_sens[i] = 0.;
	} else {
	    c_tbegtd (tabinfo->tp, tabinfo->cp_reftemp, row, &tds->ref_temp);
	    if (c_iraferr())
		return (TABLE_ERROR);
	    ntemp = c_tbagtd (tabinfo->tp, tabinfo->cp_tempsens, row,
			tds->temp_sens, 1, tds->nwl);
	    if (c_iraferr())
		return (TABLE_ERROR);
	    if (ntemp < tds->nwl) {
		c_tbtclo (tabinfo->tp);
		printf (
		"ERROR    Not all TEMPSENS values were read from TDSTAB\n");
		return (TABLE_ERROR);
	    }
	}

	return (0);
}

/* This routine closes the TDSTAB table. */

static int CloseTdsTab1 (TblInfo *tabinfo) {

	c_tbtclo (tabinfo->tp);
	if (c_iraferr())
	    return (TABLE_ERROR);

	return (0);
}

/* And this routine frees memory. */

void FreeTds1 (TdsInfo *tds) {

    /*int i;*/

        if (tds->allocated) {
            free (tds->wl);
            free (tds->temp_sens);
            tds->allocated = 0;
        }
}
