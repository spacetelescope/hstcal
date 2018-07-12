# include <stdio.h>
# include <string.h>	/* strcpy */
# include <stdlib.h>	/* calloc */

# include "c_iraf.h"
# include "xtables.h"
# include "hstio.h"

# include "stis.h"
# include "stisdef.h"
# include "hstcalerr.h"
# include "stistds.h"


typedef struct {
	IRAFPointer tp;			/* pointer to table descriptor */
	IRAFPointer cp_wl;		/* column descriptors */
	IRAFPointer cp_time;
	IRAFPointer cp_slope;
	IRAFPointer cp_intercept;
	IRAFPointer cp_nwl;
	IRAFPointer cp_nt;
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


static int OpenTdsTab (char *, TblInfo *, TdsInfo *);
static int ReadTdsTab (TblInfo *, int, TblRow *);
static int ReadTdsArray (TblInfo *, int, TdsInfo *);
static int CloseTdsTab (TblInfo *);

/* This routine gets time-dependent sensitivity info from TDSTAB (_tds)
   and saves it in the time-dependent sensitivity information structure.

   The time-dependent sensitivity table should contain the following:
	header parameters:
          REF_TIME    reference time, for COS-like TDS table
	columns:
          OPT_ELEM    Grating name.
          WAVELENGTH  Array of wavelengths (Angstroms) for correction factors.
          TIME        Times (Modified Julian Date) at endpoints of linear
                      segments (last segment doesn't have a right endpoint).
          SLOPE       Slope in correction factor within linear segments;
                      this is a 2-D array of size nwl by nt, with the index
                      over wavelength the more rapidly varying.
          INTERCEPT   Intercept for correction factor within linear segments;
                      this is a 2-D array of size nwl by nt, with the index
                      over wavelength the more rapidly varying.
          NWL         Number of wavelengths in wl and slope arrays.
          NT          Number of times (equal to the number of linear segments).
          REFTEMP     Reference temperature.
          TEMPSENS    Array (NWL) of temperature sensitivity factors.

   Rows are selected on OPT_ELEM. If a matching row is found, the array sizes
   are gotten, memory is allocated, and the arrays are read from the table.

   When done, the memory should be freed by calling FreeTds.

   Revision history:
   ----------------
   02 Jan 02  -  Adapted from GetApThr6 from cs6 (I.Busko)
   27 Dec 04  -  Also get temperature info, if present (Phil Hodge)
   19 Mar 07  -  Remove the definition of c_tbciga from this file (Phil Hodge)
   19 Sep 11  -  Add support for COS-like TDS table format; change calling
		 sequence of OpenTdsTab (Phil Hodge)
   12 Dec 11  -  In GetTds and ReadTdsArray, set status (but ignore it) from
		 the value returned by CloseTdsTab.
*/

int GetTds (char *tabname, char *opt_elem, TdsInfo *tds) {

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

	/* Open the time-dependent sensitivity table. */
	if ((status = OpenTdsTab (tabname, &tabinfo, tds)))
	    return (status);

	for (row = 1;  row <= tabinfo.nrows;  row++) {

	    if ((status = ReadTdsTab (&tabinfo, row, &tabrow)))

		return (status);

	    if (SameString (tabrow.opt_elem, opt_elem)) {

		foundit = 1;

		/* Get pedigree & descrip from the row. We need to
		   build a temporary RefTab structure just to hold
		   the table name, since it is required by the
		   RowPedigree interface.
		*/
	        InitRefTab (&tempreftab);
	        strcpy (tempreftab.name, tabname);
		if ((status = RowPedigree (&tempreftab, row,
                        tabinfo.tp, tabinfo.cp_pedigree, tabinfo.cp_descrip)))
		    return (status);
		if (tempreftab.goodPedigree == DUMMY_PEDIGREE) {
		    status = CloseTdsTab (&tabinfo);
		    return (DUMMY);
		}

		/* Read time-dependent sensitivity info into tds structure. */
		if ((status = ReadTdsArray (&tabinfo, row, tds)))
		    return (status);
	    }
	}

	if ((status = CloseTdsTab (&tabinfo)))
	    return (status);

	if (!foundit) {
	    printf ("Warning  Matching row not found in TDSTAB %s\n", tabname);
	    printf ("Warning  OPT_ELEM %s\n", opt_elem);
	    printf ("Warning  Skipping TDS correction.\n");
	    return (ROW_NOT_FOUND);
	}

	return (0);
}

/* This routine opens the time-dependent sensitivity table, finds the columns
   that we need, and gets the total number of rows in the table.

   This function reads table header keyword REF_TIME (if present) and
   looks for the INTERCEPT column.  If these are present, this function
   also sets a flag in the tds struct indicating that the table has the
   COS-like TDS format rather than the original STIS format.
*/

static int OpenTdsTab (char *tname, TblInfo *tabinfo, TdsInfo *tds) {

	tabinfo->tp = c_tbtopn (tname, IRAF_READ_ONLY, 0);
	if (c_iraferr()) {
	    printf ("ERROR    TDSTAB `%s' not found\n", tname);
	    return (OPEN_FAILED);
	}

	/* This may be changed depending on what we find in the table. */
	tds->format = COS_TDS_FORMAT;

	tabinfo->nrows = c_tbpsta (tabinfo->tp, TBL_NROWS);

	/* This table header keyword will be present if the TDS file is
	   in the COS-like format rather than the original STIS format.
	*/
	tds->ref_time = c_tbhgtd (tabinfo->tp, "REF_TIME");
	if (c_iraferr()) {
	    tds->format = ORIGINAL_TDS_FORMAT;
	    tds->ref_time = 0.;
	    clear_cvoserr();
	}

	/* Find the columns. */

	c_tbcfnd1 (tabinfo->tp, "OPT_ELEM",   &tabinfo->cp_optelem);
	c_tbcfnd1 (tabinfo->tp, "WAVELENGTH", &tabinfo->cp_wl);
	c_tbcfnd1 (tabinfo->tp, "TIME",       &tabinfo->cp_time);
	c_tbcfnd1 (tabinfo->tp, "SLOPE",      &tabinfo->cp_slope);
	c_tbcfnd1 (tabinfo->tp, "NWL",        &tabinfo->cp_nwl);
	c_tbcfnd1 (tabinfo->tp, "NT",         &tabinfo->cp_nt);
	if (tabinfo->cp_wl    == 0 ||
	    tabinfo->cp_time  == 0 ||
	    tabinfo->cp_slope == 0 ||
	    tabinfo->cp_nwl   == 0 ||
	    tabinfo->cp_nt    == 0) {
	    printf ("ERROR    Column not found in TDSTAB\n");
	    c_tbtclo (tabinfo->tp);
	    return (COLUMN_NOT_FOUND);
	}

	/* Find the columns (if present) that are used for correcting
	   for temperature dependence.
	*/
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

	/* Look for the INTERCEPT column, which we expect for a COS-like tds
	   table.
	*/
	c_tbcfnd1 (tabinfo->tp, "INTERCEPT", &tabinfo->cp_intercept);
	if (tabinfo->cp_intercept == 0)
	    tds->format = ORIGINAL_TDS_FORMAT;

	return (0);
}

/* This routine reads the column (grating name) used to select the
   correct row.
*/

static int ReadTdsTab (TblInfo *tabinfo, int row, TblRow *tabrow) {

	c_tbegtt (tabinfo->tp, tabinfo->cp_optelem, row,
                  tabrow->opt_elem, STIS_CBUF-1);
	if (c_iraferr())
	    return (TABLE_ERROR);

	return (0);
}

/* This routine reads array data from one row. */

static int ReadTdsArray (TblInfo *tabinfo, int row, TdsInfo *tds) {

	int nwl, nt, ns, ntemp, dim[2], ini, ndim, i;
	int status = 0;

	/* Find out how many elements there are in the arrays. */

	c_tbegti (tabinfo->tp, tabinfo->cp_nwl, row, &tds->nwl);
	if (c_iraferr())
	    return (TABLE_ERROR);
	c_tbegti (tabinfo->tp, tabinfo->cp_nt, row, &tds->nt);
	if (c_iraferr())
	    return (TABLE_ERROR);

	/* Allocate memory. */

	tds->wl        = (double *) calloc (tds->nwl, sizeof(double));
	tds->temp_sens = (double *) calloc (tds->nwl, sizeof(double));
	tds->time      = (double *) calloc (tds->nt,  sizeof(double));
	tds->slope     = (double **) calloc (tds->nt, sizeof(double *));
	tds->intercept = (double **) calloc (tds->nt, sizeof(double *));
	if (tds->temp_sens == NULL || tds->wl == NULL || tds->time == NULL ||
	    tds->slope == NULL || tds->intercept == NULL) {
	    status = CloseTdsTab (tabinfo);
	    return (OUT_OF_MEMORY);
	}
	for (i = 0; i < tds->nt; i++) {
	    tds->slope[i] = (double *) calloc (tds->nwl, sizeof(double));
	    tds->intercept[i] = (double *) calloc (tds->nwl, sizeof(double));
	    if (tds->slope[i] == NULL || tds->intercept[i] == NULL) {
	        c_tbtclo (tabinfo->tp);
	        return (OUT_OF_MEMORY);
	    }
	}

	tds->allocated = 1;

        /* Read 1-D arrays first. */

	nwl = c_tbagtd (tabinfo->tp, tabinfo->cp_wl, row, tds->wl, 1, tds->nwl);
	if (c_iraferr())
	    return (TABLE_ERROR);
	nt = c_tbagtd (tabinfo->tp, tabinfo->cp_time, row, tds->time,
	               1, tds->nt);
	if (c_iraferr())
	    return (TABLE_ERROR);

	/* The slope and intercept arrays are 2-dimensional, so we must get
	   the stride first.
	*/
	c_tbciga (tabinfo->tp, tabinfo->cp_slope, &ndim, dim, 2);
	ini = 1;	/* Arrays are 1-indexed in spp */
	for (i = 0; i < tds->nt; i++) {
	    ns = c_tbagtd (tabinfo->tp, tabinfo->cp_slope, row,
	                   tds->slope[i], ini, tds->nwl);
	    if (c_iraferr())
	        return (TABLE_ERROR);

	    ini += dim[0];
	}
	if (tabinfo->cp_intercept != 0) {
	    ini = 1;
	    for (i = 0; i < tds->nt; i++) {
		ns = c_tbagtd (tabinfo->tp, tabinfo->cp_intercept, row,
			       tds->intercept[i], ini, tds->nwl);
		if (c_iraferr())
		    return (TABLE_ERROR);

		ini += dim[0];
	    }
	}

	/* Check if all elements that were expected were actually read. */

	if (nwl < tds->nwl || nt < tds->nt) {
	    c_tbtclo (tabinfo->tp);
	    printf ("ERROR    Not all values were read from TDSTAB\n");
	    return (TABLE_ERROR);
	}

	/* Get the reference temperature and temperature sensitivity,
	   if they're present in this TDS table.
	*/
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

static int CloseTdsTab (TblInfo *tabinfo) {

	c_tbtclo (tabinfo->tp);
	if (c_iraferr())
	    return (TABLE_ERROR);

	return (0);
}

/* And this routine frees memory. */

void FreeTds (TdsInfo *tds) {

        int i;

        if (tds->allocated) {
            free (tds->wl);
            free (tds->temp_sens);
            free (tds->time);
            for (i = 0; i < tds->nt; free (tds->slope[i++]));
            for (i = 0; i < tds->nt; free (tds->intercept[i++]));
            free (tds->slope);
            free (tds->intercept);
            tds->allocated = 0;
        }
}
