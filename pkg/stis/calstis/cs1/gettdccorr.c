# include <stdio.h>
# include <string.h>
# include <math.h>
# include <float.h>

# include <c_iraf.h>
# include <hstio.h>
# include <xtables.h>
# include "../stis.h"
# include "calstis1.h"
# include "../stiserr.h"
# include "../stisdef.h"
# include "../stistemperature.h"

typedef struct {
	IRAFPointer tp;			/* pointer to table descriptor */
	IRAFPointer cp_mjd;		/* column descriptors */
	IRAFPointer cp_scale;
	IRAFPointer cp_norm;
	IRAFPointer cp_tmin;
	IRAFPointer cp_thermcst;
	IRAFPointer cp_pedigree;
	IRAFPointer cp_descrip;
	int nrows;			/* number of rows in table */
} TblInfo;

typedef struct {
	double mjd;
	double scale;
	double norm;
	float tmin;
	float thermcst;
} TblRow;

static double NUVFactor (double, double, double, double, double);
static int OpenTdcTab (char *, TblInfo *);
static int ReadTdcTab (TblInfo *, int, TblRow *);
static int CloseTdcTab (TblInfo *);

/* This routine computes the time-dependent correction factor for the
   NUV MAMA darks. It gets the relevant parameters from the tdc table 
   (keyword name TDCTAB).

   The tdc table should contain the following:
	header parameters:
		none needed
	columns:
		MJD:      Modified JUlian Day = JD - 2,400,000.5
		SCALE:    reference dark image scale factor
		NORM:     normalization parameter
		T_MIN:    minimum temperature
		THERMCST: thermal constant

   The table is read and then the sougth values are found by interpolation.

   If the TDCTAB keyword is not present in the header, default values for
   the tdc parameters are used.
   
   Ivo Busko, 2002 Apr 02:
	Adapted from similar code elsewhere in cs1.

   Phil Hodge, 2004 Dec 27:
	Use sts->detector_temp instead of sts->temperature.  Locally
	(in NUVFactor) convert from degrees Celsius to Kelvins.
*/

int GetTdcCorr (StisInfo1 *sts, double *factor) {

/* arguments:
StisInfo1 *sts     io: calibration switches, etc
*/
	int status;

	TblInfo tabinfo;	/* pointer to table descriptor, etc */
	TblRow tabrow;		/* values read from a table row */

	int row;		/* loop index */
	int nrows;		/* # of rows */

	double* mjd;		/* values from table are stored here   */
	double* scale;		/* to be fed into the interpolation. */
	double* norm;
	double* tmin;
	double* thermcst;

	double scale0;		/* interpolated values */
	double norm0;
	double tmin0;
	double thermcst0;
	int starti;		/* starting index for interpolation search */

	double interp1d (double, double *, double *, int, int *);

	/* Open table and find columns only if keyword exists. */

	if (strcmp (sts->tdctab.name, DEFAULT_TDC) != 0) {
	    if (status = OpenTdcTab (sts->tdctab.name, &tabinfo))
	        return (status);
	    nrows = tabinfo.nrows;
	} else 
	    nrows = 2;

	/* Alloc memory. */

	mjd      = (double *) malloc (nrows * sizeof (double));
	scale    = (double *) malloc (nrows * sizeof (double));
	norm     = (double *) malloc (nrows * sizeof (double));
	tmin     = (double *) malloc (nrows * sizeof (double));
	thermcst = (double *) malloc (nrows * sizeof (double));
	if (mjd      == NULL ||
	    scale    == NULL ||
	    norm     == NULL ||
	    tmin     == NULL ||
	    thermcst == NULL)
	    return (status = OUT_OF_MEMORY);

	if (strcmp (sts->tdctab.name, DEFAULT_TDC) != 0) {

	    /* Read rows into arrays. */

	    for (row = 1; row <= tabinfo.nrows; row++) {

	        if (status = ReadTdcTab (&tabinfo, row, &tabrow))
		    return (status);

	        if (status = RowPedigree (&sts->ccdpar, row,
                    tabinfo.tp, tabinfo.cp_pedigree, tabinfo.cp_descrip))
	            return (status);
	        if (sts->tdctab.goodPedigree == DUMMY_PEDIGREE)
	            printf ("Warning  Row %d of TDCTAB is DUMMY.\n", row);

	        mjd[row-1]      = tabrow.mjd;
	        scale[row-1]    = tabrow.scale;
	        norm[row-1]     = tabrow.norm;
	        tmin[row-1]     = (double)tabrow.tmin;
	        thermcst[row-1] = (double)tabrow.thermcst;

	        /*printf ("%d   %g  %g  %g  %g  %g\n", 
	                   row, mjd[row-1], scale[row-1],
	                   norm[row-1], tmin[row-1], thermcst[row-1]);
	        fflush (stdout);*/
	    }
	} else {

	    /* Use defaults. */

	    mjd[0]      = 0;
	    mjd[1]      = DBL_MAX;
	    scale[0]    = 9.012e20;
	    scale[1]    = 9.012e20;
	    norm[0]     = 1.00088;
	    norm[1]     = 1.00088;
	    tmin[0]     = 0.0;
	    tmin[1]     = 0.0;
	    thermcst[0] = 12710.0;
	    thermcst[1] = 12710.0;
	}

	/*printf ("%g\n", sts->expstart);
	fflush (stdout)*/;

	starti = tabinfo.nrows - 1; 
	scale0    = interp1d (sts->expstart, mjd, scale,    nrows, &starti);
	norm0     = interp1d (sts->expstart, mjd, norm,     nrows, &starti);
	tmin0     = interp1d (sts->expstart, mjd, tmin,     nrows, &starti);
	thermcst0 = interp1d (sts->expstart, mjd, thermcst, nrows, &starti);

	/*printf ("%g  %g  %g  %g\n", scale0, norm0, tmin0, thermcst0);
	printf ("%g\n", sts->detector_temp);
	fflush (stdout);*/

	*factor  = NUVFactor (sts->detector_temp, scale0, norm0, tmin0,
	                      thermcst0);

	/*printf ("%g\n", *factor);
	fflush (stdout); */

	free (mjd);
	free (scale);
	free (norm);
	free (tmin);
	free (thermcst);

	if (strcmp (sts->tdctab.name, DEFAULT_TDC) != 0) {
	    if (status = CloseTdcTab (&tabinfo))
	        return (status);
	}

	return (0);
}


/* This routine opens the tdc table, finds the columns that we
   need, and gets the total number of rows in the table
*/

static int OpenTdcTab (char *tname, TblInfo *tabinfo) {

	tabinfo->tp = c_tbtopn (tname, IRAF_READ_ONLY, 0);
	if (c_iraferr()) {
	    printf ("ERROR    TDCTAB `%s' not found.\n", tname);
	    return (OPEN_FAILED);
	}

	tabinfo->nrows = c_tbpsta (tabinfo->tp, TBL_NROWS);

	/* Find the columns. */

	c_tbcfnd1 (tabinfo->tp, "MJD", &tabinfo->cp_mjd);
	c_tbcfnd1 (tabinfo->tp, "SCALE", &tabinfo->cp_scale);
	c_tbcfnd1 (tabinfo->tp, "NORM", &tabinfo->cp_norm);
	c_tbcfnd1 (tabinfo->tp, "T_MIN", &tabinfo->cp_tmin);
	c_tbcfnd1 (tabinfo->tp, "THERMCST", &tabinfo->cp_thermcst);
	if (tabinfo->cp_mjd == 0 ||
	    tabinfo->cp_scale == 0 ||
	    tabinfo->cp_norm == 0 ||
	    tabinfo->cp_tmin == 0 ||
	    tabinfo->cp_thermcst == 0) {
	    printf ("ERROR    Column not found in TDCTAB.\n");
	    c_tbtclo (tabinfo->tp);
	    return (COLUMN_NOT_FOUND);
	}

	/* Pedigree and descrip are optional columns. */
	c_tbcfnd1 (tabinfo->tp, "PEDIGREE", &tabinfo->cp_pedigree);
	c_tbcfnd1 (tabinfo->tp, "DESCRIP", &tabinfo->cp_descrip);

	return (0);
}


/* This routine reads the relevant data from one row. */

static int ReadTdcTab (TblInfo *tabinfo, int row, TblRow *tabrow) {

	c_tbegtd (tabinfo->tp, tabinfo->cp_mjd, row, &tabrow->mjd);
	if (c_iraferr())
	    return (TABLE_ERROR);

	c_tbegtd (tabinfo->tp, tabinfo->cp_scale, row, &tabrow->scale);
	if (c_iraferr())
	    return (TABLE_ERROR);

	c_tbegtd (tabinfo->tp, tabinfo->cp_norm, row, &tabrow->norm);
	if (c_iraferr())
	    return (TABLE_ERROR);

	c_tbegtr (tabinfo->tp, tabinfo->cp_tmin, row, &tabrow->tmin);
	if (c_iraferr())
	    return (TABLE_ERROR);

	c_tbegtr (tabinfo->tp, tabinfo->cp_thermcst, row, &tabrow->thermcst);
	if (c_iraferr())
	    return (TABLE_ERROR);

	return (0);
}


/* This routine closes the tdctab table. */

static int CloseTdcTab (TblInfo *tabinfo) {

	c_tbtclo (tabinfo->tp);
	if (c_iraferr())
	    return (TABLE_ERROR);

	return (0);
}



/* The dark count rate for the NUV MAMA varies with temperature.  For this
   detector, sts->detector_temp was gotten earlier from a header keyword.  In
   this routine we will determine the scale factor by which the dark image
   should be scaled.  The factor will be one for other detectors or if
   the temperature could not be obtained from the header; those cases are
   flagged by a value of temperature less than or equal to zero.
   See the "MAMA Darks" section in Chapter 7 of the STIS Instrument Handbook.

    This algorithm was replaced in Apr 2002 with a more elaborate one in
    which the scaling parameters are variable in time. Their values are
    gotten by interpolation in the tdc table.
*/

static double NUVFactor (double temperature, double scale, double norm,
                         double tmin, double thermcst) { 

/* arguments:
double temperature  i: from header keyword (will be converted to Kelvin)
double scale        i:
double norm         i: values interpolated from tdc table
double tmin         i: 
double thermcst     i:
the function value is the factor by which the dark image should be multiplied
*/

	double factor, t_dark, tmin_k, darkrate;

	temperature += CELSIUS_TO_KELVIN;

	tmin_k = tmin + CELSIUS_TO_KELVIN;
	t_dark = (temperature > tmin_k) ? temperature : tmin_k;
	darkrate = scale * norm * exp (-thermcst / t_dark);
	factor = darkrate / 1190.3435;

	return (factor);
}
