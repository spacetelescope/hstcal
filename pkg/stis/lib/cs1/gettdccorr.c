# include <stdio.h>
# include <string.h>
# include <math.h>
# include <float.h>

# include "c_iraf.h"
# include "hstio.h"
# include "xtables.h"
# include "stis.h"
# include "calstis1.h"
# include "hstcalerr.h"
# include "stisdef.h"
# include "stistemperature.h"

# define UNDEFINED_TDC_FORMAT   0
# define POST_SM4_TDC_FORMAT    1
# define ORIG_TDC_FORMAT        2

/* TblInfo and TblRow are for the original type of TDC table. */
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

static int origTdcCorr(StisInfo1 *, double *);
static int postSM4TdcCorr(StisInfo1 *, double, double *);
static double NUVFactor(double, double, double, double, double);
static int OpenTdcTab(char *, TblInfo *);
static int ReadTdcTab(TblInfo *, int, TblRow *);
static int CloseTdcTab(TblInfo *);

/* This routine computes the time-dependent correction factor for the
   NUV MAMA darks.  It gets the relevant parameters from the tdc table
   (keyword name TDCTAB).

   Ivo Busko, 2002 Apr 02:
	Adapted from similar code elsewhere in cs1.

   Phil Hodge, 2004 Dec 27:
	Use sts->detector_temp instead of sts->temperature.  Locally
	(in NUVFactor) convert from degrees Celsius to Kelvins.

   Phil Hodge, 2012 Oct 15:
	Add median_dark to calling sequence of GetTdcCorr.  GetTdcCorr now
	just opens the table to determine which type it is and then calls
	either origTdcCorr (the original GetTdcCorr) or new function
	postSM4TdcCorr.

   Phil Hodge, 2013 Sept 13:
	In postSM4TdcCorr, change the condition for setting read_this_row,
        and add a check that read_this_row is greater than zero.

   Phil Hodge, 2013 Oct 3:
	Rename variable median_dark to mean_dark, as it now is the average,
	not the median.

   Robert Jedrzejewski, 2017 Feb 13:
        In origTdcCorr, change the first argument in the call to RowPedigree
 	from &sts->ccdpar to &sts->tdctab. 
*/

int GetTdcCorr(StisInfo1 *sts, double mean_dark, double *factor) {

/* arguments:
StisInfo1 *sts          i: calibration switches, etc
double mean_dark        i: average of good pixels of dark reference image
double *factor          o: NUV correction factor for dark image
*/
	int status;
	int tbl_type;		/* original or post-SM4 */

	IRAFPointer tp;			/* pointer to table descriptor */
	/* column descriptors */
	IRAFPointer cp_1, cp_2, cp_3, cp_4;

	tp = c_tbtopn(sts->tdctab.name, IRAF_READ_ONLY, 0);
	if (c_iraferr()) {
	    printf("ERROR    Can't open TDCTAB `%s'.\n", sts->tdctab.name);
	    return OPEN_FAILED;
	}

	/* Look for some of the columns we expect in each format TDC table. */
	c_tbcfnd1(tp, "DATE0", &cp_1);
	c_tbcfnd1(tp, "TEMP0", &cp_2);
	c_tbcfnd1(tp, "A1", &cp_3);
	c_tbcfnd1(tp, "D1", &cp_4);
	if (cp_1 != 0 && cp_2 != 0 && cp_3 != 0 && cp_4 != 0) {
	    /* columns were found */
	    tbl_type = POST_SM4_TDC_FORMAT;
	} else {
	    c_tbcfnd1(tp, "MJD", &cp_1);
	    c_tbcfnd1(tp, "SCALE", &cp_2);
	    c_tbcfnd1(tp, "NORM", &cp_3);
	    c_tbcfnd1(tp, "T_MIN", &cp_4);
	    if (cp_1 != 0 && cp_2 != 0 && cp_3 != 0 && cp_4 != 0) {
		tbl_type = ORIG_TDC_FORMAT;
	    } else {
		tbl_type = UNDEFINED_TDC_FORMAT;
	    }
	}

	c_tbtclo(tp);

	if (tbl_type == POST_SM4_TDC_FORMAT) {
	    status = postSM4TdcCorr(sts, mean_dark, factor);
	} else if (tbl_type == ORIG_TDC_FORMAT) {
	    status = origTdcCorr(sts, factor);
	} else {
	    printf("ERROR    TDCTAB `%s' is not a recognized format.\n",
		sts->tdctab.name);
	    return OPEN_FAILED;
	}

	return status;
}

/* This is the post-SM4 function for finding the NUV dark-correction factor.

   The expression for the average dark count rate (low-res pixels) is:

   count_rate = a1 * exp(-(expstart - date0) / d1) +
                a2 * exp(-(expstart - date0) / d2) +
                ta + tb * (temperature - temp0) +
                     tc * (temperature - temp0)**2

   This count_rate is then divided by mean_dark to get the correction
   factor by which the dark reference image should be multiplied to
   correct for the change in dark rate with time.

   The post-SM4 tdc table should contain the following:
	header parameters:
	    none needed
	columns:
	    DATE0	date (MJD) the parameters took effect
	    TEMP0	temperature (Celsius) to be subtracted
	    A1
	    D1
	    A2
	    D2
	    TA
	    TB
	    TC
*/

static int postSM4TdcCorr(StisInfo1 *sts, double mean_dark,
                          double *factor) {

/* arguments:
StisInfo1 *sts          i: calibration switches, etc
double mean_dark        i: mean of good pixels of dark reference image
                           (binned to 1024x1024 pixels in size)
double *factor          o: NUV correction factor for dark image
*/
        /*int status;*/

	IRAFPointer tp;			/* pointer to table descriptor */
	/* column descriptors */
	IRAFPointer cp_date0, cp_temp0,
		cp_a1, cp_d1, cp_a2, cp_d2, cp_ta, cp_tb, cp_tc;
	int nrows;			/* number of rows in table */
	/* values read from a table row */
	double date0;			/* reference date */
	double temp0;			/* reference temperature */
	double a1;
	double d1;
	double a2;
	double d2;
	double ta;
	double tb;
	double tc;
	int row;			/* loop index for row number */
	int read_this_row = -1;		/* number of row of table to read */
	double min_date0 = 0.;		/* min value of data0 */
	double best_date0 = 0.;		/* date0 in row to use */
	double d_date;			/* expstart - date0 */
	double d_temp;			/* temperature - temp0 */
	double count_rate;		/* computed count rate at EXPSTART */

	tp = c_tbtopn(sts->tdctab.name, IRAF_READ_ONLY, 0);
	if (c_iraferr()) {
	    printf("ERROR    Can't open TDCTAB `%s'.\n", sts->tdctab.name);
	    return OPEN_FAILED;
	}
	nrows = c_tbpsta(tp, TBL_NROWS);
	if (nrows < 1) {
	    printf("ERROR    TDCTAB %s is empty.\n", sts->tdctab.name);
	    c_tbtclo(tp);
	    return TABLE_ERROR;
	}

	/* Look for all of the columns we require. */
	c_tbcfnd1(tp, "DATE0", &cp_date0);
	c_tbcfnd1(tp, "TEMP0", &cp_temp0);
	c_tbcfnd1(tp, "A1", &cp_a1);
	c_tbcfnd1(tp, "D1", &cp_d1);
	c_tbcfnd1(tp, "A2", &cp_a2);
	c_tbcfnd1(tp, "D2", &cp_d2);
	c_tbcfnd1(tp, "TA", &cp_ta);
	c_tbcfnd1(tp, "TB", &cp_tb);
	c_tbcfnd1(tp, "TC", &cp_tc);
	if (cp_date0 == 0 || cp_temp0 == 0 ||
	    cp_a1 == 0 || cp_d1 == 0 ||
	    cp_a2 == 0 || cp_d2 == 0 ||
	    cp_ta == 0 || cp_tb == 0 || cp_tc == 0) {
	    printf("ERROR    Column not found in TDCTAB.\n");
	    c_tbtclo(tp);
	    return COLUMN_NOT_FOUND;
	}

	/* Find the minimum value in the DATE0 column. */
	for (row = 1;  row <= nrows;  row++) { /* row number is one-indexed */
	    c_tbegtd(tp, cp_date0, row, &date0);
	    if (c_iraferr()) {
		c_tbtclo(tp);
		return TABLE_ERROR;
	    }
	    if (row == 1)
		min_date0 = date0;
	    else if (date0 < min_date0)
		min_date0 = date0;
	}
	if (sts->expstart < min_date0) {
	    printf("ERROR    Exposure start time precedes earliest date "
	           "in TDCTAB.\n");
	    c_tbtclo(tp);
	    return TABLE_ERROR;
	}

	/* Select the row such that DATE0 is less than but closest to
	   the exposure start time.
	*/
	best_date0 = min_date0;			/* initial value */
	for (row = 1;  row <= nrows;  row++) {
	    c_tbegtd(tp, cp_date0, row, &date0);
	    if (c_iraferr()) {
		c_tbtclo(tp);
		return TABLE_ERROR;
	    }
	    if (date0 < sts->expstart && (row == 1 || date0 > best_date0)) {
		best_date0 = date0;
		read_this_row = row;
	    }
	}
	if (read_this_row < 1) {
	    c_tbtclo(tp);
	    printf("Warning:  No valid row found in TDCTAB %s\n",
	           sts->tdctab.name);
	    *factor = 1.;
	    return 0;
	}

	c_tbegtd(tp, cp_date0, read_this_row, &date0);
	if (c_iraferr()) {
	    c_tbtclo(tp);
	    return TABLE_ERROR;
	}
	d_date = sts->expstart - date0;

	c_tbegtd(tp, cp_temp0, read_this_row, &temp0);
	if (c_iraferr()) {
	    c_tbtclo(tp);
	    return TABLE_ERROR;
	}
	d_temp = sts->detector_temp - temp0;

	c_tbegtd(tp, cp_a1, read_this_row, &a1);
	if (c_iraferr()) {
	    c_tbtclo(tp);
	    return TABLE_ERROR;
	}

	c_tbegtd(tp, cp_d1, read_this_row, &d1);
	if (c_iraferr()) {
	    c_tbtclo(tp);
	    return TABLE_ERROR;
	}

	c_tbegtd(tp, cp_a2, read_this_row, &a2);
	if (c_iraferr()) {
	    c_tbtclo(tp);
	    return TABLE_ERROR;
	}

	c_tbegtd(tp, cp_d2, read_this_row, &d2);
	if (c_iraferr()) {
	    c_tbtclo(tp);
	    return TABLE_ERROR;
	}

	c_tbegtd(tp, cp_ta, read_this_row, &ta);
	if (c_iraferr()) {
	    c_tbtclo(tp);
	    return TABLE_ERROR;
	}

	c_tbegtd(tp, cp_tb, read_this_row, &tb);
	if (c_iraferr()) {
	    c_tbtclo(tp);
	    return TABLE_ERROR;
	}

	c_tbegtd(tp, cp_tc, read_this_row, &tc);
	if (c_iraferr()) {
	    c_tbtclo(tp);
	    return TABLE_ERROR;
	}

	c_tbtclo(tp);

	count_rate = a1 * exp(-d_date / d1) +
	             a2 * exp(-d_date / d2) +
	             ta + tb * d_temp +
	                  tc * d_temp * d_temp;

	if (count_rate > 100.) {
	    printf("Warning  TDC dark count rate = %.6g\n", count_rate);
	}
	*factor = count_rate / mean_dark;
	if (fabs(*factor) > 100. || fabs(*factor) < 0.01) {
	    printf("Warning  TDC correction factor = %.6g.\n", *factor);
	}

	return 0;
}

/* This is the original (pre-SM4) function for finding the NUV
   dark-correction factor.

   The tdc table should contain the following:
	header parameters:
		none needed
	columns:
		MJD:      Modified Julian Day = JD - 2,400,000.5
		SCALE:    reference dark image scale factor
		NORM:     normalization parameter
		T_MIN:    minimum temperature
		THERMCST: thermal constant

   The table is read and then the sought values are found by interpolation.
*/

static int origTdcCorr (StisInfo1 *sts, double *factor) {

/* arguments:
StisInfo1 *sts          i: calibration switches, etc
double *factor          o: NUV correction factor for dark image
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

	/* Open table and find columns. */
	if ((status = OpenTdcTab (sts->tdctab.name, &tabinfo)))
	    return status;
	nrows = tabinfo.nrows;

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

	/* Read rows into arrays. */

	for (row = 1; row <= tabinfo.nrows; row++) {

	    if ((status = ReadTdcTab (&tabinfo, row, &tabrow)))
		return status;

	    if ((status = RowPedigree (&sts->tdctab, row,
                    tabinfo.tp, tabinfo.cp_pedigree, tabinfo.cp_descrip)))
		return status;
	    if (sts->tdctab.goodPedigree == DUMMY_PEDIGREE)
		printf ("Warning  Row %d of TDCTAB is DUMMY.\n", row);

	    mjd[row-1]      = tabrow.mjd;
	    scale[row-1]    = tabrow.scale;
	    norm[row-1]     = tabrow.norm;
	    tmin[row-1]     = (double)tabrow.tmin;
	    thermcst[row-1] = (double)tabrow.thermcst;
	}

	starti = tabinfo.nrows - 1;
	scale0    = interp1d (sts->expstart, mjd, scale,    nrows, &starti);
	norm0     = interp1d (sts->expstart, mjd, norm,     nrows, &starti);
	tmin0     = interp1d (sts->expstart, mjd, tmin,     nrows, &starti);
	thermcst0 = interp1d (sts->expstart, mjd, thermcst, nrows, &starti);

	*factor  = NUVFactor (sts->detector_temp, scale0, norm0, tmin0,
	                      thermcst0);

	free (mjd);
	free (scale);
	free (norm);
	free (tmin);
	free (thermcst);

	if ((status = CloseTdcTab (&tabinfo)))
	    return status;

	return 0;
}


/* This routine opens the tdc table, finds the columns that we
   need, and gets the total number of rows in the table
*/

static int OpenTdcTab (char *tname, TblInfo *tabinfo) {

	tabinfo->tp = c_tbtopn (tname, IRAF_READ_ONLY, 0);
	if (c_iraferr()) {
	    printf ("ERROR    TDCTAB `%s' not found.\n", tname);
	    return OPEN_FAILED;
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
	    return COLUMN_NOT_FOUND;
	}

	/* Pedigree and descrip are optional columns. */
	c_tbcfnd1 (tabinfo->tp, "PEDIGREE", &tabinfo->cp_pedigree);
	c_tbcfnd1 (tabinfo->tp, "DESCRIP", &tabinfo->cp_descrip);

	return 0;
}


/* This routine reads the relevant data from one row. */

static int ReadTdcTab (TblInfo *tabinfo, int row, TblRow *tabrow) {

	c_tbegtd (tabinfo->tp, tabinfo->cp_mjd, row, &tabrow->mjd);
	if (c_iraferr())
	    return TABLE_ERROR;

	c_tbegtd (tabinfo->tp, tabinfo->cp_scale, row, &tabrow->scale);
	if (c_iraferr())
	    return TABLE_ERROR;

	c_tbegtd (tabinfo->tp, tabinfo->cp_norm, row, &tabrow->norm);
	if (c_iraferr())
	    return TABLE_ERROR;

	c_tbegtr (tabinfo->tp, tabinfo->cp_tmin, row, &tabrow->tmin);
	if (c_iraferr())
	    return TABLE_ERROR;

	c_tbegtr (tabinfo->tp, tabinfo->cp_thermcst, row, &tabrow->thermcst);
	if (c_iraferr())
	    return TABLE_ERROR;

	return 0;
}


/* This routine closes the tdctab table. */

static int CloseTdcTab (TblInfo *tabinfo) {

	c_tbtclo (tabinfo->tp);
	if (c_iraferr())
	    return TABLE_ERROR;

	return 0;
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

	return factor;
}
