# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "stis.h"
# include "calstis4.h"
# include "hstcalerr.h"

# define NPARTS1   8
# define NPARTS2  12	/* should be at least as large as NPARTS1 */
# define MAX_MAD   6.

static int QMedian (double [], short [], int, double *);
static int QCompare (const void *, const void *);
static void CheckMedians (double [], int, double, double []);
static int MedianSections (double [], short [],
	int, int, int, double [], double [], int *);
static int GetMAD (double [], double [], short [], int, double *);
static void FlagBad (double [], double [], short [], int, double);

/* The slit pattern in v extends over nominally 1024 pixels (depending on
   binning), with perhaps 90 pixels of buffer zone (flagged by qv) on either
   end.  There will typically be two occulting bars in the slit pattern,
   where the counts drop to near zero.  This routine fits a curve to the
   slit pattern, hopefully not affected by the occulting bars, then
   normalizes and inverts the pattern.  The result should be approximately
   one where the occulting bars are located and near zero everywhere else.
   The purpose of this is to give a good result when cross correlated with
   a template slit pattern that is one at the predicted location of an
   occulting bar and zero elsewhere.

   Phil Hodge, 2000 Mar 27:
	File created.

   Phil Hodge, 2000 Nov 1:
	Return NO_GOOD_DATA if all pixels are flagged as bad.
*/

int InvertSlit (double v[], short qv[], int nv, double inv[], FILE *dbg) {

/* arguments:
double v[]         i: slit pattern
short qv[]         i: data quality flags for v
int nv             i: size of v and qv arrays
double inv[]       o: inverted slit pattern
FILE *dbg          i: file handle for debug output
*/

	short *qvtemp;		/* local copy of qv */
	double *xslit;		/* array of independent var. values for v */
	double *sm_slit;	/* smoothed slit pattern */
	double global_median;	/* median of all (good) data in slit pattern */
	double mad;		/* median of absolute values of deviations */
	int first_good, last_good;
	double *xmed;		/* index for each median */
	double *med;		/* array of medians */
	double *med_ok;		/* medians after replacing negative values */
	int nmed;		/* size of xmed and med arrays */
	int i;
	int status;

	qvtemp = malloc (nv * sizeof(short));
	xslit = malloc (nv * sizeof(double));
	sm_slit = malloc (nv * sizeof(double));
	if (qvtemp == NULL || xslit == NULL || sm_slit == NULL)
	    return (OUT_OF_MEMORY);

	/* Assign values for the independent variable, which we'll use
	   for the spline fit.  Copy the data quality flags to a local
	   array so we can flag additional points, in order to exclude
	   them when taking medians.  (Flagged values in qv will typically
	   be 32767.  qvtemp is set to 1 instead, to make it easier to
	   plot in case debug output was specified.)
	*/
	for (i = 0;  i < nv;  i++) {
	    xslit[i] = (double) i;
	    if (qv[i] == 0)
		qvtemp[i] = 0;
	    else
		qvtemp[i] = 1;
	}

	/* Find the first and last values that are not flagged as bad. */
	first_good = -1;
	last_good = -1;
	for (i = 0;  i < nv;  i++) {
	    if (qv[i] == 0) {
		first_good = i;
		break;
	    }
	}
	for (i = nv-1;  i >= 0;  i--) {
	    if (qv[i] == 0) {
		last_good = i;
		break;
	    }
	}
	if (first_good < 0)
	    return (NO_GOOD_DATA);

	/* Allocate NPARTS2 elements because it's larger than NPARTS1. */
	med = malloc (NPARTS2 * sizeof (double));
	med_ok = malloc (NPARTS2 * sizeof (double));
	xmed = malloc (NPARTS2 * sizeof (double));
	if (med == NULL || med_ok == NULL || xmed == NULL)
	    return (OUT_OF_MEMORY);

	/* The following is an initial rejection cycle. */

	if ((status = MedianSections (v, qv, first_good, last_good,
                                      NPARTS1, xmed, med, &nmed)))
	    return (status);
	if ((status = splint_nr (xmed, med, nmed, xslit, sm_slit, nv)))
	    return (status);
	if ((status = GetMAD (v, sm_slit, qv, nv, &mad)))
	    return (status);
	FlagBad (v, sm_slit, qvtemp, nv, mad);

	/* Now do the fit again, ignoring pixels just flagged as outliers,
	   which should include the occulting bars themselves.
	*/

	if ((status = MedianSections (v, qvtemp, first_good, last_good,
                                      NPARTS2, xmed, med, &nmed)))
	    return (status);

	/* Replace medians that are less than or equal to zero with
	   neighboring values, or with the global median.
	*/
	if ((status = QMedian (v, qvtemp, nv, &global_median)))
	    return (status);
	CheckMedians (med, nmed, global_median, med_ok);

	if ((status = splint_nr (xmed, med_ok, nmed, xslit, sm_slit, nv)))
	    return (status);

	/* Invert the slit pattern.  Note that we test on qv here, not
	   qvtemp, since we expect that qvtemp will be flagged at the
	   occulting bars, and we need the actual v values there.
	*/
	for (i = 0;  i < nv;  i++) {
	    if (qv[i] == 0 && sm_slit[i] > 0.)
		inv[i] = (sm_slit[i] - v[i]) / sm_slit[i];
	    else
		inv[i] = 0.;
	}

	/* Write info to debug file. */
	if (dbg != NULL) {

	    fprintf (dbg,
"# (InvertSlit) pixel, median of slit illumination, corrected median:\n");
	    for (i = 0;  i < nmed;  i++) {
		fprintf (dbg, "%6.1f %.6g %.6g\n",
			xmed[i] + 1., med[i], med_ok[i]);
	    }

	    fprintf (dbg,
"# pixel, slit illumination, smoothed slit, inverted slit, dq, dq_local:\n");
	    for (i = 0;  i < nv;  i++) {
		fprintf (dbg, "%d %.6g %.6g %.6g %d %d\n",
			i+1, v[i], sm_slit[i], inv[i], qv[i], qvtemp[i]);
	    }
	    fflush (dbg);
	}

	free (qvtemp);
	free (xslit);
	free (sm_slit);
	free (med);
	free (med_ok);
	free (xmed);

	return (0);
}

/* This function returns the median of a double-precision array v,
   excluding elements flagged as bad by the qv array.
   The function value will normally be zero, but if there are no good
   data, the function value will be NO_GOOD_DATA.
*/

static int QMedian (double v[], short qv[], int nv, double *median) {

/* arguments:
double v[]      i: input array
short qv[]      i: array of data quality flags
int nv          i: size of arrays
double *median  o: median of v, excluding bad data
*/

	double *vt;		/* a copy of v, but with only good data */
	int ngood;		/* number of good values in vt */
	int i;

	if ((vt = malloc (nv * sizeof (double))) == NULL)
	    return (OUT_OF_MEMORY);

	ngood = 0;
	for (i = 0;  i < nv;  i++) {
	    if (qv[i] == 0) {
		vt[ngood] = v[i];
		ngood++;
	    }
	}
	if (ngood < 1) {
	    *median = 0.;
	    free (vt);
	    return (NO_GOOD_DATA);
	}

	if (ngood == 1) {
	    *median = vt[0];
	    free (vt);
	    return (0);
	}

	qsort (vt, ngood, sizeof (double), QCompare);

	i = ngood / 2;
	if (ngood == i * 2)
	    *median = (vt[i-1] + vt[i]) / 2.;
	else
	    *median = vt[i];

	free (vt);

	return (0);
}

/* This function is used to sort the array. */

static int QCompare (const void *vp, const void *vq) {

	const double *p = vp;
	const double *q = vq;

	if (*p > *q)
	    return (1);
	else if (*p < *q)
	    return (-1);
	else
	    return (0);
}

/* The illumination shouldn't be less than zero.  This routine copies the
   array of medians of sections (gotten by MedianSection) from med to med_ok,
   and it replaces zero or negative values in med_ok with neighboring
   values.  We start in the middle and work toward each end, since the
   illumination is likely to be good near the middle of the image.  If
   the first value (in the middle) is not positive, we will replace it
   with the global_median, since there's no good neighboring value.
*/

static void CheckMedians (double med[], int nmed, double global_median,
		double med_ok[]) {

/* arguments:
double med[]           i: array of medians in regions of slit illumination
int nmed               i: size of med and med_ok arrays
double global_median   i: median of all (good) slit illumination
double med_ok[]        o: copy of med, but with non-positive values replaced
*/

	double previous_value;		/* previous value from med array */
	int first;
	int i;

	for (i = 0;  i < nmed;  i++)
	    med_ok[i] = med[i];
	previous_value = med[nmed/2];

	first = 1;
	for (i = nmed / 2;  i < nmed;  i++) {
	    if (med[i] <= 0.) {
		if (first)
		    med_ok[i] = global_median;
		else
		    med_ok[i] = previous_value;
	    }
	    previous_value = med_ok[i];
	    first = 0;
	}

	first = 1;
	for (i = nmed / 2 - 1;  i >= 0;  i--) {
	    if (med[i] <= 0.) {
		if (first)
		    med_ok[i] = global_median;
		else
		    med_ok[i] = previous_value;
	    }
	    previous_value = med_ok[i];
	    first = 0;
	}
}

/* This routine takes the median of each of (up to) nparts sections of
   the slit illumination pattern v, returning the results in med.  If
   there's no good data in any of these sections, that section will not
   be included in med; i.e. the section locations (given by xmed) need
   not be regularly spaced.
*/

static int MedianSections (double v[], short qv[],
	int first_good, int last_good, int nparts,
	double xmed[], double med[], int *nmed) {

/* arguments:
double v[]         i: slit pattern
short qv[]         i: data quality flags for v
int first_good     i: index of first good value in v
int last_good      i: index of last good value in v
int nparts         i: number of sections into which v is to be divided
double xmed[]      o: index of middle of each section (array size is nparts)
double med[]       o: median of each section (array size is nparts)
int *nmed          o: number of elements of xmed & med that are actually used
*/

	double y;	/* median of a section */
	int nvals;	/* nominal number of elements in each section */
	int n;		/* actual number of values in current section */
	int start;	/* offset to start of current section */
	int i;
	int status;

	nvals = (last_good - first_good + 1) / nparts;
	/* Round up, so we don't lose anything on the end. */
	if (nvals * nparts < last_good-first_good+1)
	    nvals++;

	*nmed = 0;
	start = first_good;
	for (i = 0;  i < nparts;  i++) {

	    if (start >= last_good)
		break;

	    if (start + nvals - 1 <= last_good)
		n = nvals;
	    else
		n = last_good - start + 1;

	    if ((status = QMedian (v+start, qv+start, n, &y))) {

		if (status == OUT_OF_MEMORY)	/* else ignore this section */
		    return (status);

	    } else {

		med[*nmed] = y;
		xmed[*nmed] = start + (double) (n-1) / 2.;
		(*nmed)++;
	    }
	    start += n;
	}

	return (0);
}

/* This routine gets the median of the absolute values of the deviations
   of v from sm_slit.
*/

static int GetMAD (double v[], double sm_slit[], short qv[], int nv,
		double *mad) {

/* arguments:
double v[]         i: slit pattern
double sm_slit[]   i: smoothed slit pattern
short qv[]         i: data quality flags for v
int nv             i: size of v, sm_slit and qv arrays
double *mad        o: median of absolute values of deviations from fit
*/

	double *diff;		/* abs difference between data and fit */
	int i;
	int status;

	diff = malloc (nv * sizeof(double));
	if (diff == NULL)
	    return (OUT_OF_MEMORY);

	for (i = 0;  i < nv;  i++)
	    diff[i] = fabs (v[i] - sm_slit[i]);

	if ((status = QMedian (diff, qv, nv, mad)))
	    return (status);

	free (diff);

	return (0);
}

/* This routine sets qvtemp to 1 for each element of v that deviates
   too much from sm_slit.  This routine does not reset any values in
   qvtemp; if they were previously non-zero they will remain non-zero.
*/

static void FlagBad (double v[], double sm_slit[], short qvtemp[], int nv,
		double mad) {

/* arguments:
double v[]         i: slit pattern
double sm_slit[]   i: smoothed slit pattern
short qvtemp[]     io: local copy of data quality flags for v
int nv             i: size of v, sm_slit and qvtemp arrays
double mad         i: median of absolute values of deviations from fit
*/

	double diff;		/* abs difference between data and fit */
	int i;

	for (i = 0;  i < nv;  i++) {
	    diff = fabs (v[i] - sm_slit[i]);
	    if (diff > MAX_MAD * mad)
		qvtemp[i] = 1;
	}
}
