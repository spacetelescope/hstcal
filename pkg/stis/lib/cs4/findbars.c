# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include "stis.h"
# include "calstis4.h"
# include "hstcalerr.h"

# define ARCSEC_PER_DEGREE  3600.
# define MIN_BARWEIGHT      0.5		/* fraction of expected peak */
# define MAX_BARWEIGHT      1.5		/* fraction of expected peak */
# define CENTROID_CUTOFF    0.5		/* fraction of actual peak */
# define OUTLIER_CUTOFF     0.3		/* pixel; for rejecting outliers */

/* maximum deviations of individual bars from each other */
# define MAX_DIFF_WARNING   5.
# define MAX_DIFF_ERROR    10.

/* possible values of 'bad' in Centroid */
# define BAR_OK                   0
# define SEARCH_RANGE_TOO_SMALL   1
# define PEAK_VALUE_OUT_OF_RANGE  2

static void ScaleBars (double [], double [], int, double, double,
		double [], double []);
static void MakeTemplateBar (double, double, double [], int, FILE *);
static void XCorr (double [], double [], int, double [], int);
static void Centroid (double [], int, double, double *, double *, FILE *);
static int AvgShifts (double [], double [], int, double *);

/* This routine finds the shift of a long-slit spectrum in the spatial
   direction by finding the occulting bars.
   Nominally, the values in the input data are expected to be centered
   on crpix.  If the values are located at a larger pixel number, then
   the shift will be positive.

   A status of NO_GOOD_DATA will be returned if there was no disastrous
   error but the shift could not be determined because all bars were
   off the edge of the image or were flagged as bad.

   Phil Hodge, 1998 Dec 11:
	Add sts to calling sequence, for sp_range and dbg.

   Phil Hodge, 2000 Mar 27:
	Rewrite.  Remove verbose from calling sequence.

   Phil Hodge, 2002 Aug 26:
	Replace the algorithm in Centroid for computing the mean shift.

   Phil Hodge, 2005 May 25:
	In MakeTemplateBar, initialize exactly nv elements of template
	to zero; previously the first loop was from 0 to i_lower, and
	i_lower isn't necessarily less than nv.
*/

int FindBars (StisInfo4 *sts, int nbars, double barlocn[], double barwidth[],
	double *v, short *qv, int nv,
	double crpix, double cdelt,
	double *shift) {

/* arguments:
StisInfo4 *sts     i: calibration switches and info
int nbars          i: number of occulting bars
double barlocn[]   i: array giving location of each bar (arcsec)
double barwidth[]  i: array giving width of each bar (arcsec)
double v[]         i: slit pattern
short qv[]         i: data quality flags for v
int nv             i: size of v and qv arrays
double crpix       i: reference pixel in input data (zero indexed)
double cdelt       i: degrees per pixel
double *shift      o: the shift, in pixels
*/

	double *pbarlocn, *pbarwidth;	/* locations and widths in pixels */
	double *inv;		/* inverted slit pattern */
	double *template;	/* cross corr. slit with this to get shift */
	double *barshift, *barweight;	/* shift and weight for each bar */
	double *xc;		/* cross correlation of inv and template */
	int i;
	int status;
	int InvertSlit (double [], short [], int, double [], FILE *);

	*shift = 0.;		/* initial value */

	pbarlocn = malloc (nv * sizeof (double));
	pbarwidth = malloc (nv * sizeof (double));
	inv = malloc (nv * sizeof (double));
	template = malloc (nv * sizeof (double));
	barshift = malloc (nbars * sizeof (double));
	barweight = malloc (nbars * sizeof (double));
	xc = malloc (sts->sp_range * sizeof (double));
	if (pbarlocn == NULL || pbarwidth == NULL ||
		inv == NULL || template == NULL ||
		barshift == NULL || barweight == NULL ||
		xc == NULL) {
	    return (OUT_OF_MEMORY);
	}

	/* Convert bar info from arcseconds to pixels. */
	ScaleBars (barlocn, barwidth, nbars, crpix, cdelt, pbarlocn, pbarwidth);

	/* Normalize and invert the slit pattern. */
	if ((status = InvertSlit (v, qv, nv, inv, sts->dbg)))
	    return (status);

	/* Find each occulting bar. */
	for (i = 0;  i < nbars;  i++) {

	    if (pbarwidth[i] > 0.) {

		MakeTemplateBar (pbarlocn[i], pbarwidth[i],
			template, nv, sts->dbg);

		XCorr (inv, template, nv, xc, sts->sp_range);

		Centroid (xc, sts->sp_range, pbarwidth[i],
			&barshift[i], &barweight[i], sts->dbg);

	    } else {

		barshift[i] = 0.;
		barweight[i] = 0.;
	    }

	    if (sts->verbose) {
		if (barweight[i] == 0.) {
		    printf ("Warning  Bar %d could not be found.\n", i+1);
		} else {
		    printf ("         shift of bar %d is %.3f, weight = %.5g\n",
			i+1, barshift[i], barweight[i]);
		}
	    }
	}

	if ((status = AvgShifts (barshift, barweight, nbars, shift)))
	    return (status);

	free (pbarlocn);
	free (pbarwidth);
	free (inv);
	free (template);
	free (barshift);
	free (barweight);
	free (xc);

	return (0);
}

static void ScaleBars (double barlocn[], double barwidth[], int nbars,
		double crpix, double cdelt,
		double pbarlocn[], double pbarwidth[]) {

/* arguments:
double barlocn[]    i: array giving location of each bar (arcsec)
double barwidth[]   i: array giving width of each bar (arcsec)
int nbars           i: number of occulting bars
double crpix        i: pixel zero point for barlocn (reference pixel)
double cdelt        i: degrees per pixel
double pbarlocn[]   o: barlocn converted to pixels
double pbarwidth[]  o: barwidth converted to pixels
*/

	double scale;		/* arcseconds per pixel */
	int i;

	scale = cdelt * ARCSEC_PER_DEGREE;

	for (i = 0;  i < nbars;  i++) {
	    pbarlocn[i] = crpix + barlocn[i] / scale;
	    pbarwidth[i] = barwidth[i] / scale;
	}
}

/* This routine creates a template array which is zero everywhere except
   at the location of one occulting bar, where it is one.  Note that this
   template is for just one bar, not all the bars on the slit.
*/

static void MakeTemplateBar (double pbarlocn, double pbarwidth,
		double template[], int nv, FILE *dbg) {

/* arguments:
double pbarlocn    i: location of current bar (pixels)
double pbarwidth   i: width of current bar (pixels)
double template[]  o: array to be cross correlated with inverted observed slit
int nv             i: size of template array
FILE *dbg          i: file handle for debug output
*/

	double lower, upper;	/* expected locations of bar edges */
	int i_lower, i_upper;	/* nearest integers to lower & upper */
	int lower_ch, upper_ch;	/* loop limits, chopped off at boundaries */
	double fraction;	/* fraction of a pixel covered by the bar */
	int i;

	lower = pbarlocn - pbarwidth / 2.;
	upper = pbarlocn + pbarwidth / 2.;
	i_lower = NINT (lower);
	i_upper = NINT (upper);

	for (i = 0;  i < nv;  i++) {
	    template[i] = 0.;
	}

	if (i_upper > 0 && i_lower < nv - 1) {

	    lower_ch = i_lower < 0 ? 0 : i_lower + 1;
	    upper_ch = i_upper >= nv ? nv - 1 : i_upper - 1;
	    for (i = lower_ch;  i <= upper_ch;  i++) {
		template[i] = 1.;
	    }
	}

	for (i = i_upper + 1;  i < nv;  i++) {
	    template[i] = 0.;
	}

	if (i_lower >= 0 && i_lower < nv) {
	    fraction = 0.5 - lower + i_lower;
	    template[i_lower] = fraction;
	}

	if (i_upper >= 0 && i_upper < nv) {
	    fraction = upper + 0.5 - i_upper;
	    template[i_upper] = fraction;
	}

	if (dbg != NULL) {
	    fprintf (dbg, "\n");
	    fprintf (dbg,
		"# (FindBars) Expected locations of bar edges are %.6g, %.6g\n",
			lower, upper);
	    fflush (dbg);
	}
}

/* This routine takes the cross correlation between the inverted observed
   slit pattern and the template.  The location of the peak in the cross
   correlation indicates the shift of the occulting bar that is included
   in the template.  If the shift is zero, the peak in xc will be at the
   middle pixel, which is (range - 1) / 2 using zero indexing.  If the
   peak is at a larger pixel number than the middle, the shift is positive.

   The regions where the slit and template do not overlap are not
   included in the sum.  This is appropriate because the slit pattern is
   inverted, i.e. near zero within an occulting bar and near one outside
   any occulting bar.
*/

static void XCorr (double inv[], double template[], int nv,
		double xc[], int range) {

/* arguments:
double inv[]          i: normalized and inverted observed slit pattern
double template[]     i: template slit, with just one occulting bar
int nv                i: size of inv and template arrays
double xc[]           o: cross correlation of inv and template
int range             i: size of xc array (an odd number)
*/

	double sum;		/* for accumulating the sum at xc[j] */
	int i, j;
	int half;		/* (range - 1) / 2 */

	half = (range - 1) / 2;

	/* zero offset */
	sum = 0.;
	for (i = 0;  i < nv;  i++)
	    sum += inv[i] * template[i];
	xc[half] = sum;

	for (j = 1;  j <= half;  j++) {

	    /* shift the template toward the left */
	    sum = 0.;
	    for (i = 0;  i < nv-j;  i++)
		sum += inv[i] * template[i+j];
	    xc[half-j] = sum;

	    /* shift the template toward the right */
	    sum = 0.;
	    for (i = j;  i < nv;  i++)
		sum += inv[i] * template[i-j];
	    xc[half+j] = sum;
	}
}

/* This routine finds the location of the peak in the cross correlation.
   Because the inverted spectrum and the template are normalized, the
   expected height of the peak is equal to the width of the current
   occulting bar, and the width of the peak should be twice that width.
   The weight (barweight) will be set to the peak divided by the expected
   value of the peak.  If the weight is less than or greater than cutoff
   values, the weight will be set to zero so that this bar will not be
   included when the average shift is computed later.
*/

static void Centroid (double xc[], int range, double pbarwidth,
		double *barshift, double *barweight, FILE *dbg) {


/* arguments:
double xc[]        i: cross correlation of inv and template
int range          i: size of xc array (an odd number)
double pbarwidth   i: the width of the current bar (pixels)
double *barshift   o: the shift found for this bar
double *barweight  o: the weight for the shift
FILE *dbg          i: file handle for debug output
*/

	double maxval;		/* maximum value in xc */
	int imax;		/* index of maxval in xc */
	double peak;		/* location of centroid in xc */
	int middle;		/* pixel in xc corresponding to zero shift */
	int bad;		/* true if bar was probably not found */
	int i;

	double minval;		/* lower limit when finding location */
	double value;		/* a value below max in xc */
	double x_left, x_right;	/* interpolated locations at which xc=value */
	double *cent;		/* array of center values */
	double median;		/* for rejecting outliers */
	double sum, sumv;	/* for accumulating sums */
	int n;			/* number of valid elements in cent */
	int foundit;

	/* range is odd, middle is zero indexed */
	middle = range / 2;

	/* Find the maximum value and its location. */
	bad = BAR_OK;				/* initial values */
	imax = 0;
	maxval = xc[0];
	for (i = 1;  i < range;  i++) {
	    if (xc[i] > maxval) {
		imax = i;
		maxval = xc[i];
	    }
	}
	if (imax == 0 || imax == range-1)
	    bad = SEARCH_RANGE_TOO_SMALL;

	/* The peak should have a value of barwidth. */
	*barweight = maxval / pbarwidth;

	if (*barweight < MIN_BARWEIGHT || *barweight > MAX_BARWEIGHT)
	    bad = PEAK_VALUE_OUT_OF_RANGE;

	if (bad != BAR_OK) {

	    printf ("Warning  Skipping current occulting bar ... \\\n");

	    if (bad == SEARCH_RANGE_TOO_SMALL) {

		printf (
"Warning  Peak of cross correlation is at end of search range. \\\n");
		printf (
"Warning  This probably means the search range is too small; \\\n");
		printf (
"Warning  check the value of SP_RANGE in the WCPTAB.\n");

	    } else if (bad == PEAK_VALUE_OUT_OF_RANGE) {

		printf (
"Warning  Peak of cross correlation is %.6g of the expected value, \\\n",
			*barweight);
		printf (
"Warning  which is outside the allowed range %.6g to %.6g\n",
			MIN_BARWEIGHT, MAX_BARWEIGHT);
	    }

	    if (dbg != NULL) {
		fprintf (dbg, "# cross correlation:\n");
		for (i = 0;  i < range;  i++) {
		    fprintf (dbg, "%.6g", xc[i]);
		    if (i == middle)
			fprintf (dbg, " <-- nominal peak is here");
		    fprintf (dbg, "\n");
		}
		fprintf (dbg, "# skipping this bar ...\n");
		if (bad == SEARCH_RANGE_TOO_SMALL) {
		    fprintf (dbg,
	"#  peak of cross correlation is at end of search range\n");
		} else if (bad == PEAK_VALUE_OUT_OF_RANGE) {
		    fprintf (dbg,
	"#  peak of cross correlation is %.6g of the expected value,\n",
			*barweight);
		    fprintf (dbg,
	"#  which is outside the allowed range %.6g to %.6g\n",
			MIN_BARWEIGHT, MAX_BARWEIGHT);
		}
		fflush (dbg);
	    }
	    *barshift = 0.;
	    *barweight = 0.;
	    return;
	}

	cent = malloc (range * sizeof (double));
	if (cent == NULL) {
	    *barshift = 0.;
	    *barweight = 0.;
	    printf ("Warning  Out of memory in FindBars.\n");
	    return;
	}

	/* The nominal shape of the xc curve is:
		flat near the left end
		linear with a slope of +1 to the left of the peak
		linear with a slope of -1 to the right of the peak
		flat near the right end
	   A horizontal line below the peak (but above the flat regions
	   at the ends) will intersect the xc curve at two points, and
	   the midpoint between those two points is a measure of the
	   location of the bar.  This section finds a sequence of such
	   measurements at values separated by unit increments below
	   the peak; the average of these (with the appropriate offset)
	   is taken to be the shift.
	*/
	minval = maxval * CENTROID_CUTOFF;
	for (value = maxval-1., n = 0;  value > minval;  value--, n++) {

	    /* find i such that xc[i] <= value and xc[i+1] > value */
	    foundit = 0;
	    for (i = imax;  i > 0;  i--) {
		if (xc[i] <= value && xc[i+1] > value) {
		    foundit = 1;
		    break;
		}
	    }
	    if (!foundit)
		break;

	    /* interpolate to get x_left */
	    x_left = (double)i + (value - xc[i]) / (xc[i+1] - xc[i]);

	    /* find i such that xc[i] >= value and xc[i+1] < value */
	    foundit = 0;
	    for (i = imax;  i < range-1;  i++) {
		if (xc[i] >= value && xc[i+1] < value) {
		    foundit = 1;
		    break;
		}
	    }
	    if (!foundit)
		break;

	    /* interpolate to get x_right */
	    x_right = (double)i + (xc[i] - value) / (xc[i] - xc[i+1]);

	    cent[n] = (x_left + x_right) / 2.;
	}

	peak = (double) imax;			/* default value */
	if (n == 0) {
	    *barshift = 0.;
	    *barweight = 0.;
	} else if (n == 1) {
	    *barshift = cent[0] - (double) middle;
	    printf ("Warning  Only one point used for bar location.\n");
	} else {
	    /* reject outliers from cent */
	    median = MedianDouble (cent, n, 0);
	    for (i = 0;  i < n;  i++) {
		if (fabs (cent[i] - median) > OUTLIER_CUTOFF) {
		    cent[i] = -1.;		/* rejected */
		}
	    }

	    /* compute the mean and stddev of the remaining elements of cent */
	    sumv = 0.;
	    sum = 0.;
	    for (i = 0;  i < n;  i++) {
		if (cent[i] >= 0.) {
		    sum++;
		    sumv += cent[i];
		}
	    }
	    if (sum > 0.) {
		peak = sumv / sum;
		*barshift = peak - (double) middle;
	    } else {
		*barshift = 0.;
		*barweight = 0.;
	    }
	}

	free (cent);

	/* Write info to debug file. */
	if (dbg != NULL) {

	    double offset;		/* offset of peak from imax */
	    offset = peak - (double) imax;

	    fprintf (dbg,
		"# (FindBars) height of cross correlation should be %.6g\n",
			pbarwidth);
	    fprintf (dbg, "# (FindBars) weight for current bar = %.6g\n",
			*barweight);
	    fprintf (dbg, "# (FindBars) cross correlation:\n");
	    for (i = 0;  i < range;  i++) {
		fprintf (dbg, "%.6g", xc[i]);
		if (i == middle)
		    fprintf (dbg, " <-- nominal peak is here");
		if (i == imax) {
		    fprintf (dbg,
			" <-- peak is ");
		    if (offset == 0.)
			fprintf (dbg, "here");
		    else if (offset > 0.)
			fprintf (dbg, "%.6g down from here", offset);
		    else
			fprintf (dbg, "%.6g up from here", -offset);
		}
		fprintf (dbg, "\n");
	    }
	    fflush (dbg);
	}

	return;
}

static int AvgShifts (double barshift[], double barweight[], int nbars,
		double *shift) {

/* arguments:
double barshift[]    i: shift for each bar
double barweight[]   i: weight for each shift
int nbars            i: number of bars
double *shift        o: the average shift
*/

	double sumsw;		/* sum of (shift * weight) */
	double sumw;		/* sum of weights */
	int i;
	double min_shift, max_shift, min_weight;
	int wmin;		/* index of shift with minimum weight */
	int ngood;		/* number of shifts with non-zero weight */

	/* First check that the bars were actually found. */
	sumw = 0.;
	for (i = 0;  i < nbars;  i++)
	    sumw += barweight[i];
	if (sumw < MIN_BARWEIGHT)
	    return (NO_GOOD_DATA);

	ngood = 0;
	for (i = 0;  i < nbars;  i++) {

	    if (barweight[i] > MIN_BARWEIGHT) {	/* the weight could be zero */

		if (ngood == 0) {
		    min_shift = barshift[i];
		    max_shift = barshift[i];
		    wmin = i;
		    min_weight = barweight[i];
		} else {
		    if (barshift[i] < min_shift)
			min_shift = barshift[i];
		    if (barshift[i] > max_shift)
			max_shift = barshift[i];
		    if (barweight[i] < min_weight) {
			wmin = i;
			min_weight = barweight[i];
		    }
		}
		ngood++;
	    }
	}

	if (ngood < 1)
	    return (NO_GOOD_DATA);

	if (max_shift - min_shift > MAX_DIFF_WARNING)
	    printf ("Warning  Shifts of individual bars differ by %.2f\n",
		max_shift - min_shift);

	if (max_shift - min_shift > MAX_DIFF_ERROR)
	    return (NO_GOOD_DATA);

	/* Take the average of the shifts.  If the range of shifts is
	   large, exclude the one with the lowest weight.
	*/
	if (max_shift - min_shift <= MAX_DIFF_WARNING)
	    wmin = -10;		/* in the loop over i, this will never match */

	ngood = 0;
	sumw = 0.;
	sumsw = 0.;
	for (i = 0;  i < nbars;  i++) {

	    if (i == wmin) {
		printf ("Warning  bar %d excluded due to low weight\n", i+1);
		continue;		/* skip this point */
	    }

	    if (barweight[i] > MIN_BARWEIGHT) {
		sumw += barweight[i];
		sumsw += barshift[i] * barweight[i];
		ngood++;
	    }
	}

	if (ngood > 0)
	    *shift = sumsw / sumw;

	return (0);
}
