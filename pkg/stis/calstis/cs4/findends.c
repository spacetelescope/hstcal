# include <stdio.h>
# include "stis.h"
# include "calstis4.h"
# include "err.h"

# define ARCSEC_PER_DEGREE  3600.
# define HIGH_TO_LOW  (-1)
# define LOW_TO_HIGH  (1)

# define MAX_DIFF_WARNING   5.
# define MAX_DIFF_ERROR    10.

static void UpdateRange (int, double, double *, double *);
static int CheckRange (int, double, double);

/* This routine finds the shift of a spectrum in the spatial direction
   by finding the ends of the slit.
   Nominally, the values in the input data are expected to be centered
   on crpix.  If the values are located at a larger pixel number, then
   the shift will be positive.

   A status of NO_GOOD_DATA will be returned if there was no disastrous
   error, but the shift could not be determined because both ends were
   off the edge of the image or were flagged as bad.

   Phil Hodge, 1998 Dec 11:
	Add sts to calling sequence, for sp_range and dbg.
*/

int FindEnds (StisInfo4 *sts, double length,
	double *v, short *qv, int nv,
	double crpix, double cdelt, int verbose,
	double *shift) {

/* arguments:
StisInfo4 *sts     i: calibration switches and info
double length      i: width of slit in spatial direction
double *v          i: array to be cross correlated with template of aperture
short *qv          i: data quality flags for v
int nv             i: size of v and qv arrays
double crpix       i: reference pixel in input data (zero indexed)
double cdelt       i: degrees per pixel
int verbose        i: print individual shifts?
double *shift      o: the shift, in pixels
*/

	int status;

	double scale;		/* arcseconds per pixel */
	double slit_end;	/* pixel location of end of slit */
	double sum;		/* for averaging the shifts */
	int nedges;		/* number of edges included in average */
	double locn0, locn;	/* estimated and actual location of an edge */
	double c7_shift;	/* shift in units of calstis7 pixels */
	double min_shift, max_shift;	/* min & max shifts */

	int FindEdge (StisInfo4 *,
		double *, short *, int, int, double, double *);

	scale = cdelt * ARCSEC_PER_DEGREE;

	/* Get the pixel location of the lower (or left) end of the slit. */
	slit_end = crpix - length / scale / 2.;

	/* Find each edge of the slit, and average the shifts. */
	sum = 0.;
	nedges = 0;

	/* lower edge of slit */
	locn0 = slit_end;
	if ((status = FindEdge (sts, v, qv, nv, LOW_TO_HIGH, locn0, &locn))) {
	    if (status == NO_GOOD_DATA)		/* not fatal yet */
		status = 0;
	    else
		return (status);		/* serious error */
	    if (verbose)
		printf ("         shift of lower edge is undetermined\n");
	} else {
	    nedges++;
	    sum += (locn - locn0);
	    if (verbose)
		printf ("         shift of lower edge is %.2f\n",
			(locn - locn0));
	    UpdateRange (nedges, locn - locn0, &min_shift, &max_shift);
	}

	/* upper edge of slit */
	locn0 = slit_end + length / scale;
	if ((status = FindEdge (sts, v, qv, nv, HIGH_TO_LOW, locn0, &locn))) {
	    if (status == NO_GOOD_DATA)
		status = 0;
	    else
		return (status);
	    if (verbose)
		printf ("         shift of upper edge is undetermined\n");
	} else {
	    nedges++;
	    sum += (locn - locn0);
	    if (verbose)
		printf ("         shift of upper edge is %.2f\n",
			(locn - locn0));
	    UpdateRange (nedges, locn - locn0, &min_shift, &max_shift);
	}

	if (CheckRange (nedges, min_shift, max_shift)) {
	    *shift = 0.;
	    return (NO_GOOD_DATA);
	}

	c7_shift = sum / (double)nedges;

	/* Convert the shift from calstis7 pixels to pixels in raw image. */
	*shift = c7_shift;

	return (0);
}

static void UpdateRange (int nedges, double current_shift,
		double *min_shift, double *max_shift) {

	if (nedges <= 1) {

	    *min_shift = current_shift;
	    *max_shift = current_shift;

	} else {

	    if (current_shift < *min_shift)
		*min_shift = current_shift;

	    if (current_shift > *max_shift)
		*max_shift = current_shift;
	}
}

/* This routine returns one (implying shift not determined) if nedges is
   zero.  If nedges is one, zero is returned, as no further checking can
   be done.  If nedges is two, this routine checks the difference between
   the min and max shifts and prints a warning if the range is too large.
   If the range is still larger, the function will return one, indicating
   that the shift was not reliably determined.

   Note that this routine does not set or use status.
*/

static int CheckRange (int nedges, double min_shift, double max_shift) {

	if (nedges < 1)
	    return (1);		/* no data */
	else if (nedges == 1)
	    return (0);

	if (max_shift - min_shift > MAX_DIFF_WARNING)
	    printf ("Warning  Shifts of lower and upper edges differ by %.2f\n",
		max_shift - min_shift);

	if (max_shift - min_shift > MAX_DIFF_ERROR)
	    return (1);		/* difference is too large */
	else
	    return (0);
}
