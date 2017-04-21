# include <stdio.h>
# include <stdlib.h>		/* calloc */
# include "stis.h"		/* for NINT */
# include "calstis4.h"
# include "err.h"

/* The size of the cross correlation array must be at least three pixels
   because we will fit a quadratic to the three values nearest the peak.
*/
# define MIN_RANGE       3

# define ARCSEC_PER_DEGREE  3600.

/* This routine finds the shift of an echelle spectrum in the spatial
   direction by cross correlating the data with a model of the slit.
   Nominally, the values in the input data are expected to be centered
   on crpix; if the values are located at a larger pixel number, then
   the shift will be positive.

   Phil Hodge, 1998 Dec 11:
	Add sts to calling sequence, for sp_range and dbg.
*/

int XCEchelle (StisInfo4 *sts, double length,
	double *v, short *qv, int nv, short sdqflags,
	double crpix, double cdelt,
	double *shift) {

/* arguments:
StisInfo4 *sts     i: calibration switches and info
double length      i: width of slit in spatial direction
double *v          i: array to be cross correlated with template of aperture
short *qv          i: data quality flags for v
int nv             i: size of v and qv arrays
short sdqflags     i: "serious" data quality flags
double crpix       i: reference pixel in input data (zero indexed)
double cdelt       i: degrees per pixel
double *shift      o: the shift, in pixels
*/

	int status;

	double *tslit;		/* model of the slit */
	double center;		/* expected location of center of slit */
	double low, high;	/* expected locations of ends of slit */
	int ilow, ihigh;	/* portion that is fully within slit */
	int range;		/* size of xc */
	double c7_shift;	/* shift in units of calstis7 pixels */
	int i;

	int XCPeak (StisInfo4 *, double *, short *, double *,
		int, int, short, double *);

	/* Convert length from arcseconds to pixels. */
	length /= (cdelt * ARCSEC_PER_DEGREE);

	range = sts->sp_range;
	if (range > nv / 2)
	    range = nv / 2;
	if (range < MIN_RANGE)
	    range = MIN_RANGE;
	if (range / 2 * 2 == range)
	    range++;			/* must be odd */

	center = crpix;
	low = center - length / 2.;
	high = center + length / 2.;
	ilow = NINT (low);		/* nearest integer */
	ihigh = NINT (high);

	/* Off the image? */
	if (ilow < 0 || ilow >= nv || ihigh < 0 || ihigh >= nv)
	    return (NO_GOOD_DATA);

	tslit = (double *) calloc (nv, sizeof (double));
	if (tslit == NULL)
	    return (OUT_OF_MEMORY);

	/* Make a model of the slit. */

	for (i = ilow+1;  i <= ihigh-1;  i++)
	    tslit[i] = 1.;

	/* Now assign a fraction at each end. */
	if (ilow >= 0)
	    tslit[ilow] = low + 0.5 - (double)ilow;
	if (ihigh < nv)
	    tslit[ihigh] = (double)ihigh - (high - 0.5);

	/* Do the cross correlation, and find the shift. */
	if ((status = XCPeak (sts, v, qv, tslit, nv, range, sdqflags,
                              &c7_shift))) {
	    free (tslit);
	    return (status);
	}

/* ###
Convert shift in calstis7 pixels to pixels in raw image.  (how?)
To first approximation, this conversion is a NO-OP.
*/
	*shift = c7_shift;

	free (tslit);

	return (0);
}
