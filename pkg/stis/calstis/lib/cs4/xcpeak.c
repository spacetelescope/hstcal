# include <stdio.h>
# include <stdlib.h>	/* calloc */
# include "stis.h"
# include "calstis4.h"
# include "hstcalerr.h"

static int CrossCorr4 (double *, short *, double *, int, short, double *, int);

/* This routine calls a routine to do the cross correlation and then
   finds and returns the location of the peak within the xc array.

   Phil Hodge, 1998 Oct 5:
	Change status value 1051 to GENERIC_ERROR_CODE.

   Phil Hodge, 1998 Dec 11:
	Add sts to calling sequence, and print cross correlation to debug file.

   Phil Hodge, 2000 Mar 28:
	Print a newline to the debug file before printing XC array.
*/

int XCPeak (StisInfo4 *sts, double *v, short *qv, double *tspec,
		int nv, int range, short sdqflags, double *shift) {

/* arguments:
double *v, *tspec  i: arrays to be cross correlated
short *qv          i: data quality flags for v
int nv             i: size of arrays
int range          i: size of xc (must be odd)
short sdqflags     i: "serious" data quality flags
double *shift      o: pixel shift from tspec to v
*/

	int status;

	double *xc;		/* buffer to contain cross correlation */
	double maxval;		/* maximum value in xc */
	int imax;		/* index of maxval in xc */
	int i;
	double peak;		/* index of peak in xc */
	int middle;		/* index of xc corresponding to zero shift */
	int ii;			/* loop index for writing to debug file */

	int PeakQuad3 (double *, double *);

	middle = range / 2;

	xc = (double *) calloc (range, sizeof (double));
	if (xc == NULL)
	    return (OUT_OF_MEMORY);

	if ((status = CrossCorr4 (v, qv, tspec, nv, sdqflags, xc, range))) {
	    free (xc);
	    return (status);
	}

	/* Find the peak in the cross correlation.  First find max. */
	imax = 0;				/* initial values */
	maxval = xc[0];
	for (i = 1;  i < range;  i++) {
	    if (xc[i] > maxval) {
		imax = i;
		maxval = xc[i];
	    }
	}

	if (imax == 0 || imax == range-1) {

	    printf ("Warning  Peak in cross correlation is at end of range.\n");
	    peak = 0.;
	    status = NO_GOOD_DATA;

	} else {

	    if (imax == 0)
		i = 0;
	    else if (imax == range-1)
		i = range - 3;
	    else
		i = imax - 1;

	    /* Location of peak of quadratic through xc[i], xc[i+1], xc[i+2],
		relative to i+1.
	    */
	    if ((status = PeakQuad3 (xc+i, &peak))) {
		free (xc);
		return (status);
	    }

	    /* Add i+1 to get the index of the peak, then subtract the
		nominal location (middle) to get the shift.
	    */
	    *shift = peak + (i + 1) - middle;

	    status = 0;
	}

	/* Write info to debug file. */
	if (sts->dbg != NULL) {
	    fprintf (sts->dbg, "\n");
	    fprintf (sts->dbg, "# (XCPeak) cross correlation:\n");
	    for (ii = 0;  ii < range;  ii++) {
		fprintf (sts->dbg, "%.6g", xc[ii]);
		if (ii == middle)
		    fprintf (sts->dbg, " <-- nominal peak is here");
		if (ii == imax) {
		    fprintf (sts->dbg,
			" <-- peak is ");
		    if (peak == 0.)
			fprintf (sts->dbg, "here");
		    else if (peak > 0.)
			fprintf (sts->dbg, "%.6g down from here", peak);
		    else
			fprintf (sts->dbg, "%.6g up from here", -peak);
		}
		fprintf (sts->dbg, "\n");
	    }
	    fflush (sts->dbg);
	}

	free (xc);

	return (status);
}

/* This routine cross correlates x and y, putting the result in xc.
   The middle pixel of xc (i.e. element xc[(range-1)/2]) corresponds
   to zero shift between x and y.

   Element xc[0] is the sum of products starting with y shifted left
   with respect to x.  There are nv elements in the arrays to be cross
   correlated, but elements near the ends are excluded so that the total
   number of products in the sum is nv - (range-1), regardless of the
   shift.  Elements of x that are flagged as bad by qx are set to zero
   (in a copy of x) to reduce their effect on the cross correlation.
   It is an error to have no good pixels at all for any element of x.
*/

static int CrossCorr4 (double *x, short *qx, double *y, int nv,
	short sdqflags, double *xc, int range) {

/* arguments:
double *x, *y      i: arrays to be cross correlated
short *qx          i: data quality flags for x
int nv             i: size of arrays
short sdqflags     i: "serious" data quality flags
double *xc         o: cross correlation
int range          i: size of xc (must be odd)
*/

	double *tx;		/* local copy of x */
	int skip;		/* skip this many at beginning and end */
	double sum;		/* sum of products */
	int ngood;		/* number of pixels not flagged as bad */
	int i, j;

	if (range < 3 ||
	    range / 2 * 2 == range ||
	    range > nv)
	    return (GENERIC_ERROR_CODE);

	skip = (range - 1) / 2;			/* range is odd */

	/* Copy data from x to scratch, setting to zero any element
	   flagged as bad.
	*/
	if ((tx = (double *) calloc (nv, sizeof (double))) == NULL)
	    return (OUT_OF_MEMORY);
	ngood = 0;
	for (i = skip;  i < nv-skip;  i++) {
	    if ( ! (qx[i] & sdqflags) ) {
		tx[i] = x[i];
		ngood++;
	    }
	}
	if (ngood <= 0) {
	    free (tx);
	    return (NO_GOOD_DATA);
	}

	/* Do the cross correlation. */
	for (j = 0;  j < range;  j++) {
	    sum = 0.;
	    for (i = skip;  i < nv-skip;  i++)
		sum += tx[i] * y[i+skip-j];
	    xc[j] = sum;
	}

	free (tx);
	return (0);
}
