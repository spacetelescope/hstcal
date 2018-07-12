# include <stdio.h>
# include <stdlib.h>	/* calloc */
# include "stis.h"
# include "calstis4.h"
# include "hstcalerr.h"

/* This routine finds an edge (expected to be fairly sharp) in 1-D data.

   Status return of NO_GOOD_DATA means we cannot find the current edge.
   That will not be a fatal error unless we can't find any edge at all.

   I'm being lazy here and doing the cross correlation on the entire v
   array.  If this turns out to be a good algorithm, I should reduce the
   range of the cross correlation.

   Phil Hodge, 1997 Sept 11:
	Reduce the interval [first,last] if there are bad pixels,
	rather than just giving up.

   Phil Hodge, 1998 Dec 11:
	Add sts to calling sequence, for sp_range and dbg; write info
	to debug file.
*/

# define HIGH_TO_LOW  (-1)
# define LOW_TO_HIGH  (1)

int FindEdge (StisInfo4 *sts, double *v, short *qv, int nv,
		int maxmin, double locn0, double *locn) {

/* arguments:
StisInfo4 *sts   i: calibration switches and info
double *v        i: input data, collapsed along row or column
short *qv        i: data quality flags for v
int nv           i: size of v and qv arrays
int maxmin       i: +1 if edge goes from low to high; -1 if high to low
double locn0     i: expected location of edge
double *locn     o: actual location of edge
*/

	int status;

	double edge[] = {-1., 0., 1.};
	int sizedge;		/* size of edge array */
	int half;		/* half of sizedge */
	double *xc;		/* array for cross correlation with edge */
	double extreme;		/* maximum (or min) value in xc */
	int ilocn0;		/* nearest integer to locn0 */
	int ipeak;		/* index of extreme in xc */
	int i, j;
	double sum;		/* for doing cross correlation */
	double peak;		/* accurate location of peak in xc */
	int first, last;	/* range of index in v for doing cross corr */
	int i_min, i_max;	/* possibly reduced range of first, last */
	int ii;			/* loop index for writing to debug file */

	int PeakQuad3 (double *, double *);

	/* If the estimated location is off the edge of the data,
	   we can't find this edge, but it may not be a fatal error.
	*/
	ilocn0 = NINT (locn0);
	if (ilocn0 < 0 || ilocn0 >= nv) {
	    if (sts->verbose) {
		printf (
	"         Estimated location (%d) of edge is off the image.\n",
			ilocn0+1);
	    }
	    return (NO_GOOD_DATA);
	}

	sizedge = sizeof (edge) / sizeof (double);
	half = sizedge / 2;

	/* Search for edge between first and last. */
	first = (int) (locn0 + 0.5) - sts->sp_range / 2;
	last = first + sts->sp_range - 1;
	if (first < 0)
	    first = 0;
	if (last > nv - 1)
	    last = nv - 1;

	/* There must be no bad pixels near the edge that we're looking for.
	   It's not a fatal error if there are any, but we can't find the
	   current edge.
	   Start at locn0 and search in both directions for an interval
	   (including locn0) that is not flagged as bad.  If some bad
	   pixels were found, reduce [first,last] to the good interval.
	*/
	i_min = last;  i_max = first;		/* initial values */
	for (i = ilocn0;  i <= last;  i++) {
	    if (qv[i])
		break;
	    else
		i_max = i;
	}
	for (i = ilocn0;  i >= first;  i--) {
	    if (qv[i])
		break;
	    else
		i_min = i;
	}
	if (i_max - i_min < sizedge + 1) {
	    return (NO_GOOD_DATA);
	} else {
	    first = i_min;
	    last = i_max;
	}
	if (sts->dbg != NULL) {
	    fprintf (sts->dbg,
	"# (FindEdge) Searching for an edge between %d and %d inclusive\n",
			first, last);
	}

	if ((xc = (double *) calloc (nv, sizeof (double))) == NULL)
	    return (OUT_OF_MEMORY);

	/* Cross correlate v with edge mask.  (using all of v)
	   If the maximum (or minimum) is at xc[ipeak], then the
	   edge is centered at ipeak.
	*/
	for (i = half;  i <= nv-half-3;  i++) {
	    sum = 0.;
	    for (j = 0;  j < sizedge;  j++)
		sum += (edge[j] * v[i+j-half]);
	    xc[i] = sum;
	}

	/* Find the maximum (or minimum, depending on maxmin) in the
	   cross correlation.
	*/
	ipeak = first;				/* initial values */
	extreme = xc[first];
	if (maxmin == LOW_TO_HIGH) {
	    /* look for a maximum */
	    for (i = first;  i <= last;  i++) {
		if (xc[i] > extreme) {
		    ipeak = i;
		    extreme = xc[i];
		}
	    }
	} else {				/* high to low */
	    /* look for a minimum */
	    for (i = first;  i <= last;  i++) {
		if (xc[i] < extreme) {
		    ipeak = i;
		    extreme = xc[i];
		}
	    }
	}

	if (ipeak <= first || ipeak >= last) {

	    printf (
		"Warning  Edge not found, appears to be at end of range.\n");
	    peak = 0.;
	    status = NO_GOOD_DATA;

	} else {

	    if (ipeak <= first)
		i = first;
	    else if (ipeak >= last)
		i = last - 2;
	    else
		i = ipeak - 1;

	    /* Location of peak of quadratic through xc[i], xc[i+1], xc[i+2],
		relative to i+1.
	    */
	    if ((status = PeakQuad3 (xc+i, &peak)))
		return (status);

	    /* Add index of middle pixel to get location (zero indexed). */
	    *locn = peak + (double)(i) + 1.0;

	    status = 0;
	}

	/* Write info to debug file. */
	if (sts->dbg != NULL) {
	    fprintf (sts->dbg, "# (FindEdge) cross correlation, ");
	    if (maxmin == LOW_TO_HIGH)
		fprintf (sts->dbg, "looking for a maximum:\n");
	    else
		fprintf (sts->dbg, "looking for a minimum:\n");
	    for (ii = first;  ii <= last;  ii++) {
		fprintf (sts->dbg, "%.6g", xc[ii]);
		if (ii == ilocn0) {
		    fprintf (sts->dbg,
			" <-- nominal location is ");
		    if (locn0 == (double)ilocn0)
			fprintf (sts->dbg, "here");
		    else if (locn0 > (double)ilocn0)
			fprintf (sts->dbg, "%.6g down from here",
				locn0 - ilocn0);
		    else
			fprintf (sts->dbg, "%.6g up from here",
				-(locn0 - ilocn0));
		}
		if (ii == ipeak) {
		    fprintf (sts->dbg,
			" <-- edge is ");
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

	return (status);
}
