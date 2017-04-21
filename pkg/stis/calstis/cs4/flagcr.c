# include <stdio.h>
# include <stdlib.h>	/* malloc */
# include <string.h>
# include <math.h>	/* fabs, qsort */

# include "hstio.h"
# include "stis.h"
# include "calstis4.h"
# include "err.h"
# include "stisdq.h"	/* for DATAMASKED and DATAREJECT */

static int FlagLine (StisInfo4 *,
	float *, short *, int, float *, float *, float *);
static int CleanMean (float *, int, float *,
	double, double, double *, double *);
static void SaveDQ (SingleGroup *, char *, int);

/* This routine finds cosmic rays, and it flags them in the DQ extension
   with the DATAREJECT bit.  The data are assumed to have been 2-D
   rectified, so that when we examine one column at a time, we will be
   looking at homogeneous data, except that some pixels will correspond
   to points on the slit, and others will be beyond the ends of the slit
   or behind an occulting bar.

   Phil Hodge, 1998 Dec 11:
	Use parameters from WCP table; change calling sequences of
	FlagLine and CleanMean.  Include SaveDQ to write the DQ extension
	to a new FITS file if a debug file was specified; add extver to
	the calling sequence of FlagCR to support this option.

   Phil Hodge, 2002 Mar 4:
	In CleanMean, in the section for accumulating the sum of values and
	the sum of squares, assign each value to a double precision variable
	before squaring it.
	In FlagLine, when checking for dq = DATAMASKED, ignore any other
	data quality bits that are not included in sdqflags.
*/

int FlagCR (StisInfo4 *sts, SingleGroup *in, int extver) {

/* arguments:
StisInfo4 *sts     i: calibration switches and info
SingleGroup *in    io: input data (DQ array can be modified)
*/

	/* scratch space */
	float *sci;		/* SCI data */
	float *illum;		/* SCI data, through the slit */
	float *masked;		/* SCI data, missing the slit */
	float *absdiff;		/* fabs (sci - median) */
	short *dq;		/* data quality; flag CRs here */

	int nelem;		/* size of arrays */
	int i, j;

	if (sts->detector != CCD_DETECTOR)
	    return (0);

	if (sts->dispaxis == 1)
	    nelem = in->sci.data.ny;
	else if (sts->dispaxis == 2)
	    nelem = in->sci.data.nx;
	else
	    return (0);

	sci = malloc (nelem * sizeof (float));
	illum = malloc (nelem * sizeof (float));
	masked = malloc (nelem * sizeof (float));
	absdiff = malloc (nelem * sizeof (float));
	dq = malloc (nelem * sizeof (short));

	if (sci == NULL || illum == NULL || masked == NULL ||
		absdiff == NULL || dq == NULL) {
	    printf ("ERROR    (FlagCR) can't allocate memory.\n");
	    return (OUT_OF_MEMORY);
	}

	/* For each column (row), copy the data and data quality to scratch,
	   search for cosmic rays, and copy the DQ flags back into the column.
	*/
	if (sts->dispaxis == 1) {
	    for (i = 0;  i < in->sci.data.nx;  i++) {	/* for each column */
		for (j = 0;  j < nelem;  j++) {
		    sci[j] = Pix (in->sci.data, i, j);	/* note index j */
		    dq[j] = DQPix (in->dq.data, i, j);
		}
		if (FlagLine (sts, sci, dq, nelem,
				illum, masked, absdiff) > 0) {
		    for (j = 0;  j < nelem;  j++)
			DQSetPix (in->dq.data, i, j, dq[j]);
		}
	    }
	} else if (sts->dispaxis == 2) {
	    for (j = 0;  j < in->sci.data.ny;  j++) {	/* for each row */
		for (i = 0;  i < nelem;  i++) {
		    sci[i] = Pix (in->sci.data, i, j);	/* note index i */
		    dq[i] = DQPix (in->dq.data, i, j);
		}
		if (FlagLine (sts, sci, dq, nelem,
				illum, masked, absdiff) > 0) {
		    for (i = 0;  i < nelem;  i++)
			DQSetPix (in->dq.data, i, j, dq[i]);
		}
	    }
	}

	/* save DQ extension in a FITS file */
	if (sts->dbg != NULL)
	    SaveDQ (in, sts->dbgfile, extver);

	free (dq);
	free (absdiff);
	free (masked);
	free (illum);
	free (sci);

	return (0);
}

/* Find cosmic rays in current 1-D array.

   In the illuminated region (through the slit), values more than three
   sigma greater than the mean will be flagged as cosmic ray hits.  In the
   region off the slit or behind an occulting bar, however, in order to be
   flagged as a cosmic ray, the value must not only be more than three sigma
   high, but it must also be farther than one sigma greater than the mean of
   the illuminated data.  The latter test is needed because of the slop in
   the mode select mechanism; just because the pixel is flagged as masked
   doesn't mean it really is off the slit or behind a bar.

   The function value is the number of cosmic rays flagged.
*/

static int FlagLine (StisInfo4 *sts,
		float *sci, short *dq, int nelem,
		float *illum, float *masked, float *absdiff) {

/* arguments:
float sci[nelem]     i: science data for current column
short dq[nelem]      io: data quality flags; updated if cosmic ray found
int nelem            i: size of arrays
	scratch space:
float illum[nelem]   io: data through the slit
float masked[nelem]  io: data which missed the slit, or behind a bar
float absdiff[nelem] io: for absolute values of differences
*/

	/* mean and standard deviation for illuminated and masked regions */
	double i_mean, i_stddev, m_mean, m_stddev;
	int n_illum = 0;	/* number of illuminated pixels */
	int n_masked = 0;	/* number of off-slit or behind bar */
	int ncr = 0;		/* number of cosmic rays found */
	int ni, nm;		/* number of non-rejected elements */
	int i;
	short sdqflags;		/* sdqflags with DATAMASKED excluded */
	short datamasked;	/* the DATAMASKED data quality value */

	/* Locally remove the bit that flags out of bounds or behind an
	   occulting bar.
	*/
	datamasked = (short)DATAMASKED;
	sdqflags = sts->sdqflags & ~datamasked;

	/* Separate the SCI values into two arrays, depending on whether
	   the point corresponds to through the slit or off the slit.
	*/
	for (i = 0;  i < nelem;  i++) {
	    if (!(dq[i] & sdqflags)) {		/* no serious DQ flag */
		if (dq[i] & datamasked) {	/* off the slit */
		    masked[n_masked] = sci[i];
		    n_masked++;
		} else {			/* through the slit */
		    illum[n_illum] = sci[i];
		    n_illum++;
		}
	    }					/* else ignore it */
	}

	/* Find the mean and standard deviation of each array, ignorning
	   outliers.
	*/
	ni = CleanMean (illum, n_illum, absdiff,
			sts->min_mad, sts->mad_reject, &i_mean, &i_stddev);
	nm = CleanMean (masked, n_masked, absdiff,
			sts->min_mad, sts->mad_reject, &m_mean, &m_stddev);

	if (ni == 0 && nm == 0)
	    return (0);			/* no data means no cosmic rays */
	else if (ni == 0)
	    i_stddev = m_stddev;	/* for want of a better value */
	else if (nm == 0)
	    m_stddev = i_stddev;

	/* Now find cosmic rays (outliers in the positive direction)
	   and flag them with DATAREJECT in the dq array.
	*/
	for (i = 0;  i < nelem;  i++) {
	    if (dq[i] & datamasked) {
		if (sci[i] > m_mean + sts->nsigma_cr * m_stddev &&
		fabs (sci[i] - i_mean) > sts->nsigma_illum * i_stddev) {
		    dq[i] |= DATAREJECT;
		    ncr++;
		}
	    } else if (!(dq[i] & sts->sdqflags)) {
		if (sci[i] > i_mean + sts->nsigma_cr * i_stddev) {
		    dq[i] |= DATAREJECT;
		    ncr++;
		}
	    }
	}

	return (ncr);
}

/* Find the mean and standard deviation, ignoring outliers.
   The function value is the number of elements that were used
   to compute the mean and standard deviations, i.e. the number
   that were not rejected.

   If there are only two elements, the smaller of the two is taken as
   the mean.
*/

static int CleanMean (float *sci, int nelem, float *absdiff,
		double min_mad, double mad_reject,
		double *mean, double *stddev) {

/* arguments:
float sci[nelem]     io: science data for current column (will be sorted)
int nelem            i: size of arrays
float absdiff[nelem] io: scratch for absolute values of differences
double min_mad;      i: lower limit to median absolute deviation
double mad_reject;   i: rejection criterion for outliers
double mean          o: mean of sci array values
double stddev        o: standard deviation of sci array
*/

	float median;			/* median of sci */
	float mad;	/* median of abs value of deviations from median */
	double dsci;			/* = sci, but double precision */
	double sum = 0., sumsq = 0.;	/* for mean and std dev */
	double dsum;			/* = nsum */
	int nsum = 0;			/* number of good values */
	int i;
	int inplace = 1;		/* sort the array in-place */

	if (nelem < 3) {
	    if (nelem == 1) {
		*mean = sci[0];
		*stddev = 0.;
	    } else if (nelem == 2) {
		*mean = sci[0] < sci[1] ? sci[0] : sci[1];	/* min */
		*stddev = fabs (sci[0] - sci[1]);
	    } else {
		*mean = 0.;
		*stddev = 0.;
	    }
	    return (nelem);
	}

	median = MedianFloat (sci, nelem, inplace);

	/* Find the median of the absolute values of the deviations. */
	for (i = 0;  i < nelem;  i++)
	    absdiff[i] = fabs (sci[i] - median);
	mad = MedianFloat (absdiff, nelem, inplace);
	/* Set a lower limit to the median deviation. */
	mad = (mad < min_mad) ? min_mad : mad;

	/* Accumulate sums, ignoring outliers. */
	for (i = 0;  i < nelem;  i++) {
	    if (fabs (sci[i] - median) < mad_reject * mad) {
		dsci = sci[i];
		sum += dsci;
		sumsq += (dsci * dsci);
		nsum++;
	    }
	}

	if (nsum < 1) {			/* shouldn't be possible */
	    *mean = 0.;
	    *stddev = 0.;
	} else if (nsum == 1) {
	    *mean = sum;
	    *stddev = 0.;
	} else {
	    dsum = (double)nsum;
	    *mean = sum / dsum;
	    *stddev = sqrt (dsum / (dsum-1.) *
			(sumsq / dsum - (*mean * *mean)));
	}

	return (nsum);
}

/* This routine is for saving the DQ extension in a FITS file with
   the same root name as the dbgfile.  
   If an error occurs in this function, a warning will be printed,
   and it will otherwise be ignored.
*/

static void SaveDQ (SingleGroup *in, char *dbgfile, int extver) {

	char *dqname;
	int option = 0;

	if ((dqname = malloc ((STIS_LINE+1) * sizeof(char))) == NULL) {
	    printf ("Warning  Out of memory in SaveDQ (debug option).\n");
	    return;
	}
	/* Construct FITS file name. */
	strcpy (dqname, dbgfile);
	strcat (dqname, ".fits");

	/* Write the DQ extension to the debug FITS file. */
	putDQ (dqname, extver, &(in->dq), option);
	if (hstio_err()) {
	    printf (
	"Warning  SaveDQ couldn't write DQ extension (debug option); \\\n");
	    printf ("Warning  %s\n", hstio_errmsg());
	    clear_hstioerr();
	}

	free (dqname);
}
