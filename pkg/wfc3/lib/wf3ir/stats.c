# include <math.h>
# include <float.h>
# include <stdio.h>
# include <stdlib.h>
#include "hstcal.h"
# include "hstio.h"
# include "wf3.h"
# include "trlbuf.h"

/* STATS: Compute mean, median, mode, stdv, min, max of unflagged pixels
** in a SCI image. */

int stats (SingleNicmosGroup *in, int x1, int x2, int y1, int y2,
	     float low, float high, short dqmask,
	     float *mean, float *median, float *mode, float *stdv, float *min,
	     float *max) {

/* Arguments:
**	in	i: input image
**	x1, x2	i: column limits for statistics calculation
**	y1, y2	i: row limits for statistics calculation
**	low     i: lower data rejection threshold
**	high	i: upper data rejection threshold
**	dqmask	i: mask value for rejecting data
**	mean	o: mean value of unrejected pixels
**	median	o: median value of unrejected pixels
**	mode	o: modal value of unrejected pixels
**	stdv	o: standard deviation of unrejected pixels
**	min	o: minimum unrejected pixel value
**	max	o: maximum unrejected pixel value
*/

	/* Local variables */
	int i, j;		/* loop indexes */
	int npix;		/* number of unrejected pixels */
	float *arr;		/* array of good pixel values */
	double val;		/* SCI pixel value */
	double sumx, sumx2;	/* sum and sum**2 of unflagged pixel values */

	/* Function definitions */
	int sort (float *, int);
	float findMode (float *, int);
	float findMedian (float *, int);

	/* Initialize the counters and results */
	sumx    = 0;
	sumx2   = 0;
	*mean   = 0;
	*median = 0;
	*mode	= 0;
	*stdv   = 0;
	*min    = 0;
	*max    = 0;

	/* Allocate memory for the temporary array */
	arr = (float *) calloc(in->sci.data.nx*in->sci.data.ny, sizeof(float));
	if (arr == NULL) {
	    sprintf (MsgText, "Memory allocation failure in stats");
	    trlerror (MsgText);
	    return (1);
	}

	/* Loop through the requested image section, computing the
	** sum, sum of the squares, min, max and number of unflagged
	** pixels. */
	npix = 0;
	for (j=y1; j<=y2; j++) {
	     for (i=x1; i<=x2; i++) {
		  val = Pix(in->sci.data,i,j);
		  if (!(dqmask & DQPix(in->dq.data,i,j)) &&
		      val > low && val < high) {
		      arr[npix] = val;
		      sumx  += val;
		      sumx2 += val*val;
		      npix++;
		      if (npix == 1) {
			  *min = val;
			  *max = val;
		      } else {
			  if (val < *min) *min = val;
			  if (val > *max) *max = val;
		      }
		  }
	     }
	}

	if (npix == 0) {
	    free (arr);
	    return (0);
	}

	/* Compute the mean */
	*mean = sumx / npix;

	/* Sort the array of good values */
	if (sort(arr-1, npix)) {
	    free (arr);
	    return (1);
	}

	/* Find the mode */
	*mode = findMode (arr, npix);

	/* Find the median */
	*median = findMedian (arr, npix);
	
	/* Compute the standard deviation */
	if (npix > 1) {
	    *stdv = npix/(npix-1.) * (sumx2/npix - (*mean)*(*mean));
	    if (*stdv >= 0)
		*stdv = sqrt (*stdv);
	    else
		*stdv = 0.0;
	}

	/* Free local memory */
	free (arr);

	/* Successful return */
	return (0);

}
 
/* QSTATS: Compute mean, min, max, and snr of unflagged pixels
** in a SCI image. */

int qstats (SingleNicmosGroup *in, int x1, int x2, int y1, int y2,
	    short dqmask, int *numgood, float *mean, float *min, float *max,
	    float *snrmean, float *snrmin, float *snrmax, float *errmean,
	    float *errmin, float *errmax) {

/* Arguments:
**	in	i: input image
**	x1, x2	i: column limits for statistics calculation
**	y1, y2	i: row limits for statistics calculation
**	dqmask	i: mask value for rejecting data
**	numgood o: number of unrejected pixels
**	mean	o: mean value of unrejected pixels
**	min	o: minimum unrejected pixel value
**	max	o: maximum unrejected pixel value
**	snrmean o: mean sci/err pixel value
**	snrmin	o: minimum sci/err pixel value
**	snrmax	o: maximum sci/err pixel value
**	errmean	o: mean value of unrejected err image values
**	errmin	o: minimum unrejected err image value
**	errmax	o: maximum unrejected err image value
*/

	/* Local variables */
	int i, j;		/* loop indexes */
	int npix;		/* number of unrejected pixels */
	int num_bad_err;
	int dimx, dimy;
	int area;
	double value;		/* SCI pixel value */
	double valsum;
	double errval;
	double errsum;
	double snr;
	double snrsum;

	/* Initialize the counters and results */
	npix     = 0;
	num_bad_err = 0;
	valsum   = 0.;
	errsum   = 0.;
	snrsum   = 0.;
	dimx = x2 - x1 + 1;
	dimy = y2 - y1 + 1;

	/* Loop through the requested image section, computing the
	** statistical quantities. */
	npix = 0;
	for (j=y1; j<=y2; j++) {
	     for (i=x1; i<=x2; i++) {

		  if (!(dqmask & DQPix(in->dq.data,i,j))) {

		      value  = Pix(in->sci.data,i,j);
		      errval = Pix(in->err.data,i,j);

		      if (errval <= 0.) {
			  num_bad_err++;
			  continue;
		      } else {
			  snr = value / errval;
		      }

		      if (npix < 1) {
			  valsum = value;
			 *min    = (float) value;
			 *max    = (float) value;
			  errsum = errval;
			 *errmin = (float) errval;
			 *errmax = (float) errval;
			  snrsum = snr;
			 *snrmin = (float) snr;
			 *snrmax = (float) snr;
			  npix++;
		      } else {
			  valsum += value;
			  errsum += errval;
			  snrsum += snr;
			  if (value < *min) 
			      *min = value;
			  else if (value > *max)
			      *max = value;
			  if (errval < *errmin) 
			      *errmin = errval;
			  else if (errval > *errmax)
			      *errmax = errval;
			  if (snr < *snrmin) 
			      *snrmin = snr;
			  else if (snr > *snrmax)
			      *snrmax = snr;
			  npix++;
		      }
		  }
	     }
	}

	*numgood = npix;

	if (npix > 0) {
	    *mean    = (float) (valsum / (double) npix);
	    *errmean = (float) (errsum / (double) npix);
	    *snrmean = (float) (snrsum / (double) npix);
	} else {
	    area = dimy * dimx;
	    if (area == 0) {
		trlwarn ("Output image size is zero.");
	    } else if (num_bad_err > 0) {
		if (num_bad_err == area) {
		    trlwarn ("No ERR values > 0.");
		} else {
		    trlwarn (
			"All output pixels either flagged as bad or ERR <= 0.");
		}
	    } else {
		trlwarn ("All output pixels flagged as bad.");
	    }
	}

	/* Successful return */
	return (0);

}
 

# define BINF 0.10
# define MAX(a,b) (a > b ? a : b)

float findMode (float *arr, int npts) {

	int i;
	int ibinh, ipos;
	float diff, dmin;
	float mode;

    ipos=1;
    
	/* Check for trivial cases */
	if (npts == 0)
	    mode = 0.0;

        else if (npts == 1)
            mode = arr[0];

        else if (npts == 2)
            mode = 0.5 * (arr[0]+arr[1]);

	/* Compute mode */
        else {

	    /* Form histogram of values */
	    ibinh = MAX(1, (int)(BINF*(float)npts));

	    /* Calculate the reciprocal of the density and find the max value */
	    dmin = FLT_MAX;
	    for (i = ibinh; i < npts-ibinh; i++) {
		 diff = arr[i+ibinh] - arr[i-ibinh];
		 if (diff < dmin) {
		     dmin = diff;
		     ipos = i;
		 }
	    }

            mode = arr[ipos];
	}

	return (mode);

}

# undef BINF
# undef MAX

float findMedian (float *arr, int npts) {

	float median;

	/* Check for trivial cases */
	if (npts == 0)
	    median = 0.0;

	else if (npts == 1)
	    median = arr[0];

	else if ((npts % 2) == 0)
	    median = 0.5 * (arr[npts/2-1] + arr[npts/2]);
	else
	    median = arr[npts/2];
	
	return (median);
}

float findStdv (float *arr, int npts) {

	int i;
	float sum2, mean, stdv;

	if (npts <= 1)
	    stdv = 0.0;

	else {

	    sum2 = 0.0; mean = 0.0;
	    for (i = 0; i < npts; i++) {
		 mean += arr[i];
		 sum2 += arr[i]*arr[i];
	    }

	    mean = mean / npts;

	    stdv = npts/(npts-1.) * (sum2/npts - mean*mean);
	    if (stdv >= 0)
		stdv = sqrt (stdv);
	    else
		stdv = 0.0;
	}

	return (stdv);
}

