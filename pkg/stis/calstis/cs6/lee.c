# include <stdio.h>
# include <math.h>

# include "xtables.h"
# include "hstio.h"

# include "stis.h"
# include "err.h"
# include "stisdq.h"
# include "calstis6.h"

static void AddPixel (int, RowContents *, short, double *, double *, int *);
static void RemovePixel (int, RowContents *, short, double *, double *, int *);
static double LSF (double, double, double, double, int);


/*  
   Apply a Lee statistics filter to the background array and recompute 
   the net array. This is the last step in the HS scattered light 
   correction algorithm.



   Revision history:
   ----------------
   19 Jan 00  -  Implemented (I.Busko)
*/

int Lee (StisInfo6 *sts, RowContents *row_cont) {

/* arguments:
StisInfo6 *sts         i:  calibration switches and info
RowContents *row_cont  i:  output row arrays
*/
	short	mask;
	int	ipix, npix, lfilt;
	double	sum, sum2, *hold;

	mask = sts->sdqflags & DETECTORPROB;

	lfilt = sts->lfilter / 2;

	/* Alloc holding buffer. */
	if ((hold = (double *) calloc (row_cont->npts, sizeof (double))) == 
            NULL)
	    return (OUT_OF_MEMORY);

	/* Initialize sums. */
	sum  = 0.0;
	sum2 = 0.0;
	npix = 0;
	for (ipix = 0; ipix <= lfilt; ipix++)
	    AddPixel (ipix, row_cont, mask, &sum, &sum2, &npix);

	/* Compute first filtered value. */
	hold[0] = LSF (sum, sum2, row_cont->back[0], row_cont->error[ipix],
                       npix);

	/* Compute filtered values at first edge. */
	for (ipix = 1; ipix <= lfilt; ipix++) {
	    AddPixel (ipix + lfilt, row_cont, mask, &sum, &sum2, &npix);
	    hold[ipix] = LSF (sum, sum2, row_cont->back[ipix], 
                              row_cont->error[ipix], npix);
	}

	/* Compute filtered values at "internal" region. */
	for (ipix = lfilt+1; ipix <= row_cont->npts-1-lfilt; ipix++) {
	    AddPixel    (ipix + lfilt,     row_cont, mask, &sum, &sum2, &npix);
	    RemovePixel (ipix - lfilt - 1, row_cont, mask, &sum, &sum2, &npix);
	    hold[ipix] = LSF (sum, sum2, row_cont->back[ipix], 
                              row_cont->error[ipix], npix);
	}

	/* Compute filtered values at second edge. */
	for (ipix = row_cont->npts-lfilt; ipix <row_cont->npts; ipix++) {
	    RemovePixel (ipix - lfilt - 1, row_cont, mask, &sum, &sum2, &npix);
	    hold[ipix] = LSF (sum, sum2, row_cont->back[ipix], 
                              row_cont->error[ipix], npix);
	}

	/* Recompute background and net arrays. */
	for (ipix = 0; ipix < row_cont->npts; ipix++) {
	    row_cont->back[ipix] = hold[ipix];
	    row_cont->net[ipix] = row_cont->gross[ipix] - row_cont->back[ipix];
	}

	free (hold);

	return (0);
}


/* 
   Add pixel value in accumulators, taking care of masked pixels.
*/

static void AddPixel (int ipix, RowContents *row_cont, short mask, 
                      double *sum, double *sum2, int *npix) {
	double pix;

	if (!(row_cont->dq[ipix] & mask)) {
	    pix = row_cont->back[ipix];
	    *sum  += pix;
	    *sum2 += pix * pix;
	    (*npix)++;
	}
}


/* 
   Remove pixel value from accumulators, taking care of masked pixels.
*/

static void RemovePixel (int ipix, RowContents *row_cont, short mask,
                         double *sum, double *sum2, int *npix) {
	double pix;

	if (!(row_cont->dq[ipix] & mask)) {
	    pix = row_cont->back[ipix];
	    *sum  -= pix;
	    *sum2 -= pix * pix;
	    (*npix)--;
	}
}


/*
   Lee statistical filter
*/

static double LSF (double sum, double sum2, double value, double error,
                   int npix) {

	double mean, var, varn, beta;

	if (npix > 2) {
	    mean = sum / npix;
	    var  = (sum2 - sum * sum / npix) / (npix - 1);
	    varn = error * error;
	    if (var > 0.0)
	        beta = (var - varn) / var;
	    else
	        return (0.0);
	    beta = (beta > 0.0) ? beta : 0.0;
	    return (beta * value + (1.0 - beta) * mean);
	} else
	    return (0.0);
}
