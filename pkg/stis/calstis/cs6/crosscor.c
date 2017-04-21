# include <stdio.h>
# include <stdlib.h>
# include <math.h>	/* fabs */
# include <float.h>

# include "xtables.h"
# include "hstio.h"

# include "stis.h"
# include "stisdef.h"
# include "calstis6.h"
# include "err.h"
# include "stisdq.h"

/*
   Cross-correlation routine.

   This is in fact a peak-finding routine that acts in the 1-D signal
   that results from co-adding the 2-D signal along the dispersion (A1)
   direction. The scan along this direction takes into account the
   true shape of the spectrum described by its nominal A2 position
   as a function of A1 position. This description is conveyed by a
   trace list that should span the interval in A2 covered by the
   search range.

   WARNING: changes in this code must be perhaps paralleled with the
   same changes in DefineBackRegions function. This function executes
   a similar algorithm and there are lots of code duplication (IB, 01/10/00).




   Revision history:
   ----------------
   20 Feb 97  -  Implemented (I.Busko)
   10 Apr 97  -  Changes after code review (IB):
                 - fixed size of allocated accumulator memory.
   02 May 97  -  Fixed parabolic interpolation formula (IB)
   17 Jul 97  -  Switched to unweighted sum (IB)
   17 Jul 97  -  Check for maxsearch = 0 (IB)
   24 Sep 97  -  Interpolate in trace table for each co-add line (IB).
   03 Jun 98  -  Store croscorr flux values (IB)
   17 Dec 98  -  Avoid geocoronal Lya (IB)
   07 Nov 00  -  Avoid too small weights at image top and bottom lines (IB)
   07 Feb 06  -  Add an argument to Interp2D (PEH)
*/

int CrossCorr (StisInfo6 *sts, SpTrace *trc, SingleGroup *in,
               FloatHdrData *ssgx, FloatHdrData *ssgy, int maxsearch,
               int avoid1, int avoid2) {

/* arguments:
StisInfo6 *sts      io: calibration switches and info
SpTrace *trc;       i:  full list of spectrum traces
SingleGroup *in	    i:  input image
FloatHdrData ssgx;  i:  small-scale distortion in X (not used)
FloatHdrData ssgy;  i:  small-scale distortion in Y (not used)
int maxsearch       i:  maximum range for cross correlation
int avoid1, avoid2; i:  Lya region to avoid (in physical pixels)
*/

	int status;

        SpTrace *trace_y;       /* interpolated spectrum trace */
	int	ipix;		/* index of image pixel in the A1 direction */
	int	rpix;		/* index of reference pixel in the A1 dir. */
	int	jpix;		/* index of crosscor array element */
	double	*sum;		/* accumulators */
	double	*wei;
	double	y_nom;		/* nominal coordinate in the A2 direction */
	double	iy_nom;		/* above quantity in image pixel units */
	double	y_trc;		/* trace-corrected coordinate */
	double	low_end; 	/* nominal endpoints of the crosscor function*/
	double	high_end;
	float	oSci, oErr;	/* interpolated values from image array */
	short	oDQ;
	double	hold;
	double	hwei;
	int	i, imax;

	int InterpTrace6 (SpTrace **, double, SpTrace **);
	void FreeTrace6 (SpTrace **);

	if (maxsearch == 0) {
	    sts->crscroff  = 0.0;
	    sts->cc_peakf  = 0.0;
	    sts->cc_peak2f = 0.0;
	    sts->cc_peak3f = 0.0;

	    return (0);
	}

	/* Alloc and clear accumulator memory. */
	if ((sum = (double *) calloc (2*(maxsearch)+1, sizeof(double)))== NULL)
	    return (OUT_OF_MEMORY);
	if ((wei = (double *) calloc (2*(maxsearch)+1, sizeof(double)))== NULL)
	    return (OUT_OF_MEMORY);

	/* Compute nominal endpoints of crosscorr range. */
	low_end  = sts->nm_a2center - (double)maxsearch;
	high_end = sts->nm_a2center + (double)maxsearch;

	/* Loop over the crosscor range in 1 pixel steps. */
	jpix = 0;
	trace_y = NULL;

	for (y_nom = low_end; y_nom <= high_end; y_nom += 1.0) {

	    /* Interpolate in trace table. */
	    if ((status = InterpTrace6 (&trc, y_nom, &trace_y)))
	        return (status);

	    /* Loop over image pixels in the A1 direction. */
	    for (ipix = 0; ipix < in->sci.data.nx ; ipix++) {

	        /* Translate image pixel index into reference pixel index. */
	        rpix = (int)((ipix - sts->ltv[0]) / sts->ltm[0]);

	        /* Add trace curve to get actual position of trace. */
	        if (rpix > 0 && rpix < trace_y->nelem)
	            y_trc = y_nom + trace_y->a2displ[rpix];
	         else
	            y_trc = y_nom;

	        /* Translate back into image pixel index. */
	        iy_nom = y_trc * sts->ltm[1] + sts->ltv[1];

		/* Interpolate along column, checking for out of bounds. */
		Interp2D (in, sts->cc_sdqflags, (double)ipix, iy_nom, 1.0,
	                 WGT_VARIANCE, &oSci, &oErr, &oDQ);
	        hwei = 1.0;

	        /* Avoid geocoronal Lya region. */
	        if (avoid1 != avoid2) {
	            if (ipix > avoid1 && ipix < avoid2) {
	                oSci = 0.0;
	                hwei = 0.0;
	            }
	        }

	        /* Update sum and weight, discarding flagged data. */
	        if ((oDQ != DETECTORPROB) && !(oDQ & sts->cc_sdqflags)) {
		    sum[jpix] += oSci;
		    wei[jpix] += hwei;
	        }
	    }

	    /* Bump crosscor function index. */
	    jpix++;
	}

	/* Free interpolated trace structure. */
	FreeTrace6 (&trace_y);

	/* Compute weighted average and locate maximum and minimum. */
	sts->cc_peakf = - DBL_MAX;
	sts->cc_minf  =   DBL_MAX;
	imax = 0;
	for (i = 0; i < jpix; i++) {
            /* Modified on 11/7/00 to account for too small weights at
               the image top and bottom. They may combine with non-rejected
               cosmic rays and throw off the detected signal peak.
            */
	    if (wei[i] > 0.0) {
	    /*if ((sts->echelle && (wei[i] > 0.0)) ||
               ((!sts->echelle) && (wei[i] > (0.0 *
               (float)(in->sci.data.nx))))) {*/
	        sum[i] /= wei[i];
	        if (sum[i] > sts->cc_peakf) {
	            sts->cc_peakf = sum[i];
	            imax = i;
	        }
	        if (sum[i] < sts->cc_minf)
	            sts->cc_minf = sum[i];
	    } else {
	        sum[i] = 0.0;
	    }

	    /* This is used later in the parabolic fit. */
	    wei[i] = (double)i;
	}

	sts->cc_avef = 0.0;

	if (imax > 0 && imax < jpix-1) {

	    /* Compute average around peak. */
	    sts->cc_avef = (sum[imax] + sum[imax+1] + sum[imax-1]) /
	                   (wei[imax] + wei[imax+1] + wei[imax-1]);

	    /* See if peak satisfies threshold criterion. */
	    sts->cc_peakf  = sum[imax];
	    sts->cc_peak2f = sum[imax+1];
	    sts->cc_peak3f = sum[imax-1];
	    if (sts->cc_peakf  < sts->cc_thresh ||
	        sts->cc_peak2f < sts->cc_thresh ||
	        sts->cc_peak3f < sts->cc_thresh) {
	        sts->crscroff  = 0.0;
	        sts->cc_peakf  = 0.0;
	        sts->cc_peak2f = 0.0;
	        sts->cc_peak3f = 0.0;
	        free (sum);
	        free (wei);
	        return (1);
	    }

	    /* Do parabolic interpolation. It is an error if maximum is
               located at any one of the extremes of the search range,
               since we might be missing the maximum position by an unknown
               amount. We definitely need better acceptance criteria for
               deciding when the peak finding fails.
            */
	    hold = 2.0 * sum[imax] - sum[imax+1] - sum[imax-1];
	    if (hold != 0.0) {
	        sts->crscroff = wei[imax] - 0.5*(sum[imax-1] - sum[imax+1])/
                                hold;
	        sts->crscroff -= (double)maxsearch;
	        if (fabs (sts->crscroff) <= (double)maxsearch) {
	            sts->cc_peakf  = 0.0;
	            sts->cc_peak2f = 0.0;
	            sts->cc_peak3f = 0.0;
	            free (sum);
	            free (wei);
	            return (0);
	        } else {
	            sts->crscroff  = 0.0;
	            sts->cc_peakf  = 0.0;
	            sts->cc_peak2f = 0.0;
	            sts->cc_peak3f = 0.0;
	            free (sum);
	            free (wei);
	            return (1);
	        }
	    } else {
	        sts->crscroff  = 0.0;
	        sts->cc_peakf  = 0.0;
	        sts->cc_peak2f = 0.0;
	        sts->cc_peak3f = 0.0;
	        free (sum);
	        free (wei);
	        return (1);
	    }
	} else {
	    sts->crscroff  = 0.0;
	    sts->cc_peakf  = 0.0;
	    sts->cc_peak2f = 0.0;
	    sts->cc_peak3f = 0.0;
	    free (sum);
	    free (wei);
	    return (1);
	}
}
