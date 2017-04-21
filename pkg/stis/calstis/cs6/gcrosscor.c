# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <time.h>
# include <math.h>		/* fabs */

# include "xtables.h"
# include "hstio.h"

# include "stis.h"
# include "calstis6.h"
# include "err.h"
# include "stisdq.h"

static double average (double *, int *, int, double, double *);
static int OpenSGeo1 (char *, FloatHdrData *, FloatHdrData *);
static void ReportCrossCorr (StisInfo6 *, int, double, double);

/*
   Compute global crosscor offset for all orders.

   The global offset is a sigma-clipped average of the individual
   crosscor offsets in each spectral order.

   There is a lot of code duplication in this module, since it was
   basically derived from do1dx.c. This is due to the fact the global
   offset mode wasn't part of the original requirements, but was
   added much later into calstis6.




   Revision history:
   ----------------
   22 Jun 98  -  Implemented (I.Busko)
   20 Jul 98  -  Always get aperture description (because of A2 offset, IB)
   20 Nov 98  -  Report results when in verbose mode (IB)
   06 Feb 02  -  Gets reference order A2 position for MSM/blaze correction (IB)
   21 Apr 05  -  Initialize slit.gac_allocated to 0 (PEH)
*/

int GCrossCorr (StisInfo6 *sts, SingleGroup *in, int minorder,
                int maxorder, double *gcrosscor, int mref, double *yoff) {

/* arguments:
StisInfo6    *sts                 i: calibration switches and info
SingleGroup  *in                  i: input image
int          minorder, maxorder   i: minimum and maximum spectral orders
double       *gcrosscor           o: global croscor offset
int          mref                 i: reference order number
double       *yoff                o: reference order croscor offset
*/

	int status;

	XtractInfo *extract;	/* list of extraction parameter records */
	XtractInfo *extract_o;	/* one record extracted from extract */
	SpTrace *trace;		/* list of spectrum traces */
        SpTrace *trace_y;       /* spectrum trace interpolated at y */
        /* CoordInfo *coords;  */ /* list of plate scale records */
        /* CoordInfo *coord_o; */ /* one record extracted from coords */
	ApInfo slit;		/* description of slit */
	FloatHdrData ssgx;	/* small-scale distortion in X */
	FloatHdrData ssgy;	/* small-scale distortion in Y */

	double sigma;		/* residual crosscorrelation stddev. */
 	int sporder;		/* current sporder */
	int norder;		/* sporder counter */
	int maxsearch;		/* maximum range in cross correlation */
	int a2phys;		/* nominal A2CENTER in image coordinates */
	int trc_status;         /* flag for trace status */
	int success;		/* at least one order successfully extrac. ? */
	int skipping;		/* skipping this order ? */
	int ccstatus;		/* crosscor failed ? */

	int CrossCorr (StisInfo6 *, SpTrace *, SingleGroup *, FloatHdrData *,
                       FloatHdrData *, int, int, int);
	void AddOffsets6 (StisInfo6 *, ApInfo *);
	int GetTrace6 (StisInfo6 *, int, SpTrace **);
	void FreeTrace6 (SpTrace **);
	void FreeThroughput6 (ApInfo *);
/*      int GetSDC6 (StisInfo6 *, CoordInfo **, int *, int *);          */
	int GetApDes6 (StisInfo6 *, ApInfo *);
	int GetApThr6 (StisInfo6 *, ApInfo *);
	int GetGrpInfo6 (StisInfo6 *, Hdr *);
	int GetXtract (StisInfo6 *, XtractInfo **, int *, int *);
	int InterpTrace6 (SpTrace **, double, SpTrace **);
	int ReturnXtract (XtractInfo **, int, XtractInfo **);
/*      int ReturnCoord6 (CoordInfo **, int, CoordInfo **);      */
	void FreeXtract (XtractInfo **);

	/* If zero range was specified in command line, skip analysis. */
	if (sts->maxsearch != NO_RANGE) {
	    if (sts->maxsearch == 0) {
	        *gcrosscor = 0.0;
	        return (0);
	    }
	}

	/* Set flags to indicate that memory has not been allocated yet. */
	slit.allocated  = 0;
	slit.gac_allocated  = 0;
	extract     = NULL;
	extract_o   = NULL;
	trace       = NULL;
	trace_y     = NULL;

	if (sts->verbose == 1 || sts->verbose == 2)
	    printf ("         Begin offset analysis.\n");

     	/* Get aperture description. */
	if ((status = GetApDes6 (sts, &slit)))
	    return (status);

	/* Get aperture throughput info. */
	if (sts->fluxcorr == PERFORM) {
	    if ((status = GetApThr6 (sts, &slit)))
		return (status);
	}

	/* Read the small-scale distortion data into memory. */
	if (sts->sgeocorr == PERFORM) {
	    if ((status = OpenSGeo1 (sts->sdstfile.name, &ssgx, &ssgy)))
		return (status);
	}

	/* Get 1-D extraction parameters for all spectral orders.
           The extraction table is the main driver that controls the
           scanning of spectral orders in echelle data.
        */
	if ((status = GetXtract (sts, &extract, &minorder, &maxorder)))
	    return (status);
        if (sts->x1d_o == DUMMY) {
	    printf ("ERROR    DUMMY pedigree entry in extraction table.\n");
	    return (status);
	}

	sts->echelle = (maxorder > 1);

	/* Get keyword values from image extension header. */
	if ((status = GetGrpInfo6 (sts, &(in->sci.hdr))))
	    return (status);

	norder = 0;
	*yoff  = 0.0;

	/* Reset flags that control loop termination. */
	success  = 0;
	skipping = 0;

	/* This flag controls skipping by dummy reference entries. */
	sts->x1d_o = PERFORM;

	/* Do for each spectral order in input image. */
	for (sporder =  minorder;
             sporder <= maxorder;
             sporder++) {

	    /* If last order was skipped but at least one previous
               order was successfully extracted, assume orders are
               falling off from the detector and finish current IMSET.
            */
	    if (success == 1 && skipping == 1)
	        break;

	    /* Get the XtractInfo record for this order. */
	    if ((status = ReturnXtract (&extract, sporder, &extract_o)))
	        return (status);

            /* Extract the CoordInfo record for this order. This is
               used for now only to get the appropriate plate scale
               in the A2 direction (CDELT2). The plate scale is used
               to scale any eventual POSTARG2 value to reference pixel
               units. Plate scale may be a function of spectral order,
               thuw we have to read it here and not as a global parameter.
            */
	    /*
            if ((status = ReturnCoord6 (&coords, sporder, &coord_o)))
                return (status);
	    */
	    /* Get spectrum trace info for this order. If no matching
               row in trace table, silently skip to next order. Notice
               that this implies opening and closing the trace table
               every time. This shouldn't be a burden in elapsed time
               if the extraction and trace table are reasonably matched.
                If any matching trace table entry is DUMMY, skip and warn.
            */
	    trc_status = GetTrace6 (sts, sporder, &trace);
	    if (trc_status == NO_TRACE) {
	        status = 0;
	        FreeTrace6 (&trace);
	        FreeXtract (&extract_o);
	        continue;
	    } else if (trc_status == ERROR_TRACE)
	        return (trc_status);
	    else if (sts->x1d_o == DUMMY) {
	        FreeTrace6 (&trace);
	        FreeXtract (&extract_o);
	        continue;
	    }

	    /* Compute total offset = "dither" + wavecal + aperture. */
	    AddOffsets6 (sts, &slit);

	    /* Compute nominal spectrum trace position. */
	    sts->nm_a2center = trace->a2center + sts->total_offset[1];

	    /* POS TARG correction in spatial direction. The plate
               scale read from table is in arcsec per reference pixel,
               thus the POS TARG correction can be directly added to
               the nominal center position at this point.
            */
/*	    sts->nm_a2center += sts->pos_targ[1] / coord_o->cdelt2;   */

	    /* Interpolate in the trace table. */
	    if ((status = InterpTrace6 (&trace, sts->nm_a2center, &trace_y)))
	        return (status);
	    sts->nm_a2center = trace_y->a2center;

	    /* Check if nominal center is within the data array
               boundaries. Skip order if off-limits. This should take
               care of both sub-arrays and wandering off-detector
               spectral orders.
            */
	    if (sts->echelle) {	   /* added by PEH 1997 Dec 12 in do1dx */
		a2phys = (int)(sts->nm_a2center * sts->ltm[1] +
                         sts->ltv[1]);
		if (a2phys < 0 || a2phys >= in->sci.data.ny) {
	            skipping = 1;
	            FreeTrace6 (&trace);
	            FreeTrace6 (&trace_y);
	            FreeXtract (&extract_o);
	            continue;
		} else
	            skipping = 0;
	    } else
	        skipping = 0;

	    /* Compute peak-finding function. If supplied, use search range
               from command line. Only perform peak-find if the range is
               meaningful, otherwise set offset to zero. Turn off Lya
               excision region.
            */
	    if (sts->maxsearch == NO_RANGE)
	        maxsearch = extract_o->maxsearch;
	    else
	        maxsearch = sts->maxsearch;
	    if (maxsearch != NO_RANGE)
	        ccstatus = CrossCorr (sts, trace, in, &ssgx, &ssgy, maxsearch,
                                      0, 0);
	    else
	        sts->crscroff = 0.0;

	    /* Store offset data into appropriate arrays for later use. */
	    sts->cc_off[norder]   = sts->crscroff;
	    sts->cc_spord[norder] = sporder;
	    sts->cc_rej[norder]   = ccstatus;
	    norder++;

	    /* Store data for reference order for MSM/blaze correction. */

	    if (sporder == mref && ccstatus == 0)
	        *yoff = sts->crscroff;

	    /* Free memory for this spectral order. */
	    FreeTrace6 (&trace);
	    FreeTrace6 (&trace_y);
	    FreeXtract (&extract_o);
	}

	/* Compute sigma-cleaned average offset (2 iterations) */
	sigma = 0.0;
	*gcrosscor = average (sts->cc_off, sts->cc_rej, norder, 2.0, &sigma);
	*gcrosscor = average (sts->cc_off, sts->cc_rej, norder, 2.0, &sigma);

	/* Report */
	if (sts->verbose == 2)
	    ReportCrossCorr (sts, norder, *gcrosscor, sigma);

	/* Free memory for all orders. */
	FreeThroughput6 (&slit);
	FreeXtract (&extract);
	if (sts->sgeocorr == PERFORM) {
	    freeFloatHdrData (&ssgx);
	    freeFloatHdrData (&ssgy);
	}

	if (sts->verbose == 1 || sts->verbose == 2)
	    printf ("         End offset analysis.\n");

	return (0);
}




/*
    Compute average, discarding flagged data points. On output,
    the flag array is updated to include the flagged data points.
*/

static double average (double *data, int *flag, int ndata, double rejf,
                       double *sigma) {

/* arguments:
double *data          i:  the data array
int    *flag          io: the flag array.
int    ndata          i:  the array size
double  rejf          i;  the k-sigma rejection factor
*/
	double av;		/* average  */
	double sig;		/* standard deviation */
	int i, j, np;

	/* Compute sums. */
	av  = 0.0;
	sig = 0.0;
	np  = 0;
	for (i =  0;  i < ndata; i++) {
	    if (flag[i] == 0) {
	        av  += data[i];
	        sig += data[i] * data[i];
	        np++;
	    }
	}

	/* Compute stats */
	if (np >= 2) {
	    sig = (sig - av * av / np) / (np - 1);
	    if (sig > 0.0)
	        sig = sqrt (sig);
	    else
	        sig = 0.0;
	    av /= np;
	    *sigma = sig;

	    /* flag deviant data */
	    for (i =  0;  i < ndata; i++) {
	        if (fabs (data[i] - av) > (sig * rejf))
	            flag[i] = 1;
	    }

	    /* Update stats */
	    av = 0.0;
	    sig = 0.0;
	    j = 0;
	    for (i =  0;  i < ndata; i++) {
	        if (flag[i] == 0) {
	            av  += data[i];
	            sig += data[i] * data[i];
	            j++;
	        }
	    }
	    if (j >= 2) {
	        sig = (sig - av * av / j) / (j - 1);
	        if (sig > 0.0)
	            sig = sqrt (sig);
	        else
	            sig = 0.0;
	        *sigma = sig;
	    } else
	        *sigma = 0.0;

	    return (av / j);
	} else
	    return (0.0);
}


/*
    Report results of global cross correlation analysis.

*/

static void ReportCrossCorr (StisInfo6 *sts, int norder, double gcc,
                             double sigma) {

	int i, c;

	printf ("\n         Sporder Rej.  Offset\n");

	for (i = 0; i < norder; i++) {
	    if (sts->cc_rej[i])
	        c = '*';
	    else
	        c = ' ';
	    printf ("            %d    %c    %g\n",
               sts->cc_spord[i], c, sts->cc_off[i]);
	}
	printf ("         Average       %g  (%g)\n\n", gcc, sigma);
}



/*
   Open the two images in a small-scale geometric correction reference file.

   The small-scale geometric correction is NOT yet implemented in calstis6.
*/

static int OpenSGeo1 (char *sdstfile, FloatHdrData *ssgx, FloatHdrData *ssgy) {

/* arguments:
char *sdstfile          i: name of file to open
FloatHdrData *ssgx      o: small-scale distortion in X
FloatHdrData *ssgy      o: small-scale distortion in Y
*/

	initFloatHdrData (ssgx);
	initFloatHdrData (ssgy);

	getFloatHD (sdstfile, "AXIS1", 1, ssgx);
	if (hstio_err())
	    return (OPEN_FAILED);
	getFloatHD (sdstfile, "AXIS2", 1, ssgy);
	if (hstio_err())
	    return (OPEN_FAILED);

	return (0);
}


