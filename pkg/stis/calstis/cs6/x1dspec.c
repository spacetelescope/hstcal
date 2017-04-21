# include <stdio.h>
# include <math.h>
# include <stdlib.h>

# include "xtables.h"
# include "hstio.h"

# include "stis.h"
# include "err.h"
# include "stisdq.h"
# include "calstis6.h"

# define BOX_LOWER	-1	/* position of current pixel in */
# define BOX_MID	0	/* the extraction box           */
# define BOX_UPPER	1

static double BoxTilt (XtractInfo *, int);
static void AccumPix (StisInfo6 *, SingleGroup *, int, int, double, double *,
                       double *, double *, double *, double *, double *,
                       double *, short *, int *, int, int, double *, double *,
                       short *, double, double *, short, double *, int);
static void AccPixUnw (StisInfo6 *, SingleGroup *, int, int, double, double *,
                       double *, double *, double *, double *, double *,
                       short *, int *, int, int, double *, short *, double,
                       double, double *, int);
static void AccPixOpt (StisInfo6 *, SingleGroup *, int, int, double *,
                       double *, double *, double *, double *, double *,
                       double *, short *, int *, int, int, double *,
                       double, double, double, double *, short, int);
static void CleanPixel (StisInfo6 *, SingleGroup *, int, int, double *,
                        int, double, SingleGroup *, int);
static void CleanPixels (StisInfo6 *, SingleGroup *, int, int, int,
                        double *, int, double, SingleGroup *, int);
static void ComputeBack (StisInfo6 *, int, double *, double *, double *, int);
static void ReplaceBack (StisInfo6 *, double, double, double *, double *);
static void UpdateProfile (float, double *, double *,short *, int, int);

/*
   Perform the actual 1-D extraction, including background subtraction
   and optional profile generation.




   Revision history:
   ----------------
   20 Feb 97  -  Implemented (I.Busko)
   10 Apr 97  -  Changes after code review (IB):
                 - replaced literals by PI constant.
   11 Apr 97  -  OPR 33787: replaced flagged pixels by interpolated value (IB).
   11 Apr 97  -  OPR 33790: flag extracted pixels that had more than 10%
                 of their flux affected by input flagged pixels (IB).
   16 Apr 97  -  Revised OPR 33787 and 33790: do not replace flagged pixels by
                 interpolated value, just omit them and set the flag (IB)
   17 Apr 97  -  Replaced scalar bktilt by a polynomial description.
                 NO tilted spectrum box !! (IB)
   17 Jul 97  -  Background is in counts/sec (IB)
   18 Aug 97  -  Turn off data quality flag checking. This counteracts the
                 behavior prescribed in OPRs 33787 and 33790 (IB)
      Aug 97  -  Background is total in aperture, not per pixel (IB)
   13 Apr 98  -  Replace debug switch by extrloc, remove debug file, output
                 extraction info to table (IB)
   23 Jun 98  -  Profile generator (IB)
   29 Jun 98  -  Rejection ranges in profile generator (IB)
   10 Jul 98  -  Error-of-the-mean was replaced by error-of-the-sum (IB)
   20 Jul 98  -  Normalizes error by exposure time (IB)
   22 Sep 98  -  Optimal extraction (IB)
   29 Sep 98  -  Moved sdqflags = 0 to another place (IB)
   14 Oct 98  -  Add X1D_DISCARDED macro (IB)
   27 Oct 98  -  Add debug mechanism (IB)
   11 Nov 98  -  1-indexed output profile index (IB)
   13 Nov 98  -  Always alloc memory for interpolated profile (IB)
   18 Nov 98  -  1-indexed output a2center (IB)
   18 Nov 98  -  Background value and its error from command line (IB)
   11 Dec 98  -  Weigth or variance image in optimal extraction (IB)
   16 Dec 98  -  optimal extraction uses background error, not variance (IB)
   02 Dec 99  -  PIX_FILL is stored in profile for excluded pixels (IB)
   02 Dec 99  -  xxx_FILL constants moved to calstis6.h file  (IB)
   07 Dec 99  -  Rejection flag in profile builder (IB)
   13 Dec 99  -  Deviant pixel cleaning in optimal extraction (IB)
   16 Dec 99  -  Weights image updated with rejection flags (IB)
   17 Dec 99  -  Flux normalization factor in opt. extr. rejection (IB)
   17 Jan 00  -  Scattered light correction algorithm (IB)
   11 Apr 00  -  Interpolate profile so it is one pixel smaller (IB)
   25 Apr 00  -  Fixed box computation for binned data (IB)
   16 Jun 00  -  2-D array with profile Y positions (IB)
   13 Jul 00  -  Redefined extrsize as double (IB)
   31 Oct 00  -  Profile centered on centroid instead of trace (IB)
   01 Nov 00  -  Profile offset (IB)
   11 Dec 00  -  Subsampled profile (IB)
   30 Jan 01  -  Centroid output (IB)
   01 Mar 01  -  Turn on background subtraction in profile generator (IB)
   27 Feb 04  -  If background smoothing will be done, set the error
                 contribution of the background to zero (PEH)
   21 Oct 11  -  Copy ERROR to new column NET_ERROR (PEH)

*/

int X1DSpec (StisInfo6 *sts, SpTrace *trc, XtractInfo *xtr,
             double extrsize, SingleGroup *in, SingleGroup *outw,
             FloatHdrData *ssgx, FloatHdrData *ssgy,
             IntensArray *inta, RowContents *row_cont) {


/* arguments:
StisInfo6 *sts         i:  calibration switches and info
SpTrace *trc;          i:  trace parameters
XtractInfo *xtr;       i:  extraction parameters
double extrsize;       i:  extraction box size
SingleGroup *in;       i:  input image
SingleGroup *outw      o:  output weight image in optimal extr. mode
FloatHdrData ssgx;     i:  small-scale distortion in X (not used)
FloatHdrData ssgy;     i:  small-scale distortion in Y (not used)
IntensArray *inta;     i:  array with intensities for optimal extraction
RowContents *row_cont  o:  output row arrays
*/
	int	ipix;		/* index of image pixel in the A1 direction */
	int	rpix;		/* index of reference pixel in the A1 dir. */
	short	oDQ;		/* output DQ flag */
	double	sum;		/* accumulators */
	double	wei;
	double	err;
	double	sumback;
	double	sumnet;
	double	sump;
	double	sum_discarded;
	double 	norm;
	int	discarded;	/* are there discarded pixels ? */
	double	y_box;		/* center of extraction box */
	double	ilow_end; 	/* nominal endpoints of the extraction box */
 	double	ihigh_end;	/* in image pixel units*/
	int     j1, j2;         /* nearest integer endpoints (internal) */
	double  bktilt;         /* angle of background extraction boxes */
	double  s1, s2;		/* fractional pixel areas outside interval */
	double	intens;		/* intensity value for optimal extraction */
	double	weight;		/* opt. weight (goes into weight image) */
	double	*iprofile;	/* interpolated profile for opt. extraction */
	double  *profile_yabs;	/* positions corresponding to profile pixels */
        double  pix_back;       /* extracted value corrected for background */
	int	j, k, debug;
	short	pmask = 0;
	int	status;
	double	j0, centroid, centnorm, offset;
	double x1, x2;
	int ix1, ix2, ix;

	int CalcBack (StisInfo6 *, XtractInfo *, SingleGroup *,
                       FloatHdrData *, FloatHdrData *, int, double, int);

	/* Output extraction info (in image, not reference, pixels) */

	if (sts->extrloc) {
	    row_cont->a2center = sts->cc_a2center *
                                 sts->ltm[1] + sts->ltv[1] + 1.0;
	    row_cont->extrsize = extrsize * sts->ltm[1];
	    /* See if there are command-line overrides. */
	    if (sts->bksize[0] == NO_SIZE)
	        row_cont->bk1size = xtr->bksize[0] * sts->ltm[1];
            else
	        row_cont->bk1size = sts->bksize[0] * sts->ltm[1];
	    if (sts->bksize[1] == NO_SIZE)
	        row_cont->bk2size = xtr->bksize[1] * sts->ltm[1];
            else
	        row_cont->bk2size = sts->bksize[1] * sts->ltm[1];
	    if (sts->bkoffset[0] == NO_SIZE)
	        row_cont->bk1offset = xtr->bkoffset[0] * sts->ltm[1];
	    else
	        row_cont->bk1offset = sts->bkoffset[0] * sts->ltm[1];
	    if (sts->bkoffset[1] == NO_SIZE)
	        row_cont->bk2offset = xtr->bkoffset[1] * sts->ltm[1];
	    else
	        row_cont->bk2offset = sts->bkoffset[1] * sts->ltm[1];
	    if (sts->maxsearch == NO_RANGE)
	        row_cont->maxsearch = xtr->maxsearch * sts->ltm[1];
	    else
	        row_cont->maxsearch = sts->maxsearch * sts->ltm[1];
	}

	/* Alloc memory for interpolated profile. This must be always
           done, to avoid a rui exception.
        */
	if ((iprofile = (double *) malloc ((sts->profile_y + 1) *
                                  sizeof(double))) == NULL)
	    return (OUT_OF_MEMORY);

	/* If building a profile, DQ array associated with the profile
           raw array must be reset.
        */
	if (sts->do_profile) {
	for (ipix = 0; ipix < in->sci.data.nx ; ipix++)
	    for (j = 0; j < sts->profile_y; j++)
                sts->profile_dq[ipix][j] = STIS_OK;
	}

	/* Allocate memory for profile positions array. This array will
           store individual pixel (column) offsets that will be later used
           to compute the centroid. The centroid in turn is used to correct
           the final offset value that is passed to the drizzle routine.
        */
	if (sts->do_profile) {
	    if ((profile_yabs = (double *) malloc ((sts->profile_y + 1) *
                                sizeof(double))) == NULL)
	        return (OUT_OF_MEMORY);
	}

	/* Loop over physical image pixels in the A1 direction. */
	for (ipix = 0; ipix < in->sci.data.nx ; ipix++) {
            float extrlocy;            /* output extraction location */

	    /* Debug control. */
	    if (ipix == 300)
	        debug = ipix;
	    else
	        debug = 0;

	    /* Translate physical pixel index into reference pixel index. */
	    rpix = (int)((ipix - sts->ltv[0]) / sts->ltm[0]);

	    /* Compute background box tilt for current A1 reference position.*/
	    bktilt = BoxTilt (xtr, rpix);

	    /* If there is a command-line override, use it instead. */
	    if (sts->bktilt != NO_TILT)
	        bktilt = sts->bktilt;

	    /* Compute once and for all the tilted background box
               projections.
            */
	    sts->sin_bktilt = sin (bktilt / 180.0 * PI);
	    sts->cos_bktilt = cos (bktilt / 180.0 * PI);

	    /* Add trace curve to cross-correlation result
               to get actual position of trace. This is the
               extraction box center.
            */
	    if (rpix >= 0 && rpix < trc->nelem)
	        y_box = sts->cc_a2center + trc->a2displ[rpix];
	     else
	        y_box = sts->cc_a2center;

	    /* Output extraction position. */
            extrlocy = (float) (y_box * sts->ltm[1] + sts->ltv[1] + 1.0);
	    if (sts->extrloc)
	        row_cont->extrlocy[ipix] = extrlocy;

	    /* Compute nominal endpoints of spectrum extraction box.
               This must be performed in reference pixel coordinates.
            */
	    ilow_end   = (y_box * sts->ltm[1] + sts->ltv[1]) -
                         extrsize * sts->ltm[1] / 2.0;
	    ihigh_end = (y_box * sts->ltm[1] + sts->ltv[1]) +
                         extrsize * sts->ltm[1] / 2.0;

	    /* Compute background coefficients. */
	    if (sts->backcorr == PERFORM || sts->optimal) {
	        if ((status = CalcBack (sts, xtr, in, ssgx, ssgy, ipix,
                                        y_box, debug)))
	            return (status);
	    }

	    /* Get intensity value used as scale factor in optimal
               extraction algorithm.
	    */
	    if (inta->allocated)
	        intens = inta->intens[ipix];
	    else
	        intens = 1.0;

	    /* Compute nearest integer endpoints defining an
               internal interval, and fractional pixel areas
               that remain outside this interval.
            */
	    j1 = ((ilow_end  - (int)ilow_end)  < 0.5) ? (int)ilow_end :
                                                        (int)ilow_end + 1;
	    j2 = ((ihigh_end - (int)ihigh_end) < 0.5) ? (int)ihigh_end :
                                                        (int)ihigh_end + 1;
	    s1 = 0.5 - (ilow_end - j1);
	    s2 = 0.5 + ihigh_end - j2;

	    /* Offset to be added to integer pixel locations to get back
               the actual Y pixel location. This is used to build the array
               of Y locations associated with the profile data (01/18/01, IB)
            */
	    if (sts->do_profile)
	        j0 = ilow_end  - (int)ilow_end - 1.0;
	    else {
	        j0 = ilow_end  - (int)ilow_end;
                if (j0 > 0.5)
                    j0--;
	    }

	    if (sts->optimal) {
	        /* Integrate subsampled profile into actual profile values
                   to be used at current image column.
                */
	        sum = 0.0;
	        for (j = 0; j < sts->profile_y; j++) {

	            offset = j0;

	            /* Compute where in subpixel array the edges of
                       the current pixel are.
                    */
	            x1 = (j - (int)(sts->profile_y / 2) - 0.5 - offset) *
                          sts->subscale + (int)(sts->subprof_size / 2);

	            x2 = (j - (int)(sts->profile_y / 2) + 0.5 - offset) *
                          sts->subscale + (int)(sts->subprof_size / 2);

	            ix1 = (int)x1;
	            ix2 = (int)x2;

	            /* Add contributions from edge areas. */

	            iprofile[j] += sts->subprof[ipix][ix1] *
                                             (1.0 - x1 + ix1);
	            iprofile[j] += sts->subprof[ipix][ix2] *
                                             (x2 - ix2);

	            /* Add contributions from subpixels in between. */

	            if ((ix2 - ix1) > 1) {
	                for (ix = ix1 + 1; ix < ix2; ix++) {
	                    iprofile[j] += sts->subprof[ipix][ix];
	                 }
	            }
	            sum += iprofile[j];
	        }
                if (sum > 0.0)
                    for (j = 0; j <= sts->profile_y; j++)
                        iprofile[j] /= sum;

	        /* Clean up pixels based on the expected profile. */

	        CleanPixels (sts, in, ipix, j1, j2, iprofile, (int)ilow_end,
                             intens, outw, debug);
	    }

	    /* Clear accumulators. */
	    sum           = 0.0;
	    sumback       = 0.0;
	    sumnet        = 0.0;
	    sump          = 0.0;
	    sum_discarded = 0.0;
	    discarded     = 0;
	    err           = 0.0;
	    wei           = 0.0;
	    oDQ           = 0;

	    /* Accumulate the fractional pixel areas at both
               endpoints of extraction box; update profile array.
               If optimal extraction, store weights in weight image.
            */
	    if (j1 > 0 && j1 < in->sci.data.ny) {
	        /* temporary pmask is necessary to avoid segmentation
                   violation when in non-optimal mode.
                */
	        if (sts->optimal)
	            pmask = DQPix (outw->dq.data, ipix, j1);
	        AccumPix (sts, in, ipix, j1, s1, &sum, &sumback,
                          &sumnet, &sump, &err, &wei, &sum_discarded, &oDQ,
                          &discarded, BOX_LOWER, j1, iprofile,
                          sts->profile[ipix], sts->profile_dq[ipix],
                          intens, &weight, pmask, &pix_back, debug);
	        if (sts->do_profile)
	            UpdateProfile ((float)pix_back,
                                   sts->profile[ipix], profile_yabs,
                                   sts->profile_dq[ipix], 0, j1);
	        if (sts->optimal)
	            Pix (outw->sci.data, ipix, j1) = weight;
	    }

	    if (j2 > j1 && j2 > 0 && j2 < in->sci.data.ny) {
	        if (sts->optimal)
	            pmask = DQPix (outw->dq.data, ipix, j2);
	        AccumPix (sts, in, ipix, j2, s2, &sum, &sumback,
                          &sumnet, &sump, &err, &wei, &sum_discarded, &oDQ,
                          &discarded, BOX_UPPER, j1, iprofile,
                          sts->profile[ipix], sts->profile_dq[ipix],
                          intens, &weight, pmask, &pix_back, debug);

	        if (sts->do_profile)
	            UpdateProfile ((float)pix_back,
                                   sts->profile[ipix], profile_yabs,
                                   sts->profile_dq[ipix],
                                   (sts->profile_y)-1, j2);
	        if (sts->optimal)
	            Pix (outw->sci.data, ipix, j2) = weight;
	    }

	    /* Reset pixel pointers at bottom and top of extraction
               box to point to the first and last integral pixels in
               box. Perform the summation and update profile. If optimal
               extraction, store weights in weight image.
            */
	    if (j2 > (j1 + 1)) {
	        j1++;
	        j2--;
	        k = 1;
	        for (j = j1; j <= j2; j++, k++) {
	            if (j > 0 && j < in->sci.data.ny) {
	                if (sts->optimal)
	                    pmask = DQPix (outw->dq.data, ipix, j);
	                AccumPix (sts, in, ipix, j, 1.0, &sum,
                                  &sumback, &sumnet, &sump, &err, &wei,
                                  &sum_discarded, &oDQ, &discarded,
                                  BOX_MID, j1-1, iprofile,
                                  sts->profile[ipix], sts->profile_dq[ipix],
                                  intens, &weight, pmask, &pix_back, debug);

	                if (sts->do_profile)
	                    UpdateProfile ((float)pix_back,
                                           sts->profile[ipix], profile_yabs,
                                           sts->profile_dq[ipix], k, j);
	                if (sts->optimal)
	                    Pix (outw->sci.data, ipix, j) = weight;
	            }
	        }
	    }

	    /* Subtract centroid from profile, compute its distance from
               trace, and correct for fractional pixel distance.
            */
	    if (sts->do_profile) {
	        j1--;
	        j2++;
	        if (j2 > (j1 + 1)) {
	            centroid = 0.0;
	            centnorm = 0.0;
	            for (j = j1, k = 0; j <= j2; j++, k++) {
	                if (j > 0 && j < in->sci.data.ny) {
                            profile_yabs[k] -= y_box;
	                    if (sts->profile[ipix][k] != PIX_FILL) {
                                centroid += sts->profile[ipix][k] *
                                            profile_yabs[k];
                                centnorm += sts->profile[ipix][k];
	                    }
	                }
	            }
	            if (centnorm > 0.0) {
	                centroid /= centnorm;
	                sts->profile_offset[ipix] = centroid + j0 + 0.5;
	                sts->profile_centroid[ipix] = centroid;
	            }
	        }
	    }

	    /* Centroiding doesn't work with low s/n. */
/*
	    if (sts->do_profile)
	        sts->profile_offset[ipix] = j0 + 0.5;
*/
	    /* If more than 10% of the total flux was affected by input
               flagged pixels, set output's 12th bit flag.
            */
	    /*  if (discarded) {                                       */
	    /*     if (sum > 0.0) {                                    */
	    /*          if (sum_discarded / sum > 0.1)                 */
	    /*              oDQ |= X1D_DISCARDED;                      */
	    /*      } else                                             */
	    /*          oDQ |= X1D_DISCARDED;                          */
	    /*  }                                                      */

	    /* OPRs were revised but the code that fixes bad data is still
               in place. The following code is to enforce the new OPR, and
               replaces the above commented-out lines. The output spectrum
               value will still be fixed in the way described in function
               AccPixUnw, which should work quite well in most cases.
            */
	    if (discarded)
	        oDQ |= X1D_DISCARDED;

	    /* Add any background-related bit flag(s). */
	    oDQ |= sts->dqbck;

	    /* Store gross, net, background and flux normalization values. */
	    row_cont->dq[ipix] = oDQ;
	    if (wei > 0.0 && sts->exptime > 0.0) {
                float gross, back, net, error;

                gross = (float) (sum / wei);
	        if (sts->optimal) {
                    back = (float) sumback;
                    sts->profile_rejf[ipix] = err > 0.0 ? sump/err: 0.0;
	        } else {
	            back = (float) (sumback / wei);
                }
                net   = (float) (sumnet / wei);
                error = err > 0.0 ? sqrt(err / wei) : 0.0;

	        row_cont->gross[ipix] = gross / sts->exptime;
                row_cont->back[ipix]  = back / sts->exptime;
                row_cont->net[ipix]   = net / sts->exptime;
		/* at this point, both error and net_error are in counts/s */
                row_cont->error[ipix] = error / sts->exptime;
                row_cont->net_error[ipix] = error / sts->exptime;
	    } else {
	        row_cont->gross[ipix] = 0.0F;
	        row_cont->back[ipix]  = 0.0F;
	        row_cont->net[ipix]   = 0.0F;
	        row_cont->error[ipix] = 0.0F;
	        row_cont->net_error[ipix] = 0.0F;
	    }

	    /* Resample and normalize profile. */

/*  NOTE THAT THIS DISABLES THE S/N CHECKING/CUTTOF IN CPROFILE. */


	    if (sts->do_profile) {
	        if (j2 > (j1 + 1)) {
	            norm     = 0.0;
	            for (j = j1, k = 0; j <= j2; j++, k++) {
	                if (j > 0 && j < in->sci.data.ny) {
	                    norm += sts->profile[ipix][k];
	                }
	            }
	            if (norm != 0.0) {
	                for (j = j1, k = 0; j <= j2; j++, k++) {
	                    if (j > 0 && j < in->sci.data.ny)
                                sts->profile[ipix][k] /= norm;
	                }
	            }
	        }
	    }
	}

	if (sts->do_profile)
	    free (profile_yabs);
	free (iprofile);

	return (0);
}



/*  Routine that extracts a single pixel contribution to the accumulators.
    It is just a front-end to algorithm-specific routines that implement
    unweighted and optimal extraction. Background computation, identical
    in both algorithms, is computed here too.
*/

static void AccumPix (StisInfo6 *sts, SingleGroup *in, int i, int j,
                      double area, double *sum, double *sumback,
                      double *sumnet, double *sump, double *err, double *wei,
                      double *sum_discarded, short *oDQ,
                      int *discarded, int box_pos, int ilow_end,
                      double *iprofile, double *profile, short *profile_dq,
                      double intens, double *weight, short pmask,
                      double *pix_back, int debug) {

/* arguments:
StisInfo6 *sts         i:  calibration switches and info
SingleGroup *in	       i:  input image
int i, j;              i:  pixel indices in physical image units
double area;           i:  fraction of pixel area being extracted
double *sum...         io: accumulators
int *discarded;        o:  are there discarded pixels ?
int box_pos;           i:  pixel is in box extremities or internal
int ilow_end;          i:  image pixel where the box's low end sits
double *iprofile;      i:  interpolated profile for optimal extraction
double *profile;       o:  profile builder basic output
short *profile_dq;     o:  profile builder DQ output
double intens;         i:  intensity value used in optimal extraction
double weight;         o:  weight (goes into output weight image)
short pmask;           i:  profile DQ mask
double *pix_back;      o:  the extracted pixel corrected for background
int debug;             i:  debug flag / pixel index
*/
	double	back, backvar, backerr;

	/* Compute background and update its accumulator. */
	ComputeBack (sts, j, &back, &backvar, &backerr, debug);
	*sumback += back * area;

	/* If background smoothing will be done, set the error contribution
	   of the background to zero.
	*/
	if (sts->bks_mode != BKS_OFF) {
	    backvar = 0.;
	    backerr = 0.;
	}

	/* Perform algorithm-specific extraction. Notice that optimal
           extraction uses the background error, not the variance.
        */
	if (sts->optimal)
	    AccPixOpt (sts, in, i, j, sum, sumback, sumnet, sump,
                          err, wei, sum_discarded, oDQ, discarded,
                          box_pos, ilow_end, iprofile, back, backerr,
                          intens, weight, pmask, debug);
	else
	    AccPixUnw (sts, in, i, j, area, sum, sumback, sumnet,
                       err, wei, sum_discarded, oDQ, discarded,
                       box_pos, ilow_end, profile, profile_dq, back,
                       backvar, pix_back, debug);
}




/*  The following two routines are the core of calstis6. They perform the
    actual extraction, either in unweighted mode with background
    subtraction, or in optimal mode.
*/




/*  Optimal extraction.

    The algorithm currently implemented does NOT attempt to interpolate
    over bad pixels, as the unweighted algorithm does. It also does NOT
    keep track of the number of discarded pixels.
*/

static void AccPixOpt (StisInfo6 *sts, SingleGroup *in, int i, int j,
                       double *sum, double *sumback,
                       double *sumnet, double *sump, double *err, double *wei,
                       double *sum_discarded, short *oDQ, int *discarded,
                       int box_pos, int ilow_end, double *profile,
                       double back, double backerr, double intens,
                       double *weight, short pmask, int debug) {

/* arguments:
StisInfo6 *sts         i:  calibration switches and info
SingleGroup *in	       i:  input image
int i, j;              i:  pixel indices in physical image units
double area;           i:  fraction of pixel area being extracted
double *sum...         io: accumulators
int *discarded;        o:  are there discarded pixels ?
int box_pos;           i:  pixel is in box extremities or internal
int ilow_end;          i:  image pixel where the box's low end sits
double *profile;       i:  normalized profile
double back;           i:  background value
double backerr;        i:  background error (variance / sqrt(n))
double intens;         i:  normalization intensity spectrum value
double weight;         o:  weight (goes into output weight image)
short pmask;           i:  profile DQ mask
int debug;             i:  debug flag / pixel index
*/
	short  dq;
	double pix, p, var, mask;
	double back2, backerr2;

	/* Replace background by command-line overriding values.
           Note that these overriding values are applicable only in
           the case of optimal extraction, hence a separate function.
        */
	ReplaceBack (sts, back, backerr, &back2, &backerr2);

	/* Get pixel value and its "serious" DQ mask. */
	pix = Pix (in->sci.data, i, j);
	dq  = DQPix (in->dq.data, i, j) & sts->sdqflags;

	if((j - ilow_end) < sts->profile_y) {

	    p = profile[j - ilow_end];

	    /* Only accept good pixels that weren't rejected by s-clip. */
            if ((dq == 0) && (pmask == 0)) {

	        /* Compute mask, variance and weight. */
	        mask   = dq == 0 ? 1.0 : 0.0;
	        var    = sts->rn2 + backerr2;
	        var   += fabs (intens * p + back2 * sts->gain);
	        if (var > 0.0) {
	            if (sts->variance_image)
	                *weight = var;
	            else
	                *weight = mask * p / var;
	        } else
	            *weight = 0.0;

	        /* Update accumulators. */
	        *sumnet  += *weight * (pix - back2);
	        *err     += mask * p;
	        *wei     += *weight * p;
                *oDQ     |= dq;
	        *sump    += p;
	    } else {

	        (*discarded)++;
	        *weight = 0.0;

	        /* *err contains the accumulated profile values actually
                   used. Here we accumulate all the existing good profile
                   values that were left behind by the calculation above.
                   These correspond to pixels with good data quality but
                   which were rejected by the profile sigma-clip cleaning.
                   The ratio is used later to normalize the flux.
                */
/*

Turned off by Brian's request 3/23/01

	        if (pmask == X1D_DISCARDED)
	            *sump += p;
*/
	    }



/* This dumps profile values instead of weight values into output image. */
/*
*weight = p;
*/


	}
}



/*  Unweighted extraction.

    This routine updates the accumulators with the proper quantities
    depending on the value of the BACKCORR switch. Note that the error
    is computed based on the ERR array value.

    This routine also supports the profile generator, by depositing each
    pixel into the appropriate element of the 1-D profile array.
*/

static void AccPixUnw (StisInfo6 *sts, SingleGroup *in, int i, int j,
                       double area, double *sum, double *sumback,
                       double *sumnet, double *err, double *wei,
                       double *sum_discarded, short *oDQ, int *discarded,
                       int box_pos, int ilow_end, double *profile,
                       short *profile_dq, double back, double backvar,
                       double *pix_back, int debug) {

/* arguments:
StisInfo6 *sts         i:  calibration switches and info
SingleGroup *in	       i:  input image
int i, j;              i:  pixel indices in physical image units
double area;           i:  fraction of pixel area being extracted
double *sum...         io: accumulators
int *discarded;        o:  are there discarded pixels ?
int box_pos;           i:  pixel is in box extremities or internal
int ilow_end;          i:  image pixel where the box's low end sits
double *profile;       io: 1-D profile array for the i-th image column
short *profile_dq;     io: 1-D profile DQ array for the i-th image column
double back;           i:  background value
double backvar;        i:  background variance
double *pix_back;      o:  the extracted pixel corrected for background
int debug;             i:  debug flag / pixel index
*/

	int     nint;
	double  pix, perr;

	/* If current pixel is not flagged by any serious DQ flags,
           use it directly. Otherwise attempt to interpolate from
           neighboring pixels. If either one of these is also flagged,
           use just the good one. If both are flagged, give up.
           Notice that when handling pixels at either top or bottom
           extremes of the extraction box, just one pixel is available
           for interpolation or replacement. If it is bad, give up.
           The revised OPRs demand that the pixel fill value must be
           zero (so a discarded pixel won't contribute to the result).
        */
        if (!(DQPix (in->dq.data, i, j) & sts->sdqflags)) {
	    pix  = Pix (in->sci.data, i, j);
	    perr = Pix (in->err.data, i, j);
            *oDQ |= DQPix (in->dq.data, i, j);
	} else {
	    switch (box_pos) {
	    case BOX_LOWER:
                if (!(DQPix (in->dq.data, i, j+1) & sts->sdqflags)) {
	            pix  = Pix (in->sci.data, i, j+1);
	            perr = Pix (in->err.data, i, j+1);
                    *oDQ |= DQPix (in->dq.data, i, j+1);
	        } else {
	            pix  = PIX_FILL;
	            perr = ERR_FILL;
	        }
	        *sum_discarded += pix;
	        (*discarded)++;
	        break;
	    case BOX_UPPER:
                if (!(DQPix (in->dq.data, i, j-1) & sts->sdqflags)) {
	            pix  = Pix (in->sci.data, i, j-1);
	            perr = Pix (in->err.data, i, j-1);
                    *oDQ |= DQPix (in->dq.data, i, j-1);
	        } else {
	            pix  = PIX_FILL;
	            perr = ERR_FILL;
	        }
	        *sum_discarded += pix;
	        (*discarded)++;
	        break;
	    case BOX_MID:
	        nint = 0;
	        pix  = 0.0;
	        perr = 0.0;
                if (!(DQPix (in->dq.data, i, j+1) & sts->sdqflags)) {
	            pix  = Pix (in->sci.data, i, j+1);
	            perr = Pix (in->err.data, i, j+1) *
                           Pix (in->err.data, i, j+1);
                    *oDQ |= DQPix (in->dq.data, i, j+1);
	            nint++;
	        }
                if (!(DQPix (in->dq.data, i, j-1) & sts->sdqflags)) {
	            pix  += Pix (in->sci.data, i, j-1);
	            perr += Pix (in->err.data, i, j-1) *
                            Pix (in->err.data, i, j-1);
                    *oDQ |= DQPix (in->dq.data, i, j-1);
	            nint++;
	        }
	        if (nint > 0) {
	            pix /= nint;
	            perr = sqrt (perr / nint);
	        } else {
	            pix  = PIX_FILL;
	            perr = ERR_FILL;
	        }
	        *sum_discarded += pix;
	        (*discarded)++;
	        break;
	    }
	}

	/* Update remaining accumulators. */
	*sum     += pix * area;
	*sumnet  += (pix - back) * area;
	*err     += ((perr * perr) + backvar) * area;

	/* This does not accummulate, it just acts as a no-op
           in the calling routine.
        */
	*wei  = 1.0;

        /* Return value. */
        if (pix != PIX_FILL)
            *pix_back = (pix - back);
	else
            *pix_back = PIX_FILL;
}


/*  Flags deviant pixels based on the expected profile. This is used
    by the optimal extraction code.

    The main output of this function is in the DQ array of the input
    SingleGroup image. The DQ array of the weights image is also
    updated.

    The code in this routine partialy duplicates code found elsewhere
    in this file. This is a consequence of the later inclusion of this
    algorithm in the optimal extraction suite.
*/

static void CleanPixels (StisInfo6 *sts, SingleGroup *in, int ipix,
                         int j1, int j2, double *iprofile, int ilow_end,
                         double intens, SingleGroup *outw, int debug) {

/* arguments:
StisInfo6 *sts         i:  calibration switches and info
SingleGroup *in	       i:  input image
int ipix;              i:  image column
int j1, j2;            i:  extreme indices in physical image units
double *iprofile;      i:  interpolated profile for optimal extraction
int ilow_end;          i:  image pixel where the profile's low end sits
double intens;         i:  intensity value used in optimal extraction
SingleGroup *outw      o:  output weight image in optimal extr. mode
int debug;             i:  debug control
*/
	int j, j3, j4;

	/* First clean the two extreme pixels. */
	if (j1 > 0 && j1 < in->sci.data.ny)
	    CleanPixel (sts, in, ipix, j1, iprofile, ilow_end, intens,
                        outw, debug);
	if (j2 > 0 && j2 < in->sci.data.ny)
	    CleanPixel (sts, in, ipix, j2, iprofile, ilow_end, intens,
                        outw, debug);

	/* Then the remaining ones. Reset pixel pointers at bottom and
           top of extraction box to point to the first and last integral
           pixels in box.
        */
	j3 = j1;
	j4 = j2;
	if (j4 > (j3 + 1)) {
	    j3++;
	    j4--;
	    for (j = j3; j <= j4; j++)
	        CleanPixel (sts, in, ipix, j, iprofile, ilow_end, intens,
                            outw, debug);
	}
}


/*  Flags a single deviant pixel based on the expected profile. This
    is used by the optimal extraction code.

    The main output of this function is in the DQ array of the input
    SingleGroup image. The DQ array of the weights image is also
    updated.

    The code in this routine partialy duplicates code found elsewhere
    in this file. This is a consequence of the later inclusion of this
    algorithm in the optimal extraction suite.
*/

static void CleanPixel (StisInfo6 *sts, SingleGroup *in, int ipix, int j,
                        double *iprofile, int ilow_end, double intens,
                        SingleGroup *outw, int debug) {

/* arguments:
StisInfo6 *sts         i:  calibration switches and info
SingleGroup *in	       i:  input image
int ipix;              i:  image column in physical image units
int j;                 i:  image line in physical image units
double *iprofile;      i:  interpolated profile for optimal extraction
int ilow_end;          i:  image pixel where the profile's low end sits
double intens;         i:  intensity value used in optimal extraction
SingleGroup *outw      o:  output weight image in optimal extr. mode
int debug;             i:  debug control
*/
	short  dq, fullmask;
	double pix, p, var, crit;
	double	back, backvar, backerr, back2, backerr2;

	/* Compute background and eventually replace by command
           line overriding values.
        */
	ComputeBack (sts, j, &back, &backvar, &backerr, debug);
	ReplaceBack (sts, back, backerr, &back2, &backerr2);

	/* Get pixel value and its "serious" DQ mask. */
	pix = Pix (in->sci.data, ipix, j);
	dq  = DQPix (in->dq.data, ipix, j) & sts->sdqflags;

        /* Compute relevant quantities and flag pixel. Both input and
           output SingleGroups receive the flag. The output weights image
           will be used for diagnostic purposes. The input SingleGroup
           must get the flags as well, since the data will be extracted
           from it. Note that the rejection flag in the input SingleGroup
           is ORed with whatever is already set in the DQ array.
        */
	if (((j - ilow_end) < sts->profile_y) &&  !dq) {
	    p    = iprofile[j - ilow_end];
	    var  = sts->rn2 + backerr2;
	    var += fabs (intens * p + back2 * sts->gain);
	    if (var > 0.0) {
	        crit = pix - back2 - p * intens;
	        crit = crit * crit / var;
	        if (crit > (sts->sclip * sts->sclip)) {
	            fullmask = DQPix (in->dq.data, ipix, j) | X1D_DISCARDED;
	            DQSetPix (in->dq.data,   ipix, j, fullmask);
	            DQSetPix (outw->dq.data, ipix, j, X1D_DISCARDED);
	        }
	    }
	}
}



/*  Updates the profile array. */

static void UpdateProfile (float pix, double *profile, double *profile_y,
                           short *profile_dq, int index, int index_image) {

	if (pix != PIX_FILL) {
	    profile   [index] = pix;
	    profile_dq[index] = STIS_OK;
	} else {
	    profile   [index] = PIX_FILL;
	    profile_dq[index] = NO_GOOD_DATA;
	}
        profile_y[index] = index_image;
}



/*
    Compute background value at a given pixel position in the
    extraction box.
*/

static void ComputeBack (StisInfo6 *sts, int j, double *back,
                         double *backvar, double *backerr, int debug) {

	double ix[1], iy[1];

	void ComputePoly (double *, int, double *, int, double *);

	/* The background coefficients are defined for physical (image)
           pixel addresses.

           The original code was tweaked in order to support the
           scattered light correction algorithm.
        */
	if (!sts->scatter) {
	    if (sts->backcorr == PERFORM || sts->optimal) {
	        *back    = sts->bck[0]  + sts->bck[1]  * j;
	        *backvar = sts->vbck[0] + sts->vbck[1] * j * j;
	        *backerr = sts->ebck;
	    } else {
	        *back    = 0.0;
	        *backvar = 0.0;
	        *backerr = 0.0;
	    }
	} else {
	    ix[0] = (float)j;
	    ComputePoly (ix, 1, sts->bck, BACKP, iy);
	    *back    = iy[0];
	    *backvar = 0.0;
	    *backerr = 0.0;
	}
}



/*
    Replace background values by command line overrides.
*/

static void ReplaceBack (StisInfo6 *sts, double back, double backerr,
                         double *back2, double *backerr2) {

	/* If background value was supplied in command line, use it.
           If only value, but no error, was supplied, variance is
           set to zero.
        */
	if (sts->backval == NO_VALUE) {
	    *back2 = back;
	    *backerr2 = backerr * backerr;
	} else {
	    *back2 = sts->backval;
	    if (sts->backerr == NO_VALUE)
	        *backerr2 = 0.0;
	    else
	        *backerr2 = sts->backerr * sts->backerr;
	}
}



/*  Computes background tilt angle as a function of reference pixel
    position in the A1 direction. Pixel position is expressed wrt
    the physical array edge, *not* its center, and starts with 1,
    not zero. These choices depend on the definitions used to build
    the XTRACTAB reference table.
*/

static double BoxTilt (XtractInfo *xtr, int pix) {

	double angle;
	int    i;

	angle = 0.0;

	for (i = 0; i < xtr->ncoeffbk; i++)
	    angle += xtr->bktcoeff[i] * pow ((double)(pix+1), (double)i);

	return (angle);
}

