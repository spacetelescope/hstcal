# include <stdio.h>
# include <stdlib.h>
# include <math.h>	/* fabs */
# include <float.h>

# include "xtables.h"
# include "hstio.h"

# include "stis.h"
# include "stisdef.h"
# include "calstis6.h"
# include "hstcalerr.h"
# include "stisdq.h"

# define	FWHM_X		1.0	/* 2-D Gaussian kernel */
# define	FWHM_Y		2.5	/* 2-D Gaussian kernel */
# define	PDEGREE		4

static void ComputeGaussian (double, int, float *);


/*
   Re-define background regions for scatterd light correction algorithm,
   and smooth them.

   This function is an adaptation of the code in CrossCorr (whatch for code
   duplication !). It builds a similarly averaged 1-D spatial profile and
   looks for the regions at the profile wings that can be used to determine
   the background. The result is used to redefine values for the background
   extraction box sizes and offsets stored in the StisInfo6 structure.

   This function also smooths the relevant background regions in the
   input image, to be later used by the 1-D extraction code. The limits
   of the region (image lines) are output to the caller.



   Revision history:
   ----------------
   10 Jan 00  -  Implemented (I.Busko)
   07 Feb 06  -  Add an argument to Interp2D (PEH)
*/

int DefineBackRegions (StisInfo6 *sts, SpTrace *trc, CoordInfo *coord,
                       SingleGroup *in, FloatHdrData *ssgx,
                       FloatHdrData *ssgy, int *ilow_end, int *ihigh_end) {

/* arguments:
StisInfo6 *sts      io: calibration switches and info
SpTrace *trc;       i:  full list of spectrum traces
CoordInfo *coord;   i:  rectified image size information
SingleGroup *in	    i:  input image
FloatHdrData ssgx;  i:  small-scale distortion in X (not used)
FloatHdrData ssgy;  i:  small-scale distortion in Y (not used)
int *ilow_end;	    o:  pointers in image array to the subimage
int *ihigh_end;	        that maps into the convolution buffer
*/

	int status;

        SpTrace *trace_y;       /* interpolated spectrum trace */
	int	ipix;		/* index of image pixel in the A1 direction */
	int	rpix;		/* index of reference pixel in the A1 dir. */
	int	jpix;		/* index of array element */
	double	*sum;		/* accumulators */
	double	*wei;
	double	y_nom;		/* nominal coordinate in the A2 direction */
	double	iy_nom;		/* above quantity in image pixel units */
	double	y_trc;		/* trace-corrected coordinate */
	double	low_end; 	/* nominal endpoints of profile function*/
	double	high_end;
	float	oSci, oErr;	/* interpolated values from image array */
	short	oDQ;
	double	hwei;
	double	*px, *py, *pw;	/* polynomial fittting buffers */
	double	*coeff;
	int	ndata;		/* # data points in polynomial fit */
	float	**buffers;	/* 2-D convolution buffer (sci) */
	float	**buffere;	/* 2-D convolution buffer (err) */
	int	bufx, bufy;	/* 2-D buffer dimensions */
	int	boff;		/* padding offset in buffer */
	float	*kernel;	/* 1-D convolution kernel */
	int	ksize;		/* kernel size */
	float	ymax;		/* used to find maxima */
	int	imax;		/* index of peak value in array */
	int	size;		/* size of region in A2 direction */
	int	i1, i2;		/* used to define the central 40% of profile */
	int	l1, l2;		/* internal end points of profile wings */
	double	sz1, sz2;	/* computed background boxes */
	double	offst1, offst2;
	int	i, j, k, l;	/* generic indices indices */
	int	npix;

	int InterpTrace6 (SpTrace **, double, SpTrace **);
	void FreeTrace6 (SpTrace **);
	int FitPoly (double *, double *, double *, int, int, double *);
	void ComputePoly (double *, int, double *, int, double *);

	npix = in->sci.data.nx;

	/* FIRST STEP: FIND BACKGROUND REGIONS. */

	/* Get size of 2-D region. */
	size = coord->npix;

	/* Alloc and clear accumulator memory. */
	if ((sum = (double *) calloc (size+1, sizeof(double)))== NULL)
	    return (OUT_OF_MEMORY);
	if ((wei = (double *) calloc (size+1, sizeof(double)))== NULL)
	    return (OUT_OF_MEMORY);

	/* Compute endpoints of spatial profile. Use the crosscor
           corrected value; this assumes that cross correlation
           is executed before this function is called.
         */
	low_end  = sts->cc_a2center - ((double)size)/2.;
	high_end = sts->cc_a2center + ((double)size)/2.;

	/* Loop over the spatial profile range in 1 pixel steps. */
	jpix = 0;
	trace_y = NULL;

	for (y_nom = low_end; y_nom <= high_end; y_nom += 1.0) {

	    /* Interpolate in trace table. */

	    if ((status = InterpTrace6 (&trc, y_nom, &trace_y)))
	        return (status);

	    /* Loop over image pixels in the A1 direction. */
	    for (ipix = 0; ipix < npix ; ipix++) {

	        /* Translate image pixel index into reference pixel index. */
	        rpix = (int)((ipix - sts->ltv[0]) / sts->ltm[0]);

	        /* Add trace curve to get actual position of trace. */
	        if (rpix > 0 && rpix < trace_y->nelem)
	            y_trc = y_nom + trace_y->a2displ[rpix];
	         else
	            y_trc = y_nom;

	        /* Translate back into image pixel index. */
	        iy_nom = y_trc * sts->ltm[1] + sts->ltv[1];

	        /* Store for subsequent use in 2-D filtering. */
	        if (y_nom == low_end && ipix == 0)
	            *ilow_end = iy_nom;

		/* Interpolate along column, checking for out of bounds. */
		Interp2D (in, sts->cc_sdqflags, (double)ipix, iy_nom, 1.0,
	                 WGT_VARIANCE, &oSci, &oErr, &oDQ);
	        hwei = 1.0;   /* this weight is kept in place in case we
                                 decide to weight with the error array. This
                                 was implemented in the original crosscor
                                 routine and later turned off. */

	        /* Update sum and weight, discarding flagged data. */

	        if ((oDQ != DETECTORPROB) && !(oDQ & sts->cc_sdqflags)) {
	            sum[jpix] += oSci;
	            wei[jpix] += hwei;
	        }

	    }

	    /* Bump profile index. */
	    jpix++;
	}

	/* Set image pointers used in 2-D convolution. Add one extra
           line at each extremity to give some additional room, but avoid
           addressing out of bounds.
        */
	*ihigh_end = *ilow_end + jpix + 1;
	(*ilow_end)--;
	if (*ilow_end < 0)
	    *ilow_end = 0;
	if (*ihigh_end >= in->sci.data.ny)
	    *ihigh_end = in->sci.data.ny - 1;

	/* Free interpolated trace structure. */
	FreeTrace6 (&trace_y);

	/* Locate internal end points of profile's wings. These are found
           by disregarding the central 40% of the spatial profile.
        */
	i1 = (int)((double)jpix * 0.3);
	i2 = (int)((double)jpix * 0.7);

	/* Fit a 4th degree polynomial to the profile wings. The independent
           variable is the relative pixel address within the profile array,
           not the physical image coordinate.
        */
	if ((px = (double *) calloc (jpix, sizeof(double)))== NULL)
	    return (OUT_OF_MEMORY);
	if ((py = (double *) calloc (jpix, sizeof(double)))== NULL)
	    return (OUT_OF_MEMORY);
	if ((pw = (double *) calloc (jpix, sizeof(double)))== NULL)
	    return (OUT_OF_MEMORY);
	if ((coeff = (double *) calloc (PDEGREE+3, sizeof(double)))== NULL)
	    return (OUT_OF_MEMORY);

	ndata = 0;
	for (i = 0; i < i1; i++) {
	    px[ndata] = (double)i;
	    py[ndata] = sum[i];
	    pw[ndata] = 1.0;
	    ndata++;
	}
	for (i = i2; i < jpix; i++) {
	    px[ndata] = (double)i;
	    py[ndata] = sum[i];
	    pw[ndata] = 1.0;
	    ndata++;
	}

	FitPoly (px, py, pw, ndata, PDEGREE, coeff);

	/* Now we subtract this polynomial from the raw data, including
           the central region that was previously skipped.
        */
	for (i = 0; i < jpix; i++)
	    px[i] = (float)i;

	ComputePoly (px, jpix, coeff, PDEGREE, py);

	for (i = 0; i < jpix; i++)
	    py[i] = sum[i] - py[i];

	/* Look for the peak in the profile, and then find the two
           points where the profile falls to 0.035 of the peak.
        */
	imax = - jpix;
        ymax = - FLT_MAX;
        for (i = 0; i < jpix; i++) {
            if (py[i] > ymax) {
                ymax = py[i];
	        imax = i;
	    }
        }
	if (imax < 0)
	    return (ERROR_RETURN);
	ymax *= 0.035;
        for (l1 = 1; l1 < imax; l1++) {
            if (py[l1] > ymax) {
	        l1--;
	        break;
	    }
	}
        for (l2 = jpix-1; l2 > imax; l2--) {
            if (py[l2] > ymax) {
	        l2++;
	        break;
	    }
	}

	free (coeff);
	free (pw);
	free (py);
	free (px);
	free (wei);
	free (sum);

	/* Finally, compute new background extraction boxes. */
	sz1 = l1;
	sz2 = jpix - l2;
	offst1 = -((double)imax        - sz1 / 2.);
	offst2 =   (double)(jpix-imax) - sz2 / 2.;

	/* They must be physically meaningful. These criteria are arbitrary ! */
	if (sz1 < 4 || sz2 < 4 || offst1 > -4 || offst2 < 4) {
	    printf (
  "Warning  No reliable signal. Using default background extraction boxes.\n");
	    sz1 = sts->bksize[0];
	    sz2 = sts->bksize[1];
	    offst1 = sts->bkoffset[0];
	    offst2 = sts->bkoffset[1];
	}

	sts->bksize[0]   = sz1;
	sts->bksize[1]   = sz2;
	sts->bkoffset[0] = offst1;
	sts->bkoffset[1] = offst2;

	/* SECOND STEP: FILTER WITH 2-D GAUSSIAN. */

	/* We convolve first along the A1 direction, then along A2.
           A temporary 2-D buffer holds the affected image lines with
           enough zero padding to cancel wraparound effects. Both sci
           and err are processed since the Lee filter needs the correct
           noise model for the filtered data.
        */

	/* Buffer dimensions.  */
	bufx = npix;
	bufy = *ihigh_end - *ilow_end + 1 + 2 * (int)FWHM_Y;
	boff = (int)FWHM_Y;

	/* Alloc zeroed 2-D buffer. */
	if ((buffers = (float **) malloc (bufy * sizeof (float *))) == NULL)
	    return (OUT_OF_MEMORY);
	if ((buffere = (float **) malloc (bufy * sizeof (float *))) == NULL)
	    return (OUT_OF_MEMORY);
	for (i = 0; i < bufy; i++) {
	    if ((buffers[i] = (float *) calloc (bufx, sizeof (float))) == NULL)
	        return (OUT_OF_MEMORY);
	    if ((buffere[i] = (float *) calloc (bufx, sizeof (float))) == NULL)
	        return (OUT_OF_MEMORY);
	}

	/* Compute 1-D kernel in A1 direction. */
	ksize = 2 * (int)FWHM_X + 1;
	if ((kernel = (float *) malloc (ksize * sizeof (float))) == NULL)
	    return (OUT_OF_MEMORY);
	ComputeGaussian (FWHM_X, ksize, kernel);

	/* Convolve in the A1 direction, storing result in buffer. No
           padding exists in the original data array, so we must explictly
           check for out-of-bounds in the A1 direction. The buffer is
           filled such that zero padding results in the A2 direction.
        */
	jpix = boff;
	for (j = *ilow_end; j < *ihigh_end; j++) {
	    for (i = 0; i < npix; i++) {
	        for (k = 0; k < ksize; k++) {
	            l = i - ksize/2 + k;
	            /* We do not normalize for rejected pixels yet. */
	            if (l >=0 && l < npix) {
	                if (!(DQPix (in->dq.data, l, j) & sts->sdqflags))
	                    buffers[jpix][i] += Pix(in->sci.data, l, j) *
                                                kernel[k];
	                    buffere[jpix][i] += Pix(in->err.data, l, j) *
                                                Pix(in->err.data, l, j) *
                                                kernel[k];
	            }
	        }
	    }
	    jpix++;
	}
	free (kernel);

	/* Compute 1-D kernel in A2 direction. */
	ksize = 2 * (int)FWHM_Y + 1;
	if ((kernel = (float *) malloc (ksize * sizeof (float))) == NULL)
	    return (OUT_OF_MEMORY);
	ComputeGaussian (FWHM_Y, ksize, kernel);

	/* Convolve in the A2 direction, with the result stored back in
           the original image array. Now the data source (buffer) is
           zero-padded so there is no need to check for out-of-bounds.
        */
	for (i = 0; i < npix; i++) {
	    jpix = boff;
	    for (j = *ilow_end; j < *ihigh_end; j++) {
	        Pix(in->sci.data, i, j) = 0.0;
	        Pix(in->err.data, i, j) = 0.0;
	        for (k = 0; k < ksize; k++) {
	            l = jpix - ksize/2 + k;
	            Pix(in->sci.data, i, j) += buffers[l][i] * kernel[k];
	            Pix(in->err.data, i, j) += buffere[l][i] * kernel[k];
	        }
	        jpix++;
	    }
	}
	free (kernel);
	for (i = 0; i < npix; i++) {
	    for (j = *ilow_end; j < *ihigh_end; j++) {
	        if (Pix(in->err.data, i, j) > 0.0)
	            Pix(in->err.data, i, j) = sqrt (Pix(in->err.data, i, j));
	    }
	}

	for (i = 0; i < bufy; free (buffers[i++]));
	for (i = 0; i < bufy; free (buffere[i++]));
	free (buffers);
	free (buffere);

	return (STIS_OK);
}



/*
    Function to compute normalized 1-D Gaussian.
*/

static void ComputeGaussian (double fwhm, int size, float *array) {

/* arguments:
double	fwhm;			i: FWHM of Gaussian
int	size;			i: size of output array
float	*array;			o: array with Gaussian values
*/
	int	i;
	float	amax;

	for (i = 0; i < size; i++)
	    array[i] = exp (-2.70927 *
                       pow((double)((i - size/2) / fwhm), 2.));
	amax = 0.0;
	for (i = 0; i < size; amax += array[i++]);
	for (i = 0; i < size; array[i++] /= amax);
}
