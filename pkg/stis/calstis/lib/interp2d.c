# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "hstio.h"

# include "stis.h"
# include "stisdq.h"

/* This routine interpolates in a 2-D array to determine the value at one
   pixel.  The values for science data, error and data quality flags are
   returned.

   The location of the output pixel as mapped back into the input image
   (these are given by ix and iy) is used to determine what interpolation
   to use.  If the output pixel maps back very close to the input pixel
   (see NEAR_CENTER, above), that pixel's value is copied to output,
   scaled by the jacobian.  If the output pixel maps back to very close to
   a single column or row, linear interpolation is used on the two pixels
   nearest the mapping of the output pixel center.  Otherwise, bilinear
   interpolation is used on the four pixels nearest the center.

   The output data quality is obtained by ORing the data quality flags
   of the input pixels (one, two, or four pixels, as described above)
   that are closest to the mapped image of the center of the output pixel.

   Phil Hodge, 1997 Oct 15:
	Flag out of bounds with DETECTORPROB instead of DATAMASKED.

   Phil Hodge, 1998 Dec 17:
	Remove the code to set sdqflags to zero locally.  If sdqflags
	is non-zero, therefore, data quality will be checked and used.

   Phil Hodge, 1999 Feb 15:
	Remove the sections of code having to do with checking whether
	[ix,iy] is close to the center of a pixel; different interpolation
	was used in that case, and that resulted in discontinuities in the
	results.
	Change the interpolation of errors from:
		sqrt (sum ( w[i] * err[i]^2 )) / sum (w[i])
	to
		sqrt (sum ( (w[i] * err[i])^2 )) / sum (w[i])

   Phil Hodge, 2006 Feb 20:
	Add the err_algorithm argument, and use it to specify which
	algorithm will be used to interpolate the error array:

	    err_algorithm = WGT_VARIANCE (weight the variance):
		sqrt (sum ( w[i] * err[i]^2 ))

	    err_algorithm = WGT_ERROR (weight the error):
		sqrt (sum ( (w[i] * err[i])^2 ))
*/

void Interp2D (SingleGroup *in, short sdqflags,
	double ix, double iy, double jacobian,
	int err_algorithm,
	float *oSci, float *oErr, short *oDQ) {

/* arguments:
SingleGroup *in   i: input data
short sdqflags    i: which data quality flags are serious
double ix, iy     i: pixel location in input
double jacobian   i: factor for conserving flux
float *oSci       o: interpolated value in science array
float *oErr       o: interpolated error value
short *oDQ        o: interpolated data quality value
*/

	float p, q, r, s;	/* 1-D weights */
	float w0, w1, w2, w3;	/* weight for each pixel */
	float e0, e1, e2, e3;	/* error * weight for each pixel */
	float sumw;		/* sum of weights */
	float value;		/* interpolated value (sci or err) */
	int inx, iny;		/* size of input image */
	int ii, ij;		/* a pixel location in input */
	int iix, iiy; 		/* nearest integers to ix & iy */
	int ngood;		/* number of neighbors (1-4) flagged as good */

	inx = in->sci.data.nx;
	iny = in->sci.data.ny;

	/* Round off to nearest integer, and check for out of bounds. */
	iix = NINT(ix);
	iiy = NINT(iy);
	if (iix < 0 || iix > inx-1 ||
	    iiy < 0 || iiy > iny-1) {
	    /* The current point is outside the input image. */
	    *oSci = 0.;			/* harmless values */
	    *oErr = 0.;
	    *oDQ = DETECTORPROB;
	    return;
	}

	/* Use bilinear interpolation on the four pixels centered
	   near (ix,iy).  These pixels are:
		[ii, ij+1]  [ii+1, ij+1]
		[ii, ij  ]  [ii+1, ij  ]
	*/

	ii = (int)ix;
	if (ii < 0)
	    ii = 0;
	if (ii > inx-2)
	    ii = inx-2;
	q = ix - ii;		/* weights for X direction */
	p = 1.0F - q;

	ij = (int)iy;
	if (ij < 0)
	    ij = 0;
	if (ij > iny-2)
	    ij = iny-2;
	s = iy - ij;		/* weights for Y direction */
	r = 1.0F - s;

	/* Assign a weight for each of the four pixels. */

	ngood = 0;
	if (DQPix (in->dq.data, ii, ij) & sdqflags) {
	    w0 = 0.0F;
	} else {
	    w0 = p * r;		/* lower left pixel */
	    ngood++;
	}
	if (DQPix (in->dq.data, ii+1, ij) & sdqflags) {
	    w1 = 0.0F;
	} else {
	    w1 = q * r;		/* lower right pixel */
	    ngood++;
	}
	if (DQPix (in->dq.data, ii, ij+1) & sdqflags) {
	    w2 = 0.0F;
	} else {
	    w2 = p * s;		/* upper left pixel */
	    ngood++;
	}
	if (DQPix (in->dq.data, ii+1, ij+1) & sdqflags) {
	    w3 = 0.0F;
	} else {
	    w3 = q * s;		/* upper right pixel */
	    ngood++;
	}

	sumw = w0 + w1 + w2 + w3;

	if (ngood == 0 || sumw <= 0.) {

	    /* All four are bad, or the sum of weights is zero,
		so just copy values from the nearest pixel.
	    */
	    *oSci = Pix (in->sci.data, iix, iiy) * jacobian;
	    *oErr = Pix (in->err.data, iix, iiy) * sqrt (jacobian);

	} else {

	    if (ngood < 4) {
		/* Normalize the weights so their sum is one. */
		w0 /= sumw;
		w1 /= sumw;
		w2 /= sumw;
		w3 /= sumw;
	    }

	    /* Interpolate. */

	    value = w0 * Pix (in->sci.data, ii, ij) +
		    w1 * Pix (in->sci.data, ii+1, ij) +
		    w2 * Pix (in->sci.data, ii, ij+1) +
		    w3 * Pix (in->sci.data, ii+1, ij+1);
	    *oSci = jacobian * value;

	    e0 = Pix (in->err.data, ii, ij);
	    e1 = Pix (in->err.data, ii+1, ij);
	    e2 = Pix (in->err.data, ii, ij+1);
	    e3 = Pix (in->err.data, ii+1, ij+1);
	    if (err_algorithm == WGT_VARIANCE) {
		value = (w0 * e0 * e0) + (w1 * e1 * e1) +
		        (w2 * e2 * e2) + (w3 * e3 * e3);
	    } else if (err_algorithm == WGT_ERROR) {
		value = (w0 * w0 * e0 * e0) + (w1 * w1 * e1 * e1) +
		        (w2 * w2 * e2 * e2) + (w3 * w3 * e3 * e3);
	    } else {
		value = -1.;	/* cannot happen */
	    }
	    if (jacobian * value > 0.)
		*oErr = sqrt (jacobian * value);
	    else
		*oErr = 0.;
	}

	*oDQ = DQPix (in->dq.data, ii, ij) |
	       DQPix (in->dq.data, ii+1, ij) |
	       DQPix (in->dq.data, ii, ij+1) |
	       DQPix (in->dq.data, ii+1, ij+1);
}
