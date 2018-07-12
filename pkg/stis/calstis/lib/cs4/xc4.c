/* This file contains the following:
	ReadRealImage
	XCWavecal
internal:
	FindOffset
*/

# include <stdio.h>

# include "hstio.h"
# include "stis.h"
# include "hstcalerr.h"

static void FindOffset (CmplxArray *, int, int, double *, double *);

/* ReadRealImage allocates memory for z->data, copies the SCI extension
   of a SingleGroup into the real part of z->data, and sets the imaginary
   part of z->data to zero.

   Phil Hodge, 2000 Jan 5
   Phil Hodge, 2000 July 21:
	Extract generic functions, moving them to lib/xcfft.c.
   Phil Hodge, 2004 July 23:
	Replace TOO_FAR = 64 with too_far = z->nx / 2, because for wavecals
	taken with 6X0.2 the cross correlation can be broad.
   Phil Hodge, 2011 Jan 5:
	Change XCWavecal so it leaves cwave unchanged, putting the cross
	correlation in a new argument, crosscorr.
*/

int ReadRealImage (SingleGroup *in, int x_sect[], int y_sect[],
		CmplxArray *z) {

/* arguments:
SingleGroup *in   i: input data
int x_sect[]      i: first and last pixels in the X axis
int y_sect[]      i: first and last pixels in the Y axis
CmplxArray *z     o: complex data, copied from input
*/
	int nx, ny;
	int i, j;

	nx = in->sci.data.nx;
	ny = in->sci.data.ny;

	if (AllocCmplxArray (z, nx, ny))
	    return (OUT_OF_MEMORY);

	/* Copy to the complex array.  The data portion of z was initialized
	   to zero when z was allocated.
	*/
	for (j = y_sect[0];  j <= y_sect[1];  j++) {
	    for (i = x_sect[0];  i <= x_sect[1];  i++) {
		RPIX2D (z, i, j) = Pix (in->sci.data, i, j);
	    }
	}

	return (0);
}

/* This multiplies the complex conjugate of clamp by cwave (to take
   the cross correlation between them), then takes the inverse Fourier
   transform of the product, leaving the resulting cross correlation
   in crosscorr.
*/

int XCWavecal (CmplxArray *clamp, CmplxArray *cwave, CmplxArray *crosscorr,
		double *shiftx, double *shifty) {

/* arguments:
   On input, clamp and cwave should contain the Fourier transforms of
   the 2-D image data.
   On output, crosscorr will contain the cross correlation of the input
   cwave and clamp.  Neither clamp nor cwave will be modified.

CmplxArray *clamp        i: FT of template wavecal, created from lamp spectrum
CmplxArray *cwave        i: FT of wavecal data (SCI extension)
CmplxArray *crosscorr    o: cross correlation of clamp and cwave
double *shiftx, *shifty  o: shift from clamp to cwave, in each axis
*/

	int i, j;		/* loop indexes */
	int nx, ny;		/* image size */
	float max;		/* maximum value in cross correlation */
	int lmaxx, lmaxy;	/* location of maximum */
	double offsetx, offsety;	/* offset of peak from lmaxx, lmaxy */
	int status;

	nx = cwave->nx;
	ny = cwave->ny;

	/* Multiply the complex conjugate of the first by the second. */
	for (j = 0;  j < ny;  j++) {
	    for (i = 0;  i < nx;  i++) {
		RPIX2D (crosscorr, i, j) =
			RPIX2D (clamp, i, j) * RPIX2D (cwave, i, j) +
			IPIX2D (clamp, i, j) * IPIX2D (cwave, i, j);
		IPIX2D (crosscorr, i, j) =
			RPIX2D (clamp, i, j) * IPIX2D (cwave, i, j) -
			IPIX2D (clamp, i, j) * RPIX2D (cwave, i, j);
	    }
	}

	/* Take the inverse Fourier transform of the product. */
	if ((status = ifft2d (crosscorr)))
	    return (status);

	/* Find the peak. */
	max = RPIX2D (crosscorr, 0, 0);		/* initial values */
	lmaxx = 0;
	lmaxy = 0;
	for (j = 0;  j < ny;  j++) {
	    for (i = 0;  i < nx;  i++) {
		if (max < RPIX2D (crosscorr, i, j)) {
		    max = RPIX2D (crosscorr, i, j);
		    lmaxx = i;
		    lmaxy = j;
		}
	    }
	}

	FindOffset (crosscorr, lmaxx, lmaxy, &offsetx, &offsety);

	if (lmaxx < nx / 2)
	    *shiftx = lmaxx + offsetx;
	else
	    *shiftx = lmaxx - nx + offsetx;

	if (lmaxy < ny / 2)
	    *shifty = lmaxy + offsety;
	else
	    *shifty = lmaxy - ny + offsety;

	return (0);
}

static void FindOffset (CmplxArray *z, int lmaxx, int lmaxy,
		double *offsetx, double *offsety) {

	float bkg;		/* estimate of background level */
	float cutoff;		/* include if value is greater than this */
	float value;		/* value at a pixel */
	double sumv, sumxv, sumyv;
	int i, j, ii, jj;
	int i1, i2, j1, j2;	/* range of pixels for getting bkg */
	int xdir, ydir;		/* +1 or -1, indicates search direction */
	int xrange, yrange;	/* over which to increment sums */
	int done;
	int too_far;		/* for check on width of cross correlation */

	/* default values */
	*offsetx = 0.;
	*offsety = 0.;

	/* half the size of the image is a reasonable upper limit */
	too_far = z->nx / 2;

	/* Get an estimate of the background level. */
	j1 = z->ny/2 - z->ny/8;
	j2 = z->ny/2 + z->ny/8;
	i1 = z->nx/2 - z->nx/8;
	i2 = z->nx/2 + z->nx/8;
	sumv = 0.;
	for (j = j1;  j < j2;  j++)
	    for (i = i1;  i < i2;  i++)
		sumv += RPIX2D (z, i, j);
	bkg = sumv / ((j2 - j1) * (i2 - i1));

	/* cutoff = bkg + (peak - bkg) / 3. */
	cutoff = bkg + (RPIX2D (z, lmaxx, lmaxy) - bkg) / 3.;

	/* Start at the peak, and move in the X direction toward the
	   middle of the image until the value is below the cutoff.
	   This gives the range of pixels which we will use for
	   determining the location.
	*/
	if (lmaxx < z->nx / 2)
	    xdir = +1;
	else
	    xdir = -1;
	done = 0;

	for (i = lmaxx;  !done;  i += xdir) {
	    if (xdir * (i - lmaxx) > too_far || i < 0 || i >= z->nx) {
		printf ("Warning:  peak in cross correlation is too broad\n");
		return;
	    }
	    if (RPIX2D (z, i, lmaxy) < cutoff) {
		xrange = xdir * (i - lmaxx);
		done = 1;
	    }
	}
	if (lmaxy < z->ny / 2)
	    ydir = +1;
	else
	    ydir = -1;
	done = 0;
	for (j = lmaxy;  !done;  j += ydir) {
	    if (ydir * (j - lmaxy) > too_far || j < 0 || j >= z->ny) {
		printf ("Warning:  peak in cross correlation is too broad\n");
		return;
	    }
	    if (RPIX2D (z, lmaxx, j) < cutoff) {
		yrange = ydir * (j - lmaxy);
		done = 1;
	    }
	}

	/* Take sums over location of peak plus or minus the range. */

	sumv = 0.;
	sumxv = 0.;
	sumyv = 0.;
	for (j = lmaxy - yrange;  j <= lmaxy + yrange;  j++) {

	    /* periodic boundary condition */
	    if (j < 0)
		jj = j + z->ny;
	    else if (j >= z->ny)
		jj = j - z->ny;
	    else
		jj = j;

	    for (i = lmaxx - xrange;  i <= lmaxx + xrange;  i++) {

		if (i < 0)
		    ii = i + z->nx;
		else if (i >= z->nx)
		    ii = i - z->nx;
		else
		    ii = i;

		value = RPIX2D (z, ii, jj);
		if (value > cutoff) {
		    sumv += value;
		    sumxv += ((i - lmaxx) * value);
		    sumyv += ((j - lmaxy) * value);
		}
	    }
	}

	if (sumv != 0.) {
	    *offsetx = sumxv / sumv;
	    *offsety = sumyv / sumv;
	}
}
