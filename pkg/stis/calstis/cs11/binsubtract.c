/* This file contains the following routines:
	LimitsInWav
	BinSubtract
	BinSameBin
	BinSumBin
*/
# include <stdio.h>
# include <math.h>	/* fabs, sqrt */

# include "c_iraf.h"
# include "hstio.h"

# include "stis.h"
# include "calstis11.h"
# include "err.h"
# include "stisdef.h"

/* Add to (or subtract from) corner locations before rounding to nearest
   integer, to avoid roundoff error and to resolve NINT(0.5) ambiguity.

   Phil Hodge, 1998 Oct 5:
	Change status value 1031 to GENERIC_ERROR_CODE.
*/
# define SMALL_AMOUNT  0.001

static void LimitsInWav (SingleGroup *, SingleGroup *,
		double *, double *, double *, double *,
		int *, int *, int *, double *);
static int BinSameBin (SingleGroup *, SingleGroup *, double,
		int *, int *, int *);
static int BinSumBin (SingleGroup *, SingleGroup *, double,
		int *, int *, int *, double *);

/* This routine does the actual subtraction of the science image
   from the wavecal.  If the science image is binned differently
   or is a different subimage from the wavecal, this is taken into
   account.

   The ERR extension is not modified.  If the science image is binned
   differently from the wavecal, the DQ extension of the wavecal will be
   updated.

	mapping from wavecal pixel to science data pixel
	xS = (xW - ltvW) * ltmS / ltmW + ltvS
*/

int BinSubtract (SingleGroup *wav, SingleGroup *sci, double ratio,
		int verbose) {

/* arguments:
SingleGroup *wav    io: input is wavecal; output is wavecal minus science data
SingleGroup *sci    i: science data
double ratio        i: ratio of exposure times, wavecal to science
int verbose         i: print info about first pixel and binning?
*/

	int status;

	double ltmW[2], ltvW[2];	/* for wavecal */
	double ltmS[2], ltvS[2];	/* for science data */
	/* X & Y coordinates of the lower left corner and upper right corner
	   of the science image mapped into the wavecal, but truncated at
	   the borders of the wavecal in the case that the science image
	   extends outside the wavecal.
	*/
	int lower_left[2], upper_right[2];
	int sci_start[2];	/* first sci pixel corresp. to lower_left */
	double bin[2];		/* number of science pixels per wavecal pixel */

	/* Get the linear transformation from reference coordinates
	   to image pixel coordinates for the wavecal and science data.
	*/
	if ((status = GetLT0 (&wav->sci.hdr, ltmW, ltvW))) /* zero indexed LTV */
	    return (status);
	if ((status = GetLT0 (&sci->sci.hdr, ltmS, ltvS)))
	    return (status);

	if (ltmW[0] <= 0. || ltmS[0] <= 0. ||
	    ltmW[1] <= 0. || ltmS[1] <= 0.) {
	    printf ("ERROR    Can't handle LTM less than or equal to zero.\n");
	    return (GENERIC_ERROR_CODE);
	}

	/* Get the location of the science image within the wavecal image,
	   in units of pixel coordinates in the wavecal.  Also check that
	   the science image covers the wavecal image, and print a warning
	   if it does not.
	*/
	LimitsInWav (wav, sci, ltmW, ltvW, ltmS, ltvS,
		lower_left, upper_right, sci_start, bin);

	if (verbose) {
	    if (lower_left[0] > 0 || lower_left[1] > 0 ||
		upper_right[0] < wav->sci.data.nx - 1 ||
		upper_right[1] < wav->sci.data.ny - 1) {
		printf (
		"         Science image covers this subset of wavecal: \\\n");
		printf (
		"         x:  %d to %d,  y:  %d to %d\n",
			lower_left[0]+1, upper_right[0]+1,
			lower_left[1]+1, upper_right[1]+1);
	    } else {
		printf ("         Science image covers entire wavecal.\n");
	    }
	    if (sci_start[0] > 0 || sci_start[1] > 0) {
		printf ("         Starting pixel in science image = %d, %d\n",
			sci_start[0]+1, sci_start[1]+1);
	    }
	    printf (
"         Number of science image pixels per wavecal pixel = %.6g, %.6g\n",
			bin[0], bin[1]);
	}

	if (ltmW[0] == ltmS[0] && ltmW[1] == ltmS[1]) {

	    /* Same binning, possibly different subset. */
	    if ((status = BinSameBin (wav, sci, ratio,
                                      lower_left, upper_right, sci_start)))
		return (status);

	} else {

	    /* different binning */
	    if ((status = BinSumBin (wav, sci, ratio,
                                     lower_left, upper_right, sci_start, bin)))
		return (status);
	}

	return (0);
}

/* This routine maps the lower left and upper right corner pixels of the
   science image into the pixel coordinate system of the wavecal, truncating
   at the edges of the wavecal.  The start and end pixels in the wavecal
   corresponding to the corners of the science image are returned; these
   can be used as loop limits in the wavecal when subtracting the science
   image.  A warning is printed if the science image does not cover the
   wavecal.

   After finding lower_left (in the wavecal), the lower left corner of
   that pixel is mapped back into the science image and rounded to an
   integer to find the starting pixel in the science image.  The ratio
   of bin sizes is also computed, to be used to determine how many
   science pixels to sum over when for one wavecal pixel, or vice versa.

   Using zero indexing with N pixels, the pixel numbers run from 0 to N-1.
   In the IRAF convention, the integer value is the center of the pixel,
   so the outer edges of the image are at -0.5 and N-0.5.

	xW = (xS - ltvS) * ltmW / ltmS + ltvW
*/

static void LimitsInWav (SingleGroup *wav, SingleGroup *sci,
	double *ltmW, double *ltvW, double *ltmS, double *ltvS,
	int *lower_left, int *upper_right, int *sci_start, double *bin) {

/* arguments:
SingleGroup *wav         i: wavecal
SingleGroup *sci         i: science data
double ltmW[2], ltvW[2]  i: linear transformation for wavecal to ref coords
double ltmS[2], ltvS[2]  i: linear transformation for science to ref coords
int lower_left[2]        o: lower left corner of science image in wavecal
int upper_right[2]       o: upper right corner of science image in wavecal
int sci_start[2]         o: first sci pixel corresponding to lower_left
double bin[2]            o: number of sci pixels per wavecal pixel
*/

	double pixel;		/* a corner in science image */
	double ll[2], ur[2];	/* lower left and upper right corners */

	/* left edge of first science image pixel */
	pixel = -0.5;
	ll[0] = (pixel - ltvS[0]) * ltmW[0] / ltmS[0] + ltvW[0];

	/* lower edge of first science image pixel */
	pixel = -0.5;
	ll[1] = (pixel - ltvS[1]) * ltmW[1] / ltmS[1] + ltvW[1];

	/* right edge of last science image pixel */
	pixel = sci->sci.data.nx - 0.5;
	ur[0] = (pixel - ltvS[0]) * ltmW[0] / ltmS[0] + ltvW[0];

	/* upper edge of last science image pixel */
	pixel = sci->sci.data.ny - 0.5;
	ur[1] = (pixel - ltvS[1]) * ltmW[1] / ltmS[1] + ltvW[1];

	/* Strink ranges by a small amount to allow for roundoff. */
	ll[0] += SMALL_AMOUNT;
	ll[1] += SMALL_AMOUNT;
	ur[0] -= SMALL_AMOUNT;
	ur[1] -= SMALL_AMOUNT;

	/* Round off to nearest integer, and truncate at image edges. */
	lower_left[0] = NINT(ll[0]);
	lower_left[1] = NINT(ll[1]);
	upper_right[0] = NINT(ur[0]);
	upper_right[1] = NINT(ur[1]);
	if (lower_left[0] < 0)
	    lower_left[0] = 0;
	if (lower_left[1] < 0)
	    lower_left[1] = 0;
	if (upper_right[0] > wav->sci.data.nx - 1)
	    upper_right[0] = wav->sci.data.nx - 1;
	if (upper_right[1] > wav->sci.data.ny - 1)
	    upper_right[1] = wav->sci.data.ny - 1;

	/* Is science image too small? */
	if (lower_left[0] > 0 || upper_right[0] < wav->sci.data.nx - 1 ||
	    lower_left[1] > 0 || upper_right[1] < wav->sci.data.ny - 1) {

	    printf (
	"Warning  The science image does not cover the wavecal image.\n");
	}

	/* Now that we've found the starting pixel in the wavecal, map
	   the corner of that one pixel back into the science image.
	*/

	/* left edge */
	pixel = lower_left[0] - 0.5 + SMALL_AMOUNT;
	ll[0] = (pixel - ltvW[0]) * ltmS[0] / ltmW[0] + ltvS[0];

	/* lower edge */
	pixel = lower_left[1] - 0.5 + SMALL_AMOUNT;
	ll[1] = (pixel - ltvW[1]) * ltmS[1] / ltmW[1] + ltvS[1];

	sci_start[0] = NINT (ll[0]);
	sci_start[1] = NINT (ll[1]);

	bin[0] = ltmS[0] / ltmW[0];
	bin[1] = ltmS[1] / ltmW[1];
}

/* This routine is used for the case that the science data array and the
   wavecal are binned the same.  The zero points don't have to be the
   same, however.  Only the SCI extension values are modified.
*/

static int BinSameBin (SingleGroup *wav, SingleGroup *sci, double ratio,
	int *lower_left, int *upper_right, int *sci_start) {

/* arguments:
SingleGroup *wav         io: wavecal
SingleGroup *sci         i: science data
double ratio             i: ratio of exposure times
int lower_left[2]        i: lower left corner of science image in wavecal
int upper_right[2]       i: upper right corner of science image in wavecal
int sci_start[2]         i: first sci pixel corresponding to lower_left
*/

	int off[2];		/* offsets of science from wavecal */
	int i, j;		/* loop indices in wavecal */

	/* Simple expression because the binning is the same. */
	off[0] = sci_start[0] - lower_left[0];
	off[1] = sci_start[1] - lower_left[1];

	/* SCI extension (don't modify ERR or DQ extension) */
	for (j = lower_left[1];  j <= upper_right[1];  j++) {
	    for (i = lower_left[0];  i <= upper_right[0];  i++) {
		Pix (wav->sci.data, i, j) = Pix (wav->sci.data, i, j) -
		ratio * Pix (sci->sci.data, i+off[0], j+off[1]);
	    }
	}

	return (0);
}

static int BinSumBin (SingleGroup *wav, SingleGroup *sci, double ratio,
	int *lower_left, int *upper_right, int *sci_start, double *bin) {

	double x0, y0;		/* starting pixel in sci before roundoff */
	int i0, j0;		/* actual starting pixel in science */
	int m, n;		/* loop indices in wavecal */
	int i, j;		/* loop indices in science image */
	int range[2];		/* range of sci pixels to sum over */
	double sum;		/* for summing pixel values */
	short dqsum;		/* combined data quality */
	float combined_ratio;	/* includes exposure time and binning */

	/* If the science image is more coarsely binned than the wavecal,
	   we need to take only a fraction of the science pixel value.
	   If it's the wavecal that's more coarsely binned, however,
	   we already take that into consideration when we sum the pixel
	   values.
	   Also set range = bin, with 1 as a lower limit.
	*/
	combined_ratio = ratio;
	for (i = 0;  i < 2;  i++) {
	    if (bin[i] < 1.) {
		combined_ratio *= bin[i];
		range[i] = 1;
	    } else {
		range[i] = NINT (bin[i]);
	    }
	}

	/* SCI and DQ extensions (don't modify ERR) */

	j0 = sci_start[1];
	y0 = (double) j0 + SMALL_AMOUNT;
	/* for each line of wavecal ... */
	for (n = lower_left[1];  n <= upper_right[1];  n++) {
	    i0 = sci_start[0];
	    x0 = (double) i0 + SMALL_AMOUNT;
	    /* for each sample in current wavecal line ... */
	    for (m = lower_left[0];  m <= upper_right[0];  m++) {
		sum = 0.;
		dqsum = 0;
		for (j = j0;  j < j0+range[1];  j++) {
		    for (i = i0;  i < i0+range[0];  i++) {
			sum += Pix (sci->sci.data, i, j);
			dqsum |= DQPix (wav->dq.data, m, n);
		    }
		}
		Pix (wav->sci.data, m, n) = Pix (wav->sci.data, m, n) -
			combined_ratio * sum;
		DQSetPix (wav->dq.data, m, n, dqsum);
		x0 += bin[0];
		i0 = (int) x0;
	    }
	    y0 += bin[1];
	    j0 = (int) y0;
	}

	return (0);
}
