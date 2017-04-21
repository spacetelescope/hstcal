# include <stdio.h>

# include "c_iraf.h"
# include "hstio.h"

# include "stis.h"
# include "calstis1.h"
# include "err.h"		/* for REF_TOO_SMALL */
# include "stisdef.h"

/* This routine finds the bin factors and corner location to use
   when calling bin2d to extract and bin the appropriate subset of
   a reference image to match a science image.
   If the science image has zero offset and is the same size and
   binning as the reference image, same_size will be set to true;
   otherwise, the values of rx, ry, x0, y0 will be assigned.
   Normally the science image (x) will be binned the same or more
   than the reference image (y); rx and ry will be the bin size of
   the science image divided by the bin size of the reference image,

   If the binning of the reference image is greater than the binning
   of the science image, then status = REF_TOO_SMALL will be returned.
   In this case, the ratios (rx and ry) of the bin sizes will be the
   reference image (y) bin size divided by the science image (x) bin size.
   This is not necessarily an error.

   The high_res flag will be set to true if the bin size of the
   reference image is one in the dispersion direction.  This is
   really needed only for MAMA data where the dispersion is high
   enough that Doppler convolution must be applied to the reference
   data.

   Phil Hodge, 1998 Sept 24:
	rx, ry, x0, y0 were not assigned for the case that the two images
	are the same size and offset; same_size was not assigned for the
	case that the reference image is binned more than the science image.
*/

int FindBin (StisInfo1 *sts, SingleGroup *x, SingleGroup *y,
	int *same_size, int *high_res,
	int *rx, int *ry, int *x0, int *y0) {

/* arguments:
StisInfo1 *sts    i: calibration switches and info
SingleGroup *x    i: science image
SingleGroup *y    i: reference image
int *same_size    o: true if zero offset and same size and binning
int *high_res     o: true if bin size is one in dispersion direction
int *rx, *ry      o: ratio of bin sizes
int *x0, *y0      o: location of start of subimage in ref image
*/

	int status;

	int sci_bin[2];			/* bin size of science image */
	int sci_corner[2];		/* science image corner location */
	int ref_bin[2];			/* bin size of reference image */
	int ref_corner[2];		/* ref image corner location */
	int rsize;			/* 1 for CCD, 2 for MAMA */
	int cshift[2];			/* shift of sci relative to ref */
	int ratiox, ratioy;		/* local variables for rx, ry */
	int xzero, yzero;		/* local variables for x0, y0 */

	/* Get bin sizes of science and reference images from headers. */
	rsize = (sts->detector == CCD_DETECTOR) ? 1 : 2;
	if ((status = GetCorner (&x->sci.hdr, rsize, sci_bin, sci_corner)))
	    return (status);
	if ((status = GetCorner (&y->sci.hdr, rsize, ref_bin, ref_corner)))
	    return (status);

	/* High-res is really only relevant for the MAMA detectors. */
	*high_res = (sts->dispaxis == 1 && ref_bin[0] == 1) ||
                    (sts->dispaxis == 2 && ref_bin[1] == 1);

	if (sci_corner[0] == ref_corner[0] &&
	    sci_corner[1] == ref_corner[1] &&
	    sci_bin[0] == ref_bin[0] &&
	    sci_bin[1] == ref_bin[1] &&
	    x->sci.data.nx == y->sci.data.nx &&
	    x->sci.data.ny == y->sci.data.ny) {

	    /* We can use the reference image directly, without binning
		and without extracting a subset.
	    */
	    *same_size = 1;
	    *rx = 1;
	    *ry = 1;
	    *x0 = 0;
	    *y0 = 0;

	} else if (ref_bin[0] > sci_bin[0] ||
		   ref_bin[1] > sci_bin[1]) {

	    /* Reference image is binned more than the science image. */
	    *same_size = 0;

	    *rx = ref_bin[0] / sci_bin[0];
	    *ry = ref_bin[1] / sci_bin[1];
	    *x0 = (sci_corner[0] - ref_corner[0]) / ref_bin[0];
	    *y0 = (sci_corner[1] - ref_corner[1]) / ref_bin[1];

	    return (REF_TOO_SMALL);

	} else {

	    /* We must bin, extract subset, or both. */
	    *same_size = 0;

	    /* ratio of bin sizes */
	    ratiox = sci_bin[0] / ref_bin[0];
	    ratioy = sci_bin[1] / ref_bin[1];
	    if (ratiox * ref_bin[0] != sci_bin[0] ||
		ratioy * ref_bin[1] != sci_bin[1])
		return (SIZE_MISMATCH);

	    /* cshift is the offset in units of unbinned (or low-res)
		pixels.  Divide by ref_bin to convert to units of pixels
		in the reference image.
	    */
	    cshift[0] = sci_corner[0] - ref_corner[0];
	    cshift[1] = sci_corner[1] - ref_corner[1];
	    xzero = cshift[0] / ref_bin[0];
	    yzero = cshift[1] / ref_bin[1];
	    if (xzero * ref_bin[0] != cshift[0] ||
		yzero * ref_bin[1] != cshift[1]) {
		printf (
		"Warning  Subimage offset not divisible by bin size.\n");
	    }
	    *rx = ratiox;
	    *ry = ratioy;
	    *x0 = xzero;
	    *y0 = yzero;
	}

	return (0);
}
