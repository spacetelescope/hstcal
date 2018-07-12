# include <stdio.h>

# include "stis.h"
# include "calstis7.h"
# include "stissizes.h"

/* This routine changes the image size and world coordinate system values
   that describe the output image, if the input image is binned (if CCD)
   or high-res (if MAMA), i.e. if the ltm values differ from one.  The
   parameters that are changed have been read from either the spectroscopic
   distortion table (SDCTAB) or the image distortion table (IDCTAB).
   The output image size and pixel spacing will be scaled by the input
   pixel size, and the reference pixel location will be adjusted.  The
   CRVALi (coordinate values at the reference pixel) will not be changed.

   Note that plate_scale is assigned in sts.

   NOTE:  need to swap axes if dispaxis = 2.

   Phil Hodge, 1997 Sept 12:
	Include input image size in calling sequence;
	reduce output size if input is a subarray.

   Phil Hodge, 2000 Oct 20:
	In the section for modifying the output image size depending on
	subarray, the test for echelle data was obscure, and it didn't work
	for imaging mode.  Change it to a test on sts->first_order.
*/

int ScaleWCS (StisInfo7 *sts, int nx, int ny, CoordInfo *coord_o) {

/* arguments:
StisInfo7 *sts      io: calibration switches and info (plate_scale updated)
int nx, ny           i: size of input image
CoordInfo *coord_o  io: size and coordinate info for output
*/

	double temp[2];
	int fullbinned[2];	/* binned size of a full frame image */
	double plate_scale;	/* before scaling according to binning */

	/* Reduce the output image size by a factor if the input is binned. */
	temp[0] = (double)coord_o->npix[0] * sts->ltm[0];
	temp[1] = (double)coord_o->npix[1] * sts->ltm[1];
	coord_o->npix[0] = NINT (temp[0]);
	coord_o->npix[1] = NINT (temp[1]);

	/* Reduce the output image size if the input is a subarray. */
	if (sts->detector == CCD_DETECTOR) {
	    temp[0] = CCD_NPIX_X * sts->ltm[0];
	    temp[1] = CCD_NPIX_Y * sts->ltm[1];
	} else {
	    temp[0] = MAMA_NPIX_X * sts->ltm[0];
	    temp[1] = MAMA_NPIX_Y * sts->ltm[1];
	}
	fullbinned[0] = NINT (temp[0]);
	fullbinned[1] = NINT (temp[1]);
	if (nx < fullbinned[0])
	    coord_o->npix[0] -= (fullbinned[0] - nx);
	/* The second condition was included because for echelle data the
	   output image size is supposed to be much smaller than the input
	   image size.
	*/
	if (ny < fullbinned[1] && sts->first_order)
	    coord_o->npix[1] -= (fullbinned[1] - ny);

	/* arcseconds per reference pixel */
	if (sts->dispaxis == 1)
	    plate_scale = coord_o->cdelt[1];	/* CDELT2 */
	else if (sts->dispaxis == 2)
	    plate_scale = coord_o->cdelt[0];	/* CDELT1 */

	/* scale the plate scale to image pixels */
	sts->plate_scale[0] = plate_scale / sts->ltm[0];
	sts->plate_scale[1] = plate_scale / sts->ltm[1];

	coord_o->cdelt[0] /= sts->ltm[0];
	coord_o->cdelt[1] /= sts->ltm[1];

	/* convert crpix to image pixel coordinates */
	coord_o->crpix[0] = coord_o->crpix[0] * sts->ltm[0] + sts->ltv[0];
	coord_o->crpix[1] = coord_o->crpix[1] * sts->ltm[1] + sts->ltv[1];

	return (0);
}
