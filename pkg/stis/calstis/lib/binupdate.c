# include <stdio.h>

# include "stis.h"
# include "err.h"

/* This routine modifies some of the coordinate parameters in place
   to account for a linear transformation, a change in scale and
   a shift of the origin.  The parameters that are modified are
   crpix, the cd matrix, the diagonal elements of ltm, and ltv.

   The bin factor must be positive.  The offset may be positive, negative,
   or zero.  When an output image is a binned subset of an input image
   (like blkavg), we will have bin > 1, and offset >= 0.  When an output
   image is expanded (like blkrep), we will have bin < 1, and offset can
   be <= 0.

   Note that offset is in units of the smaller pixels.  That is,
   when binning (bin > 1), offset is in units of the input pixels;
   when unbinning (bin < 1), offset is in units of the output pixels.

   The CD matrix is mapped to the 1-D array cd in the following order:
	CD1_1, CD1_2, CD2_1, CD2_2

   Phil Hodge, 1998 Oct 5:
	Change status value 1991 to INTERNAL_ERROR.
*/

int BinUpdate (double *block, double *offset,
	double *ltm, double *ltv, double *cd, double *crpix) {

/* arguments:
double block[2]         i: bin factor
double offset[2]        i: offset of binned image
double ltm[2], ltv[2]   io: transformation from reference to image pixels
double cd[4]            io: CD matrix
double crpix[2]         io: reference pixel
*/

	/* If we're actually unbinning (blkrep rather than blkavg), then
	   the offset that we use to compute the new parameters needs
	   to be multiplied by the bin factor (which is less than one)
	   so it will be in units of the input image rather than in units
	   of the image with smaller pixels.  So we need a local variable
	   for a copy of offset, possibly scaled by bin.
	*/
	double off[2];

	if (block[0] == 0. || block[1] == 0.) {
	    printf ("ERROR    (binupdate) block size of zero\n");
	    return (INTERNAL_ERROR);
	}

	if (block[0] > 1.)
	    off[0] = offset[0];
	else
	    off[0] = offset[0] * block[0];
	cd[0] *= block[0];		/* CD1_1 */
	cd[2] *= block[0];		/* CD2_1 */
	ltm[0] /= block[0];
	ltv[0] = (ltv[0] - offset[0] + (block[0] - 1.) / 2.) / block[0];
	crpix[0] = (crpix[0] - offset[0] + (block[0] - 1.) / 2.) / block[0];

	if (block[1] > 1.)
	    off[1] = offset[1];
	else
	    off[1] = offset[1] * block[1];
	cd[1] *= block[1];		/* CD1_2 */
	cd[3] *= block[1];		/* CD2_2 */
	ltm[1] /= block[1];
	ltv[1] = (ltv[1] - offset[1] + (block[1] - 1.) / 2.) / block[1];
	crpix[1] = (crpix[1] - offset[1] + (block[1] - 1.) / 2.) / block[1];

	return (0);
}
