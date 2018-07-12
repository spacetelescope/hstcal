# include <stdio.h>
# include "stis.h"
# include "calstis1.h"
# include "hstcalerr.h"
# include "stissizes.h"	/* detector sizes, overscan sizes */

# define FULL_FRAME_READOUT  1
# define SUBARRAY_READOUT    2
# define BINNED_READOUT      3

/* This routine assigns the widths of the overscan region on each edge
   of the image by assuming particular region sizes, depending on the
   binning and which amplifier was used for readout.

	    parallel readout direction, amp A
	    |
	    |
	    v

	C               D
	  -------------
	 |             |   <--- serial readout direction, amp A
	 |             |
	 |             |
	 |             |
	 |             |
	 |             |
	  -------------
	A               B

   Phil Hodge, 1997 Oct 28:
	Remove ox1, ox2, oy1, oy2 from calling sequence;
	compute vx & vy and add them to the calling sequence.

   Phil Hodge, 1998 Oct 5:
	Change status value 1231 to GENERIC_ERROR_CODE.
*/

int FindOverscan (char *ccdamp, int nx, int ny, int *bin,
		int *trimx1, int *trimx2, int *trimy1, int *trimy2,
		int *biassect, int *vx, int *vy) {

/* arguments:
char *ccdamp          i: amplifier that was used (A-D)
int nx, ny            i: lengths of first and second image axes
int bin[2]            i: bin factor for each axis
int *trimx1, *trimx2  o: amount to trim from each end of a line
int *trimy1, *trimy2  o: amount to trim from each end of a column
int biassect[2]       o: section to use for finding bias level
int vx[2], vy[2]      o: range of pixel numbers for virtual overscan region
*/

	int readout;			/* full-frame, subarray, or binned */
	int bad = 0;			/* flag for bad image size */

	/* sizes of overscan regions, if amp A was used */
	int A_left, A_right, A_bottom, A_top;
	int A_bsect[2];			/* biassect if amp A used */

	/* sizes of overscan regions for amp actually used */
	int left, right, bottom, top;


	/* The overscan sizes depend on whether the image is full-frame,
	   subarray, or binned.  These are mutually exclusive.
	*/
	if (bin[0] == 1 && bin[1] == 1) {
	    if (ny <= CCD_NPIX_Y)
		readout = SUBARRAY_READOUT;
	    else
		readout = FULL_FRAME_READOUT;
	} else {
	    readout = BINNED_READOUT;
	}

	/* Assign overscan sizes to local variables, depending on readout.
	   In the next section we'll account for which amp was used.
	*/

	if (readout == FULL_FRAME_READOUT) {

	    A_left = OVERSCAN_LEFT;
	    A_right = OVERSCAN_RIGHT;
	    A_top = OVERSCAN_TOP;
	    A_bottom = 0;
	    A_bsect[0] = BIASSECT_1L;
	    A_bsect[1] = BIASSECT_1R;

	} else if (readout == SUBARRAY_READOUT) {

	    A_left = OVERSCAN_LEFT - SUBARRAY_BORDER;
	    A_right = OVERSCAN_RIGHT - SUBARRAY_BORDER;
	    A_top = 0;
	    A_bottom = 0;
	    A_bsect[0] = BIASSECT_1L;
	    A_bsect[1] = BIASSECT_1R - SUBARRAY_BORDER;

	} else if (readout == BINNED_READOUT) {

	    /* Add one to OVERSCAN_LEFT to include one mixed pixel,
		overscan and illuminated.  Subtract one on each end
		of CCD_NPIX_X to exclude these mixed pixels.
	    */
	    if (bin[0] == 1) {
		/* nx = 19 + 1024 + 1 + 10, with no mixed pixel */
		A_left = OVERSCAN_LEFT;
		A_right = nx - CCD_NPIX_X - A_left;
	    } else {
		A_left = (OVERSCAN_LEFT + 1) / bin[0];
		A_right = nx - (CCD_NPIX_X / bin[0] - 1) - A_left;
	    }

	    A_top = ny - CCD_NPIX_Y / bin[1];
	    A_bottom = 0;

	    if (bin[0] == 1) {
		A_bsect[0] = BIASSECT_1L_B;	/* "binned" by one */
		A_bsect[1] = BIASSECT_1R_B;
	    } else if (bin[0] == 2) {
		A_bsect[0] = BIASSECT_2L;
		A_bsect[1] = BIASSECT_2R;
	    } else if (bin[0] == 4) {
		A_bsect[0] = BIASSECT_4L;
		A_bsect[1] = BIASSECT_4R;
	    } else {
		A_bsect[0] = 0;	/* invalid */
		A_bsect[1] = -1;
	    }
	}

	/* Now check which amp was used for readout, and rearrange the
	   locations of the overscan regions accordingly.
	*/
	if (ccdamp[0] == 'A') {
	    left = A_left;
	    right = A_right;
	    bottom = A_bottom;
	    top = A_top;
	    biassect[0] = A_bsect[0];
	    biassect[1] = A_bsect[1];
	} else if (ccdamp[0] == 'B') {
	    left = A_right;
	    right = A_left;
	    bottom = A_bottom;
	    top = A_top;
	    biassect[0] = (nx - A_bsect[1] - 1);
	    biassect[1] = (nx - A_bsect[0] - 1);
	} else if (ccdamp[0] == 'C') {
	    left = A_left;
	    right = A_right;
	    bottom = A_top;
	    top = A_bottom;
	    biassect[0] = A_bsect[0];
	    biassect[1] = A_bsect[1];
	} else if (ccdamp[0] == 'D') {
	    left = A_right;
	    right = A_left;
	    bottom = A_top;
	    top = A_bottom;
	    biassect[0] = (nx - A_bsect[1] - 1);
	    biassect[1] = (nx - A_bsect[0] - 1);
	}

	/* Sanity check on image size. */
	if (readout == FULL_FRAME_READOUT) {
	    if ((nx - left - right) < CCD_NPIX_X ||
		(ny - top - bottom) < CCD_NPIX_Y) {
		bad = 1;
	    }
	} else if (readout == SUBARRAY_READOUT) {
	    if ((nx - left - right) < CCD_NPIX_X)
		bad = 1;
	}
	if (left < 0 || right < 0 || bottom < 0 || top < 0)
	    bad = 1;
	if (bad)		/* image is too small to do blevcorr */
	    return (GENERIC_ERROR_CODE);

	*trimx1 = left;
	*trimx2 = right;
	*trimy1 = bottom;
	*trimy2 = top;

	/* This is the region within the virtual overscan to use for
	   finding the drift (slope) of the bias level along the lines.
	   Currently we use the entire virtual overscan region that is
	   over the illuminated portion of the chip.
	*/
	vx[0] = *trimx1;
	vx[1] = nx - 1 - *trimx2;
	if (ccdamp[0] == 'A' || ccdamp[0] == 'B') {
	    vy[0] = ny - *trimy2;
	    vy[1] = ny - 1;
	} else if (ccdamp[0] == 'C' || ccdamp[0] == 'D') {
	    vy[0] = 0;
	    vy[1] = *trimy1 - 1;
	}

	return (0);
}
