# include <stdio.h>
# include "hstio.h"
# include "acs.h"
# include "acsinfo.h"
# include "hstcalerr.h"		/* for REF_TOO_SMALL */

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
*/

/*
** Developed to support line-by-line operations within CALACS.
** Based on FindBin but operates on 1 line of reference data at a time.
*/
int FindLine (SingleGroup *x, SingleGroupLine *y,
	int *same_size, int *rx, int *ry, int *x0, int *y0) {

/* arguments:
SingleGroup *x    i: science image
SingleGroupLine *y    i: line from reference image
int *same_size    o: true if zero offset and same size and binning
int *rx, *ry      o: ratio of bin sizes
int *x0, *y0      o: location of start of subimage in ref image
*/

	extern int status;

	int sci_bin[2];			/* bin size of science image */
	int sci_corner[2];		/* science image corner location */
	int ref_bin[2];			/* bin size of reference image */
	int ref_corner[2];		/* ref image corner location */
	int rsize;
	int cshift[2];			/* shift of sci relative to ref */
	int ratiox, ratioy;		/* local variables for rx, ry */
	int xzero, yzero;		/* local variables for x0, y0 */
	int GetCorner (Hdr *, int, int *, int *);

	/* Get bin sizes of science and reference images from headers. */
	rsize = 1;
	if (GetCorner (&x->sci.hdr, rsize, sci_bin, sci_corner))
	    return (status);
	if (GetCorner (&y->sci.hdr, rsize, ref_bin, ref_corner))
	    return (status);

	if (sci_corner[0] == ref_corner[0] &&
	    sci_corner[1] == ref_corner[1] &&
	    sci_bin[0] == ref_bin[0] &&
	    sci_bin[1] == ref_bin[1] &&
	    x->sci.data.nx == y->sci.tot_nx) {

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

	    *rx = ref_bin[0] / sci_bin[0];
	    *ry = ref_bin[1] / sci_bin[1];
	    *x0 = (sci_corner[0] - ref_corner[0]) / ref_bin[0];
	    *y0 = (sci_corner[1] - ref_corner[1]) / ref_bin[1];
	    *same_size = 0;

	} else {
        /* For subarray input images, whether they are binned or not... */
	    *same_size = 0;

	    /* ratio of bin sizes */
	    ratiox = sci_bin[0] / ref_bin[0];
	    ratioy = sci_bin[1] / ref_bin[1];
	    if (ratiox * ref_bin[0] != sci_bin[0] ||
		ratioy * ref_bin[1] != sci_bin[1])
		return (status = SIZE_MISMATCH);

	    /* cshift is the offset in units of unbinned
		pixels.  Divide by ref_bin to convert to units of pixels
		in the reference image.
	    */
	    cshift[0] = sci_corner[0] - ref_corner[0];
	    cshift[1] = sci_corner[1] - ref_corner[1];
	    xzero = cshift[0] / ref_bin[0];
	    yzero = cshift[1] / ref_bin[1];
	    if (xzero * ref_bin[0] != cshift[0] ||
		yzero * ref_bin[1] != cshift[1]) {
			trlwarn ("Subimage offset not divisible by bin size.");
	    }
	    *rx = ratiox;
	    *ry = ratioy;
	    *x0 = xzero;
	    *y0 = yzero;
	}

	return (status);
}


/* This function contains the basic algorithms and logic for determining
    the size and position of the input image relative to the reference
    image.  Callable by both FindLineHdr and FindLine.

    *** It still needs to be cleaned up ***
    WJH 30 July 1999

    *** Commented to disable "unused function" warning from compiler ***
    PLL 20 July 2015
*/
#if false
static int getBin (int *sci_corner, int *ref_corner, int *sci_bin, int
*ref_bin, int *same_size, int *rx, int *ry, int *x0, int *y0) {
    extern int status;
	int cshift[2];			/* shift of sci relative to ref */
	int ratiox, ratioy;		/* local variables for rx, ry */
	int xzero, yzero;		/* local variables for x0, y0 */


	if (sci_corner[0] == ref_corner[0] &&
	    sci_corner[1] == ref_corner[1] &&
	    sci_bin[0] == ref_bin[0] &&
	    sci_bin[1] == ref_bin[1] ){
/*	    x->sci.data.nx == y->sci.tot_nx) {  Need to add check on axis sizes  */

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

	    *rx = ref_bin[0] / sci_bin[0];
	    *ry = ref_bin[1] / sci_bin[1];
	    *x0 = (sci_corner[0] - ref_corner[0]) / ref_bin[0];
	    *y0 = (sci_corner[1] - ref_corner[1]) / ref_bin[1];
	    *same_size = 0;

	} else {

	    /* We must bin, extract subset, or both. */
	    *same_size = 0;

	    /* ratio of bin sizes */
	    ratiox = sci_bin[0] / ref_bin[0];
	    ratioy = sci_bin[1] / ref_bin[1];
	    if (ratiox * ref_bin[0] != sci_bin[0] ||
		ratioy * ref_bin[1] != sci_bin[1])
		return (status = SIZE_MISMATCH);

	    /* cshift is the offset in units of unbinned
		pixels.  Divide by ref_bin to convert to units of pixels
		in the reference image.
	    */
	    cshift[0] = sci_corner[0] - ref_corner[0];
	    cshift[1] = sci_corner[1] - ref_corner[1];
	    xzero = cshift[0] / ref_bin[0];
	    yzero = cshift[1] / ref_bin[1];
	    if (xzero * ref_bin[0] != cshift[0] ||
		yzero * ref_bin[1] != cshift[1]) {
			trlwarn ("Subimage offset not divisible by bin size.");
	    }
	    *rx = ratiox;
	    *ry = ratioy;
	    *x0 = xzero;
	    *y0 = yzero;
	}

        return 0;
}
#endif


/*
** Developed to support line-by-line operations within CALACS.
** Based on FindBin but operates on 1 line of reference data at a time.
*/
int FindLineHdr (Hdr *scihdr, Hdr *refhdr, int dimx, int refx,
	int *same_size, int *rx, int *ry, int *x0, int *y0) {

/* arguments:
Hdr *scihdr 	  i: science image header
Hdr *refhdr 	  i: reference image header
int dimx		  i: X dimension of input image
int refx		  i: X dimension of reference image
int *same_size    o: true if zero offset and same size and binning
int *rx, *ry      o: ratio of bin sizes
int *x0, *y0      o: location of start of subimage in ref image
*/

	extern int status;

	int sci_bin[2];			/* bin size of science image */
	int sci_corner[2];		/* science image corner location */
	int ref_bin[2];			/* bin size of reference image */
	int ref_corner[2];		/* ref image corner location */
	int rsize;
	int cshift[2];			/* shift of sci relative to ref */
	int ratiox, ratioy;		/* local variables for rx, ry */
	int xzero, yzero;		/* local variables for x0, y0 */
	int GetCorner (Hdr *, int, int *, int *);

	/* Get bin sizes of science and reference images from headers. */
	rsize = 1;
	if (GetCorner (scihdr, rsize, sci_bin, sci_corner))
	    return (status);
	if (GetCorner (refhdr, rsize, ref_bin, ref_corner))
	    return (status);

	if (sci_corner[0] == ref_corner[0] &&
	    sci_corner[1] == ref_corner[1] &&
	    sci_bin[0] == ref_bin[0] &&
	    sci_bin[1] == ref_bin[1] &&
	    dimx == refx) {

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

	    *rx = ref_bin[0] / sci_bin[0];
	    *ry = ref_bin[1] / sci_bin[1];
	    *x0 = (sci_corner[0] - ref_corner[0]) / ref_bin[0];
	    *y0 = (sci_corner[1] - ref_corner[1]) / ref_bin[1];
	    *same_size = 0;

	} else {

	    /* We must bin, extract subset, or both. */
	    *same_size = 0;

	    /* ratio of bin sizes */
	    ratiox = sci_bin[0] / ref_bin[0];
	    ratioy = sci_bin[1] / ref_bin[1];
	    if (ratiox * ref_bin[0] != sci_bin[0] ||
		ratioy * ref_bin[1] != sci_bin[1])
		return (status = SIZE_MISMATCH);

	    /* cshift is the offset in units of unbinned
		pixels.  Divide by ref_bin to convert to units of pixels
		in the reference image.
	    */
	    cshift[0] = sci_corner[0] - ref_corner[0];
	    cshift[1] = sci_corner[1] - ref_corner[1];
	    xzero = cshift[0] / ref_bin[0];
	    yzero = cshift[1] / ref_bin[1];
	    if (xzero * ref_bin[0] != cshift[0] ||
		yzero * ref_bin[1] != cshift[1]) {
			trlwarn ("Subimage offset not divisible by bin size.");
	    }
	    *rx = ratiox;
	    *ry = ratioy;
	    *x0 = xzero;
	    *y0 = yzero;
	}

	return (status);
}

