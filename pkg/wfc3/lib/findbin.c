# include <stdio.h>
#include "hstcal.h"
# include "hstio.h"
# include "trlbuf.h"		/* for trlwarn */
# include "hstcalerr.h"		/* for REF_TOO_SMALL */

# include <string.h>
# include "wf3.h"

/* This routine finds the bin factors and corner location to use
   when calling bin2d to extract and bin the appropriate subset of
   a reference image to match a science image.
   
   If the science image has zero offset and is the same size and
   binning as the reference image, same_size will be set to true;
   otherwise, the values of rx, ry, x0, y0 will be assigned.
   Normally the science image (x) will be binned the same as
   the reference image (y).

   If the binning of the reference image is greater than the binning
   of the science image, then the ratios (rx and ry) of the bin sizes
   will be the reference image (y) bin size divided by the science
   image (x) bin size. This is not necessarily an error.

   Revision history:

   Howard Bushouse, 2000 Dec 7:
   Added error checks for WFC3 in the case where the science image
   is binned more than the reference image.

   H. Bushouse, 2002 Mar 18:
   Added the FindBinIR routine (copy of FindBin) to accomodate use
   with WFC3 IR images.
 */

/* 
 ** Developed to support line-by-line operations within CALACS.  
 ** Based on calstis FindBin but operates on 1 line of reference data
 ** at a time.
 */
int FindLine (SingleGroup *x, SingleGroupLine *y,
        int *same_size, int *rx, int *ry, int *x0, int *y0) {

    /* arguments:
       SingleGroup *x		i: science image
       SingleGroupLine *y	i: line from reference image
       int *same_size		o: true if zero offset and same size and binning
       int *rx, *ry		o: ratio of bin sizes
       int *x0, *y0		o: location of start of subimage in ref image
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
    xzero=0;
    yzero=0;
    
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
           and without extracting a subset.  */
        *same_size = 1;
        *rx = 1;
        *ry = 1;
        *x0 = 0;
        *y0 = 0;

    } else if (ref_bin[0] < sci_bin[0] ||
            ref_bin[1] < sci_bin[1]) {

        /* Science image is binned more than reference image: This
         ** is not allowed for WFC3. */
        *same_size = 0;

        *rx = sci_bin[0] / ref_bin[0];
        *ry = sci_bin[1] / ref_bin[1];
        *x0 = 0;
        *y0 = 0;
        trlwarn("Science image binned more than reference image");
    } else if (ref_bin[0] > sci_bin[0] ||
            ref_bin[1] > sci_bin[1]) {

        /* Reference image is binned more than science image. For WFC3
         ** this could happen for the shutter-shading correction file,
         ** and the low-order flat field. */
        *same_size = 0;

        *rx = ref_bin[0] / sci_bin[0];
        *ry = ref_bin[1] / sci_bin[1];
        *x0 = (sci_corner[0] - ref_corner[0]) / ref_bin[0];
        *y0 = (sci_corner[1] - ref_corner[1]) / ref_bin[1];
        trlwarn("Reference image binned more than science image");
    } else {

        /* For subarray input images, whether they are binned or not... */
        *same_size = 0;
        

        /* ratio of bin sizes */
        ratiox = sci_bin[0] / ref_bin[0];
        ratioy = sci_bin[1] / ref_bin[1];
        if (ratiox * ref_bin[0] != sci_bin[0] ||
                ratioy * ref_bin[1] != sci_bin[1]){
            status = SIZE_MISMATCH;    
            return (status);
        }
        
        /* cshift is the offset in units of unbinned
           pixels.  Divide by ref_bin to convert to units of pixels
           in the reference image. sci corner or ref_corner can be negative if in physical overscan
         */
        cshift[0] = sci_corner[0] - ref_corner[0];
        cshift[1] = sci_corner[1] - ref_corner[1];
        
        xzero = cshift[0] / ref_bin[0];
        yzero = cshift[1] / ref_bin[1];
        if (xzero * ref_bin[0] != cshift[0] ||
                yzero * ref_bin[1] != cshift[1]) {
            trlwarn ("Subimage offset not divisible by bin size.");
        }
        sprintf(MsgText,"Subimage locations rx=%d, ry=%d, x0=%d, y0=%d",ratiox,ratioy,xzero,yzero);
        trlmessage(MsgText);
        *rx = ratiox;
        *ry = ratioy;
        *x0 = xzero;
        *y0 = yzero;
    }

    return (status);
}

/* 
 ** Developed to support line-by-line operations within CALACS.  
 ** Based on FindBin but operates on 1 line of reference data at a time.
 */
int FindLineHdr (Hdr *scihdr, Hdr *refhdr, int dimx, int refx,
        int *same_size, int *rx, int *ry, int *x0, int *y0) {

    /* arguments:
       Hdr *scihdr	i: science image header
       Hdr *refhdr	i: reference image header
       int dimx	i: X dimension of input image
       int refx	i: X dimension of reference image
       int *same_size	o: true if zero offset and same size and binning
       int *rx, *ry	o: ratio of bin sizes
       int *x0, *y0	o: location of start of subimage in ref image
     */

    extern int status;

    int sci_bin[2];			/* bin size of science image */
    int sci_corner[2];		/* science image corner location */
    int ref_bin[2];			/* bin size of reference image */
    int ref_corner[2];		/* ref image corner location */
    int rsize = 1;
    int cshift[2];			/* shift of sci relative to ref */
    int ratiox, ratioy;		/* local variables for rx, ry */
    int xzero, yzero;		/* local variables for x0, y0 */
    int GetCorner (Hdr *, int, int *, int *);

    /* Get bin sizes of science and reference images from headers. */
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
         ** and without extracting a subset.  */
        *same_size = 1;
        *rx = 1;
        *ry = 1;
        *x0 = 0;
        *y0 = 0;

    } else if (ref_bin[0] < sci_bin[0] ||
            ref_bin[1] < sci_bin[1]) {

        /* Science image is binned more than reference image: This
         ** is not allowed for WFC3. */
        *same_size = 0;

        *rx = sci_bin[0] / ref_bin[0];
        *ry = sci_bin[1] / ref_bin[1];
        status = SIZE_MISMATCH;
        return (status);

    } else if (ref_bin[0] > sci_bin[0] ||
            ref_bin[1] > sci_bin[1]) {

        /* Reference image is binned more than the science image. */
        *same_size = 0;

        *rx = ref_bin[0] / sci_bin[0];
        *ry = ref_bin[1] / sci_bin[1];
        *x0 = (sci_corner[0] - ref_corner[0]) / ref_bin[0];
        *y0 = (sci_corner[1] - ref_corner[1]) / ref_bin[1];

    } else {

        /* We must extract subset. */
        *same_size = 0;

        /* ratio of bin sizes */
        ratiox = sci_bin[0] / ref_bin[0];
        ratioy = sci_bin[1] / ref_bin[1];
        if (ratiox * ref_bin[0] != sci_bin[0] ||
                ratioy * ref_bin[1] != sci_bin[1]){
            status = SIZE_MISMATCH;
            return (status);
        }

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


int FindBin (SingleGroup *x, SingleGroup *y, int *same_size,
        int *rx, int *ry, int *x0, int *y0) {

    /* arguments:
       SingleGroup *x    i: science image
       SingleGroup *y    i: reference image
       int *same_size    o: true if zero offset and same size and binning
       int *rx, *ry      o: ratio of bin sizes
       int *x0, *y0      o: location of start of subimage in ref image
     */

    extern int status;

    int sci_bin[2];			/* bin size of science image */
    int sci_corner[2];		/* science image corner location */
    int ref_bin[2];			/* bin size of reference image */
    int ref_corner[2];		/* ref image corner location */
    int rsize = 1;
    int cshift[2];			/* shift of sci relative to ref */
    int ratiox, ratioy;		/* local variables for rx, ry */
    int xzero, yzero;		/* local variables for x0, y0 */
    int GetCorner (Hdr *, int, int *, int *);
    
    xzero=0;
    yzero=0;
    
    /* Get bin sizes of science and reference images from headers. */
    if (GetCorner (&x->sci.hdr, rsize, sci_bin, sci_corner))
        return (status);
    if (GetCorner (&y->sci.hdr, rsize, ref_bin, ref_corner))
        return (status);

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
        status = REF_TOO_SMALL;
        return (status);

    } else {

        /* We must bin, extract subset, or both. */
        *same_size = 0;

        /* ratio of bin sizes */
        ratiox = sci_bin[0] / ref_bin[0];
        ratioy = sci_bin[1] / ref_bin[1];
        if (ratiox * ref_bin[0] != sci_bin[0] ||
                ratioy * ref_bin[1] != sci_bin[1]){
            status = SIZE_MISMATCH;
            return (status);
        }

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
            trlwarn ("Subimage offset not divisible by bin size.");
        }
        *rx = ratiox;
        *ry = ratioy;
        *x0 = xzero;
        *y0 = yzero;
    }

    return (status);
}


int FindBinIR (SingleNicmosGroup *x, SingleNicmosGroup *y, int *same_size,
        int *rx, int *ry, int *x0, int *y0) {

    /* arguments:
       SingleGroup *x    i: science image
       SingleGroup *y    i: reference image
       int *same_size    o: true if zero offset and same size and binning
       int *rx, *ry      o: ratio of bin sizes
       int *x0, *y0      o: location of start of subimage in ref image
     */

    extern int status;

    int sci_bin[2];			/* bin size of science image */
    int sci_corner[2];		/* science image corner location */
    int ref_bin[2];			/* bin size of reference image */
    int ref_corner[2];		/* ref image corner location */
    int rsize = 1;
    int cshift[2];			/* shift of sci relative to ref */
    int ratiox, ratioy;		/* local variables for rx, ry */
    int xzero, yzero;		/* local variables for x0, y0 */
    int GetCorner (Hdr *, int, int *, int *);

    /* Get bin sizes of science and reference images from headers. */
    if (GetCorner (&x->sci.hdr, rsize, sci_bin, sci_corner))
        return (status);
    if (GetCorner (&y->sci.hdr, rsize, ref_bin, ref_corner))
        return (status);

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
        status = REF_TOO_SMALL;
        return (status);

    } else {

        /* We must bin, extract subset, or both. */
        *same_size = 0;

        /* ratio of bin sizes */
        ratiox = sci_bin[0] / ref_bin[0];
        ratioy = sci_bin[1] / ref_bin[1];
        if (ratiox * ref_bin[0] != sci_bin[0] ||
                ratioy * ref_bin[1] != sci_bin[1]) {
            status = SIZE_MISMATCH;  
            return (status);
            
        }

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
            trlwarn ("Subimage offset not divisible by bin size.");
        }
        *rx = ratiox;
        *ry = ratioy;
        *x0 = xzero;
        *y0 = yzero;
    }

    return (status);
}
