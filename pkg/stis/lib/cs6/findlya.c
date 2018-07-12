# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <math.h>
# include <float.h>

# include "xtables.h"
# include "hstio.h"

# include "stis.h"
# include "stisdq.h"
# include "calstis6.h"

# define    SAFETY_MARGIN 15  /* safety around the slit projected image */

static int FindPixel (double, double *, int);

/*
   Find the region in the image that is contaminated by geocoronal Lya.
   This region must be avoided by the cross-correlation procedure. This
   routine computes the avoidance window limits, in physical (image)
   units.

   Originally this routine was put in place to support the cross correlation
   procedure. Later, it was reused by the background smoothing algorithm,
   which also has to avoid smoothing over the contaminated regions.

   For grating data, the Lya approximate position is found using the CRVAL1,
   CRPIX1 and CD1_1 values read from the SCI extension header. The slit
   width in the dispersion direction is computed (approximately) from the
   slit width in arcsec stored in the slit structure, corrected by the
   plate scale in the cross-dispersion direction and the ratio of the two
   LTM values (the plate scale in the dispersion direction is in A/pixel).

   For prism data, the Lya approximate position is found by actually
   searching in a wavelength array provided by the caller.

   The Lya position is then refined by cross-correlating a rectangular
   function with same width as the slit, with the image data collapsed
   in cross-dispersion direction. The collapsing uses the DQ flag value
   applicable for the cross-correlation procedure. The brightest peak
   in the cross-correlation result is accepted as the final Lya position,
   except if it differs from the original position by more than 50 pixels.
   In that case the second brightest peak is used, and the procedure is
   repeated iteratively. If the procedure fails down to the 5th brightest
   peak, the original Lya position is used.

   The Lya avoidance window is set to be 30 pixels wider than the projected
   slit width.

   On Aug 2002 we added the second geocoronal line at 1300 A.




   Revision history:
   ----------------
   17 Dec 98  -  Implemented (I.Busko)
   08 Aug 02  -  Add line at 1300 A (IB)
   07 Nov 02  -  Fix problem with Ly alpha detection in prism data (IB)

*/

void FindLya (StisInfo6 *sts, SingleGroup *in, ApInfo *slit,
              double *wave, int size,
              int *avoid1a, int *avoid2a,
              int *avoid1b, int *avoid2b) {

/* arguments:
StisInfo6 *sts         i: calibration switches and info
SingleGroup *in	       i: input image
ApInfo *slit           i: slit width taken from here
double *wave           i: array with wavelengths (only prism data use this)
int size               i: size of wavelength array
int avoid1a, avoid2a;  o: avoidance region limits (in physical pixels)
int avoid1b, avoid2b;
*/

	float *image, *crosscor, norm, pval;
	int centera, centerb, width, pindex, offset;
	int i, j, image_index, box_index;

	/* Initialize to default. */
	*avoid1a = 0;
	*avoid2a = 0;
	*avoid1b = 0;
	*avoid2b = 0;

	/* Compute approximate centers and projected slit width.

	   We could use here a single piece of code that would handle
           both the grating and prism cases. This would require changes
           in existing, tested code, though, which is a bad idea in
           principle since this is production software.

           Note that cd[1] for the prism is negative. Also, the computed
           slit width is approximate since the cross-terms in the CD matrix
           are not zero.
        */
	if (strcmp (sts->opt_elem, "PRISM") == 0) {

	    centera = FindPixel (1216., wave, size);
	    centerb = FindPixel (1300., wave, size);
	    width   = (int) (slit->width[0] / -(sts->cd[1]) / 3600.) *
                            (wave[centera+1] - wave[centera]);

	} else {
	    centera = (int)((1216. - sts->crval[0]) / sts->cd[0] +
                      sts->crpix[0]);
	    centerb = (int)((1300. - sts->crval[0]) / sts->cd[0] +
                      sts->crpix[0]);
	    width   = (int)(slit->width[0] / sts->cd[1] / 3600. /
                      sts->ltm[1] * sts->ltm[0]);
	}

	/* If lines not present, return immediately. */
	if ((centerb + width/2) < 0 ||
	    (centera - width/2) > in->sci.data.nx)
	    return;

	/* Build collapsed version of image data. */
	image = (float *) calloc (in->sci.data.nx, sizeof (float));
	for (j = 0; j < in->sci.data.ny; j++) {
	    for (i = 0; i < in->sci.data.nx; i++) {
                if (!(DQPix (in->dq.data, i, j) & sts->cc_sdqflags))
	            image[i] += Pix (in->sci.data, i, j);
	    }
	}

	/* Cross-correlate. */
	crosscor = (float *) calloc (in->sci.data.nx, sizeof (float));
	for (i = 0; i < in->sci.data.nx; i++) {
	    norm = 0.0;
	    for (box_index = 0; box_index < width; box_index++) {
	        image_index = box_index - width/2 + i;
                if (image_index > -1 && image_index < in->sci.data.nx) {
	            crosscor[i] += image[image_index];
	            norm        += 1.0;
	        }
	    }
	    if (norm > 0.0)
	        crosscor[i] /= norm;
	}

	free (image);

	/* Find Lya peak. */
	for (j = 1; j <= 5; j++) {
	    pval = - FLT_MAX;
	    for (i = 0; i < in->sci.data.nx; i++) {
	        if (crosscor[i] > pval) {
	            pval = crosscor[i];
	            pindex = i;
	        }
	    }
	    if (abs (pindex - centera) > 50) {
	        crosscor[pindex] = - FLT_MAX;
	        continue;
	    } else {

	        /* Found acceptable peak. */
	        *avoid1a = pindex - width/2 - SAFETY_MARGIN;
	        *avoid2a = pindex + width/2 + SAFETY_MARGIN;

	        free (crosscor);

	        /* Avoidance region of 1300 line is found by simple offset. */

	        offset = centerb - centera;
	        *avoid1b = *avoid1a + offset;
	        *avoid2b = *avoid2a + offset;

	        return;
	    }
	}

	/* Didn't find acceptable peak; use approximate center then. */

	*avoid1a = centera - width/2 - SAFETY_MARGIN;
	*avoid2a = centera + width/2 + SAFETY_MARGIN;

	/* Avoidance region of 1300 line is found by simple offset. */

	offset = centerb - centera;
	*avoid1b = *avoid1a + offset;
	*avoid2b = *avoid2a + offset;

	free (crosscor);
}


/* Given a wavelength, find the corresponding pixel number by binary
   search. If the supplied wavelength is out of bounds, return -1.
   It is assumed here that the wavelengths increase with pixel number.

   This is a modifed version of a function in calstis7 by Phil Hodge.
*/

static int FindPixel (double wl, double *wave, int size) {

	int x_low, x_high;	/* pixel numbers at ends of test range */
	int x_test;		/* pixel number at middle of test range */
	double wl_test;		/* wavelength at x_test */
	double wl_last;		/* wavelength at last iteration */

	/* reference pixel range */
	x_low = 0;
	x_high = 1023;

	if (wl < wave[x_low] || wl > MAX_PRISM_WAVELENGTH)
	    return (-1);

	wl_test = -1.0;
	wl_last = 0.0;
	while (wl_last != wl_test) {

	    wl_last = wl_test;

	    x_test = (x_low + x_high) / 2;
	    wl_test = wave[x_test];

	    if (wl < wl_test)
		x_high = x_test;
	    else
		x_low = x_test;
	}

	return (x_low - 2);
}
