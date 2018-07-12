# include <stdio.h>

# include "c_iraf.h"
# include "hstio.h"
# include "stis.h"
# include "calstis1.h"
# include "hstcalerr.h"
# include "stisdef.h"

/* This routine applies the shutter shading correction:

	corrected = uncalibrated * EXPTIME / (EXPTIME + SHADFILE)

   The value of EXPTIME will (locally) be divided by the number of
   exposures that were combined by cosmic-ray rejection to make the
   current image.
*/

int doShad (StisInfo1 *sts, SingleGroup *x) {

/* arguments:
StisInfo1 *sts     i: calibration switches, etc
SingleGroup *x    io: image to be calibrated; written to in-place
*/

	int status;

	SingleGroup y, z;	/* y and z are scratch space */
	float exptime;		/* exposure time / ncombine */
	int extver = 1;		/* get this imset from shutter shading image */
	int rx, ry;		/* for binning dark down to size of x */
	int x0, y0;		/* offsets of sci image */
	int same_size;		/* true if no binning of ref image required */
	int high_res;		/* true if high-res pixels in dispersion dir */
	int avg = 1;		/* bin2d should average within each bin */

	int FindBin (StisInfo1 *, SingleGroup *, SingleGroup *,
		int *, int *, int *, int *, int *, int *);

	if (sts->shadcorr != PERFORM)
	    return (0);

	if (sts->exptime <= 0.) {
	    printf (
		"Warning  EXPTIME must be positive if SHADCORR = PERFORM.\n");
	    sts->shadcorr = IGNORED;
	    return (0);
	}

	initSingleGroup (&y);

	/* Correct for the number of exposures that were combined. */
	exptime = sts->exptime / sts->ncombine;

	/* Get the shutter shading image data. */
	getSingleGroup (sts->shad.name, extver, &y);
	if (hstio_err())
	    return (OPEN_FAILED);

	/* Compare binning of science image and reference image;
	   get the same_size flag, and get info about binning and offset
	   for use by bin2d.
	*/
	if ((status = FindBin (sts, x, &y,
                               &same_size, &high_res, &rx, &ry, &x0, &y0)))
	    return (status);

	/* Apply the shutter shading image to the input image. */

	if (same_size) {

	    /* No binning required. */

	    if ((status = multk2d (&y, 1. / exptime)))
		return (status);
	    if ((status = addk2d (&y, 1.)))
		return (status);
	    if ((status = div2d (x, &y))) {
		printf ("ERROR    (doShad) size mismatch.\n");
		return (status);
	    }
	    freeSingleGroup (&y);

	} else {

	    /* Bin the reference image down to the actual size of x. */

	    initSingleGroup (&z);
	    allocSingleGroup (&z, x->sci.data.nx, x->sci.data.ny, True);
	    if ((status = bin2d (&y, x0, y0, rx, ry, avg, &z))) {
		printf ("ERROR    (doShad) size mismatch.\n");
		return (status);
	    }
	    freeSingleGroup (&y);			/* done with y */
	    if ((status = multk2d (&z, 1. / exptime)))
		return (status);
	    if ((status = addk2d (&z, 1.)))
		return (status);
	    if ((status = div2d (x, &z)))
		return (status);
	    freeSingleGroup (&z);			/* done with z */
	}

	return (0);
}
