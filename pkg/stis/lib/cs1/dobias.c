# include <stdio.h>

# include "c_iraf.h"
# include "hstio.h"
# include "stis.h"
# include "calstis1.h"
# include "hstcalerr.h"
# include "stisdef.h"

/* This routine subtracts the bias image from x (in-place).
   For science data (as opposed to wavecals), it will normally be the
   case that two or more images will have been combined for cosmic-ray
   rejection.  Before subtracting the bias image, it must be scaled by
   the number of images that were combined, since they will have been
   added together rather than averaged.
   Note that the bias image is assumed to have already been scaled by
   the gain.

   Phil Hodge, 1998 July 30:
	If we're processing a wavecal, it's no longer a fatal error for
	the bias image to be binned more than the wavecal.  In this case
	biascorr will be reset to SKIPPED, and a warning will be printed.
*/

int doBias (StisInfo1 *sts, SingleGroup *x) {

/* arguments:
StisInfo1 *sts     i: calibration switches, etc
SingleGroup *x    io: image to be calibrated; written to in-place
*/

	int status;

	SingleGroup y, z;	/* y and z are scratch space */
	int extver = 1;		/* get this imset from bias image */
	int rx, ry;		/* for binning bias image down to size of x */
	int x0, y0;		/* offsets of sci image */
	int same_size;		/* true if no binning of ref image required */
	int high_res;		/* returned by findBin and ignored */
	int avg = 0;		/* bin2d should sum within each bin */

	int FindBin (StisInfo1 *, SingleGroup *, SingleGroup *,
		int *, int *, int *, int *, int *, int *);

	if (sts->biascorr != PERFORM)
	    return (0);

	initSingleGroup (&y);

	/* Get the bias image data. */
	getSingleGroup (sts->bias.name, extver, &y);
	if (hstio_err())
	    return (OPEN_FAILED);

	/* Compare binning of science image and reference image;
	   get the same_size flag, and get info about binning and offset
	   for use by bin2d.
	*/
	if ((status = FindBin (sts, x, &y,
                &same_size, &high_res, &rx, &ry, &x0, &y0))) {
	    if (status == REF_TOO_SMALL && sts->wavecal) {
		printf (
	"Warning  BIASFILE is binned more than the input wavecal, \\\n");
		printf ("Warning  so BIASCORR will not be performed.\n");
		sts->biascorr = SKIPPED;	/* reset BIASCORR switch */
		freeSingleGroup (&y);
		return (0);
	    } else {
		return (status);
	    }
	}

	/* Subtract the bias image from x. */

	if (same_size) {

	    /* No binning required. */

	    /* scale bias by number of images that were combined */
	    if ((status = multk2d (&y, (float)(sts->ncombine))))
		return (status);
	    if ((status = sub2d (x, &y))) {
		printf ("ERROR    (biascorr) size mismatch.\n");
		return (status);
	    }
	    freeSingleGroup (&y);

	} else {

	    /* First bin the bias image down to the actual size of x. */

	    initSingleGroup (&z);
	    allocSingleGroup (&z, x->sci.data.nx, x->sci.data.ny, True);
	    if ((status = bin2d (&y, x0, y0, rx, ry, avg, &z))) {
		printf ("ERROR    (biascorr) size mismatch.\n");
		return (status);
	    }
	    freeSingleGroup (&y);			/* done with y */

	    /* scale the bias */
	    if ((status = multk2d (&z, (float)(sts->ncombine))))
		return (status);
	    if ((status = sub2d (x, &z)))
		return (status);
	    freeSingleGroup (&z);			/* done with z */

	}

	return (0);
}
