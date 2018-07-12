# include <stdio.h>
# include <stdlib.h>		/* for calloc */
# include <string.h>

# include "c_iraf.h"
# include "hstio.h"
# include "stis.h"
# include "calstis1.h"
# include "hstcalerr.h"
# include "stisdef.h"

/* This routine divides x in-place by the flat fields.
   There are up to three flat fields.  They are read into SingleGroups,
   multiplied together (leaving the result each time in y), and the
   product of the flat fields is divided into the input x, with the
   final quotient left in x.

   Phil Hodge, 1998 Mar 13:
	Change calling sequence of DoppConv.

   Phil Hodge, 1998 Aug 6:
	Add doppmag to calling sequence of DoppConv.

   Phil Hodge, 1998 Sept 24:
	Add missing return statements after calls to unbin2d.
*/

int doFlat (StisInfo1 *sts, SingleGroup *x) {

/* arguments:
StisInfo1 *sts     i: calibration switches, etc
SingleGroup *x    io: image to be calibrated; written to in-place
*/

	int status;

	SingleGroup w, y, z;	/* scratch space */
	float *ds;		/* Doppler smearing array */
	int nds;		/* size of ds */
	int d0;			/* index in ds of Doppler = 0 */
	int border;		/* = doppmag, unless obsmode = time-tag */
	int extver = 1;		/* get this group from flat field images */
	int rx, ry;		/* for binning dark down to size of x */
	int x0, y0;		/* offsets of sci image */
	int same_size;		/* true if no binning of ref image required */
	int high_res;		/* true if high-res pixels in dispersion dir */
	int avg = 1;		/* bin2d should average within each bin */
	int nx, ny;		/* how large to make y for lfltfile */
	int dummy;

	int FindBin (StisInfo1 *, SingleGroup *, SingleGroup *,
		int *, int *, int *, int *, int *, int *);
	int MakeDopp (double, double, double, double, double, int,
		float *, int *, int *);
	int DoppConv (SingleGroup *, int, float *, int, int);

	initSingleGroup (&w);
	initSingleGroup (&y);
	initSingleGroup (&z);

	/* pixel-to-pixel flat */
	if (sts->pfltcorr == PERFORM) {
	    getSingleGroup (sts->pflt.name, extver, &y);
	    if (hstio_err())
		return (OPEN_FAILED);
	}

	/* delta flat */
	if (sts->dfltcorr == PERFORM) {
	    if (sts->pfltcorr == PERFORM) {
		getSingleGroup (sts->dflt.name, extver, &z);
		if (y.sci.data.nx != z.sci.data.nx ||
		    y.sci.data.ny != z.sci.data.ny) {
		    printf (
	"ERROR    Pixel-to-pixel flat and delta flat are not the same size.\n");
		    return (SIZE_MISMATCH);
		}
		if ((status = mult2d (&y, &z)))		/* y is the product */
		    return (status);
		freeSingleGroup (&z);
	    } else {
		getSingleGroup (sts->dflt.name, extver, &y);
	    }
	    if (hstio_err())
		return (OPEN_FAILED);
	}

	/* low-order flat */
	if (sts->lfltcorr == PERFORM) {

	    /* Get lflt into a scratch area because lflt is smaller than y. */

	    if (sts->pfltcorr == PERFORM || sts->dfltcorr == PERFORM) {

		/* This is the normal case; we already have a product in y. */

		getSingleGroup (sts->lflt.name, extver, &w);
		if (hstio_err())
		    return (OPEN_FAILED);

		allocSingleGroup (&z, y.sci.data.nx, y.sci.data.ny, True);
		if (hstio_err())
		    return (ALLOCATION_PROBLEM);

		/* Resample w to z by linear interpolation. */
		if ((status = unbin2d (&w, &z))) /* unbin w --> z */
		    return (status);
		freeSingleGroup (&w);		/* we won't need w again */

		if ((status = mult2d (&y, &z))) /* y is the product */
		    return (status);
		freeSingleGroup (&z);

	    } else {

		/* We have neither a pixel-to-pixel flat nor a delta flat. */

		getSingleGroup (sts->lflt.name, extver, &z);

		/* figure out how much to expand the low-order flat. */
		FindBin (sts, x, &z, &dummy, &dummy, &rx, &ry, &dummy, &dummy);
		status = 0;		/* ignore status = REF_TOO_SMALL */

		/* Create y.  We don't need to assign any initial values
		   because y is strictly output from unbin2d.
		*/
		nx = rx * z.sci.data.nx;
		ny = ry * z.sci.data.ny;
		allocSingleGroup (&y, nx, ny, True);
		if (hstio_err())
		    return (ALLOCATION_PROBLEM);

		if ((status = unbin2d (&z, &y))) /* unbin z --> y */
		    return (status);
		freeSingleGroup (&z);
	    }
	}

	/* Now y contains the product of (up to) three flats. */

	/* Compare binning of science image and product of flat fields;
	   get same_size and high_res flags, and get info about binning
	   and offset for use by bin2d.
	*/
	if ((status = FindBin (sts, x, &y,
                               &same_size, &high_res, &rx, &ry, &x0, &y0)))
	    return (status);

	/* Do we need to do Doppler convolution? */
	if (sts->doppcorr == PERFORM) {

	    if (!high_res) {
		printf (
	"ERROR    Doppler convolution (DOPPCORR) was specified, \\\n");
		printf (
	"ERROR    but the flat fields are binned to low-res pixels.\n");
		return (SIZE_MISMATCH);
	    }

	    /* Allocate space for the Doppler smearing array, making it
		larger than we will need.  The actual size nds will be
		updated by MakeDopp.
	    */
	    nds = 2 * (sts->doppmag + 1) + 21;
	    ds = (float *) calloc (nds, sizeof (float));

	    if ((status = MakeDopp (sts->doppzero, sts->doppmag, sts->orbitper,
                                    sts->expstart, sts->exptime, sts->dispsign,
                                    ds, &nds, &d0)))
		return (status);

	    /* Convolve y with the Doppler smearing function. */
	    if (strcmp (sts->obsmode, "TIME-TAG") == 0)
		border = 0;
	    else
		border = NINT(sts->doppmag);
	    if ((status = DoppConv (&y, border, ds, nds, d0)))
		return (status);

	    free (ds);
	}

	/* Now we've got the complete flat field in y, convolved with
	   the Doppler smearing function if necessary.  Divide x by y.
	*/

	if (same_size) {

	    /* No binning required. */

	    if ((status = div2d (x, &y))) {
		printf ("ERROR    (doFlat) size mismatch\n");
		return (status);
	    }
	    freeSingleGroup (&y);

	} else {

	    /* Bin the flat field down to the actual size of x. */

	    allocSingleGroup (&z, x->sci.data.nx, x->sci.data.ny, True);
	    if ((status = bin2d (&y, x0, y0, rx, ry, avg, &z))) {
		printf ("ERROR    (doFlat) size mismatch\n");
		return (status);
	    }
	    freeSingleGroup (&y);		/* done with y */
	    if ((status = div2d (x, &z)))
		return (status);
	    freeSingleGroup (&z);		/* done with z */
	}

	return (0);
}
