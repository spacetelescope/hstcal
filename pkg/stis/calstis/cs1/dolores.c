# include <stdio.h>

# include "c_iraf.h"
# include "hstio.h"
# include "stis.h"
# include "calstis1.h"
# include "hstcalerr.h"
# include "stisdef.h"

/* This routine bins the input image so that BINX = 2 and BINY = 2.
   The values are added, not averaged.
   On input, it may be that one of BINX or BINY is already 2; that's OK.
   If both BINX and BINY are already 2, however, this routine sets done to
   false and does NOT reallocate x.
   The coordinate parameters, including LTVi and LTMi_i, will be updated
   in the extension headers.
   NOTE:  The input x will be freed and reallocated if we actually do
   change the binning.

   Phil Hodge, 1998 Oct 5:
	Change status value 1001 to GENERIC_ERROR_CODE.

   Phil Hodge, 1999 Nov 2:
	Second argument is now SingleGroup **x instead of SingleGroup *x.
	Remove SingleGroup *y from calling sequence.
*/

int doLoRes (StisInfo1 *sts, SingleGroup **x, int *done) {

/* arguments:
StisInfo1 *sts    i: calibration switches, etc
SingleGroup **x  io: input image to be binned; output binned image
int *done         o: true if input image has been binned
*/

	SingleGroup *in;	/* local variable, for cleaner notation */
	SingleGroup *out;	/* binned copy of x */
	int status;

	/* Factors by which we compress X and Y.  These are passed to bin2d,
	   and they are also used to keep track of which axes were binned.
	*/
	int mx, my;

	/* This is a flag to bin2d meaning sum the values, don't average. */
	int avg = 0;

	*done = 0;				/* initial value */

	in = *x;		/* just for convenience of notation */

	if (sts->bin[0] == 2 && sts->bin[1] == 2)
	    return (0);

	if (sts->bin[0] < 1 || sts->bin[1] < 1) {
	    printf ("ERROR    (doLoRes) bin size is less than one.\n");
	    return (GENERIC_ERROR_CODE);
	}
	if (sts->bin[0] > 2 || sts->bin[1] > 2) {
	    printf ("ERROR    (doLoRes) bin size is greater than two.\n");
	    return (GENERIC_ERROR_CODE);
	}

	mx = (sts->bin[0] == 1) ? 2 : 1;
	my = (sts->bin[1] == 1) ? 2 : 1;

	/* Create a new data set of the appropriate (smaller) size. */
	if ((out = calloc (1, sizeof (SingleGroup))) == NULL)
	    return (OUT_OF_MEMORY);
	initSingleGroup (out);
	allocSingleGroup (out, (in->sci.data.nx) / mx, (in->sci.data.ny) / my, True);
	if (hstio_err())
	    return (OUT_OF_MEMORY);

	/* Bin the data, putting the result in out. */
	if ((status = bin2d (in, 0, 0, mx, my, avg, out)))
	    return (status);

	/* Free x, and copy out to x. */
	freeSingleGroup (*x);
	free (*x);
	*x = out;

	/* Update values related to binning in the sts struct. */
	sts->bin[0] = 2;
	sts->bin[1] = 2;

	*done = 1;			/* yes, we binned the image */
	return (0);
}
