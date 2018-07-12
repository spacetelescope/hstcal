/* This file contains:
	doNonLin
	FindGlobRate
	FindRate
	ExpandDQ
*/

# include <stdio.h>
# include <math.h>		/* for sqrt and exp, in FindRate */

# include "c_iraf.h"
# include "hstio.h"
# include "stis.h"
# include "calstis1.h"
# include "hstcalerr.h"
# include "stisdq.h"		/* for SATPIXEL */
# include "stisdef.h"

static double FindGlobRate (SingleGroup *, double);
static double FindRate (double, double);
static void ExpandDQ (SingleGroup *, float, int *, int, int);

/* This routine corrects or flags global nonlinearity, and/or it flags
   excessive local nonlinearity.

   If EXPTIME is greater than zero, the global count rate will be computed
   by summing all the counts and dividing by the exposure time, and the
   GLOBRATE keyword will then be updated with this value.  Note that this
   means that doNonLin must be called before the data values have been
   modified (e.g. by dark or flat field correction).

   If lflgcorr is PERFORM, we will check for and flag excessive global and
   local nonlinearity.  If glincorr is PERFORM, we will check for excessive
   global nonlinearity (by comparing the value of GLOBRATE and the
   GLOBAL_LIMIT from the table) and correct the nonlinearity if it's not
   too large.

   If the global linearity limit is exceeded, the gsat argument will be
   set to one, and the keyword GLOBLIM in the SCI extension header will be
   set to "EXCEEDED".  If the global linearity limit is not exceeded, gsat
   will be zero, and GLOBLIM will be "NOT-EXCEEDED".

   The count rate limit for local linearity is adjusted for bin size
   (by scaling by the pixel area if high-res is used in X and/or Y),
   but it is not known at the time of writing whether this is the
   appropriate correction.

   Phil Hodge, 2000 Mar 3:
	Compute the global count rate instead of using the value of the
	GLOBRATE keyword, and then update that keyword in the output.
*/

int doNonLin (StisInfo1 *sts, SingleGroup *x, int *gsat, int *lsat) {

/* arguments:
StisInfo1 *sts    i: calibration switches, etc
SingleGroup *x   io: image to be calibrated; written to in-place
int *gsat         o: > 0 if global saturation limit exceeded
int *lsat         o: > 0 if locally saturated pixels found
*/

	int status;

	double local_limit;	/* max local rate * exposure time */
	double ratio;		/* ratio of true to observed count rate */
	int i, j;

	/* Assign initial values. */
	*gsat = 0;
	*lsat = 0;
	Put_KeyS (&x->sci.hdr, "GLOBLIM", "NOT-EXCEEDED",
			    "global count rate exceeded?");

	if (sts->glincorr == PERFORM || sts->lflgcorr == PERFORM) {

	    /* Find the global count rate, and update the header keyword. */
	    if (sts->exptime > 0) {
		sts->globrate = FindGlobRate (x, sts->exptime);
		if ((status = Put_KeyD (&x->sci.hdr, "GLOBRATE", sts->globrate,
                                        "global count rate")))
		    return (status);
	    }

	    if (sts->globrate > sts->global_limit) {

		*gsat = 1;			/* yes, they are saturated */

		/* Global nonlinearity is excessive; set GLOBLIM. */

		if ((status = Put_KeyS (&x->sci.hdr, "GLOBLIM", "EXCEEDED",
                                        "global count rate exceeded")))
		    return (status);

	    } else if (sts->glincorr == PERFORM) {

		/* Correct global nonlinearity. */

		ratio = FindRate (sts->globrate, sts->tau);
		if ((status = multk2d (x, ratio)))
		    return (status);
	    }
	}

	/* Flag local nonlinearity if it's excessive. */
	if (sts->lflgcorr == PERFORM) {

	    local_limit = sts->local_limit * sts->exptime;
	    /* Adjust for binning.  (Is this OK?) */
	    local_limit *= (float)(sts->bin[0] * sts->bin[1]) / 4.;

	    /* If a value is nonlinear, flag it as well as pixels near it. */
	    for (j = 0;  j < x->sci.data.ny;  j++) {
		for (i = 0;  i < x->sci.data.nx;  i++) {
		    if (Pix(x->sci.data,i,j) > local_limit) {
			ExpandDQ (x, sts->expand, sts->bin, i, j);
			(*lsat)++;
		    }
		}
	    }
	}

	return (0);
}

/* This finds the count rate over the entire detector (including pixels
   flagged as bad) and over the entire exposure.  The result is returned
   as the function value.
*/

static double FindGlobRate (SingleGroup *x, double exptime) {

/* arguments:
SingleGroup *x    i: image to be calibrated
double exptime    i: exposure time
*/

	double sum;		/* the sum of all pixel values in the image */
	int i, j;

	sum = 0.;
	for (j = 0;  j < x->sci.data.ny;  j++) {
	    for (i = 0;  i < x->sci.data.nx;  i++) {
		sum += Pix(x->sci.data, i, j);
	    }
	}

	return (sum / exptime);
}

/* This routine takes the observed count rate and computes a correction
   factor which is the true count rate divided by the observed count rate.
   The expression for the observed rate y as a function of the true
   rate x is assumed to be y = x * exp (-tau * x).
*/

static double FindRate (double y, double tau) {

/* arguments:
double y         i: the observed count rate (counts / sec)
double tau       i: time constant (seconds)
function value   o: the ratio of true to observed count rate
*/

	double x0, x1, x2, x3;	/* successive estimates of true rate */
	double dx, dy;

	if (y <= 0.)
	    return (1.);

	/*
	       dy/dx = (1 - tau*x) * exp (-tau*x)
	*/
	x0 = y;				/* y is first estimate for x */

	dy = y - (x0 * exp (-tau * x0));
	dx = dy / ((1. - tau * x0) * exp (-tau * x0));
	x1 = x0 + dx;

	dy = y - (x1 * exp (-tau * x1));
	dx = dy / ((1. - tau * x1) * exp (-tau * x1));
	x2 = x1 + dx;

	dy = y - (x2 * exp (-tau * x2));
	dx = dy / ((1. - tau * x2) * exp (-tau * x2));
	x3 = x2 + dx;

	return (x3 / y);
}

/* This routine assigns (ORs) a data quality flag to a circular or
   elliptical region around the specified pixel.
*/

static void ExpandDQ (SingleGroup *x, float expand, int *bin, int i0, int j0) {

/* arguments:
SingleGroup *x      io: image to be flagged; written to in-place
float expand        i: expand by this many high-res pixels
int bin[2]          i: bin size (one or two) for each axis
int i0, j0          i: center pixel of region to flag
*/

	int ilow, ihigh, jlow, jhigh;	/* region around (i0,j0) */
	short dq;			/* a data quality value */
	int binx2, biny2;		/* squares of bin[0], bin[1] */
	int iexpand;			/* expand rounded up */
	float radius2;			/* expand squared */
	int dx, dy;			/* offsets from (i0,j0) */
	int i, j;

	if (i0 < 0 || i0 >= x->dq.data.nx || j0 < 0 || j0 >= x->dq.data.ny) {
	    printf ("Warning  (ExpandDQ) (%d,%d) is out of range.\n", i0, j0);
	    return;
	}

	binx2 = bin[0] * bin[0];
	biny2 = bin[1] * bin[1];
	iexpand = (int)(expand + 1);
	radius2 = expand * expand;

	/* Define the region within which we'll loop. */
	if (bin[0] == 1) {
	    ilow  = i0 - iexpand;
	    ihigh = i0 + iexpand;
	} else {
	    ilow  = i0 - (iexpand+1) / bin[0];
	    ihigh = i0 + (iexpand+1) / bin[0];
	}
	if (bin[1] == 1) {
	    jlow  = j0 - iexpand;
	    jhigh = j0 + iexpand;
	} else {
	    jlow  = j0 - (iexpand+1) / bin[1];
	    jhigh = j0 + (iexpand+1) / bin[1];
	}

	/* Truncate the region at the edges of the image. */
	if (ilow < 0)
	    ilow = 0;
	if (ihigh >= x->dq.data.nx)
	    ihigh = x->dq.data.nx - 1;
	if (jlow < 0)
	    jlow = 0;
	if (jhigh >= x->dq.data.ny)
	    jhigh = x->dq.data.ny - 1;

	/* Flag pixels within expand high-res pixels of (i0,j0). */
	for (j = jlow;  j <= jhigh;  j++) {
	    dy = j - j0;
	    for (i = ilow;  i <= ihigh;  i++) {
		dx = i - i0;
		if ((float)(dx * dx * binx2 + dy * dy * biny2) <= radius2) {
		    dq = SATPIXEL | DQPix(x->dq.data,i,j);
		    DQSetPix (x->dq.data, i, j, dq);
		}
	    }
	}
}
