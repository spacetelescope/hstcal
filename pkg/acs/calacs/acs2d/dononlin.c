/* This file contains:
	doNonLin
	FindRate
	ExpandDQ
*/

# include <stdio.h>
# include <math.h>		/* for sqrt and exp, in FindRate */

#include "hstcal.h"
# include "hstio.h"
# include "acs.h"
# include "acsinfo.h"
# include "hstcalerr.h"
# include "acsdq.h"		/* for SATPIXEL */

static double FindRate (double, double);
static void ExpandDQ (SingleGroup *, float, int, int);

/* This routine corrects or flags global nonlinearity, and/or it flags
   excessive local nonlinearity.
   If lflgcorr is PERFORM, we will check for and flag excessive global and
   local nonlinearity.  If glincorr is PERFORM, we will check for excessive
   global nonlinearity (by comparing the value of keyword GLOBRATE and the
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
   
   Warren Hack, 1998 June 11:
   	Initial ACS version.  No code changes necessary, other than 
		renaming structure to 'ACSInfo'.
*/

int doNonLin (ACSInfo *acs2d, SingleGroup *x, int *gsat, int *lsat) {

/* arguments:
ACSInfo *acs2d    i: calibration switches, etc
SingleGroup *x   io: image to be calibrated; written to in-place
int *gsat         o: > 0 if global saturation limit exceeded
int *lsat         o: > 0 if locally saturated pixels found
*/

	extern int status;

	double local_limit;	/* max local rate * exposure time */
	double ratio;		/* ratio of true to observed count rate */
	int i, j;

	int multk2d (SingleGroup *, float);
	int PutKeyStr (Hdr *, char *, char *, char *);

	/* Assign initial values. */
	*gsat = 0;
	*lsat = 0;
	PutKeyStr (&x->sci.hdr, "GLOBLIM", "NOT-EXCEEDED",
			    "global count rate exceeded?");

	if (acs2d->glincorr == PERFORM || acs2d->lflgcorr == PERFORM) {

	    if (acs2d->globrate > acs2d->global_limit) {

		*gsat = 1;			/* yes, they are saturated */

		/* Global nonlinearity is excessive; set GLOBLIM. */

		if (PutKeyStr (&x->sci.hdr, "GLOBLIM", "EXCEEDED",
				    "global count rate exceeded"))
		    return (status);

	    } else if (acs2d->glincorr == PERFORM) {

		/* Correct global nonlinearity. */

		ratio = FindRate (acs2d->globrate, acs2d->tau);
		if (multk2d (x, ratio))
		    return (status);
	    }
	}

	/* Flag local nonlinearity if it's excessive. */
	if (acs2d->lflgcorr == PERFORM) {

	    local_limit = acs2d->local_limit * acs2d->exptime;
        
	    /* If a value is nonlinear, flag it as well as pixels near it. */
	    for (j = 0;  j < x->sci.data.ny;  j++) {
		for (i = 0;  i < x->sci.data.nx;  i++) {
		    if (Pix(x->sci.data,i,j) > local_limit) {
			ExpandDQ (x, acs2d->expand, i, j);
			(*lsat)++;
		    }
		}
	    }
	}

	return (status);
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
   elliptical region around the specified pixel.  For every flagged
   pixel, this function will mark those pixels that are +/- EXPAND
   number of pixels from it.  Therefore, if EXPAND = 1, all pixels
   within 1 pixel from the flagged pixel will be marked.
   
   Removed all references to binning.  WJH 27 July 1999
*/

static void ExpandDQ (SingleGroup *x, float expand, int i0, int j0) {

/* arguments:
SingleGroup *x      io: image to be flagged; written to in-place
float expand        i: expand by this many pixels
int i0, j0          i: center pixel of region to flag
*/

	int ilow, ihigh, jlow, jhigh;	/* region around (i0,j0) */
	short dq;			            /* a data quality value */
	int iexpand;			        /* expand rounded up */
	float radius2;			        /* expand squared */
	int dx, dy;			            /* offsets from (i0,j0) */
	int i, j;

	if (i0 < 0 || i0 >= x->dq.data.nx || j0 < 0 || j0 >= x->dq.data.ny) {
	    sprintf (MsgText, "(ExpandDQ) (%d,%d) is out of range.", i0, j0);
	    trlwarn (MsgText);
		return;
	}

	iexpand = (int)(expand + 0.5);
	radius2 = expand * expand;

	/* Define the region within which we'll loop. */
    ilow  = i0 - iexpand;
    ihigh = i0 + iexpand;

    jlow  = j0 - iexpand;
    jhigh = j0 + iexpand;

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
		if ((float)(dx * dx + dy * dy) <= radius2) {
		    dq = SATPIXEL | DQPix(x->dq.data,i,j);
		    DQSetPix (x->dq.data, i, j, dq);
		}
	    }
	}
}
