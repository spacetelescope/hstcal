# include <stdio.h>

# include "hstio.h"
# include "stis.h"
# include "calstis7.h"
# include "stisdq.h"

/* We'll flag beyond the slit ends if the slit is shorter than this
   fraction of the image height.
*/
# define SIGNIFICANT_FRACTION  (0.5)

/* This routine sets data quality flag DATAMASKED for pixels that are
   beyond the ends of the slit or are behind an occulting bar.
   If dispaxis is not one, this routine returns without doing anything;
   this is not an error.

   Phil Hodge, 1998 July 16:
	Expand the region that is not flagged for the slit length;
	expand by 0.5 pixel on each end, and also don't flag bottom
	and top pixel of slit.

   Phil Hodge, 2000 Apr 18:
	When flagging the bars, chop off bottom & top at the image edges.
	Return if dispaxis is not one (i.e. delete the section for dispaxis=2).

   Phil Hodge, 2000 Aug 4:
	The section for flagging beyond the aperture is now only done for
	slits that are significantly smaller than the image.  The ends of
	the long slits are typically outside the image anyway, and for
	certain apertures (e.g. 52X0.1F1) we have no way of knowing where
	the center of the slit is located within the image (and therefore
	where the endpoints are).  The code assumes that crpix2 is the
	center of the slit, which should still be valid for short slits.
*/

int DataMasked (StisInfo7 *sts, CoordInfo *coords, ApInfo *slit,
	SingleGroup *out) {

/* arguments:
StisInfo7 *sts       i: calibration switches and info
CoordInfo *coords    i: coordinate parameters
ApInfo *slit         i: description of slit (for setting DQ)
SingleGroup *out     io: output data
*/

	double center;		/* location of center of slit */
	double ends[2];		/* start and end locations of slit */
	double bar[2];		/* start and end of a bar */
	int bottom, top;	/* loop limits, dispaxis = 1 */
	int i, j;		/* loop indexes for pixel number */
	int k;			/* loop index for occulting bar */
	short dq;		/* data quality flag */

	/* Values from ApInfo struct, converted from arcsec to pixels. */
	double length;		/* size in axis perpendicular to dispersion */
	double barlocn[MAX_BARS], barwidth[MAX_BARS];
	double scale;		/* arcsec per pixel along the slit */

	if (sts->dispaxis != 1)
	    return (0);

	scale = coords->cdelt[1];		/* cdelt2, arcsec/pixel */
	length = slit->width[1] / scale;

	for (k = 0;  k < slit->nbars;  k++) {
	    barlocn[k] = slit->barlocn[k] / scale;
	    barwidth[k] = slit->barwidth[k] / scale;
	}

	center = coords->crpix[1];

	/* Flag regions beyond the slit ends, if the slit is significantly
	   shorter than the image.
	*/
	if (length < SIGNIFICANT_FRACTION * out->sci.data.ny) {

	    ends[0] = center - length / 2. - 0.5;
	    ends[1] = center + length / 2. + 0.5;
	    bottom = NINT(ends[0]);
	    top = NINT(ends[1]);

	    if (bottom > out->dq.data.ny)
		bottom = out->dq.data.ny;
	    if (top < -1)
		top = -1;

	    /* Flag the regions off the bottom and top of the slit. */
	    for (j = 0;  j < bottom;  j++) {
		for (i = 0;  i < out->dq.data.nx;  i++) {
		    dq = DQPix (out->dq.data, i, j) | DATAMASKED;
		    DQSetPix (out->dq.data, i, j, dq);
		}
	    }
	    for (j = top+1;  j < out->dq.data.ny;  j++) {
		for (i = 0;  i < out->dq.data.nx;  i++) {
		    dq = DQPix (out->dq.data, i, j) | DATAMASKED;
		    DQSetPix (out->dq.data, i, j, dq);
		}
	    }
	}

	/* Flag the bars. */
	for (k = 0;  k < slit->nbars;  k++) {

	    bar[0] = center + barlocn[k] - barwidth[k] / 2.;
	    bar[1] = center + barlocn[k] + barwidth[k] / 2.;
	    bottom = NINT(bar[0]);
	    top = NINT(bar[1]);

	    if (top >= 0 && bottom < out->dq.data.ny) {

		if (bottom < 0)
		    bottom = 0;
		if (top >= out->dq.data.ny)
		    top = out->dq.data.ny - 1;

		for (j = bottom;  j <= top;  j++) {
		    for (i = 0;  i < out->dq.data.nx;  i++) {
			dq = DQPix (out->dq.data, i, j) | DATAMASKED;
			DQSetPix (out->dq.data, i, j, dq);
		    }
		}
	    }
	}

	return (0);
}
