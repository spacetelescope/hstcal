# include <stdio.h>
# include <stdlib.h>		/* for strtod */
# include <string.h>
# include <math.h>		/* for fabs */

# include "hstio.h"
# include "stis.h"
# include "calstis1.h"
# include "stissizes.h"	/* for CCD_NPIX_X and CCD_NPIX_Y */
# include "stisdq.h"

# define SHORT_SLIT    5.
# define LONG_SLIT    45.

# define STIS_APERTURE_BORDER    10	/* unit = reference (detector) pixels */
# define ARCSEC_PER_DEGREE  3600.

/* This routine flags regions beyond the boundaries of the aperture with
   DETECTORPROB, in the case of observations with the CCD.

   If the input image was not taken with the CCD detector, or if TARGNAME
   is BIAS or DARK, this function returns without flagging anything.

   The aperture size is read from the APER_FOV keyword, which is expected
   to be a string of the form "28X50", giving the height by the width in
   arcseconds.  If either of these numbers cannot be read, a warning will
   be printed, and no region will be flagged.

   If the aperture is outside the image (e.g. if the image is a subarray
   above or below the aperture), a warning will be printed, and no region
   will be flagged.

   If the height of the aperture is less than SHORT_SLIT or longer than
   LONG_SLIT (defined above), no region will be flagged.

   Sections to the left and right of the aperture are only flagged for
   imaging observations.  This is because a spectrum would likely extend
   beyond the aperture to the left or right.

   Phil Hodge, 1998 June 9:
	Function created.

   Phil Hodge, 1998 July 23:
	Check the sts->bias_or_dark flag, and don't flag anything if set.
*/

void FlagFilter (StisInfo1 *sts, ShortTwoDArray *dq,
		int nx, int ny, double ri_m[], double ri_v[]) {

/* arguments:
StisInfo1 *sts         i: calibration switches, etc
ShortTwoDArray *dq    io: data quality array
int nx, ny             i: size of input image
double ri_m[], ri_v[]  i: linear transformation from reference to image coords
*/

	char *wx, *dummy;
	double xcenter, ycenter;	/* center of detector */
	double xwidth, ywidth;		/* width of aperture in X and Y */
	double cdelt1, cdelt2;	/* arcsec per ref pixel in X and Y */
	double xlowR, xhighR;	/* aperture limits in X in ref coords */
	double ylowR, yhighR;	/* aperture limits in Y in ref coords */
	double xlowI, xhighI;	/* aperture limits in X in image coords */
	double ylowI, yhighI;	/* aperture limits in Y in image coords */
	int xlow, xhigh;	/* aperture limits rounded to int */
	int ylow, yhigh;	/* aperture limits rounded to int */
	int det_xlow, det_xhigh;	/* detector limits in image coords */
	int det_ylow, det_yhigh;
	int i, j;		/* loop indexes */
	int peculiar = 0;	/* true if image is outside aperture */
	short dq_flag;

	/* This is the second line of several warning messages. */
	char noflag[] =
	{"Warning  regions beyond the aperture will not be flagged.\n"};

	if (sts->detector != CCD_DETECTOR)
	    return;

	if (sts->bias_or_dark)
	    return;

	/* Extract the aperture size from the aper_fov keyword value.
	   If the narrower of the two dimensions is too small or too large,
	   we don't need to flag anything.
	*/
	ywidth = strtod (sts->aper_fov, &wx);	/* arcseconds */
	if (ywidth == 0.) {
	    printf ("Warning  Can't interpret APER_FOV = `%s'; \\\n",
		sts->aper_fov);
	    printf ("%s", noflag);
	    return;
	}
	wx++;				/* skip over the "X" */
	xwidth = strtod (wx, &dummy);
	if (xwidth == 0.) {
	    printf ("Warning  Can't interpret APER_FOV = `%s'; \\\n",
		sts->aper_fov);
	    printf ("%s", noflag);
	    return;
	}
	if (ywidth < SHORT_SLIT || ywidth > LONG_SLIT)
	    return;

	/* This is a comment on the notation.

	   Because of overscan, binning, and subarray, we need to map the
	   limits of the illuminated portion of the detector to image pixel
	   coordinates, and map the limits of the aperture (filter) to image
	   pixels.  We'll flag the region outside the aperture but within the
	   illuminated portion of the detector.  The variables are used as
	   follows:

	   The lower left and upper right corners of the illuminated portion
	   of the detector in image pixel coordinates are:
		(det_xlow, det_ylow) and (det_xhigh, det_yhigh).

	   The lower left and upper right corners of the aperture in image
	   pixel coordinates are:
		(xlow, ylow) and (xhigh, yhigh).

	   These are in image pixel coordinates.  Because of overscan and
	   subarray and limited aperture size, any of the above eight variables
	   can be either inside or outside the image.  If they are outside,
	   we'll truncate them at the image borders.
	*/

	/* Find the edges of the detector in the image.  We need this info
	   so we can avoid flagging within the overscan regions.
	   Note that xlowI, etc, are just used as scratch here.
	*/
	xlowI  = ri_v[0];			/* first pixel, zero indexed */
	xhighI = (CCD_NPIX_X - 1.) * ri_m[0] + ri_v[0];	/* last pixel */
	ylowI  = ri_v[1];
	yhighI = (CCD_NPIX_Y - 1.) * ri_m[1] + ri_v[1];
	det_xlow  = NINT (xlowI);
	det_xhigh = NINT (xhighI);
	det_ylow  = NINT (ylowI);
	det_yhigh = NINT (yhighI);

	xcenter = CCD_NPIX_X / 2.;		/* center of detector */
	ycenter = CCD_NPIX_Y / 2.;

	/* Compute the lower and upper limits of the aperture and convert
	   from arcseconds to reference (detector) pixels.  Note that the
	   range is widened by some amount to allow for variations in MSM
	   and slit wheel positioning.
	*/
	cdelt1 = fabs (sts->cdelt[0] * ri_m[0]) * ARCSEC_PER_DEGREE;
	cdelt2 = fabs (sts->cdelt[1] * ri_m[1]) * ARCSEC_PER_DEGREE;
	xlowR  = xcenter - xwidth / cdelt1 / 2. - STIS_APERTURE_BORDER;
	xhighR = xcenter + xwidth / cdelt1 / 2. + STIS_APERTURE_BORDER;
	ylowR  = ycenter - ywidth / cdelt2 / 2. - STIS_APERTURE_BORDER;
	yhighR = ycenter + ywidth / cdelt2 / 2. + STIS_APERTURE_BORDER;

	/* Convert the above limits to image pixel coordinates. */
	xlowI  = xlowR  * ri_m[0] + ri_v[0];
	xhighI = xhighR * ri_m[0] + ri_v[0];
	ylowI  = ylowR  * ri_m[1] + ri_v[1];
	yhighI = yhighR * ri_m[1] + ri_v[1];

	xlow  = NINT (xlowI);
	xhigh = NINT (xhighI);
	ylow  = NINT (ylowI);
	yhigh = NINT (yhighI);

	if (strcmp (sts->obstype, "IMAGING") == 0) {
	    if (xlow > nx - 1)
		peculiar = 1;
	    if (xhigh < 0)
		peculiar = 1;
	}
	if (ylow > ny - 1)
	    peculiar = 1;
	if (yhigh < 0)
	    peculiar = 1;

	if (peculiar) {
	    printf (
	"Warning  The input image appears to be outside the aperture; \\\n");
	    printf ("%s", noflag);
	    return;
	}

	/* Truncate det_xlow, det_xhigh, det_ylow, det_yhigh at the
	   image borders.
	*/
	if (det_xlow < 0)
	    det_xlow = 0;
	if (det_xhigh > nx - 1)
	    det_xhigh = nx - 1;
	if (det_ylow < 0)
	    det_ylow = 0;
	if (det_yhigh > ny - 1)
	    det_yhigh = ny - 1;

	/* Truncate xlow, xhigh, ylow, yhigh at the image borders. */
	if (xlow < 0)
	    xlow = 0;
	if (xhigh > nx - 1)
	    xhigh = nx - 1;
	if (ylow < 0)
	    ylow = 0;
	if (yhigh > ny - 1)
	    yhigh = ny - 1;

	/* Assign data quality outside the region passed by the aperture. */

	/* Flag entire lines below the aperture. */
	for (j = det_ylow;  j < ylow;  j++) {
	    for (i = det_xlow;  i <= det_xhigh;  i++) {
		dq_flag = DETECTORPROB | PDQPix (dq, i, j);
		PDQSetPix (dq, i, j, dq_flag);
	    }
	}

	/* Flag sections to the left and right of the aperture, but only
	   for imaging observations.
	*/
	if (strcmp (sts->obstype, "IMAGING") == 0) {
	    for (j = ylow;  j <= yhigh;  j++) {
		for (i = det_xlow;  i < xlow;  i++) {	/* left columns */
		    dq_flag = DETECTORPROB | PDQPix (dq, i, j);
		    PDQSetPix (dq, i, j, dq_flag);
		}
		for (i = xhigh+1;  i <= det_xhigh;  i++) {	/* right */
		    dq_flag = DETECTORPROB | PDQPix (dq, i, j);
		    PDQSetPix (dq, i, j, dq_flag);
		}
	    }
	}

	/* Flag entire lines above the aperture. */
	for (j = yhigh+1;  j <= det_yhigh;  j++) {
	    for (i = det_xlow;  i <= det_xhigh;  i++) {
		dq_flag = DETECTORPROB | PDQPix (dq, i, j);
		PDQSetPix (dq, i, j, dq_flag);
	    }
	}
}
