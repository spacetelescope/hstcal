# include <stdio.h>
# include <stdlib.h>	/* calloc */

# include "hstio.h"
# include "stis.h"
# include "calstis4.h"
# include "err.h"

# define ARCSEC_PER_DEGREE  3600.

/* This routine determines the shift of the aperture in the dispersion
   direction.
   The data are summed along a column
   to make a 1-D array, and the shift within that array is found.
   When summing a column or row, not all the data are included.  A
   certain fraction is taken, centered within the column or row.
   For each pixel in the summed 1-D array, if ALL pixels in a column
   that would have been included in the sum are flagged as bad, the
   flag for the 1-D data quality array will be set to sdqflags;
   otherwise, the flag will be set to zero.

   A pointer to the 1-D spectrum is returned, to be used as a weight
   by SpatialShift.  The median is subtracted from the spectrum, then
   negative values are set to zero.  The memory allocated for specweight
   should be freed by WaveCal.

   Phil Hodge, 1998 Dec 11:
	Use parameters from WCP table; add sts to calling sequence
	of XCWave.  Delete section for dispaxis=2.  Write median to
	the debug file.

   Phil Hodge, 2000 July 28:
	Include a check that nwl is greater than zero.

   Phil Hodge, 2000 Oct 18:
	Removal of DATAMASKED from sdqflags has been moved to GetGrpInfo7.

   Phil Hodge, 2001 Feb 23:
	Add disp to calling sequence; pass info to XCWave.
*/

int WaveShift (StisInfo4 *sts, ApInfo *slit, DispRelation *disp,
		LampInfo *lamp, SingleGroup *in,
		double **specweight, double *shift) {

/* arguments:
StisInfo4 *sts      i: calibration switches and info
ApInfo *slit        i: description of slit
DispRelation *disp  i: dispersion relation (currently only used for prism)
LampInfo *lamp      i: spectrum of calibration lamp
SingleGroup *in     i: input data
double **specweight o: summed 1-D spectrum, negative values clipped
double *shift       o: the shift, in pixels
*/

	int status;

	double crpix, crval, cdelt;	/* coord param for spatial axis */
	double slitwidth;	/* width of slit in pixels */
	double *v;		/* 1-D data */
	short *qv;		/* data quality for v */
	double median;		/* median of v */
	int nwl;		/* length of dispersion axis */
	int ngood;		/* number included in sum */
	int i, j;		/* loop indexes */
	int ifirst, ilast;	/* loop limits for summing in first axis */
	int jfirst, jlast;	/* loop limits in second axis */
	short flagval;		/* data quality flag for a pixel */

	int XCWave (StisInfo4 *, double *, double *, int,
		double *, short *, int, short,
		double,
		DispRelation *,
		double, double, double, double *);

	/* image section to use */
	ifirst = sts->wl_sect1[0];
	ilast  = sts->wl_sect1[1];
	jfirst = sts->wl_sect2[0];
	jlast  = sts->wl_sect2[1];

	/* Get info on coordinates.
	   We need the slit width (in the dispersion direction) in pixels,
	   but we have it in arcseconds.  cdelt in the cross dispersion
	   axis gives degrees per pixel, but we also have to consider the
	   possibility that the two axes are binned differently.
	*/
	crpix = sts->crpix[0];
	crval = sts->crval[0];
	cdelt = sts->cdelt[0];
	slitwidth = slit->width[0] /
	(sts->cdelt[1] * (sts->scale[0] / sts->scale[1]) * ARCSEC_PER_DEGREE);

	/* nwl is the full width of the image, not just ilast-ifirst+1,
	   for convenience in indexing.
	*/
	nwl = in->sci.data.nx;

	if (nwl <= 0) {
	    printf ("Warning  No data for shift in dispersion direction.\n");
	    *shift = UNDEFINED_SHIFT;
	    return (0);
	}
	v = calloc (nwl, sizeof(double));
	qv = calloc (nwl, sizeof(short));
	if (v == NULL || qv == NULL)
	    return (OUT_OF_MEMORY);

	/* Copy pointer to output for use by SpatialShift. */
	*specweight = v;

	/* Average the data values in the direction perpendicular to
	   the dispersion to make a 1-D array.  The data quality flag
	   for the sum will be zero (good) unless all pixels included
	   in the sum are flagged as bad.
	*/

	for (i = ifirst;  i <= ilast;  i++) {	/* loop over columns */

	    ngood = 0;
	    for (j = jfirst;  j <= jlast;  j++) {	/* sum current column */
		flagval = DQPix (in->dq.data, i, j);
		if ( ! (flagval & sts->sdqflags) ) {
		    v[i] += Pix (in->sci.data, i, j);
		    ngood++;
		}
	    }
	    if (ngood > 0)
		v[i] /= (double) ngood;
	    else
		qv[i] = sts->sdqflags;	/* this column is all bad */
	}
	for (i = 0;  i < ifirst;  i++)
	    qv[i] = sts->sdqflags;
	for (i = ilast + 1;  i < nwl;  i++)
	    qv[i] = sts->sdqflags;

	/* Find the median of the spectrum, subtract the median from the
	   spectrum, then replace negative values with zero.
	*/
	median = MedianDouble (v, nwl, 0);	/* 0 ==> DON'T sort in-place */
	for (i = 0;  i < nwl;  i++)
	    v[i] = (v[i] < median) ? 0. : (v[i] - median);

	/* Write info to the debug file. */
	if (sts->dbg != NULL) {
	    fprintf (sts->dbg,
	"# (WaveShift) %.6g has been subtracted from the observed spectrum,\n",
		median);
	    fprintf (sts->dbg, "# and values below zero truncated to zero.\n");
	}

	/* Now we have the data in v, its flags in qv, and the 1-D
	   coordinate parameters crpix, crval, and cdelt.  Find the shift.
	*/
	if ((status = XCWave (sts, lamp->wl, lamp->flux, lamp->nelem,
                              v, qv, nwl, sts->sdqflags,
                              slitwidth,
                              disp, crpix, crval, cdelt, shift)))
	    return (status);

	free (qv);

	return (0);
}
