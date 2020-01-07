# include <stdio.h>
# include <stdlib.h>	/* calloc */

# include "hstio.h"
# include "stis.h"
# include "calstis4.h"
# include "hstcalerr.h"

static int CollapsePrism (SpTrace *, SingleGroup *,
		int, int, int, int,
		double [], double [], short,
		double [], short[]);
static void LinearInterp (SingleGroup *, int, double, double *, short *);

/* This routine determines the shift of the aperture in the direction
   parallel to the slit.
   The data are averaged along a row
   to make a 1-D array, using specweight as weight, and the shift within
   that array is found.
   When summing a row or column, not all the data are included.  A
   certain fraction is taken, centered within the row or column.
   For each pixel in the summed 1-D array, if ALL pixels in a row
   that would have been included in the sum are flagged as bad, the
   flag for the 1-D data quality array will be set to sdqflags;
   otherwise, the flag will be set to zero.

   A status of NO_GOOD_DATA will be returned if there was no disastrous
   error, but the shift could not be determined.

   Phil Hodge, 1998 Dec 11:
	Use parameters from WCP table; add sts to calling sequences of
	several functions.  Delete section for dispaxis=2.

   Phil Hodge, 2000 Mar 28:
	Remove verbose from the calling sequence of FindBars.
	Don't print slit illumination to the debug file if the slit is
	a long slit (because this is done in FindBars).

   Phil Hodge, 2000 May 2:
	Print the DQ value to the debug file for every pixel, not just
	when DQ is non-zero.

   Phil Hodge, 2000 July 28:
	For long-slit aperture, include a check that there is at least one
	occulting bar.  Check that nv is greater than zero.

   Phil Hodge, 2000 Oct 18:
	Removal of DATAMASKED from sdqflags has been moved to GetGrpInfo7.

   Phil Hodge, 2001 Feb 23:
	Add trace to calling sequence; this is for prism data.
*/

int SpatialShift (StisInfo4 *sts, ApInfo *slit, SpTrace *trace,
		SingleGroup *in, double *specweight, double *shift) {

/* arguments:
StisInfo4 *sts      i: calibration switches and info
ApInfo *slit        i: description of slit
SpTrace *trace      i: list of spectral traces
SingleGroup *in     i: input data
double specweight[] i: the observed spectrum, used as weights
double *shift       o: the shift, in pixels
*/

	int status;

	double crpix, cdelt;	/* coord param for spatial axis */
	double length;		/* slit length */
	double sumw;		/* sum of weights */
	double *v;		/* 1-D data */
	short *qv;		/* data quality for v */
	int nv;			/* length of cross dispersion axis */
	int slittype;		/* code for type of slit */
	int i, j;		/* loop indexes */
	int ifirst, ilast;	/* loop limits for summing in first axis */
	int jfirst, jlast;	/* loop limits in second axis */
	short flagval;		/* data quality flag for a pixel */

	void WhichSlit (char *, int, int *);
	int XCEchelle (StisInfo4 *, double,
		double *, short *, int, short,
		double, double, double *);
	int FindEnds (StisInfo4 *, double,
		double *, short *, int,
		double, double, int, double *);
	int FindBars (StisInfo4 *, int, double *, double*,
		double *, short *, int,
		double, double, double *);

	/* image section to use */
	ifirst = sts->sp_sect1[0];
	ilast  = sts->sp_sect1[1];
	jfirst = sts->sp_sect2[0];
	jlast  = sts->sp_sect2[1];

	/* nv is the full height of the image, not just jlast-jfirst+1,
	   for convenience in indexing.
	*/
	nv = in->sci.data.ny;

	crpix = sts->crpix[1];
	cdelt = sts->cdelt[1];
	length = (double) slit->width[1];

	if (nv <= 0) {
	    printf ("Warning  No data for spatial shift.\n");
	    *shift = UNDEFINED_SHIFT;
	    return (0);
	}
	v = (double *) calloc (nv, sizeof(double));
	qv = (short *) calloc (nv, sizeof(short));
	if (v == NULL || qv == NULL)
	    return (OUT_OF_MEMORY);

	/* Average the data values in the dispersion direction to make
	   a 1-D array.  The data quality flag for the sum will be zero
	   (good) unless all pixels included in the sum are flagged as bad.
	*/
	if (sts->disp_type == PRISM_DISP) {

	    if ((status = CollapsePrism (trace, in,
			jfirst, jlast, ifirst, ilast,
                        sts->ltm, sts->ltv, sts->sdqflags,
                        v, qv)))
		return (status);

	} else {		/* already 2-D rectified data */

	    for (j = jfirst;  j <= jlast;  j++) {	/* loop over rows */

		sumw = 0.;
		for (i = ifirst;  i <= ilast;  i++) {	/* sum current row */
		    flagval = DQPix (in->dq.data, i, j);
		    if ( ! (flagval & sts->sdqflags) ) {
			v[j] += Pix (in->sci.data, i, j) * specweight[i];
			sumw += specweight[i];
		    }
		}
		if (sumw > 0.)
		    v[j] /= sumw;
		else
		    qv[j] = sts->sdqflags;	/* this row is all bad */
	    }
	}
	for (j = 0;  j < jfirst;  j++)
	    qv[j] = sts->sdqflags;
	for (j = jlast + 1;  j < nv;  j++)
	    qv[j] = sts->sdqflags;

	/* Now we have the data in v, its flags in qv, and the 1-D
	   coordinate parameters crpix and cdelt.  Find the shift.
	*/

	WhichSlit (sts->aperture, sts->dispaxis, &slittype);

	/* Write info to debug file. */
	if (sts->dbg != NULL) {

	    fprintf (sts->dbg, "\n");
	    fprintf (sts->dbg, "# (SpatialShift) Slit type is:  ");
	    if (slittype == SHORT_ECHELLE_SLIT)
		fprintf (sts->dbg, "short echelle slit\n");
	    else if (slittype == MEDIUM_ECHELLE_SLIT)
		fprintf (sts->dbg, "medium echelle slit\n");
	    else if (slittype == LONG_SLIT)
		fprintf (sts->dbg, "long slit\n");

	    if (slittype != LONG_SLIT) {
		fprintf (sts->dbg,
		"# (SpatialShift) pixel, slit illumination, DQ:\n");
		for (i = 0;  i < nv;  i++)
		    fprintf (sts->dbg, "%d %.6g %d\n", i + 1, v[i], qv[i]);
	    }
	}

	if (slittype == SHORT_ECHELLE_SLIT) {
	    status = XCEchelle (sts, length,
			v, qv, nv, sts->sdqflags, crpix, cdelt, shift);
	} else if (slittype == MEDIUM_ECHELLE_SLIT) {
	    status = FindEnds (sts, length,
			v, qv, nv, crpix, cdelt, sts->verbose, shift);
	} else if (slittype == LONG_SLIT) {
	    if (slit->nbars > 0) {
		status = FindBars (sts,
			slit->nbars, slit->barlocn, slit->barwidth,
			v, qv, nv, crpix, cdelt, shift);
	    } else {
		printf (
	"Warning  Aperture `%s' has no occulting bars, and without them \\\n",
			sts->aperture);
		printf (
	"         we can't find the shift in the spatial direction.\n");
		*shift = UNDEFINED_SHIFT;
		status = 0;		/* not a fatal error */
	    }
	} else {
	    printf ("Warning  Aperture `%s' is not supported for a wavecal.\n",
		sts->aperture);
	    status = NO_GOOD_DATA;
	}

	free (qv);
	free (v);

	return (status);
}

static int CollapsePrism (SpTrace *trace, SingleGroup *in,
		int jfirst, int jlast, int ifirst, int ilast,
		double ltm[], double ltv[], short sdqflags,
		double v[], short qv[]) {

	SpTrace *trace_y;	/* spectrum trace interpolated at y */
	double ydispl;		/* offset of spectrum in y direction */
	int i, j;		/* loop indexes */
	double x, y;		/* i & j in reference pixels */
	int i_r;		/* nearest int to x */
	double y_im;		/* j + trace, image pixel coords */
	int ngood;		/* number included in sum */
	double value;		/* value to add to v array */
	short flagval;
	int status;
	/* true if image is already in reference pixel coordinates */
	int in_ref_coords;

	int InterpTrace4 (SpTrace **, double, SpTrace **);
	void FreeTrace4 (SpTrace **);

	trace_y = NULL;

	in_ref_coords = (ltm[0] == 1. && ltm[1] == 1. &&
			 ltv[0] == 0. && ltv[1] == 0.);

	for (j = jfirst;  j <= jlast;  j++) {	/* loop over rows */

	    ngood = 0;

	    if (in_ref_coords)
		y = (double)j;
	    else
		y = ((double)j - ltv[1]) / ltm[1];

	    /* Interpolate to get the spectrum trace at j. */
	    if ((status = InterpTrace4 (&trace, y, &trace_y)))
		return (status);

	    for (i = ifirst;  i <= ilast;  i++) {	/* sum current row */

		if (in_ref_coords) {
		    ydispl = trace_y->a2displ[i];
		} else {
		    /* convert i to ref pixels, then ydispl to image pixels */
		    x = ((double)i - ltv[0]) / ltm[0];
		    i_r = NINT (x);
		    ydispl = trace_y->a2displ[i_r];
		    ydispl *= ltm[1];
		}

		y_im = (double)j + ydispl;

		/* interpolate at location of trace */
		LinearInterp (in, i, y_im, &value, &flagval);

		if ( ! (flagval & sdqflags) ) {
		    v[j] += value;
		    ngood++;
		}
	    }
	    if (ngood == 0)
		qv[j] = sdqflags;
	    else
		qv[j] = 0;
	}

	FreeTrace4 (&trace_y);

	return (0);
}

static void LinearInterp (SingleGroup *in, int i, double y_im,
		double *value, short *flagval) {

	double p, q;		/* weights */
	int j;

	j = (int)y_im;

	if (j < 0 || j >= in->sci.data.ny - 1) {
	    *value = 0.;
	    *flagval = 1;
	    return;
	}

	p = y_im - (double)j;
	q = 1. - p;

	/* interpolate the value */
	*value = q * Pix (in->sci.data, i, j) +
		 p * Pix (in->sci.data, i, j+1);

	/* take data quality flag from nearer pixel */
	if (p > 0.5)
	    *flagval = DQPix (in->dq.data, i, j);
	else
	    *flagval = DQPix (in->dq.data, i, j+1);
}
