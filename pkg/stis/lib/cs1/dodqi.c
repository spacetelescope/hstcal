/* This file contains:
	doDQI
	OpenBpixTab
	ReadBpixTab
	CloseBpixTab
	DQINormal
	DQIHigh
	FirstLast
*/

# include <stdio.h>
# include <string.h>

# include "c_iraf.h"
# include "hstio.h"
# include "xtables.h"
# include "stis.h"
# include "calstis1.h"
# include "hstcalerr.h"
# include "stisdq.h"
# include "stisdef.h"

typedef struct {
	IRAFPointer tp;			/* pointer to table descriptor */
	IRAFPointer cp_opt_elem;	/* grating name */
	IRAFPointer cp_xstart, cp_ystart; /* starting pixel for bad values */
	IRAFPointer cp_length;		/* this many bad values */
	IRAFPointer cp_axis;		/* repeat along X or Y axis (1 or 2) */
	IRAFPointer cp_flag;		/* this is the data quality value */
	int axlen1, axlen2;		/* DQ array size specified in header */
	int nrows;			/* number of rows in table */
} TblInfo;

typedef struct {
	/* values read from a table row */
	char opt_elem[STIS_CBUF+1];
	int xstart, ystart;
	int length;
	int axis;
	short flag;
} TblRow;

static int OpenBpixTab (char *, TblInfo *);
static int ReadBpixTab (TblInfo *, int, TblRow *);
static int CloseBpixTab (TblInfo *);
static void FirstLast (double *, double *, int *, int *,
		int *, int *, int *, int *);
static void DQINormal (ShortTwoDArray *, double *, TblRow *);
static void DQIHigh (ShortTwoDArray *, double *, TblRow *, int, int);

/* This routine ORs the data quality array in the input SingleGroup x
   with the pixels flagged in the data quality initialization table.
   If Doppler convolution is specified (doppcorr == perform), then the
   values from the table will be convolved with the Doppler smearing
   function.  Note that the values in the input data quality array
   (from x) will NOT be convolved.

   The data quality initialization (BPIXTAB) table should contain the
   following integer header parameters:
		SIZAXIS1, SIZAXIS2:  full size of data quality array;
			these are expected to be 1062 x 1044 for CCD data
			and 1024 x 1024 for MAMA data
   and the following integer columns:
		PIX1, PIX2:  starting pixel (one indexed) of a range
			of pixels to be assigned an initial value
		LENGTH:  number of pixels to be assigned the value
		AXIS:  1 --> X axis, 2 --> Y axis
		VALUE:  value to be ORed with data quality array

   In addition, rows will be selected by comparing the value of:
		OPT_ELEM:  grating (or mirror) name (string)
   with the value from the primary header of the science file.

   Each row of the table specifies a set of LENGTH pixels to be flagged
   in either the X or Y direction with the value VALUE.  The value assigned
   for each pixel will be the OR of VALUE and any previous value at that
   pixel.  For MAMA data, the units for PIX1, PIX2 and LENGTH will be
   low-res pixels.

   There are three potentially different coordinate systems that are
   relevant here.  The coordinates in the bpixtab are in units of reference
   pixels.  The image for which we are setting data quality flags may have
   a different pixel size and/or a different origin, due to binning,
   subimage, and CCD overscan.  Finally, we may temporarily set the data
   quality flags in a full-size scratch array, and then later bin and/or
   subset down to the actual image; this is necessary if the image is
   binned differently from the reference pixel size or if Doppler shift
   correction was applied on-board.

   Because of the different coordinate systems, we need separate mappings
   between systems.  If we use a scratch array, we will need a mapping
   from reference coordinates to scratch, and then we'll use a mapping
   from the scratch coordinates to the image coordinates.  If we don't
   use a scratch array, we'll map directly from reference coordinates
   to image coordinates.  The linear mappings are described by a matrix
   (actually just the diagonal part) and a vector:

	ri_m, ri_v:  mapping from reference coords to image
	rs_m, rs_v:  mapping from reference coords to scratch
	si_m, si_v:  mapping from scratch array to image

   ri_m & ri_v are gotten from the image header keywords (and adjusted
   for the fact that we use zero-indexed coordinates):

	ri_m[0] = LTM1_1
	ri_m[1] = LTM2_2
	ri_v[0] = LTV1 + (LTM1_1 - 1)
	ri_v[1] = LTV2 + (LTM2_2 - 1)

   rs_m & rs_v are particularly simple.  If the scratch array is high-res
   MAMA, the mapping is as follows (if not, the mapping is the identity):

	rs_m[0] = 2.
	rs_m[1] = 2.
	rs_v[0] = 0.5  (rs_v are +0.5 because coords are zero indexed)
	rs_v[1] = 0.5

   If Xr, Xs, Xi are the pixel coordinate of the same point in reference,
   scratch, and image coordinates respectively, we get the si mapping
   as follows:

	Xi = Xr * ri_m + ri_v
	Xs = Xr * rs_m + rs_v,  so Xr = (Xs - rs_v) / rs_m
   substituting for Xr,
	Xi = [(Xs - rs_v) / rs_m] * ri_m + ri_v
	   = Xs * (ri_m / rs_m) + (ri_v - rs_v * ri_m / rs_m)
   therefore
	si_m = (ri_m / rs_m)
	si_v = ri_v - rs_v * (ri_m / rs_m) = ri_v - rs_v * si_m

   Phil Hodge, 1998 Mar 13:
	Change the sections for Doppler convolution.

   Phil Hodge, 1998 Apr 3:
	In DQIHigh, if the data quality flag is DETECTORPROB, don't shift
	away from the left or right edge due to Doppler shift.

   Phil Hodge, 1998 June 8:
	Call FlagFilter.

   Phil Hodge, 1998 Sept 29:
	Move the call to FlagFilter up to just following the section for
	flagging saturated pixels.  In OpenBpixTab, delete the section
	that closes the table and returns if the table contains no rows.

   Phil Hodge, 1998 Oct 16:
	Also read the OPT_ELEM column, and select rows based on that value.
	Consistency change in DQINormal:  "ylow" to "ystart" in the following:
	    yhigh = ylow + tabrow->length - 1;

   Phil Hodge, 2007 June 28:
	Rewrite FirstLast().  The previous version was assigning incorrect
	values for first, last, and sfirst.
*/

int doDQI (StisInfo1 *sts, SingleGroup *x) {

/* arguments:
StisInfo1 *sts    i: calibration switches, etc
SingleGroup *x    io: image to be calibrated; DQ array written to in-place
*/

	int status;

	TblInfo tabinfo;	/* pointer to table descriptor, etc */
	TblRow tabrow;		/* values read from a table row */

	ShortTwoDArray ydq;		/* scratch space */

	/* mappings from one coordinate system to another */
	double ri_m[2], ri_v[2];	/* reference to image */
	double rs_m[2], rs_v[2];	/* reference to scratch */
	double si_m[2], si_v[2];	/* scratch to image */

	/* for copying from scratch array (only copy overlap region): */
	int first[2], last[2];	/* corners of overlap region in image coords */
	int sfirst[2];		/* lower left corner of overlap in scratch */
	int rbin[2];		/* bin size of image relative to ref bin size */

	int snpix[2];		/* size of scratch array */
	int npix[2];		/* size of current image */

	float *ds;		/* Doppler smearing array */
	int nds, d0;		/* size of ds and index in ds of zero point */
	int k, kmin, kmax;	/* loop index; range of indexes in ds */
	int doppmin, doppmax;	/* Doppler offsets relative to d0 */

	int in_place;		/* true if same bin size and no Doppler */
	int high_res;		/* true if Doppler or either axis is high-res */
	int i, j, i0, j0;	/* indexes for scratch array ydq */
	int m, n;		/* indexes for data quality array in x */
	short sum_dq;		/* for binning data quality array */

	int row;		/* loop index for row number */

	void FlagFilter (StisInfo1 *, ShortTwoDArray *,
		int, int, double *, double *);
	int MakeDopp (double, double, double, double, double, int,
		float *, int *, int *);

	/* We could still flag saturation even if the bpixtab was dummy. */
	if (sts->dqicorr != PERFORM && sts->dqicorr != DUMMY)
	    return (0);

	/* For the CCD, check for and flag saturation. */
	if (sts->detector == CCD_DETECTOR) {
	    for (j = 0;  j < x->sci.data.ny;  j++) {
		for (i = 0;  i < x->sci.data.nx;  i++) {
		    if ((int) Pix (x->sci.data, i, j) > sts->saturate) {
			sum_dq = DQPix (x->dq.data, i, j) | SATPIXEL;
			DQSetPix (x->dq.data, i, j, sum_dq);	/* saturated */
		    }
		}
	    }
	}

	/* Get the linear transformation between reference and input image. */
	if ((status = GetLT0 (&x->sci.hdr, ri_m, ri_v))) /* zero indexed LTV */
	    return (status);

	/* Flag regions beyond the bounderies of the aperture, for CCD data. */
	if (sts->detector == CCD_DETECTOR) {
	    FlagFilter (sts, &x->dq.data, x->dq.data.nx, x->dq.data.ny,
			ri_m, ri_v);
	}

	/* There might not be any bad pixel table.  If not, quit now. */
	if (sts->bpix.exists == EXISTS_NO || sts->dqicorr != PERFORM)
	    return (0);

	initShortData (&ydq);

	/* In some cases we can set the data quality flags directly in
	   the DQ array, but in other cases we must create a scratch
	   array and copy back to the original.  Either the original or
	   the scratch may be in high-res mode.
	*/
	if (sts->detector == CCD_DETECTOR) {

	    if (sts->bin[0] == 1 && sts->bin[1] == 1)
		in_place = 1;			/* no binning */
	    else
		in_place = 0;
	    high_res = 0;

	} else {				/* MAMA */

	    if (sts->doppcorr == PERFORM) {

		high_res = 1;

		if (sts->bin[0] == 1 && sts->bin[1] == 1)
		    in_place = 1;		/* high-res in both axes */
		else
		    in_place = 0;

	    } else {				/* no Doppler convolution */

		if (sts->bin[0] == 2 && sts->bin[1] == 2) {
		    high_res = 0;		/* both axes low-res */
		    in_place = 1;
		} else if (sts->bin[0] == 1 && sts->bin[1] == 1) {
		    high_res = 1;		/* both axes high-res */
		    in_place = 1;
		} else {
		    high_res = 1;		/* low-res in one axis */
		    in_place = 0;
		}
	    }
	}

	/* Get the other linear transformations (ri_m & ri_v were gotten
	   earlier, just after checking for saturation.)
	*/
	if (!in_place) {
	    if (high_res) {
		/* DQ array is binned finer than reference coords */
		rs_m[0] = 2.;
		rs_m[1] = 2.;
		rs_v[0] = 0.5;
		rs_v[1] = 0.5;
		/* assumes rs_m = 2, rs_v = 0.5 */
		si_m[0] = ri_m[0] * 0.5;
		si_m[1] = ri_m[1] * 0.5;
		si_v[0] = ri_v[0] - ri_m[0] * 0.25;
		si_v[1] = ri_v[1] - ri_m[1] * 0.25;
	    } else {
		/* scratch is in reference coords */
		rs_m[0] = 1.;
		rs_m[1] = 1.;
		rs_v[0] = 0.;
		rs_v[1] = 0.;
		/* assumes rs_m = 1, rs_v = 0 */
		si_m[0] = ri_m[0];
		si_m[1] = ri_m[1];
		si_v[0] = ri_v[0];
		si_v[1] = ri_v[1];
	    }
	}

	if (sts->doppcorr == PERFORM) {
	    /* Compute the Doppler smearing array, if we need it.  We need
		the size (nds) and zero point (d0), not the array itself.
	    */
	    nds = 2 * (sts->doppmag + 1) + 21;	/* reassigned by makeDopp */
	    ds = (float *) calloc (nds, sizeof (float));
	    if ((status = MakeDopp (sts->doppzero, sts->doppmag, sts->orbitper,
                                    sts->expstart, sts->exptime, sts->dispsign,
                                    ds, &nds, &d0)))
		return (status);
	    /* Find the range of non-zero elements in ds. */
	    kmin = nds - 1;		/* initial values */
	    kmax = 0;
	    for (k = 0;  k < nds;  k++) {
		if (ds[k] > 0.) {	/* there will be no negative values */
		    if (k < kmin)
			kmin = k;
		    if (k > kmax)
			kmax = k;
		}
	    }
	    /* It's the indexes relative to d0 that are important. */
	    doppmin = kmin - d0;
	    doppmax = kmax - d0;
	    free (ds);
	} else {
	    doppmin = 0;
	    doppmax = 0;
	}

	/* Open the data quality initialization table, find columns, etc. */
	if ((status = OpenBpixTab (sts->bpix.name, &tabinfo)))
	    return (status);

	/* Size of scratch image */
	if (high_res) {
	    snpix[0] = 2 * tabinfo.axlen1;
	    snpix[1] = 2 * tabinfo.axlen2;
	} else {
	    snpix[0] = tabinfo.axlen1;
	    snpix[1] = tabinfo.axlen2;
	}

	/* size of current image */
	npix[0] = x->dq.data.nx;
	npix[1] = x->dq.data.ny;

	if (!in_place) {
	    /* Allocate space for a scratch array. */
	    allocShortData (&ydq, snpix[0], snpix[1], True);
	    if (hstio_err()) {
		printf (
		"ERROR    (doDQI) couldn't allocate data quality array.\n");
		return (OUT_OF_MEMORY);
	    }
	    for (j = 0;  j < snpix[1];  j++)
		for (i = 0;  i < snpix[0];  i++)
		    DQSetPix (ydq, i, j, 0);		/* initially OK */
	}

	/* Read each row of the table, and fill in data quality values. */

	for (row = 1;  row <= tabinfo.nrows;  row++) {

	    if ((status = ReadBpixTab (&tabinfo, row, &tabrow))) {
		printf ("ERROR    Error reading BPIXTAB.\n");
		return (status);
	    }

	    if (!SameString (tabrow.opt_elem, sts->opt_elem))
		continue;

	    if (tabrow.xstart < 0 || tabrow.xstart >= tabinfo.axlen1 ||
		tabrow.ystart < 0 || tabrow.ystart >= tabinfo.axlen2) {
		printf (
	"Warning  Starting pixel (%d,%d) in BPIXTAB is out of range.\n",
			tabrow.xstart+1, tabrow.ystart+1);
		continue;			/* ignore this row */
	    }

	    /* Assign the flag value to all relevant pixels. */
	    if (in_place) {
		if (high_res)
		    DQIHigh (&x->dq.data, ri_v, &tabrow, doppmin, doppmax);
		else
		    DQINormal (&x->dq.data, ri_v, &tabrow);
	    } else {				/* use scratch array */
		if (high_res)
		    DQIHigh (&ydq, rs_v, &tabrow, doppmin, doppmax);
		else
		    DQINormal (&ydq, rs_v, &tabrow);
	    }
	}

	if ((status = CloseBpixTab (&tabinfo)))	/* done with the table */
	    return (status);

	if (!in_place) {

	    /* Get corners of region of overlap between image and
		scratch array.
	    */
	    FirstLast (si_m, si_v, snpix, npix, rbin, first, last, sfirst);

	    /* We have been writing to a scratch array ydq.  Now copy or
		bin the values down to the actual size of x.
	    */
	    j0 = sfirst[1];
	    for (n = first[1];  n <= last[1];  n++) {
		i0 = sfirst[0];
		for (m = first[0];  m <= last[0];  m++) {
		    sum_dq = DQPix (x->dq.data, m, n);
		    for (j = j0;  j < j0+rbin[1];  j++)
			for (i = i0;  i < i0+rbin[0];  i++)
			    sum_dq |= DQPix (ydq, i, j);
		    DQSetPix (x->dq.data, m, n, sum_dq);
		    i0 += rbin[0];
		}
		j0 += rbin[1];
	    }

	    freeShortData (&ydq);			/* done with ydq */
	}

	return (0);
}

static int OpenBpixTab (char *tname, TblInfo *tabinfo) {

	tabinfo->tp = c_tbtopn (tname, IRAF_READ_ONLY, 0);
	if (c_iraferr())
	    return (OPEN_FAILED);

	tabinfo->nrows = c_tbpsta (tabinfo->tp, TBL_NROWS);

	/* Find the columns. */
	c_tbcfnd1 (tabinfo->tp, "PIX1", &tabinfo->cp_xstart);
	c_tbcfnd1 (tabinfo->tp, "PIX2", &tabinfo->cp_ystart);
	c_tbcfnd1 (tabinfo->tp, "LENGTH", &tabinfo->cp_length);
	c_tbcfnd1 (tabinfo->tp, "AXIS", &tabinfo->cp_axis);
	c_tbcfnd1 (tabinfo->tp, "VALUE", &tabinfo->cp_flag);
	if (tabinfo->cp_xstart == 0 ||
	    tabinfo->cp_ystart == 0 ||
	    tabinfo->cp_length == 0 ||
	    tabinfo->cp_axis == 0 ||
	    tabinfo->cp_flag == 0) {
	    c_tbtclo (tabinfo->tp);
	    printf ("ERROR    Column not found in BPIXTAB.\n");
	    return (COLUMN_NOT_FOUND);
	}

	/* This column is optional, for backward compatibility. */
	c_tbcfnd1 (tabinfo->tp, "OPT_ELEM", &tabinfo->cp_opt_elem);

	/* Find out how large a full-size data quality array should be. */
	tabinfo->axlen1 = c_tbhgti (tabinfo->tp, "SIZAXIS1");
	if (c_iraferr()) {
	    c_tbtclo (tabinfo->tp);
	    printf (
		"ERROR    Couldn't get SIZAXIS1 from BPIXTAB header.\n");
	    return (TABLE_ERROR);
	}
	tabinfo->axlen2 = c_tbhgti (tabinfo->tp, "SIZAXIS2");
	if (c_iraferr()) {
	    c_tbtclo (tabinfo->tp);
	    printf ("ERROR    Couldn't get SIZAXIS2 from BPIXTAB header.\n");
	    return (TABLE_ERROR);
	}

	return (0);
}


/* This routine reads the relevant data from one row.  The starting
   pixel location, number of pixels and axis, and the flag value are read.
   The starting pixel (xstart, ystart) will be decremented by one to
   convert to zero-indexing.
*/

static int ReadBpixTab (TblInfo *tabinfo, int row, TblRow *tabrow) {

	/* OPT_ELEM is an optional column. */
	if (tabinfo->cp_opt_elem == 0) {
	    strcpy (tabrow->opt_elem, "ANY");
	} else {
	    c_tbegtt (tabinfo->tp, tabinfo->cp_opt_elem, row,
			tabrow->opt_elem, STIS_CBUF);
	    if (c_iraferr())
		return (TABLE_ERROR);
	}

	c_tbegti (tabinfo->tp, tabinfo->cp_xstart, row, &tabrow->xstart);
	if (c_iraferr())
	    return (TABLE_ERROR);
	tabrow->xstart--;
	c_tbegti (tabinfo->tp, tabinfo->cp_ystart, row, &tabrow->ystart);
	if (c_iraferr())
	    return (TABLE_ERROR);
	tabrow->ystart--;
	c_tbegti (tabinfo->tp, tabinfo->cp_length, row, &tabrow->length);
	if (c_iraferr())
	    return (TABLE_ERROR);
	c_tbegti (tabinfo->tp, tabinfo->cp_axis, row, &tabrow->axis);
	if (c_iraferr())
	    return (TABLE_ERROR);
	c_tbegts (tabinfo->tp, tabinfo->cp_flag, row, &tabrow->flag);
	if (c_iraferr())
	    return (TABLE_ERROR);

	if (tabrow->axis !=1 && tabrow->axis != 2) {
	    printf ("ERROR    Axis = %d in BPIXTAB, but it must be 1 or 2.\n",
			tabrow->axis);
	    c_tbtclo (tabinfo->tp);
	    return (TABLE_ERROR);
	}
	if (tabrow->length <= 0) {
	    printf (
	"ERROR    Length = %d in BPIXTAB, but it must be positive.\n",
		    tabrow->length);
	    c_tbtclo (tabinfo->tp);
	    return (TABLE_ERROR);
	}

	return (0);
}

/* This routine closes the bpixtab table. */

static int CloseBpixTab (TblInfo *tabinfo) {

	c_tbtclo (tabinfo->tp);
	if (c_iraferr())
	    return (TABLE_ERROR);

	return (0);
}


/* This routine assigns data quality values as specified in one row of
   the DQI table.  The current value is ORed with any previous value.
*/

static void DQINormal (ShortTwoDArray *ydq, double *ltv, TblRow *tabrow) {

/* arguments:
ShortTwoDArray *ydq  io: data quality array
double ltv[2]        i: offset of array within detector
TblRow *tabrow       i: data quality info read from one row
*/

	int xstart, ystart;	/* from tabrow, but scaled and shifted */
	int xlow, xhigh;	/* limits for loop on i */
	int ylow, yhigh;	/* limits for loop on j */
	int i, j;		/* indexes for scratch array ydq */
	short sum_dq;		/* for binning data quality array */

	/* Take account of possible subarray or overscan.  We use ltv,
	   but we assume ltm = 1.
	*/
	xstart = tabrow->xstart + ltv[0];
	ystart = tabrow->ystart + ltv[1];

	if (tabrow->axis == 1) {

	    /* The range for x is the repeat count. */
	    xlow = xstart;
	    xhigh = xstart + tabrow->length - 1;

	    /* out of range? */
	    if (xhigh < 0 || xlow >= ydq->nx ||
		ystart < 0 || ystart >= ydq->ny)
		return;

	    if (xlow < 0)
		xlow = 0;
	    if (xhigh >= ydq->nx)
		xhigh = ydq->nx - 1;

	    j = ystart;
	    for (i = xlow;  i <= xhigh;  i++) {
		sum_dq = tabrow->flag | PDQPix (ydq, i, j);
		PDQSetPix (ydq, i, j, sum_dq);
	    }

	} else if (tabrow->axis == 2) {

	    /* The range for y is the repeat count. */
	    ylow = ystart;
	    yhigh = ystart + tabrow->length - 1;

	    /* out of range? */
	    if (xstart < 0 || xstart >= ydq->nx ||
		yhigh < 0 || ylow >= ydq->ny)
		return;

	    if (ylow < 0)
		ylow = 0;
	    if (yhigh >= ydq->ny)
		yhigh = ydq->ny - 1;

	    i = xstart;
	    for (j = ylow;  j <= yhigh;  j++) {
		sum_dq = tabrow->flag | PDQPix (ydq, i, j);
		PDQSetPix (ydq, i, j, sum_dq);
	    }
	}
}

/* This routine assigns data quality values as specified in one row of
   the DQI table, for the case where the data quality array is binned
   in MAMA high-res pixels while the tabular values are low-res pixels.
   The range of pixels to be flagged may be extended due to the Doppler
   shift.  The current value is ORed with any previous value.
*/

static void DQIHigh (ShortTwoDArray *ydq, double *ltv,
		TblRow *tabrow, int doppmin, int doppmax) {

/* arguments:
ShortTwoDArray *ydq   io: data quality array
double ltv[2]         i: vector part of mapping from reference coords
TblRow *tabrow        i: data quality info read from one row
int doppmin, doppmax  i: offsets for Doppler shift
*/

	double temp;		/* xstart or ystart */
	int xstart, ystart;	/* from tabrow, but scaled and shifted */
	int xlength, ylength;	/* repeat count in high-res pixels */
	int xlow, xhigh;	/* limits for loop on i */
	int ylow, yhigh;	/* limits for loop on j */
	int nx, ny;		/* size of data quality array */
	int i, j;		/* indexes for scratch array ydq */
	short sum_dq;		/* for binning data quality array */

	/* xstart, ystart is the starting pixel in high-res coords,
	   obtained by mapping the lower left corner assuming ltm = 2;
	   i.e. (Xs - 0.5) = (Xr - 0.5) * 2 + ltv
	    -->  Xs = Xr * 2 + ltv - 0.5
	*/
	temp = (double)tabrow->xstart * 2. + ltv[0] - 0.5;
	xstart = NINT (temp);
	temp = (double)tabrow->ystart * 2. + ltv[1] - 0.5;
	ystart = NINT (temp);

	/* The repeat count is either two (one low-res pixel) or
	   twice the number read from the dqi table.
	*/
	if (tabrow->axis == 1) {
	    xlength = tabrow->length * 2;
	    ylength = 2;
	} else if (tabrow->axis == 2) {
	    xlength = 2;
	    ylength = tabrow->length * 2;
	}

	nx = ydq->nx;
	ny = ydq->ny;

	/* doesn't include Doppler yet */
	xlow  = xstart;
	xhigh = xstart + xlength - 1;

	/* Now include Doppler convolution. */
	if (tabrow->flag & DETECTORPROB) {
	    /* We're flagging the edge of the detector, so don't shift
		the flagged region away from the edge.
	    */
	    if (xlow != 0 && doppmin < 0)
		xlow += doppmin;
	    if (xhigh != nx-1 && doppmax > 0)
		xhigh += doppmax;
	} else {
	    xlow += doppmin;
	    xhigh += doppmax;
	}

	/* no Doppler */
	ylow  = ystart;
	yhigh = ystart + ylength - 1;

	/* Out of range? */
	if (xhigh < 0 || xlow >= nx || yhigh < 0 || ylow >= ny)
	    return;

	if (xlow < 0)
	    xlow = 0;
	if (xhigh >= nx)
	    xhigh = nx - 1;

	if (ylow < 0)
	    ylow = 0;
	if (yhigh >= ny)
	    yhigh = ny - 1;

	for (j = ylow;  j <= yhigh;  j++) {
	    for (i = xlow;  i <= xhigh;  i++) {
		sum_dq = tabrow->flag | PDQPix (ydq, i, j);
		PDQSetPix (ydq, i, j, sum_dq);
	    }
	}
}

# define TOLERANCE (0.001)

static void FirstLast (double *ltm, double *ltv, int *snpix, int *npix,
		int *rbin, int *first, int *last, int *sfirst) {

/* arguments:
double ltm[2], ltv[2]   i: linear transformation from scratch to image coords
				(adjusted for use with zero indexed coords)
int snpix[2]            i: size of scratch image
int npix[2]             i: size of actual image
int rbin[2]             o: number of pixels in scratch for one image pixel
int first[2], last[2]   o: corners of overlap region, in image coords
int sfirst[2]           o: lower left corner of overlap, in scratch array
*/

	double scr;		/* pixel coordinate in scratch array */
	int i, k;

	rbin[0] = NINT (1. / ltm[0]);
	rbin[1] = NINT (1. / ltm[1]);

	for (k = 0;  k < 2;  k++) {	/* axis number */
	    /* Search for the first pixel in the image array that maps to a
		point that is completely within the scratch array.  The left
		(lower) edge of pixel i is (i - 0.5).  Map (i - 0.5) to the
		scratch array, and if that point is within the scratch array,
		i is the first fully illuminated pixel.
	    */
	    for (i = 0;  i < npix[k];  i++) {
		scr = (i - 0.5 - ltv[k]) / ltm[k];
		/* -0.5 is left (lower) edge of first pixel in scratch */
		if (scr+TOLERANCE >= -0.5) {
		    first[k] = i;
		    sfirst[k] = NINT (scr+0.5);
		    break;
		}
	    }
	    /* Now look for the last fully illuminated pixel, using the
		right (upper) edge of the image pixels.
	    */
	    for (i = npix[k]-1;  i > 0;  i--) {
		scr = (i + 0.5 - ltv[k]) / ltm[k];
		/* compare scr with right (upper) edge of last pixel in the
		   scratch array
		*/
		if (scr-TOLERANCE <= snpix[k]-1. + 0.5) {
		    last[k] = i;
		    break;
		}
	    }
	}
}
