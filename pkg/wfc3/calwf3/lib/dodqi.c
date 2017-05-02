/* This file contains:
	doDQI
	OpenBpixTab
	ReadBpixTab
	CloseBpixTab
	DQINormal
	FirstLast
*/

# include <stdio.h>
# include <string.h>

# include "xtables.h"
# include "hstio.h"
# include "wf3.h"
# include "wf3info.h"
# include "hstcalerr.h"
# include "wf3dq.h"

# define MIN(a,b) (a < b ? a : b)

typedef struct {
	IRAFPointer tp;			/* pointer to table descriptor */
	IRAFPointer cp_xstart, cp_ystart; /* starting pixel for bad values */
	IRAFPointer cp_length;		/* this many bad values */
	IRAFPointer cp_axis;		/* repeat along X or Y axis (1 or 2) */
	IRAFPointer cp_flag;		/* this is the data quality value */
	IRAFPointer cp_amp;		/* selection columns */
	IRAFPointer cp_ccdchip;
	IRAFPointer cp_ccdgain;
	int axlen1, axlen2;		/* DQ array size specified in header */
	int nrows;			/* number of rows in table */
} TblInfo;

typedef struct {
	/* values read from a table row */
	int xstart, ystart;
	int length;
	int axis;
	short flag;
	/* values used for selecting row to be used */
	char ccdamp[SZ_CBUF+1];
	float ccdgain;
	int ccdchip;
} TblRow;

static void FirstLast (double *, double *, int *, int *, int *, int *,
		       int *, int *);

/* This routine ORs the data quality array in the input SingleGroup x
   with the pixels flagged in the data quality initialization table.

   The data quality initialization (BPIXTAB) table should contain the
   following integer header parameters:
		SIZAXIS1, SIZAXIS2:  full size of data quality array;
			these are expected to be 4096 x 2048 for ACS WFC
			and WFC3 UVIS data
			and 1024 x 1024 for ACS MAMA and HRC
			and WFC3 IR data
   and the following integer columns:
		PIX1, PIX2:  starting pixel (one indexed) of a range
			of pixels to be assigned an initial value
		LENGTH:  number of pixels to be assigned the value
		AXIS:  1 --> X axis, 2 --> Y axis
		VALUE:  value to be ORed with data quality array
   Each row of the table specifies a set of LENGTH pixels to be flagged
   in either the X or Y direction with the value VALUE.  The value assigned
   for each pixel will be the OR of VALUE and any previous value at that
   pixel.  For ACS MAMA data, the units for PIX1, PIX2 and LENGTH will be
   low-res pixels.

   There are three potentially different coordinate systems that are
   relevant here.  The coordinates in the bpixtab are in units of reference
   pixels.  The image for which we are setting data quality flags may have
   a different pixel size and/or a different origin, due to binning,
   subimage, and CCD overscan.  The linear mappings are described by a matrix
   (actually just the diagonal part) and a vector:

	ri_m, ri_v:  mapping from reference coords to image

   ri_m & ri_v are gotten from the image header keywords (and adjusted
   for the fact that we use zero-indexed coordinates):

	ri_m[0] = LTM1_1
	ri_m[1] = LTM2_2
	ri_v[0] = LTV1 + (LTM1_1 - 1)
	ri_v[1] = LTV2 + (LTM2_2 - 1)
  
   Warren Hack, 1998 June 2:
    Revised for ACS data... Based on original CALSTIS code by Phil Hodge.

   Howard Bushouse, 2001 Dec 11:
    Revised for WFC3 data - Reinstated CALSTIS code required to handle
    binned science data (had been removed for CALACS).

   H.Bushouse, 2002 June 20:
    Updated to apply only those rows that have CCDAMP, CCDGAIN, and
    (most importantly) CCDCHIP values corresponding to the input science data
    (following CALACS changes).

   H.Bushouse, 2003 Oct 16:
    Updated to use floating-point gain value for WFC3.

   H.Bushouse, 2004 Feb 20:
    Modified to make the check against matching CCDAMP and CCDGAIN values
    in the BPIXTAB optional, so that they do not have to be included if
    they are deemed to be unnecessary. Also modified to apply new A-TO-D
    saturated DQ flag (following CALACS change).

   H.Bushouse, 2006 July 17:
    Added new ToWF3RawCoords routine to adjust BPIXTAB pixel coords for 
    presence of serial virtual overscan in WFC3 raw images. Also fixed
    loop limits in assigning scratch image flags to binned science image
    flags.

   H.Bushouse, 2007 Feb 20:
    Modified to allow for wildcard values in BPIXTAB Amp, Gain, and
    Chip columns (following CALACS change).

   H.Bushouse, 2007 Dec 21:
    Updated to use new rewrite of FirstLast routine provided by P. Hodge
    from calstis.
*/

int doDQI (WF3Info *wf3, SingleGroup *x) {

/* arguments:
WF3Info *wf3    i: calibration switches, etc
SingleGroup *x    io: image to be calibrated; DQ array written to in-place
*/

	extern int status;

	TblInfo tabinfo;	/* pointer to table descriptor, etc */
	TblRow tabrow;		/* values read from a table row */

	DQHdrData ydq;		/* scratch space */

	/* mappings from one coordinate system to another */
	double ri_m[2], ri_v[2];	/* reference to image */
	double rs_m[2], rs_v[2];	/* reference to scratch */
	double si_m[2], si_v[2];	/* scratch to image */

	/* for copying from scratch array (only copy overlap region) */
	int first[2], last[2];	/* corners of overlap region in image coords */
	int sfirst[2];		/* lower left corner of overlap in scratch */
	int rbin[2];		/* bin size of image relative to ref bin size */

	int snpix[2];		/* size of scratch array */
	int npix[2];		/* size of current image */

	int in_place;		/* true if same bin size */
	int i, j, i0, j0;	/* indexes for scratch array ydq */
	int m, n;		/* indexes for data quality array in x */
	short sum_dq;		/* for binning data quality array */
	float sat;

	int row;		/* loop index for row number */
	int dimx, dimy;
	int nrows;		/* number of rows applied to DQ array */
	int sameamp, samegain, samechip;

	int GetLT0 (Hdr *, double *, double *);
	int SameInt (int, int);
	int SameFlt (float, float);
	int SameString (char *, char *);

	int OpenBpixTab (char *, TblInfo *);
	int ReadBpixTab (TblInfo *, int, TblRow *);
	int CloseBpixTab (TblInfo *);

	void DQINormal (DQHdrData *, double *, TblRow *);
	void ToWF3RawCoords (WF3Info *, double *, TblRow *);

	/* We could still flag saturation even if the bpixtab was dummy. */
	if (wf3->dqicorr != PERFORM && wf3->dqicorr != DUMMY)
	    return (status);

	/* For the CCD, check for and flag saturation. */
	sat = wf3->saturate;
	if (wf3->detector != IR_DETECTOR) {
            dimx = x->sci.data.nx;
            dimy = x->sci.data.ny;

	    for (j = 0;  j < dimy;  j++) {
		 for (i = 0;  i < dimx;  i++) {
		      /* Flag a-to-d saturated pixels with 2048 dq bit */
		      if (Pix (x->sci.data, i, j) > ATOD_SATURATE) {
			  sum_dq = DQPix (x->dq.data, i, j) | ATODSAT;
			  DQSetPix (x->dq.data, i, j, sum_dq); /* atod sat */
		      }
		      /* Flag full-well or atod saturated pixels with 256 bit */
		      if (Pix (x->sci.data, i, j) > sat || 
			  Pix (x->sci.data, i, j) > ATOD_SATURATE) {
			  sum_dq = DQPix (x->dq.data, i, j) | SATPIXEL;
			  DQSetPix (x->dq.data, i, j, sum_dq);	/* saturated */
		      }
		 }
	    }
	}

	/* Get the linear transformations. */
	if (GetLT0 (&x->sci.hdr, ri_m, ri_v))		/* zero indexed LTV */
	    return (status);

	/* There might not be any bad pixel table.  If not, quit now. */
	if (wf3->bpix.exists == EXISTS_NO || wf3->dqicorr != PERFORM)
	    return (status);

	/* In some cases (when the science data are not binned) we can
	** set the DQ flags directly in the DQ array, but in other cases
	** (when the science data are binned) we must create a scratch
	** array and copy back to the original. */

	if (wf3->bin[0] == 1 && wf3->bin[1] == 1)
	    in_place = 1;				/* no binning */
	else
	    in_place = 0;

	/* Get the other linear transformations (ri_m and ri_v were gotten
	** earlier, just after checking for saturation). */

	if (!in_place) {

	    /* scratch array will be in reference table coords */
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

	/* Open the data quality initialization table, find columns, etc. */
	if (OpenBpixTab (wf3->bpix.name, &tabinfo))
	    return (status);

	/* Size of scratch image */
	snpix[0] = tabinfo.axlen1;
	snpix[1] = tabinfo.axlen2;

	/* Size of current image */
	npix[0] = x->dq.data.nx;
	npix[1] = x->dq.data.ny;

	/* Allocate space for scratch array */
	initShortHdrData (&ydq);
	if (!in_place) {
	    allocShortHdrData (&ydq, snpix[0], snpix[1], True);
	    if (hstio_err()) {
		trlerror ("doDQI couldn't allocate data quality array.");
		return (status = OUT_OF_MEMORY);
	    }
	    for (j=0; j < snpix[1]; j++)
		 for (i=0; i < snpix[0]; i++)
		      DQSetPix (ydq.data, i, j, 0);	/* initially OK */
	}

	/* Read each row of the table, and fill in data quality values. */
	nrows = 0;
	for (row = 1;  row <= tabinfo.nrows;  row++) {

	    if (ReadBpixTab (&tabinfo, row, &tabrow)) {
		trlerror ("Error reading BPIXTAB.");
		return (status);
	    }

	    /* If CCDAMP column does not exist or it has a wildcard value,
	    ** always return a match. */
	    if (tabinfo.cp_amp == 0 || SameString(tabrow.ccdamp, "N/A")) {
		sameamp = 1;
	    } else {
		sameamp = SameString (tabrow.ccdamp, wf3->ccdamp);
	    }

	    /* If CCDGAIN column does not exist or it has a wildcard value,
	    ** always return a match. */
	    if (tabinfo.cp_ccdgain == 0 || tabrow.ccdgain == -999) {
		samegain = 1;
	    } else {
		samegain = SameFlt (tabrow.ccdgain, wf3->ccdgain);
	    }

	    /* If CCDCHIP column does not exist or it has a wildcard value,
	    ** always return a match. */
	    if (tabinfo.cp_ccdchip == 0 || tabrow.ccdchip == -999) {
		samechip = 1;
	    } else {
		samechip = SameInt (tabrow.ccdchip, wf3->chip);
	    }

	    /* Check for a match with selection criteria */
	    if (sameamp && samegain && samechip) {

		/* Adjust BPIXTAB pixel coords for presence of serial
		   virtual overscan pixels in WFC3 raw images */
		ToWF3RawCoords (wf3, ri_m, &tabrow);

	        /* Assign the flag value to all relevant pixels. */
	        if (in_place)
		    DQINormal (&x->dq, ri_v, &tabrow);
	        else
	            DQINormal (&ydq, rs_v, &tabrow);
		nrows += 1;
	    }
	}

	if (CloseBpixTab (&tabinfo))		/* done with the table */
	    return (status);

	if (nrows == 0) {
	    sprintf (MsgText, "No rows from BPIXTAB applied to DQ array.");
	    trlwarn (MsgText);

	}

	/* Copy scratch contents into input DQ array */
	if (!in_place) {

	    /* Get corners of region of overlap between image
	       and scratch array */
	    FirstLast (si_m, si_v, snpix, npix, rbin, first, last, sfirst);

	    /* We have been writing to a scratch array ydq. Now copy
	       or bin the values down to the actual size of image x */

	    j0 = sfirst[1];
	    for (n = first[1]; n <= last[1]; n++) {
		 i0 = sfirst[0];
		 for (m = first[0]; m <= last[0]; m++) {
		      sum_dq = DQPix (x->dq.data, m, n);
		      for (j = j0; j < MIN(j0+rbin[1], ydq.data.ny); j++) {
			   for (i = i0; i < MIN(i0+rbin[0], ydq.data.nx); i++) {
				if (i >= 0 && j >= 0)
				    sum_dq |= DQPix (ydq.data, i, j);
		      }}
		      DQSetPix (x->dq.data, m, n, sum_dq);
		      i0 += rbin[0];
		 }
		 j0 += rbin[1];
	    }

	    freeShortHdrData (&ydq);		/* done with ydq */
	}

	return (status);
}

int OpenBpixTab (char *tname, TblInfo *tabinfo) {

	extern int status;

	tabinfo->tp = c_tbtopn (tname, IRAF_READ_ONLY, 0);
	if (c_iraferr())
	    return (status = OPEN_FAILED);

	tabinfo->nrows = c_tbpsta (tabinfo->tp, TBL_NROWS);
	if (tabinfo->nrows < 1) {
	    return (status);			/* nothing to do */
	}

	/* Find the columns. */
	c_tbcfnd1 (tabinfo->tp, "CCDAMP", &tabinfo->cp_amp);
	c_tbcfnd1 (tabinfo->tp, "CCDCHIP", &tabinfo->cp_ccdchip);
	c_tbcfnd1 (tabinfo->tp, "CCDGAIN", &tabinfo->cp_ccdgain);
	c_tbcfnd1 (tabinfo->tp, "PIX1", &tabinfo->cp_xstart);
	c_tbcfnd1 (tabinfo->tp, "PIX2", &tabinfo->cp_ystart);
	c_tbcfnd1 (tabinfo->tp, "LENGTH", &tabinfo->cp_length);
	c_tbcfnd1 (tabinfo->tp, "AXIS", &tabinfo->cp_axis);
	c_tbcfnd1 (tabinfo->tp, "VALUE", &tabinfo->cp_flag);
	if (tabinfo->cp_xstart == 0 ||
	    tabinfo->cp_ystart == 0 ||
	    tabinfo->cp_length == 0 ||
	    tabinfo->cp_axis == 0 ||
	    tabinfo->cp_ccdchip == 0 ||
	    tabinfo->cp_flag == 0) {
	    c_tbtclo (tabinfo->tp);
	    trlerror ("Column not found in BPIXTAB.");
	    return (status = COLUMN_NOT_FOUND);
	}

	/* Find out how large a full-size data quality array should be. */
	tabinfo->axlen1 = c_tbhgti (tabinfo->tp, "SIZAXIS1");
	if (c_iraferr ()) {
	    c_tbtclo (tabinfo->tp);
	    trlerror ("Couldn't get SIZAXIS1 from BPIXTAB header.");
	    return (status = TABLE_ERROR);
	}
	tabinfo->axlen2 = c_tbhgti (tabinfo->tp, "SIZAXIS2");
	if (c_iraferr ()) {
	    c_tbtclo (tabinfo->tp);
	    trlerror ("Couldn't get SIZAXIS2 from BPIXTAB header.");
	    return (status = TABLE_ERROR);
	}

	return (status);
}


/* This routine reads the relevant data from one row.  The starting
   pixel location, number of pixels and axis, and the flag value are read.
   The starting pixel (xstart, ystart) will be decremented by one to
   convert to zero-indexing.
*/

int ReadBpixTab (TblInfo *tabinfo, int row, TblRow *tabrow) {

	extern int status;

	/* read in values from selection columns */
	/* If AMP column exists, read it, otherwise it is not necessary */
	if (tabinfo->cp_amp != 0) {
	    c_tbegtt (tabinfo->tp, tabinfo->cp_amp, row, tabrow->ccdamp,
		      SZ_CBUF);
	    if (c_iraferr())
	        return (status = TABLE_ERROR);
	}
	c_tbegti (tabinfo->tp, tabinfo->cp_ccdchip, row, &tabrow->ccdchip);
	if (c_iraferr())
	    return (status = TABLE_ERROR);
	/* If GAIN column exists, read it, otherwise it is not necessary */
	if (tabinfo->cp_ccdgain != 0) {
	    c_tbegtr (tabinfo->tp, tabinfo->cp_ccdgain, row, &tabrow->ccdgain);
	    if (c_iraferr())
	        return (status = TABLE_ERROR);
	}
	
	/* read in bad pixel value description columns */
	c_tbegti (tabinfo->tp, tabinfo->cp_xstart, row, &tabrow->xstart);
	if (c_iraferr())
	    return (status = TABLE_ERROR);
	tabrow->xstart--;
	c_tbegti (tabinfo->tp, tabinfo->cp_ystart, row, &tabrow->ystart);
	if (c_iraferr())
	    return (status = TABLE_ERROR);
	tabrow->ystart--;
	c_tbegti (tabinfo->tp, tabinfo->cp_length, row, &tabrow->length);
	if (c_iraferr())
	    return (status = TABLE_ERROR);
	c_tbegti (tabinfo->tp, tabinfo->cp_axis, row, &tabrow->axis);
	if (c_iraferr())
	    return (status = TABLE_ERROR);
	c_tbegts (tabinfo->tp, tabinfo->cp_flag, row, &tabrow->flag);
	if (c_iraferr())
	    return (status = TABLE_ERROR);

	if (tabrow->axis !=1 && tabrow->axis != 2) {
	    sprintf (MsgText, "Axis = %d in BPIXTAB, but it must be 1 or 2.",
		     tabrow->axis);
	    trlerror (MsgText);
	    c_tbtclo (tabinfo->tp);
	    return (status = TABLE_ERROR);
	}
	if (tabrow->length <= 0) {
	    sprintf (MsgText,"Length = %d in BPIXTAB, but it must be positive.",
		     tabrow->length);
	    trlerror (MsgText);
	    c_tbtclo (tabinfo->tp);
	    return (status = TABLE_ERROR);
	}

	return (status);
}

/* This routine closes the bpixtab table. */

int CloseBpixTab (TblInfo *tabinfo) {

	extern int status;

	c_tbtclo (tabinfo->tp);
	if (c_iraferr())
	    return (status = TABLE_ERROR);

	return (status);
}

/* This routine adds an offset to the DQ pixel coords read from the
** BPIXTAB to account for the presence of serial virtual overscan columns
** in some modes of WFC3 raw images.
*/

void ToWF3RawCoords (WF3Info *wf3, double *ltm, TblRow *tabrow) {

/* arguments:
WF3Info *wf3	     i: WFC3 calibration info
double ltm[2]        i: scale factor of science data
TblRow *tabrow       io: data quality info read from one row of BPIXTAB
*/

	int rbin; 	/* binning factor */

	int streq_ic (char *, char *);	/* string comparison */

	rbin = NINT (1. / ltm[0]);

	/* The correction is only needed for four-amp readouts, because they
	** are the only ones to have virtual overscan columns in the middle
	** of an image. All other two- and one-amp modes have the virtual
	** columns at either end of the image and are already accounted for
	** in the LTV values in the raw image header. */
	if (streq_ic (wf3->ccdamp, "ABCD")) {

	    /* Only apply the offset if the flagged pixel is located to
	    ** to the right (i.e. higher pixel index) of the serial virtual
	    ** overscan columns. This applies to all pixels in the domain
	    ** of the second readout amp for the chip (i.e. amp B for AB,
	    ** and amp D for CD. */
	    if (tabrow->xstart >= wf3->ampx * rbin) {

		/* Add an offset equal to the number of virtual overscan
		** columns that occur in the middle of an unbinned raw image.
		** We get the number of colums from the trim information in
		** the OSCNTAB reference table, multiplied back up to
		** unbinned space. */
		tabrow->xstart += (wf3->trimx[2] + wf3->trimx[3]) * rbin;

		/* Raw images with a binning factor of 2 use a smaller
		** trim value, therefore we need to add an extra offset. */
		if (rbin == 2)
		    tabrow->xstart += 2;
	    }
	}
}

/* This routine assigns data quality values as specified in one row of
   the DQI table.  The current value is ORed with any previous value.
*/

void DQINormal (DQHdrData *ydq, double *ltv, TblRow *tabrow) {

/* arguments:
DQHdrData *ydq	    io: data quality array
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
	    if (xhigh < 0 || xlow >= ydq->data.nx ||
		ystart < 0 || ystart >= ydq->data.ny)
		    return;

	    if (xlow < 0)
	    	xlow = 0;
	    if (xhigh >= ydq->data.nx)
    		xhigh = ydq->data.nx - 1;

	    j = ystart;
	    for (i = xlow;  i <= xhigh;  i++) {
		    sum_dq = tabrow->flag | PDQPix (&ydq->data, i, j);
		    PDQSetPix (&ydq->data, i, j, sum_dq);
	    }

	} else if (tabrow->axis == 2) {

	    /* The range for y is the repeat count. */
	    ylow = ystart;
	    yhigh = ystart + tabrow->length - 1;

	    /* out of range? */
	    if (xstart < 0 || xstart >= ydq->data.nx ||
		yhigh < 0 || ylow >= ydq->data.ny)
		    return;

	    if (ylow < 0)
		    ylow = 0;
	    if (yhigh >= ydq->data.ny)
		    yhigh = ydq->data.ny - 1;

	    i = xstart;
	    for (j = ylow;  j <= yhigh;  j++) {
		    sum_dq = tabrow->flag | PDQPix (&ydq->data, i, j);
		    PDQSetPix (&ydq->data, i, j, sum_dq);
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

