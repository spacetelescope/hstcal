/* This file contains:
 doDQI
 OpenBpixTab
 ReadBpixTab
 CloseBpixTab
 DQINormal
 */

# include <stdio.h>
# include <string.h>

#include "hstcal.h"
# include "xtables.h"
# include "hstio.h"
# include "acs.h"
# include "acsinfo.h"
# include "hstcalerr.h"
# include "acsdq.h"

typedef struct {
    IRAFPointer tp;			/* pointer to table descriptor */
    IRAFPointer cp_xstart, cp_ystart; /* starting pixel for bad values */
    IRAFPointer cp_length;		/* this many bad values */
    IRAFPointer cp_axis;		/* repeat along X or Y axis (1 or 2) */
    IRAFPointer cp_flag;		/* this is the data quality value */
    IRAFPointer cp_amp;		/* selection columns */
    IRAFPointer cp_ccdchip;
    IRAFPointer cp_ccdgain;
    int intgain;
    int nrows;			/* number of rows in table */
} TblInfo;

typedef struct {
    /* values read from a table row */
    int xstart, ystart;
    int length;
    int axis;
    short flag;
    /* values used for selecting row to be used.*/
    char ccdamp[ACS_CBUF];
    /*
     For tables with int ccdgain, set intgain=1 and populate ccdgaini.
     For tables with float ccdgain, set intgain=0 and populate ccdgainf.
     Based on value of intgain, use either ccdgaini or ccdgainf in code.
     */
    int ccdgaini;
    float ccdgainf;
    int ccdchip;
} TblRow;


static int OpenBpixTab (char *, TblInfo *);
static int ReadBpixTab (TblInfo *, int, TblRow *);
static int CloseBpixTab (TblInfo *);

static void DQINormal (DQHdrData *, double *, TblRow *);

/* This routine ORs the data quality array in the input SingleGroup x
 with the pixels flagged in the data quality initialization table.

 The data quality initialization (BPIXTAB) table should contain the
 following integer header parameters:
 SIZAXIS1, SIZAXIS2:  full size of data quality array;
 these are expected to be 4096 x 2048 for WFC data
 and 1024 x 1024 for MAMA and HRC data
 and the following integer columns:
 PIX1, PIX2:  starting pixel (one indexed) of a range
 of pixels to be assigned an initial value
 LENGTH:  number of pixels to be assigned the value
 AXIS:  1 --> X axis, 2 --> Y axis
 VALUE:  value to be ORed with data quality array
 Each row of the table specifies a set of LENGTH pixels to be flagged
 in either the X or Y direction with the value VALUE.  The value assigned
 for each pixel will be the OR of VALUE and any previous value at that
 pixel.  For MAMA data, the units for PIX1, PIX2 and LENGTH will be
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

 2001 December 3, WJH:
 Updated to apply only those rows which have CCDAMP,
 CCDGAIN, and (most importantly) CCDCHIP values corresponding to
 the input science data.

 2003 November 5, WJH:
 Modified to optionally check against CCDAMP and CCDGAIN values.
 This would support tables without these columns to apply to all
 configurations, yet still work with current tables with all columns.

 2020 April 29, MDD:
 Modified to accommodate the A-to-D saturation threshold now being defined
 in the CCDTAB. The threshold is no longer a constant and has different
 behavior pre- and post-SM4.

 2020 May 14, MDD:
 Removed flagging of full-well saturated pixels based upon a science
 pixel value being greater than a defined scalar value.  Use of a
 full-well saturation image supersedes the functionality previously
 in this routine.
 */

int doDQI (ACSInfo *acs, SingleGroup *x) {

    /* arguments:
     ACSInfo *acs    i: calibration switches, etc
     SingleGroup *x    io: image to be calibrated; DQ array written to in-place
     */

    extern int status;

    TblInfo tabinfo;			/* pointer to table descriptor, etc */
    TblRow tabrow;				/* values read from a table row */

    /* mappings from one coordinate system to another */
    double ri_m[2], ri_v[2];	/* reference to image */
    int npix_x, npix_y;			/* size of current image */

    int i, j;					/* indexes for scratch array ydq */
    short sum_dq;				/* for binning data quality array */
    int atod_sat;

    int row;					/* loop index for row number */
    int dimx, dimy;
    int nrows;   				/* number of rows applied to DQ array */
    short dq_fill = 64;         /* default fill value when no rows are applied*/
    int xpos,ypos;
    int sameamp, samegain, samechip;

    int GetLT0 (Hdr *, double *, double *);
    int SameInt (int, int);
    int SameFlt (float, float);
    int SameString (char *, char *);

    /* We could still flag saturation even if the bpixtab was dummy. */
    if (acs->dqicorr != PERFORM && acs->dqicorr != DUMMY)
        return (status);

    /* For the CCD, check for and flag saturation. */
    atod_sat = acs->atod_saturate;

    if (acs->detector != MAMA_DETECTOR) {
        dimx = x->sci.data.nx;
        dimy = x->sci.data.ny;
        for (j = 0;  j < dimy;  j++) {
            for (i = 0;  i < dimx;  i++) {
                /* Flag a-to-d saturated pixels with 2048 dq bit */
                if (Pix (x->sci.data, i, j) >= atod_sat) {
                    sum_dq = DQPix (x->dq.data, i, j) | ATODSAT;
                    DQSetPix (x->dq.data, i, j, sum_dq);	/* a-to-d saturated */
                }
            }
        }
    }

    /* Get the linear transformations. */
    if (GetLT0 (&x->sci.hdr, ri_m, ri_v))		/* zero indexed LTV */
        return (status);

    /* size of current image */
    npix_x = x->dq.data.nx;
    npix_y = x->dq.data.ny;

    /* Flag regions beyond the bounderies of the aperture, for
    CCD imaging type observations.

    if (acs->detector != MAMA_DETECTOR) {
    FlagFilter (acs, &x->dq.data, npix_x, npix_y, ri_m, ri_v);
    }
    */

    /* There might not be any bad pixel table.  If not, quit now. */
    if (acs->bpix.exists == EXISTS_NO || acs->dqicorr != PERFORM)
        return (status);

    /* Open the data quality initialization table, find columns, etc. */
    if (OpenBpixTab (acs->bpix.name, &tabinfo))
        return (status);

    nrows = 0;
    /* Read each row of the table, and fill in data quality values. */

    for (row = 1;  row <= tabinfo.nrows;  row++) {

        if (ReadBpixTab (&tabinfo, row, &tabrow)) {
            trlerror ("Error reading BPIXTAB.");
            return (status);
        }

        /* If CCDAMP column does not exist, always return a match.*/
        if (tabinfo.cp_amp == 0 || SameString(tabrow.ccdamp,"N/A")) {
            sameamp = 1;
        } else {
            sameamp = SameString (tabrow.ccdamp, acs->ccdamp);
        }

        /* If CCDGAIN column does not exist, always return a match. */
        if (tabinfo.cp_ccdgain == 0) {
            /* No ccdgain column at all, set samegain to 1(yes) */
            samegain = 1;
        } else {
            /* We have a ccdgain column, check the value */
            if(tabinfo.intgain == 1) {
                if (tabrow.ccdgaini == -999) {
                    samegain = 1;
                } else {
                    samegain = SameInt (tabrow.ccdgaini, (int)acs->ccdgain);
                }
            } else {
                samegain = SameFlt(tabrow.ccdgainf,-999.0);
                if (samegain == 0) {
                    samegain = SameFlt (tabrow.ccdgainf, acs->ccdgain);
                }
            }
        }
        /* If CCDCHIP column does not exist, always return a match. */
        if (tabinfo.cp_ccdchip == 0 || tabrow.ccdchip == -999) {
            samechip = 1;
        } else {
            samechip = SameInt (tabrow.ccdchip, acs->chip);
        }

        /* Check whether the row matches the conditions. */
        if (sameamp && samegain && samechip) {
            /* Assign the flag value to all relevant pixels. */
            DQINormal (&x->dq, ri_v, &tabrow);
            nrows += 1;
        }
    }

    if (CloseBpixTab (&tabinfo))		/* done with the table */
        return (status);

    if (nrows == 0) {
        sprintf(MsgText,"No rows from BPIXTAB applied to DQ array.");
        trlwarn(MsgText);
        /* This code will mark the first pixel with a value of 64
         to prevent CALACS from crashing when no pixels are marked bad. */
        sprintf(MsgText,"Inserting single-pixel DQ place-holder at (1,1).");
        trlwarn(MsgText);
        xpos = (int)ri_m[0];
        ypos = (int)ri_m[1];
        sum_dq = DQPix (x->dq.data, xpos, ypos) | dq_fill;
        DQSetPix (x->dq.data, xpos, ypos, sum_dq);

    }

    return (status);
}

static int OpenBpixTab (char *tname, TblInfo *tabinfo) {

    extern int status;
    int colnum, datatype, lendata, lenfmt;
    char *colname;
    char *colunits;
    char *colfmt;

    if ((colname = calloc (SZ_COLNAME+1, sizeof(char))) == NULL) {
        trlerror ("Out of memory.\n");
        return (OUT_OF_MEMORY);
    }
    if ((colunits = calloc (ACS_CBUF+1, sizeof(char))) == NULL) {
        trlerror ("Out of memory.\n");
        return (OUT_OF_MEMORY);
    }
    if ((colfmt = calloc (ACS_CBUF+1, sizeof(char))) == NULL) {
        trlerror ("Out of memory.\n");
        return (OUT_OF_MEMORY);
    }

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
    /* get info on ccdgain column to determine whether we
     have int or float values to read in.
     */

    c_tbcinf(tabinfo->cp_ccdgain, &colnum, colname, colunits, colfmt, &datatype, &lendata, &lenfmt);
    if (datatype == IRAF_INT) {
        tabinfo->intgain = 1;
    } else {
        tabinfo->intgain = 0;
    }
    free(colname);
    free(colunits);
    free(colfmt);

    return (status);
}


/* This routine reads the relevant data from one row.  The starting
 pixel location, number of pixels and axis, and the flag value are read.
 The starting pixel (xstart, ystart) will be decremented by one to
 convert to zero-indexing.
 */

static int ReadBpixTab (TblInfo *tabinfo, int row, TblRow *tabrow) {

    extern int status;

    /* read in values from selection columns...*/
    /* If this column exists, read it, otherwise it is not necessary */
    if (tabinfo->cp_amp != 0) {
        c_tbegtt (tabinfo->tp, tabinfo->cp_amp, row,
                  tabrow->ccdamp, ACS_CBUF-1);
        if (c_iraferr())
            return (status = TABLE_ERROR);
    }
    /* If this column exists, read it, otherwise it is not necessary */
    if (tabinfo->cp_ccdgain != 0) {
        if (tabinfo->intgain == 1) {
            c_tbegti (tabinfo->tp, tabinfo->cp_ccdgain, row, &tabrow->ccdgaini);
        } else {
            c_tbegtr (tabinfo->tp, tabinfo->cp_ccdgain, row, &tabrow->ccdgainf);
        }
        if (c_iraferr())
            return (status = TABLE_ERROR);
    }

    c_tbegti (tabinfo->tp, tabinfo->cp_ccdchip, row, &tabrow->ccdchip);
    if (c_iraferr())
        return (status = TABLE_ERROR);

    /* read in bad pixel value description columns*/
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
        sprintf (MsgText, "Axis = %d in BPIXTAB, but it must be 1 or 2.", tabrow->axis);
        trlerror (MsgText);
        c_tbtclo (tabinfo->tp);
        return (status = TABLE_ERROR);
    }
    if (tabrow->length <= 0) {
        sprintf (MsgText,"Length = %d in BPIXTAB, but it must be positive.", tabrow->length);
        trlerror (MsgText);
        c_tbtclo (tabinfo->tp);
        return (status = TABLE_ERROR);
    }

    return (status);
}

/* This routine closes the bpixtab table. */

static int CloseBpixTab (TblInfo *tabinfo) {

    extern int status;

    c_tbtclo (tabinfo->tp);
    if (c_iraferr())
        return (status = TABLE_ERROR);

    return (status);
}


/* This routine assigns data quality values as specified in one row of
 the DQI table.  The current value is ORed with any previous value.
 */

static void DQINormal (DQHdrData *ydq, double *ltv, TblRow *tabrow) {

    /* arguments:
     DQHdrData *ydq  	io: data quality array
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
        if (xhigh < 0 || xlow >= ydq->data.nx || ystart < 0 || ystart >= ydq->data.ny)
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
        if (xstart < 0 || xstart >= ydq->data.nx || yhigh < 0 || ylow >= ydq->data.ny)
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
