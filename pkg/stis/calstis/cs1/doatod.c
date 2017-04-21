# include <stdio.h>
# include <stdlib.h>		/* for malloc */
# include <math.h>		/* for fabs */

# include "hstio.h"
# include "c_iraf.h"
# include "xtables.h"
# include "stis.h"
# include "calstis1.h"
# include "err.h"
# include "stisdq.h"		/* for SATPIXEL */
# include "stisdef.h"

typedef struct {
	IRAFPointer tp;			/* pointer to table descriptor */
	IRAFPointer cp_amp;		/* pointers to column descriptors */
	IRAFPointer cp_gain;
	IRAFPointer cp_key;
	IRAFPointer cp_keyval;
	IRAFPointer cp_nelem;
	IRAFPointer cp_atod;
	IRAFPointer cp_pedigree;
	IRAFPointer cp_descrip;
	int nrows;			/* number of rows in table */
} TblInfo;

/* This is information used to identify the correct table row. */
typedef struct {
	char ccdamp[STIS_CBUF+1]; 	/* CCD amplifier read out (A,B,C,D) */
	int ccdgain;			/* commanded gain of the CCD */
	char ref_key[STIS_CBUF+1];	/* name of image header keyword */
	double ref_key_value;		/* compare with ref_key from image */
} TblRow;

/* This is the array of corrected values and the size of that array. */
typedef struct {
	int nelem;			/* length of the correction array */
	float *atod;			/* the correction array */
} TblArray;

static int OpenAtoDTab (char *, TblInfo *);
static int ReadAtoDTab (TblInfo *, int, TblRow *);
static int ReadAtoDArray (TblInfo *, int, TblArray *);
static int CloseAtoDTab (TblInfo *);

/* This routine gets the analog-to-digital correction from a table
   and applies it to each element of the SCI data extension.

   The A-to-D table should contain the following:
	header parameters:
		none needed
	columns:
		CCDAMP:  identifies which amp was used (string A-D)
		CCDGAIN:  commanded gain of the CCD (int)
		REF_KEY:  image header keyword name (char array)
		REF_KEY_VALUE:  for comparison with image header value
		NELEM:  number of elements in the A-to-D array (int)
		ATOD:  the array of A-to-D corrections (float)

   The table is read to find the row for which CCDAMP and CCDGAIN match
   the amplifier and gain (previously read from the image header).
   Then that row is read to get the REF_KEY string and REF_KEY_VALUE.
   The REF_KEY string is the name of a keyword in the image primary
   header.  The value of that keyword is read from the image header,
   and the absolute value of the difference between that value and the
   value in the table is taken.  The row for which that difference is
   minimum is the row from which the correction array is read.  The
   array length is read from that row as well, and the length may differ
   from row to row, although the maximum length is fixed.

   For each pixel in the SCI extension, the input pixel value (an integer)
   is used as an index into the correction array, and the output pixel
   value is the ATOD correction array value for that index.

   Values that on input are greater than or equal to the number of elements
   in the ATOD array will be set to the value of the last element of that
   array, and they will be flagged as saturated.

   Phil Hodge, 1998 Sept 28:
	When flagging a saturated pixel, OR the current data quality value
	with SATPIXEL, instead of just assigning SATPIXEL to the dq array.

   Phil Hodge, 1998 Oct 16:
	Change status value 1010 to GENERIC_ERROR_CODE.
	In doAtoD, call GetKeyD outside of the if block.

   Phil Hodge, 2000 Jan 13:
	Add one to ccdamp and ref_key buffer sizes.
*/

int doAtoD (StisInfo1 *sts, SingleGroup *x) {

/* arguments:
StisInfo1 *sts     i: calibration switches, etc
SingleGroup *x    io: image to be calibrated; written to in-place
*/

	int status;

	TblInfo tabinfo;	/* pointer to table descriptor, etc */
	TblRow tabrow;		/* values read from a table row */
	TblArray tabarray;	/* correction array read from table row */

	int foundit;		/* row found in table? */
	int row;		/* loop index for row number */
	int row_min;		/* row with closest temperature (min dt) */
	double ref_key_value;	/* value gotten from image header */
	double dt, dt_min;	/* for finding temperature in table */
	int ival;		/* input science data value from x */
	int i, j;
	short dq;		/* a data quality value */

	int no_default = 0;	/* missing keyword is fatal error */

	if (sts->atodcorr != PERFORM)
	    return (0);

	if (sts->ncombine > 1) {
	    printf (
"ERROR    NCOMBINE is already > 1 before ATODCORR has been performed.\n");
	    return (GENERIC_ERROR_CODE);
	}

	/* Open the A-to-D table. */
	if ((status = OpenAtoDTab (sts->atod.name, &tabinfo)))
	    return (status);

	/* Find the row with value closest to the temperature. */

	foundit = 0;
	for (row = 1;  row <= tabinfo.nrows;  row++) {
	    if ((status = ReadAtoDTab (&tabinfo, row, &tabrow)))
		return (status);
	    if (SameString (tabrow.ccdamp, sts->ccdamp) &&
		SameInt (tabrow.ccdgain, sts->ccdgain)) {
		/* Get value from header. */
		if ((status = Get_KeyD (&x->sci.hdr, tabrow.ref_key, no_default,
                                        0., &ref_key_value)))
		    return (status);
		dt = fabs (ref_key_value - tabrow.ref_key_value);
		if (!foundit) {
		    foundit = 1;	/* assign initial values */
		    dt_min = dt;
		    row_min = row;
		} else if (dt < dt_min) {	/* update dt_min */
		    dt_min = dt;
		    row_min = row;
		}
	    }
	}

	if (!foundit) {
	    printf (
	"ERROR    CCD amp %s, gain %d, not found in ATODTAB `%s'.\n",
		sts->ccdamp, sts->ccdgain, sts->atod.name);
	    CloseAtoDTab (&tabinfo);
	    return (TABLE_ERROR);
	}

	/* Get pedigree & descrip from the row. */
	if ((status = RowPedigree (&sts->atod, row_min,
                tabinfo.tp, tabinfo.cp_pedigree, tabinfo.cp_descrip)))
	    return (status);
	if (sts->atod.goodPedigree == DUMMY_PEDIGREE) {
	    sts->atodcorr = DUMMY;
	    CloseAtoDTab (&tabinfo);
	    return (0);
	}

	/* Reread the appropriate row to get the correction array. */
	if ((status = ReadAtoDArray (&tabinfo, row_min, &tabarray)))
	    return (status);

	/* Apply this correction to each pixel in the image.  At this
	   stage the values should still be integers, so assigning
	   a value to an integer (ival) should not result in truncation.
	*/
	for (j = 0;  j < x->sci.data.ny;  j++) {
	    for (i = 0;  i < x->sci.data.nx;  i++) {
		ival = (int) Pix (x->sci.data, i, j);
		if (ival >= tabarray.nelem) {
		    Pix (x->sci.data, i, j) = tabarray.atod[tabarray.nelem-1];
		    dq = DQPix (x->dq.data, i, j) | SATPIXEL;
		    DQSetPix (x->dq.data, i, j, dq);	/* saturated */
		} else if (ival >= 0) {
		    Pix (x->sci.data, i, j) = tabarray.atod[ival];
		}		/* else if ival < 0, no change */
	    }
	}

	free (tabarray.atod);
	if ((status = CloseAtoDTab (&tabinfo)))
	    return (status);

	return (0);
}

/* This routine opens the analog to digital correction table, finds the
   columns that we need, and gets the total number of rows in the table.
   The columns are CCDAMP, CCDGAIN, REF_KEY, REF_KEY_VALUE, NELEM, and
   ATOD.
*/

static int OpenAtoDTab (char *tname, TblInfo *tabinfo) {

	tabinfo->tp = c_tbtopn (tname, IRAF_READ_ONLY, 0);
	if (c_iraferr()) {
	    printf ("ERROR    ATODTAB `%s' not found.\n", tname);
	    return (OPEN_FAILED);
	}

	tabinfo->nrows = c_tbpsta (tabinfo->tp, TBL_NROWS);

	/* Find the columns. */
	c_tbcfnd1 (tabinfo->tp, "CCDAMP", &tabinfo->cp_amp);
	c_tbcfnd1 (tabinfo->tp, "CCDGAIN", &tabinfo->cp_gain);
	c_tbcfnd1 (tabinfo->tp, "REF_KEY", &tabinfo->cp_key);
	c_tbcfnd1 (tabinfo->tp, "REF_KEY_VALUE", &tabinfo->cp_keyval);
	c_tbcfnd1 (tabinfo->tp, "NELEM", &tabinfo->cp_nelem);
	c_tbcfnd1 (tabinfo->tp, "ATOD", &tabinfo->cp_atod);
	if (tabinfo->cp_amp == 0 ||
	    tabinfo->cp_gain == 0 ||
	    tabinfo->cp_key == 0 ||
	    tabinfo->cp_keyval == 0 ||
	    tabinfo->cp_nelem == 0 ||
	    tabinfo->cp_atod == 0) {
	    printf ("ERROR    Column not found in ATODTAB.\n");
	    c_tbtclo (tabinfo->tp);
	    return (COLUMN_NOT_FOUND);
	}

	/* Pedigree and descrip are optional columns. */
	c_tbcfnd1 (tabinfo->tp, "PEDIGREE", &tabinfo->cp_pedigree);
	c_tbcfnd1 (tabinfo->tp, "DESCRIP", &tabinfo->cp_descrip);

	return (0);
}

/* This routine reads the relevant data from one row.  The amplifier
   number, gain, and reference key value are gotten.
*/

static int ReadAtoDTab (TblInfo *tabinfo, int row, TblRow *tabrow) {

	c_tbegtt (tabinfo->tp, tabinfo->cp_amp, row,
			tabrow->ccdamp, STIS_CBUF);
	if (c_iraferr())
	    return (TABLE_ERROR);

	c_tbegti (tabinfo->tp, tabinfo->cp_gain, row, &tabrow->ccdgain);
	if (c_iraferr())
	    return (TABLE_ERROR);

	c_tbegtt (tabinfo->tp, tabinfo->cp_key, row,
			tabrow->ref_key, STIS_CBUF);
	if (c_iraferr())
	    return (TABLE_ERROR);

	c_tbegtd (tabinfo->tp, tabinfo->cp_keyval, row, &tabrow->ref_key_value);
	if (c_iraferr())
	    return (TABLE_ERROR);

	return (0);
}

/* This routine reads the array data from one row.  The number of elements
   in the array is gotten, the array is allocated, and the correction data
   are read into the array.
*/

static int ReadAtoDArray (TblInfo *tabinfo, int row, TblArray *tabarray) {

	int nret;		/* number of elements actually read */

	/* Find out how many elements there are in the ATOD array,
	   and allocate space for the array to be read from the table.
	*/
	c_tbegti (tabinfo->tp, tabinfo->cp_nelem, row, &tabarray->nelem);
	if (c_iraferr())
	    return (TABLE_ERROR);

	tabarray->atod = (float *) malloc (tabarray->nelem * sizeof(float));
	if (tabarray->atod == NULL)
	    return (OUT_OF_MEMORY);

	nret = c_tbagtr (tabinfo->tp, tabinfo->cp_atod, row,
			tabarray->atod, 1, tabarray->nelem);
	if (c_iraferr())
	    return (TABLE_ERROR);

	if (nret < tabarray->nelem) {
	    printf (
"ERROR    CORRECTION array in row %d of ATODTAB is too short.\n", row);
	    free (tabarray->atod);
	    return (TABLE_ERROR);
	}

	return (0);
}

/* This routine closes the atodtab table. */

static int CloseAtoDTab (TblInfo *tabinfo) {

	c_tbtclo (tabinfo->tp);
	if (c_iraferr())
	    return (TABLE_ERROR);

	return (0);
}
