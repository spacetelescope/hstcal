# include <stdio.h>
# include <stdlib.h>		/* for malloc */
# include <math.h>		/* for fabs */
#include "hstcal.h"
# include "xtables.h"

# include "hstio.h"
# include "wf3.h"
# include "wf3info.h"
# include "hstcalerr.h"
# include "wf3dq.h"		/* for SATPIXEL */

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
	char ccdamp[SZ_CBUF+1]; 	/* CCD amplifier read out (A,B,C,D) */
	float ccdgain;			/* commanded gain of the CCD */
	char ref_key[SZ_CBUF+1];	/* name of image header keyword */
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
		CCDGAIN:  commanded gain of the CCD (float)
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


   Original CALSTIS version slightly modified for ACS use.
   ACS version slightly modified for WFC3 use.

   H.Bushouse 2003 Oct 16: Modified to use floating-point gain values
                           for WFC3.
*/

int doAtoD (WF3Info *wf3, SingleGroup *x) {

/* arguments:
WF3Info *wf3	 i: calibration switches, etc
SingleGroup *x	io: image to be calibrated; written to in-place
*/

	extern int status;

	TblInfo tabinfo;	/* pointer to table descriptor, etc */
	TblRow tabrow;		/* values read from a table row */
	TblArray tabarray;	/* correction array read from table row */

	int foundit;		/* row found in table? */
	int row;		/* loop index for row number */
	int row_min;		/* row with matching keyword */
	double ref_key_value;	/* value gotten from image header */
	double dt, dt_min;	/* for finding desired row in table */
	int ival;		/* input science data value from x */
	int i, j;
	short dqval;

	int GetKeyDbl (Hdr *, char *, int, double, double *);
	int RowPedigree (RefTab *, int, IRAFPointer, IRAFPointer, IRAFPointer);
	int SameFlt (float, float);
	int SameString (char *, char *);

    row_min=0;
    dt_min=0.0f;

	if (wf3->atodcorr != PERFORM)
	    return (status);

	if (wf3->ncombine > 1) {
	    trlerror("NCOMBINE is already > 1 before ATODCORR has been performed.");
	    return (status = 1010);
	}

	/* Open the A-to-D table. */
	if (OpenAtoDTab (wf3->atod.name, &tabinfo))
	    return (status);

	/* Find the row with value closest to the temperature. */

	foundit = 0;
	for (row = 1;  row <= tabinfo.nrows;  row++) {
	    if (ReadAtoDTab (&tabinfo, row, &tabrow))
		return (status);
	    if (SameString (tabrow.ccdamp, wf3->ccdamp) &&
		SameFlt (tabrow.ccdgain, wf3->ccdgain)) {
		if (GetKeyDbl (x->globalhdr, tabrow.ref_key, NO_DEFAULT,
			       0., &ref_key_value))
		    return (status);
		if (!foundit) {
		    /* assign initial values */
		    foundit = 1;
		    row_min = row;
		    dt_min = fabs (ref_key_value - tabrow.ref_key_value);
		} else {
		    /* Get value from image, and update dt_min. */
		    dt = fabs (ref_key_value - tabrow.ref_key_value);
		    if (dt < dt_min) {
			dt_min = dt;
			row_min = row;
		    }
		}
	    }
	}

	if (!foundit) {
	    trlerror("CCD amp %s, gain %g, not found in ATODTAB `%s'.", wf3->ccdamp, wf3->ccdgain, wf3->atod.name);
	    CloseAtoDTab (&tabinfo);
	    return (status = TABLE_ERROR);
	}

	/* Get pedigree & descrip from the row. */
	if (RowPedigree (&wf3->atod, row_min,
		tabinfo.tp, tabinfo.cp_pedigree, tabinfo.cp_descrip))
	    return (status);
	if (wf3->atod.goodPedigree == DUMMY_PEDIGREE) {
	    wf3->atodcorr = DUMMY;
	    CloseAtoDTab (&tabinfo);
	    return (status);
	}

	/* Reread the appropriate row to get the correction array. */
	if (ReadAtoDArray (&tabinfo, row_min, &tabarray))
	    return (status);

	/* Apply this correction to each pixel in the image.  At this
	   stage the values should still be integers, so assigning
	   a value to an integer (ival) should not result in truncation.
	*/
	dqval = 0;
	for (j = 0;  j < x->sci.data.ny;  j++) {
	    for (i = 0;  i < x->sci.data.nx;  i++) {
		ival = (int) Pix (x->sci.data, i, j);
		if (ival >= tabarray.nelem) {
		    Pix (x->sci.data, i, j) = tabarray.atod[tabarray.nelem-1];
			dqval = SATPIXEL | DQPix (x->dq.data, i, j);
		    DQSetPix (x->dq.data, i, j, dqval);	/* saturated */
		} else if (ival >= 0) {
		    Pix (x->sci.data, i, j) = tabarray.atod[ival];
		}		/* else if ival < 0, no change */
	    }
	}

	free (tabarray.atod);
	if (CloseAtoDTab (&tabinfo))
	    return (status);

	return (status);
}

/* This routine opens the analog to digital correction table, finds the
   columns that we need, and gets the total number of rows in the table.
   The columns are CCDAMP, CCDGAIN, REF_KEY, REF_KEY_VALUE, NELEM, and
   ATOD.
*/

static int OpenAtoDTab (char *tname, TblInfo *tabinfo) {

	extern int status;

	tabinfo->tp = c_tbtopn (tname, IRAF_READ_ONLY, 0);
	if (c_iraferr()) {
	    trlerror("ATODTAB `%s' not found.", tname);
	    return (status = OPEN_FAILED);
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
	    trlerror("Column not found in ATODTAB.");
	    c_tbtclo (tabinfo->tp);
	    return (status = COLUMN_NOT_FOUND);
	}

	/* Pedigree and descrip are optional columns. */
	c_tbcfnd1 (tabinfo->tp, "PEDIGREE", &tabinfo->cp_pedigree);
	c_tbcfnd1 (tabinfo->tp, "DESCRIP", &tabinfo->cp_descrip);

	return (status);
}

/* This routine reads the relevant data from one row.  The amplifier
   number, gain, and reference key value are gotten.
*/

static int ReadAtoDTab (TblInfo *tabinfo, int row, TblRow *tabrow) {

	extern int status;

	c_tbegtt (tabinfo->tp, tabinfo->cp_amp, row,
			tabrow->ccdamp, SZ_CBUF);
	if (c_iraferr())
	    return (status = TABLE_ERROR);

	c_tbegtr (tabinfo->tp, tabinfo->cp_gain, row, &tabrow->ccdgain);
	if (c_iraferr())
	    return (status = TABLE_ERROR);

	c_tbegtt (tabinfo->tp, tabinfo->cp_key, row,
			tabrow->ref_key, SZ_CBUF);
	if (c_iraferr())
	    return (status = TABLE_ERROR);

	c_tbegtd (tabinfo->tp, tabinfo->cp_keyval, row, &tabrow->ref_key_value);
	if (c_iraferr())
	    return (status = TABLE_ERROR);

	return (status);
}

/* This routine reads the array data from one row.  The number of elements
   in the array is gotten, the array is allocated, and the correction data
   are read into the array.
*/

static int ReadAtoDArray (TblInfo *tabinfo, int row, TblArray *tabarray) {

	extern int status;
	int nret;		/* number of elements actually read */

	/* Find out how many elements there are in the ATOD array,
	   and allocate space for the array to be read from the table.
	*/
	c_tbegti (tabinfo->tp, tabinfo->cp_nelem, row, &tabarray->nelem);
	if (c_iraferr())
	    return (status = TABLE_ERROR);

	tabarray->atod = (float *) malloc (tabarray->nelem * sizeof(float));
	if (tabarray->atod == NULL)
	    return (status = OUT_OF_MEMORY);

	nret = c_tbagtr (tabinfo->tp, tabinfo->cp_atod, row,
			tabarray->atod, 1, tabarray->nelem);
	if ( (status = c_iraferr()) )
	    return (status = TABLE_ERROR);

	if (nret < tabarray->nelem) {
	    trlerror("CORRECTION array in row %d of ATODTAB is too short.", row);
	    free (tabarray->atod);
	    return (status = TABLE_ERROR);
	}

	return (status);
}

/* This routine closes the atodtab table. */

static int CloseAtoDTab (TblInfo *tabinfo) {

	extern int status;

	c_tbtclo (tabinfo->tp);
	if (c_iraferr())
	    return (status = TABLE_ERROR);

	return (status);
}
