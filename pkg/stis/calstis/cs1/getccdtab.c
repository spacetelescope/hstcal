# include <stdio.h>

# include "c_iraf.h"
# include "hstio.h"
# include "xtables.h"
# include "stis.h"
# include "calstis1.h"
# include "err.h"
# include "stisdef.h"

# define DEFAULT_BLEV_CLIP 50.

typedef struct {
	IRAFPointer tp;			/* pointer to table descriptor */
	IRAFPointer cp_amp;		/* column descriptors */
	IRAFPointer cp_ccdgain;
	IRAFPointer cp_ccdoffset;
	IRAFPointer cp_bin1;
	IRAFPointer cp_bin2;
	IRAFPointer cp_atodgain;
	IRAFPointer cp_bias;
	IRAFPointer cp_readnoise;
	IRAFPointer cp_saturate;
	IRAFPointer cp_blev_clip;
	IRAFPointer cp_pedigree;
	IRAFPointer cp_descrip;
	int nrows;			/* number of rows in table */
} TblInfo;

typedef struct {
	char ccdamp[STIS_CBUF+1];
	int ccdgain;
	int bin[2];
	int ccdoffset;
	float atodgain;
	float ccdbias;
	float readnoise;
	float saturate;
	float blev_clip;
} TblRow;

static int OpenCCDTab (char *, TblInfo *);
static int ReadCCDTab (TblInfo *, int, TblRow *);
static int CloseCCDTab (TblInfo *);

/* This routine gets the gain, bias and readnoise for the CCD from
   the CCD parameters table (keyword name CCDTAB).

   The CCD parameters table should contain the following:
	header parameters:
		none needed
	columns:
		CCDAMP:  identifies which amp was used (string A-D)
		CCDGAIN:  commanded gain of the CCD (int)
		CCDOFFST:  commanded bias (int)
		BINAXIS1, BINAXIS2:  commanded bin sizes (int)
		ATODGAIN:  actual gain of the CCD (double)
		CCDBIAS:  typical value of bias (double)
		READNSE:  typical value of readout noise (double)
		SATURATE:  CCD saturation threshold
		BLEV_CLIP:  sigma for clipping in virtual overscan region

   The table is read to find the row for which the CCDAMP matches the
   expected amplifier name (a single letter, from the image header) and
   for which the commanded CCD gain, commanded CCD bias, and bin sizes
   match the values read from the image header.  Then that row is read
   to get the actual gain, bias, readnoise, and saturation level.

   Phil Hodge, 2000 Jan 13:
	Add one to ccdamp buffer size.

   Phil Hodge, 2004 July 28:
	Also get BLEV_CLIP.
*/

int GetCCDTab (StisInfo1 *sts) {

/* arguments:
StisInfo1 *sts     io: calibration switches, etc
*/

	int status;

	TblInfo tabinfo;	/* pointer to table descriptor, etc */
	TblRow tabrow;		/* values read from a table row */

	int row;		/* loop index */
	int foundit;		/* has the correct row been found? */

	/* Open the CCD parameters table and find columns. */
	if ((status = OpenCCDTab (sts->ccdpar.name, &tabinfo)))
	    return (status);

	/* Check each row for a match with ccdamp, ccdgain, ccdoffst,
	   binaxis1, and binaxis2, and get info from the matching row.
	*/

	foundit = 0;
	for (row = 1;  row <= tabinfo.nrows;  row++) {

	    /* Read the current row into tabrow. */
	    if ((status = ReadCCDTab (&tabinfo, row, &tabrow)))
		return (status);

	    if (SameString (tabrow.ccdamp, sts->ccdamp) &&
		SameInt (tabrow.ccdgain, sts->ccdgain) &&
		SameInt (tabrow.ccdoffset, sts->ccdoffset) &&
		SameInt (tabrow.bin[0], sts->binaxis[0]) &&
		SameInt (tabrow.bin[1], sts->binaxis[1])) {

		foundit = 1;
		if ((status = RowPedigree (&sts->ccdpar, row,
                        tabinfo.tp, tabinfo.cp_pedigree, tabinfo.cp_descrip)))
		    return (status);
		if (sts->ccdpar.goodPedigree == DUMMY_PEDIGREE)
		    printf ("Warning  Row %d of CCDTAB is DUMMY.\n", row);
		sts->atodgain = tabrow.atodgain;
		sts->ccdbias = tabrow.ccdbias;
		sts->readnoise = tabrow.readnoise;
		sts->saturate = tabrow.saturate;
		sts->blev_clip = tabrow.blev_clip;
		break;
	    }
	}

	if (!foundit) {
	    printf ("ERROR    Matching row not found in CCDTAB `%s'.\n",
			sts->ccdpar.name);
	    printf (
		"ERROR    CCDAMP %s, CCDGAIN %d, CCDOFFST %d, BINAXIS %d,%d.\n",
		sts->ccdamp, sts->ccdgain, sts->ccdoffset,
		sts->binaxis[0], sts->binaxis[1]);
	    CloseCCDTab (&tabinfo);
	    return (TABLE_ERROR);
	}

	if ((status = CloseCCDTab (&tabinfo)))  /* close the table */
	    return (status);

	return (0);
}

/* This routine opens the CCD parameters table, finds the columns that we
   need, and gets the total number of rows in the table.  The columns are
   CCDAMP, CCDGAIN, CCDBIAS, and READNSE.
*/

static int OpenCCDTab (char *tname, TblInfo *tabinfo) {

	tabinfo->tp = c_tbtopn (tname, IRAF_READ_ONLY, 0);
	if (c_iraferr()) {
	    printf ("ERROR    CCDTAB `%s' not found.\n", tname);
	    return (OPEN_FAILED);
	}

	tabinfo->nrows = c_tbpsta (tabinfo->tp, TBL_NROWS);

	/* Find the columns. */
	c_tbcfnd1 (tabinfo->tp, "CCDAMP", &tabinfo->cp_amp);
	c_tbcfnd1 (tabinfo->tp, "CCDGAIN", &tabinfo->cp_ccdgain);
	c_tbcfnd1 (tabinfo->tp, "BINAXIS1", &tabinfo->cp_bin1);
	c_tbcfnd1 (tabinfo->tp, "BINAXIS2", &tabinfo->cp_bin2);
	c_tbcfnd1 (tabinfo->tp, "CCDOFFST", &tabinfo->cp_ccdoffset);
	c_tbcfnd1 (tabinfo->tp, "ATODGAIN", &tabinfo->cp_atodgain);
	c_tbcfnd1 (tabinfo->tp, "CCDBIAS", &tabinfo->cp_bias);
	c_tbcfnd1 (tabinfo->tp, "READNSE", &tabinfo->cp_readnoise);
	c_tbcfnd1 (tabinfo->tp, "SATURATE", &tabinfo->cp_saturate);
	c_tbcfnd1 (tabinfo->tp, "BLEV_CLIP", &tabinfo->cp_blev_clip);
	if (tabinfo->cp_amp == 0 ||
	    tabinfo->cp_ccdgain == 0 ||
	    tabinfo->cp_bin1 == 0 ||
	    tabinfo->cp_bin2 == 0 ||
	    tabinfo->cp_ccdoffset == 0 ||
	    tabinfo->cp_atodgain == 0 ||
	    tabinfo->cp_bias == 0 ||
	    tabinfo->cp_readnoise == 0 ||
	    tabinfo->cp_saturate == 0) {
	    printf ("ERROR    Column not found in CCDTAB.\n");
	    c_tbtclo (tabinfo->tp);
	    return (COLUMN_NOT_FOUND);
	}

	/* BLEV_CLIP is an optional column. */

	/* Pedigree and descrip are optional columns. */
	c_tbcfnd1 (tabinfo->tp, "PEDIGREE", &tabinfo->cp_pedigree);
	c_tbcfnd1 (tabinfo->tp, "DESCRIP", &tabinfo->cp_descrip);

	return (0);
}

/* This routine reads the relevant data from one row.  The amplifier
   number, CCD gain, bias, and readnoise are gotten.
*/

static int ReadCCDTab (TblInfo *tabinfo, int row, TblRow *tabrow) {

	c_tbegtt (tabinfo->tp, tabinfo->cp_amp, row,
			tabrow->ccdamp, STIS_CBUF);
	if (c_iraferr())
	    return (TABLE_ERROR);

	c_tbegti (tabinfo->tp, tabinfo->cp_ccdgain, row, &tabrow->ccdgain);
	if (c_iraferr())
	    return (TABLE_ERROR);

	c_tbegti (tabinfo->tp, tabinfo->cp_ccdoffset, row, &tabrow->ccdoffset);
	if (c_iraferr())
	    return (TABLE_ERROR);

	c_tbegti (tabinfo->tp, tabinfo->cp_bin1, row, &tabrow->bin[0]);
	if (c_iraferr())
	    return (TABLE_ERROR);
	c_tbegti (tabinfo->tp, tabinfo->cp_bin2, row, &tabrow->bin[1]);
	if (c_iraferr())
	    return (TABLE_ERROR);

	c_tbegtr (tabinfo->tp, tabinfo->cp_atodgain, row, &tabrow->atodgain);
	if (c_iraferr())
	    return (TABLE_ERROR);

	c_tbegtr (tabinfo->tp, tabinfo->cp_bias, row, &tabrow->ccdbias);
	if (c_iraferr())
	    return (TABLE_ERROR);

	c_tbegtr (tabinfo->tp, tabinfo->cp_readnoise, row, &tabrow->readnoise);
	if (c_iraferr())
	    return (TABLE_ERROR);

	c_tbegtr (tabinfo->tp, tabinfo->cp_saturate, row, &tabrow->saturate);
	if (c_iraferr())
	    return (TABLE_ERROR);

	if (tabinfo->cp_blev_clip == 0) {
	    tabrow->blev_clip = DEFAULT_BLEV_CLIP;
	} else {
	    c_tbegtr (tabinfo->tp, tabinfo->cp_blev_clip, row,
			&tabrow->blev_clip);
	    if (c_iraferr())
		return (TABLE_ERROR);
	}

	return (0);
}

/* This routine closes the ccdtab table. */

static int CloseCCDTab (TblInfo *tabinfo) {

	c_tbtclo (tabinfo->tp);
	if (c_iraferr())
	    return (TABLE_ERROR);

	return (0);
}
