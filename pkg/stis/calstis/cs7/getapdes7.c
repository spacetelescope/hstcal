# include <stdio.h>

# include "c_iraf.h"
# include "hstio.h"
# include "xtables.h"

# include "stis.h"
# include "calstis7.h"
# include "err.h"
# include "stisdef.h"

typedef struct {
	IRAFPointer tp;			/* pointer to table descriptor */
	IRAFPointer cp_aperture;	/* column descriptors */
	IRAFPointer cp_width[2];
        IRAFPointer cp_angle;
	IRAFPointer cp_offset[2];
	IRAFPointer cp_nbars;
	IRAFPointer cp_barlocn[3];
	IRAFPointer cp_barwidth[3];
	IRAFPointer cp_pedigree;
	IRAFPointer cp_descrip;
	int nrows;			/* number of rows in table */
} TblInfo;

typedef struct {
	char aperture[STIS_CBUF+1];
} TblRow;

static int OpenApTab (char *, TblInfo *);
static int ReadApTab (TblInfo *, int, TblRow *);
static int ReadApArray (TblInfo *, int, ApInfo *);
static int CloseApTab (TblInfo *);

/* This routine gets aperture information from APDESTAB (_apd).

   The aperture description table should contain the following:
	header parameters:
		none needed
	columns:
		APERTURE:  aperture name (string)
		WIDTH1:  width along X axis (float)
		WIDTH2:  width along Y axis (float)
                ANGLE:   tilt angle of aperture
		OFFSET1, OFFSET2:  offset from nominal position (float)
		NBARS:  number of occulting bars (int)
		BAR1LOCN, BAR2LOCN, BAR3LOCN:  arcsec from end (float)
		BAR1WIDTH, BAR2WIDTH, BAR3WIDTH:  widths of bars (float)

   The table is read to find the row for which the value of APERTURE
   is the same as in the input image header.  For that row, the
   information about the slit is read in.  The width, offset, bar
   locations and bar widths are all in arcseconds.

   Phil Hodge, 1998 Mar 17:
	Remove the section for getting REF_APER from the table header.

   Phil Hodge, 2000 Jan 13:
	Add one to aperture buffer size.

   Paul Barrett, 2004 Feb 11:
        Add code to read tilt angle from new column of aperture table.
*/

int GetApDes7 (StisInfo7 *sts, ApInfo *slit) {

/* arguments:
StisInfo7 *sts   i: calibration switches and info
ApInfo *slit     o: description of slit
*/

	int status;

	TblInfo tabinfo;	/* pointer to table descriptor, etc */
	TblRow tabrow;		/* values read from a table row */

	int row;		/* number of rows, and loop index */
	int k;			/* loop index for occulting bars */
	int foundit = 0;	/* true if aperture found in table */

	/* Open the aperture description table. */
	if ((status = OpenApTab (sts->apdestab.name, &tabinfo)))
	    return (status);

	for (k = 0;  k < MAX_BARS;  k++) {	/* initial values */
	    slit->barlocn[k] = 0.;
	    slit->barwidth[k] = 0.;
	}

	/* Check each row for a match with aperture. */

	for (row = 1;  row <= tabinfo.nrows;  row++) {

	    if ((status = ReadApTab (&tabinfo, row, &tabrow)))
		return (status);

	    if (SameString (tabrow.aperture, sts->aperture)) {

		foundit = 1;

		/* Get pedigree & descrip from the row. */
		if ((status = RowPedigree (&sts->apdestab, row,
                        tabinfo.tp, tabinfo.cp_pedigree, tabinfo.cp_descrip)))
		    return (status);
		if (sts->apdestab.goodPedigree == DUMMY_PEDIGREE) {
		    printf ("Warning  DUMMY pedigree in row %d of %s.\n",
			row, sts->apdestab.name);
		    sts->x2dcorr_o = DUMMY;
		    CloseApTab (&tabinfo);
		    return (0);
		}

		/* Read aperture info into slit structure. */
		if ((status = ReadApArray (&tabinfo, row, slit)))
		    return (status);

		break;
	    }
	}

	if ((status = CloseApTab (&tabinfo)))
	    return (status);

	if (!foundit) {
	    printf ("Warning  APERTURE %s not found in APDESTAB %s.\n",
		sts->aperture, sts->apdestab.name);
	    sts->x2dcorr_o = OMIT;
	}

	return (0);
}

/* This routine opens the aperture description table, finds the columns
   that we need, and gets the total number of rows in the table.
*/

static int OpenApTab (char *tname, TblInfo *tabinfo) {

	tabinfo->tp = c_tbtopn (tname, IRAF_READ_ONLY, 0);
	if (c_iraferr()) {
	    printf ("ERROR    APDESTAB `%s' not found.\n", tname);
	    return (OPEN_FAILED);
	}

	tabinfo->nrows = c_tbpsta (tabinfo->tp, TBL_NROWS);

	/* Find the columns. */
	c_tbcfnd1 (tabinfo->tp, "APERTURE", &tabinfo->cp_aperture);
	c_tbcfnd1 (tabinfo->tp, "WIDTH1", &tabinfo->cp_width[0]);
	c_tbcfnd1 (tabinfo->tp, "WIDTH2", &tabinfo->cp_width[1]);
	c_tbcfnd1 (tabinfo->tp, "OFFSET1", &tabinfo->cp_offset[0]);
	c_tbcfnd1 (tabinfo->tp, "OFFSET2", &tabinfo->cp_offset[1]);
	c_tbcfnd1 (tabinfo->tp, "NBARS", &tabinfo->cp_nbars);
	c_tbcfnd1 (tabinfo->tp, "BAR1LOCN", &tabinfo->cp_barlocn[0]);
	c_tbcfnd1 (tabinfo->tp, "BAR1WIDTH", &tabinfo->cp_barwidth[0]);
	c_tbcfnd1 (tabinfo->tp, "BAR2LOCN", &tabinfo->cp_barlocn[1]);
	c_tbcfnd1 (tabinfo->tp, "BAR2WIDTH", &tabinfo->cp_barwidth[1]);
	c_tbcfnd1 (tabinfo->tp, "BAR3LOCN", &tabinfo->cp_barlocn[2]);
	c_tbcfnd1 (tabinfo->tp, "BAR3WIDTH", &tabinfo->cp_barwidth[2]);
	if (tabinfo->cp_aperture == 0 ||
	    tabinfo->cp_width[0] == 0 || tabinfo->cp_width[1] == 0 ||
	    tabinfo->cp_offset[0] == 0 || tabinfo->cp_offset[1] == 0 ||
	    tabinfo->cp_nbars == 0 ||
	    tabinfo->cp_barlocn[0] == 0 || tabinfo->cp_barwidth[0] == 0 ||
	    tabinfo->cp_barlocn[1] == 0 || tabinfo->cp_barwidth[1] == 0 ||
	    tabinfo->cp_barlocn[2] == 0 || tabinfo->cp_barwidth[2] == 0) {
	    printf ("ERROR    Column not found in APDESTAB.\n");
	    c_tbtclo (tabinfo->tp);
	    return (COLUMN_NOT_FOUND);
	}

	/* Pedigree, descrip, and angle are optional columns. */
	c_tbcfnd1 (tabinfo->tp, "PEDIGREE", &tabinfo->cp_pedigree);
	c_tbcfnd1 (tabinfo->tp, "DESCRIP", &tabinfo->cp_descrip);
	c_tbcfnd1 (tabinfo->tp, "ANGLE",  &tabinfo->cp_angle);

	return (0);
}

/* This routine reads the column (aperture name) used to select the
   correct row.
*/

static int ReadApTab (TblInfo *tabinfo, int row, TblRow *tabrow) {

	c_tbegtt (tabinfo->tp, tabinfo->cp_aperture, row,
			tabrow->aperture, STIS_CBUF);
	if (c_iraferr())
	    return (TABLE_ERROR);

	return (0);
}

/* Read info from the current row.
   Not all the values we need are arrays; the name ReadApArray was chosed
   just for consistency with other similar routines.
*/

static int ReadApArray (TblInfo *tabinfo, int row, ApInfo *slit) {

        if (tabinfo->cp_angle == 0) {
            slit->angle = REF_ANGLE;
        } else {
            c_tbegtd (tabinfo->tp, tabinfo->cp_angle   , row, &slit->angle);
        }
	c_tbegtd (tabinfo->tp, tabinfo->cp_width[0], row, &slit->width[0]);
	c_tbegtd (tabinfo->tp, tabinfo->cp_width[1], row, &slit->width[1]);
	c_tbegtd (tabinfo->tp, tabinfo->cp_offset[0], row, &slit->ap_offset[0]);
	c_tbegtd (tabinfo->tp, tabinfo->cp_offset[1], row, &slit->ap_offset[1]);
	c_tbegti (tabinfo->tp, tabinfo->cp_nbars, row, &slit->nbars);
	if (slit->nbars > 0) {
	    c_tbegtd (tabinfo->tp, tabinfo->cp_barlocn[0], row,
			&slit->barlocn[0]);
	    c_tbegtd (tabinfo->tp, tabinfo->cp_barwidth[0], row,
			&slit->barwidth[0]);
	    if (slit->nbars > 1) {
		c_tbegtd (tabinfo->tp, tabinfo->cp_barlocn[1], row,
			    &slit->barlocn[1]);
		c_tbegtd (tabinfo->tp, tabinfo->cp_barwidth[1], row,
			    &slit->barwidth[1]);
		if (slit->nbars > 2) {
		    c_tbegtd (tabinfo->tp, tabinfo->cp_barlocn[2], row,
				&slit->barlocn[2]);
		    c_tbegtd (tabinfo->tp, tabinfo->cp_barwidth[2], row,
				&slit->barwidth[2]);
		}
	    }
	}
	if (c_iraferr())
	    return (TABLE_ERROR);

	return (0);
}

/* This routine closes the apdestab table. */

static int CloseApTab (TblInfo *tabinfo) {

	c_tbtclo (tabinfo->tp);
	if (c_iraferr())
	    return (TABLE_ERROR);

	return (0);
}
