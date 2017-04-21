# include <stdio.h>

# include "xtables.h"
# include "hstio.h"

# include "stis.h"
# include "stisdef.h"
# include "calstis6.h"
# include "err.h"


typedef struct {
	IRAFPointer tp;			/* pointer to table descriptor */
	IRAFPointer cp_aperture;	/* column descriptors */
	IRAFPointer cp_offset[2];
	IRAFPointer cp_width[2];
	IRAFPointer cp_pedigree;
	IRAFPointer cp_descrip;
	int nrows;			/* number of rows in table */
} TblInfo;

typedef struct {
	char aperture[STIS_CBUF];
} TblRow;


static int OpenApTab (char *, TblInfo *);
static int ReadApTab (TblInfo *, int, TblRow *);
static int ReadApArray (TblInfo *, int, ApInfo *);
static int CloseApTab (TblInfo *);


/* This routine gets aperture information from APDESTAB (_apd).

   The aperture description table should contain the following:
		APERTURE:  aperture name (string)
		OFFSET1, OFFSET2:  offset from nominal position (float)
		WIDTH1:  width along X axis (float)
		WIDTH2:  width along Y axis (float)

   The table is read to find the row for which the value of APERTURE
   is the same as in the input image header. For that row, the
   information about the slit is read in. The offsets and widths are
   in arcseconds.

   WIDTH1 is used only to find the projected slit width when
   computing the geocoronal Lya avoidance rgion.




   Revision history:
   ----------------
   20 Feb 97  -  Borrowed from calstis7 (I.Busko)
   20 Feb 97  -  Removed input from several columns (IB)
   24 Feb 97  -  Rename routine to avoid conflict with cs7 (IB)
   10 Apr 97  -  Changes after code review (IB):
                 - replaced literal by ROW_NOT_FOUND constant
   08 May 97  -  Conform to new _trl standard (IB)
   17 Dec 98  -  Read WIDTH values (IB)
   23 Feb 00  -  Write WIDTH values into StisInfo6 structure (IB)
*/

int GetApDes6 (StisInfo6 *sts, ApInfo *slit) {

/* arguments:
StisInfo6 *sts   i: calibration switches and info
ApInfo *slit     o: description of slit
*/

	int status;

	TblInfo tabinfo;	/* pointer to table descriptor, etc */
	TblRow tabrow;		/* values read from a table row */

	int row;		/* number of rows, and loop index */
	int foundit = 0;	/* true if aperture found in table */

	/* Open the aperture description table. */
	if ((status = OpenApTab (sts->apdestab.name, &tabinfo)))
	    return (status);

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
		    sts->dispcorr = DUMMY;
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
	    printf ("ERROR    APERTURE %s not found in APDESTAB %s\n",
		sts->aperture, sts->apdestab.name);
	    return (ROW_NOT_FOUND);
	}

	/* used by the IDT algorithm */
	sts->ap_xsize = slit->width[0];
	sts->ap_ysize = slit->width[1];

	return (0);
}



/* This routine opens the aperture description table, finds the columns
   that we need, and gets the total number of rows in the table.
*/

static int OpenApTab (char *tname, TblInfo *tabinfo) {

	tabinfo->tp = c_tbtopn (tname, IRAF_READ_ONLY, 0);
	if (c_iraferr()) {
	    printf ("ERROR    APDESTAB `%s' not found\n", tname);
	    return (OPEN_FAILED);
	}

	tabinfo->nrows = c_tbpsta (tabinfo->tp, TBL_NROWS);

	/* Find the columns. */
	c_tbcfnd1 (tabinfo->tp, "APERTURE", &tabinfo->cp_aperture);
	c_tbcfnd1 (tabinfo->tp, "OFFSET1",  &tabinfo->cp_offset[0]);
	c_tbcfnd1 (tabinfo->tp, "OFFSET2",  &tabinfo->cp_offset[1]);
	c_tbcfnd1 (tabinfo->tp, "WIDTH1", &tabinfo->cp_width[0]);
	c_tbcfnd1 (tabinfo->tp, "WIDTH2", &tabinfo->cp_width[1]);
	if (tabinfo->cp_aperture  == 0 ||
	    tabinfo->cp_width[0]  == 0 ||
            tabinfo->cp_width[1]  == 0 ||
	    tabinfo->cp_offset[0] == 0 ||
            tabinfo->cp_offset[1] == 0) {
	    printf ("ERROR    Column not found in APDESTAB\n");
	    c_tbtclo (tabinfo->tp);
	    return (COLUMN_NOT_FOUND);
	}

	/* Pedigree and descrip are optional columns. */
	c_tbcfnd1 (tabinfo->tp, "PEDIGREE", &tabinfo->cp_pedigree);
	c_tbcfnd1 (tabinfo->tp, "DESCRIP", &tabinfo->cp_descrip);

	return (0);
}



/* This routine reads the column (aperture name) used to select the
   correct row.
*/

static int ReadApTab (TblInfo *tabinfo, int row, TblRow *tabrow) {

	c_tbegtt (tabinfo->tp, tabinfo->cp_aperture, row,
			tabrow->aperture, STIS_CBUF-1);
	if (c_iraferr())
	    return (TABLE_ERROR);

	return (0);
}



/* Read info from the current row. */

static int ReadApArray (TblInfo *tabinfo, int row, ApInfo *slit) {

	c_tbegtr (tabinfo->tp, tabinfo->cp_offset[0], row, &slit->ap_offset[0]);
	c_tbegtr (tabinfo->tp, tabinfo->cp_offset[1], row, &slit->ap_offset[1]);
	c_tbegtd (tabinfo->tp, tabinfo->cp_width[0], row, &slit->width[0]);
	c_tbegtd (tabinfo->tp, tabinfo->cp_width[1], row, &slit->width[1]);
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
