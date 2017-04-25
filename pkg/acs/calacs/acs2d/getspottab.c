# include <stdio.h>

# include "hstio.h"
# include "xtables.h"
# include <time.h>

# include "acs.h"
# include "hstcalerr.h"

typedef struct {
	IRAFPointer tp;			/* pointer to table descriptor */
	IRAFPointer cp_date;		/* column descriptors */
	IRAFPointer cp_shiftx;
	IRAFPointer cp_shifty;
	int nrows;			/* number of rows in table */
} TblInfo;

typedef struct {
	char date[ACS_CBUF];	/* date for shift */
    time_t dtime;
	float shiftx;
	float shifty;
} TblRow;

static int OpenSpotTab (char *, TblInfo *);
static int ReadSpotTab (TblInfo *, int, TblRow *);
static int CloseSpotTab (TblInfo *);

/* This routine gets information from the SPOTFLAT shift table.

   The SPOTFLAT shift table should contain the following:
   floating point header parameter:
		none needed
   and the following columns:
		DATE:  (char *) date in form of 'dd/mm/yy'
        SHIFTX: (float) shift of spots in X direction (pixels)
        SHIFTY: (float) shift of spots in Y direction (pixels)
        
   The table is read to find the row that comes closest to the observation
   date, and that row is read to obtain the shifts in X and Y.
   
   Input date will be in form of seconds since 1970 as produced by mktime().
*/

int GetSpotTab (char *spottab, time_t date, float *shiftx, float *shifty) {

/* arguments:
ACS2dInfo *acs2d     io: calibration switches, etc
*/

	extern int status;

	TblInfo tabinfo;	/* pointer to table descriptor, etc */
	TblRow tabrow;		/* values read from a table row */

	int row;			/* loop index */
    int min_row;
    time_t delta_date, delta;

	/* Open the SPOT shift table and find columns. */
	if (OpenSpotTab (spottab, &tabinfo))
	    return (status);

	/* Check each row for a match with detector, and get the info
	   from the matching row.
	*/
    delta_date = 9999999999;
    min_row = tabinfo.nrows;
	for (row = 1;  row <= tabinfo.nrows;  row++) {

	    /* Read the current row into tabrow. */
	    if (ReadSpotTab (&tabinfo, row, &tabrow))
		return (status);
        
        delta = abs(date - tabrow.dtime);
        
	    if (delta < delta_date) {
            min_row = row;
            delta_date = delta;
        }            
    }
    
    /* Now that we have determined which row has the date closest
        to the observation date, read in that row as the final result.
    */  		
    ReadSpotTab (&tabinfo, min_row, &tabrow);
    
    *shiftx = tabrow.shiftx;
    *shifty = tabrow.shifty;

	if (CloseSpotTab (&tabinfo))		/* close the table */
	    return (status);

	return (status);
}

/* This routine opens the MAMA linearity table, finds the columns that we
   need, gets one header keyword, and gets the total number of rows in the
   table.  The columns are DETECTOR, GLOBAL_LIMIT, LOCAL_LIMIT, and TAU.
*/

static int OpenSpotTab (char *tname, TblInfo *tabinfo) {

	extern int status;

	tabinfo->tp = c_tbtopn (tname, IRAF_READ_ONLY, 0);
	if (c_iraferr()) {
	    sprintf (MsgText, "SPOTTAB `%s' not found.", tname);
	    trlerror (MsgText);
		return (status = OPEN_FAILED);
	}

	tabinfo->nrows = c_tbpsta (tabinfo->tp, TBL_NROWS);
    
	/* Find the columns. */
	c_tbcfnd1 (tabinfo->tp, "DATE", &tabinfo->cp_date);
	c_tbcfnd1 (tabinfo->tp, "SHIFTX", &tabinfo->cp_shiftx);
	c_tbcfnd1 (tabinfo->tp, "SHIFTY", &tabinfo->cp_shifty);
	if (tabinfo->cp_date == 0 ||
	    tabinfo->cp_shiftx == 0 ||
	    tabinfo->cp_shifty == 0 ) {
	    trlerror ("Column not found in SPOTTAB.");
	    c_tbtclo (tabinfo->tp);
	    return (status = COLUMN_NOT_FOUND);
	}

	return (status);
}

/* This routine reads the relevant data from one row.  The date,
   shiftx, shifty are read.
*/

static int ReadSpotTab (TblInfo *tabinfo, int row, TblRow *tabrow) {

	extern int status;
    int parseTabDate(char *, time_t *);

	c_tbegtt (tabinfo->tp, tabinfo->cp_date, row,
			tabrow->date, ACS_CBUF-1);
	if (c_iraferr())
	    return (status = TABLE_ERROR);

    /* Convert the date read in from the row into
        time in seconds for comparison with the DATE-OBS value */
    status = parseTabDate(tabrow->date, &tabrow->dtime);
	c_tbegtr (tabinfo->tp, tabinfo->cp_shiftx, row, &tabrow->shiftx);
	if (c_iraferr())
	    return (status = TABLE_ERROR);

	c_tbegtr (tabinfo->tp, tabinfo->cp_shifty, row, &tabrow->shifty);
	if (c_iraferr())
	    return (status = TABLE_ERROR);

	return (status);
}

/* This routine closes the mlintab table. */

static int CloseSpotTab (TblInfo *tabinfo) {

	extern int status;

	c_tbtclo (tabinfo->tp);
	if (c_iraferr())
	    return (status = TABLE_ERROR);

	return (status);
}
