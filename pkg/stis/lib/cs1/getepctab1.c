# include <stdio.h>
# include <math.h>

# include "c_iraf.h"
# include "hstio.h"
# include "xtables.h"
# include "stis.h"
# include "calstis1.h"
# include "hstcalerr.h"
# include "stisdef.h"

typedef struct {
    IRAFPointer tp;		/* pointer to table descriptor */
    IRAFPointer cp_mjd;		/* column descriptors */
    IRAFPointer cp_temp;
    IRAFPointer cp_pedigree;
    IRAFPointer cp_descrip;
    int nrows;			/* number of rows in table */
} TblInfo;

typedef struct {
    double mjd;
    double temp;
} TblRow;

static int OpenEPCTab(char *, TblInfo *);
static int ReadEPCTab(TblInfo *, int, TblRow *);
static int CloseEPCTab(TblInfo *);

/*
  This routine gets the mjd and temperature from the EPC parameters
  calibration table (keyword name EPCTAB).

  The EPC parameters table should contain the following:
	header parameters:
		none needed
	columns:
		TIME:     gives the Modified Julian Day of the value.
		VALUE:    gives the temperature of the CCD housing.

  All rows are read from the table. The data are filtered for
  discrepant temperature values using a simple filter algorithm.  The
  MJD and temperature arrays are then saved for later use to
  accurately determine the time-averaged temperature of the side 2
  CCD during dark subtraction.

  Paul Barrett, 2003 Sep 18:
	Copied from cs1/getcddtab.c and modified.

  Paul Barrett, 2004 Jun 25:
        Modified for engineering EPC table data.

*/

int GetEPCTab (StisInfo1 *sts, float sigma) {

/* arguments:
StisInfo1 *sts     io: calibration switches, etc
*/

    int status;

    TblInfo tabinfo;	/* pointer to table descriptor, etc */
    TblRow  tabrow;	/* values read from a table row */

    int row;		/* loop index */
    double *mjd;        /* MJD field of table */
    double *temp;       /* temperature field of table */

    sts->epc_rows = 0;

    if (sts->epctab.name == NULL)
        return (0);

    /*  Open the EPC parameters table and find columns. */
    if ((status = OpenEPCTab(sts->epctab.name, &tabinfo)))
        return (status);

    /*  Alloc memory  */

    mjd  = (double *) malloc(tabinfo.nrows * sizeof(double));
    temp = (double *) malloc(tabinfo.nrows * sizeof(double));
    if (mjd == NULL || temp == NULL)
        return (status = OUT_OF_MEMORY);

    sts->epc_mjd  = (double *) malloc(tabinfo.nrows * sizeof(double));
    sts->epc_temp = (double *) malloc(tabinfo.nrows * sizeof(double));
    if (sts->epc_mjd == NULL || sts->epc_temp == NULL)
        return (status = OUT_OF_MEMORY);

    /*  Read each row.  */

    for (row = 1; row <= tabinfo.nrows; row++) {

        /* Read the current row into tabrow. */
        if ((status = ReadEPCTab(&tabinfo, row, &tabrow)))
            return (status);

        mjd[row-1]  = tabrow.mjd;
        temp[row-1] = tabrow.temp;
    }
    /*
     *  Simple filtering procedure for discrepant temperatures.
     */
    for (row = 0; row < tabinfo.nrows; row++) {
        double d1, d2, p1, p2;
        /*  Calculate difference between adjacent temperatures. */
        d1 = temp[row] - ((row < 1) ? temp[row]: temp[row-1]);
        d2 = temp[row] - ((row > tabinfo.nrows-2) ? temp[row]: temp[row+1]);
        p1 = exp(-0.5*d1*d1/sigma/sigma);
        p2 = exp(-0.5*d2*d2/sigma/sigma);
        /*  Check for invalid temperature range or discrepant
         *  temperature readings.
         */
        if (!((temp[row] < 15. || temp[row] > 25.) ||
              (row == 0 && p2 < 0.05) ||
              (row == tabinfo.nrows-1 && p1 < 0.05) ||
              (p1 < 0.05 && p2 < 0.05))) {
            sts->epc_mjd[sts->epc_rows]  = mjd[row];
            sts->epc_temp[sts->epc_rows] = temp[row];
            sts->epc_rows++;
        }
    }

    if ((status = CloseEPCTab(&tabinfo))) /* close the table */
        return (status);

    free(temp);
    free(mjd);

    return (0);
}

/* This routine opens the EPC parameters table, finds the columns that we
   need, and gets the total number of rows in the table.  The columns are
   TIME and VALUE.
*/

static int OpenEPCTab(char *tname, TblInfo *tabinfo) {

    tabinfo->tp = c_tbtopn(tname, IRAF_READ_ONLY, 0);
    if (c_iraferr()) {
        printf("Warning  EPCTAB `%s' not found.\n", tname);
        return (OPEN_FAILED);
    }

    tabinfo->nrows = c_tbpsta(tabinfo->tp, TBL_NROWS);

    /* Find the columns. */
    c_tbcfnd1(tabinfo->tp, "TIME",  &tabinfo->cp_mjd);
    c_tbcfnd1(tabinfo->tp, "VALUE", &tabinfo->cp_temp);

    if (tabinfo->cp_mjd == 0 || tabinfo->cp_temp == 0) {
        printf ("Warning  Column not found in EPCTAB.\n");
        c_tbtclo(tabinfo->tp);
        return (COLUMN_NOT_FOUND);
    }

    return (0);
}

/* This routine reads the relevant data (TIME and VALUE) from one row.
 */

static int ReadEPCTab(TblInfo *tabinfo, int row, TblRow *tabrow) {

    c_tbegtd(tabinfo->tp, tabinfo->cp_mjd, row, &tabrow->mjd);
    if (c_iraferr())
        return (TABLE_ERROR);

    c_tbegtd(tabinfo->tp, tabinfo->cp_temp, row, &tabrow->temp);
    if (c_iraferr())
        return (TABLE_ERROR);

    return (0);
}

/* This routine closes the epctab table. */

static int CloseEPCTab(TblInfo *tabinfo) {

    c_tbtclo (tabinfo->tp);
    if (c_iraferr())
        return (TABLE_ERROR);

    return (0);
}
