# include <stdio.h>
# include <stdlib.h>	/* malloc */

# include "xtables.h"
# include "hstio.h"

# include "stis.h"
# include "stisdef.h"
# include "calstis6.h"
# include "hstcalerr.h"


typedef struct {
	IRAFPointer tp;			/* pointer to table descriptor */
	IRAFPointer cp_sporder;		/* column descriptors */

	IRAFPointer cp_npts;
	IRAFPointer cp_nptsoff;
	IRAFPointer cp_subscale;
	IRAFPointer cp_minw;
	IRAFPointer cp_maxw;
	IRAFPointer cp_minp;
	IRAFPointer cp_maxp;
	IRAFPointer cp_sn;
	IRAFPointer cp_profoff;
	IRAFPointer cp_prof;

	IRAFPointer cp_pedigree;
	IRAFPointer cp_descrip;
	int nrows;			/* number of rows in table */
} TblInfo;

typedef struct {
	int sporder;
} TblRow;


static int OpenProfileTab (char *, TblInfo *);
static int ReadProfileTab (TblInfo *, int, TblRow *);
static int ReadProfileArray (TblInfo *, int, ProfileArray **, double*);
static int CloseProfileTab (TblInfo *);


/* This routine gets profile information from OPROFTAB and saves the
   info in the profile information structure. This information is used
   by the optimal extraction algorithm.

   This routine supports the profile file format curently created by
   calstis6 when running in "profile generator" mode. This format
   assumes that each table row stores a cross-dispersion profile
   associated with a single spectral order and a single wavelength
   interval inside that order.

   The intensity info table should contain the following:
	header parameters:
		none needed
	columns:
	  SPORDER:     spectral order number (short)
	  NPTS:        number of elements in profile array (short)
	  NPTS_OFFSET: number of elements in offsets array (short)
	  SUB_FACTOR:  subsampling factor (double)
	  MIN_WAVE:    minimum wavelength for this profile (double)
	  MAX_WAVE:    maximum wavelength for this profile (double)
	  MIN_PIX:     minimum pixel for this profile (short)
	  MAX_PIX:     maximum pixel for this profile (short)
	  S_N:         signal-to-noise (float)
	  PROF:        subsampled profile values (double)
	  PROF_OFFSET: offset values for current bin (double)

   The table is read to find all rows for which the values of SPORDER are
   the same as in the function parameter. There can be several such rows,
   each with a profile corresponding to a different wavelength/pixel
   range. All these rows are read into memory, pointed to by ProfileArray.
   The SPTRCTAB table need not be sorted.

   When done, memory should be freed by calling FreeProfileArray.

   For now, the same subsampling factor value is assumed to be stored
   in all rows of the profile table.



   Revision history:
   ----------------
   17 Sep 98  -  Implemented (I.Busko)
   04 Dec 00  -  Subsampled profiles (IB)
*/

int GetProfile (StisInfo6 *sts, int sporder, ProfileArray **profa) {

/* arguments:
StisInfo6 *sts        i: calibration switches and info
int sporder;          i: spectral order number
ProfileArray *profa;  o: list with profile arrays for current sporder
*/

	int status;

	TblInfo tabinfo;	/* pointer to table descriptor, etc */
	TblRow tabrow;		/* values read from a table row */

	int row;		/* loop index */

	int CheckProfile (ProfileArray **);
	void FreeProfileArray (ProfileArray **);

	/* Open the profile table. */
	if ((status = OpenProfileTab (sts->pftab.name, &tabinfo)))
	    return (status);

	for (row = 1;  row <= tabinfo.nrows;  row++) {

	    if ((status = ReadProfileTab (&tabinfo, row, &tabrow)))
		return (status);

	    /* Check for matching sporder. */
	    if (SameInt (tabrow.sporder, sporder)) {

	        /* Get pedigree & descrip from the row. */
	        if ((status = RowPedigree (&sts->pftab, row,
                        tabinfo.tp, tabinfo.cp_pedigree, tabinfo.cp_descrip)))
	            return (status);

	        /* Read profiles into structure. */
	        if ((status = ReadProfileArray (&tabinfo, row, profa,
                                                &(sts->subscale))))
	            return (status);
	    }
	}

	/* Sanity check. */
	if ((status = CheckProfile (profa))) {
	    FreeProfileArray (profa);
	    if (status < 0) {
		printf ("ERROR    Matching row not found in OPROFAB %s\n",
				sts->pftab.name);
		printf ("ERROR    SPORDER %d\n", sporder);
		return (ROW_NOT_FOUND);
	    } else {
		return (status);
	    }
	}

	if ((status = CloseProfileTab (&tabinfo)))
	    return (status);

	return (0);
}



/* This routine opens the profile table table, finds the columns
   that we need, and gets the total number of rows in the table.
*/

static int OpenProfileTab (char *tname, TblInfo *tabinfo) {

	tabinfo->tp = c_tbtopn (tname, IRAF_READ_ONLY, 0);
	if (c_iraferr()) {
	    printf ("ERROR    OPROFTAB `%s' not found\n", tname);
	    return (OPEN_FAILED);
	}

	tabinfo->nrows = c_tbpsta (tabinfo->tp, TBL_NROWS);

	/* Find the columns. */
	c_tbcfnd1 (tabinfo->tp, "SPORDER", &tabinfo->cp_sporder);
	c_tbcfnd1 (tabinfo->tp, "NPTS", &tabinfo->cp_npts);
	c_tbcfnd1 (tabinfo->tp, "NPTS_OFFSET", &tabinfo->cp_nptsoff);
	c_tbcfnd1 (tabinfo->tp, "SUB_FACTOR", &tabinfo->cp_subscale);
	c_tbcfnd1 (tabinfo->tp, "MIN_WAVE", &tabinfo->cp_minw);
	c_tbcfnd1 (tabinfo->tp, "MAX_WAVE", &tabinfo->cp_maxw);
	c_tbcfnd1 (tabinfo->tp, "MIN_PIX", &tabinfo->cp_minp);
	c_tbcfnd1 (tabinfo->tp, "MAX_PIX", &tabinfo->cp_maxp);
	c_tbcfnd1 (tabinfo->tp, "S_N", &tabinfo->cp_sn);
	c_tbcfnd1 (tabinfo->tp, "PROF_OFFSET", &tabinfo->cp_profoff);
	c_tbcfnd1 (tabinfo->tp, "PROF", &tabinfo->cp_prof);
	if (tabinfo->cp_sporder == 0 || tabinfo->cp_npts == 0 ||
	    tabinfo->cp_minw    == 0 || tabinfo->cp_maxw == 0 ||
	    tabinfo->cp_minp    == 0 || tabinfo->cp_maxp == 0 ||
	    tabinfo->cp_profoff == 0 || tabinfo->cp_prof == 0 ||
	    tabinfo->cp_sn      == 0 || tabinfo->cp_nptsoff == 0) {
	    c_tbtclo (tabinfo->tp);
	    printf ("ERROR    Column not found in OPROFTAB\n");
	    return (COLUMN_NOT_FOUND);
	}

	/* Pedigree and descrip are optional columns. */
	c_tbcfnd1 (tabinfo->tp, "PEDIGREE", &tabinfo->cp_pedigree);
	c_tbcfnd1 (tabinfo->tp, "DESCRIP", &tabinfo->cp_descrip);

	return (0);
}


/* This routine reads the SPORDER column used to select the sought rows. */

static int ReadProfileTab (TblInfo *tabinfo, int row, TblRow *tabrow) {

	c_tbegti (tabinfo->tp, tabinfo->cp_sporder, row, &tabrow->sporder);
	if (c_iraferr())
	    return (TABLE_ERROR);

	return (0);
}



/* This routine reads the data from one row into the profile structure.
   Several scalar column values and two arrays are gotten. Memory for
   these arrays is allocated; FreeProfileArray takes care of freeing it.
*/

static int ReadProfileArray (TblInfo *tabinfo, int row, ProfileArray **profa,
                             double *subscale) {
	int status;

	int npts1, npts2;
	ProfileArray *newp;
	int NewProfile (ProfileArray **, ProfileArray *);

	if ((newp = (ProfileArray *) malloc (sizeof (ProfileArray))) == NULL) {
	    printf ("ERROR    Can't allocate memory.\n");
	    return (OUT_OF_MEMORY);
	}
	newp->next = NULL;

	/* Get profile scalar data. */
	c_tbegti (tabinfo->tp, tabinfo->cp_npts,     row, &newp->npts);
	c_tbegti (tabinfo->tp, tabinfo->cp_nptsoff,  row, &newp->nptsoff);
	c_tbegtd (tabinfo->tp, tabinfo->cp_subscale, row, subscale);
	c_tbegtd (tabinfo->tp, tabinfo->cp_minw,     row, &newp->minw);
	c_tbegtd (tabinfo->tp, tabinfo->cp_maxw,     row, &newp->maxw);
	c_tbegti (tabinfo->tp, tabinfo->cp_minp,     row, &newp->minp);
	c_tbegti (tabinfo->tp, tabinfo->cp_maxp,     row, &newp->maxp);
	c_tbegtr (tabinfo->tp, tabinfo->cp_sn,       row, &newp->sn);

	/* Alloc array memory. */
	newp->profoff = (double *) malloc (newp->nptsoff * sizeof (double));
	newp->prof    = (double *) malloc (newp->npts * sizeof (double));
	if (newp->profoff == NULL || newp->prof == NULL) {
	    printf ("ERROR    Can't allocate memory.\n");
	    return (OUT_OF_MEMORY);
	}

	/* Get profile array data. */
	npts1 = c_tbagtd (tabinfo->tp, tabinfo->cp_profoff, row,
			 newp->profoff, 1, newp->nptsoff);
	if (c_iraferr())
	    return (TABLE_ERROR);
	npts2 = c_tbagtd (tabinfo->tp, tabinfo->cp_prof, row,
			 newp->prof, 1, newp->npts);
	if (c_iraferr())
	    return (TABLE_ERROR);

	if (npts1 != newp->nptsoff) {
	    c_tbtclo (tabinfo->tp);
	    printf ("ERROR    Inconsistent array info in OPROFTAB\n");
	    return (TABLE_ERROR);
	}

	/* Insert newp into the profile list. */
	if ((status = NewProfile (profa, newp)))
	    return (status);

	free (newp->profoff);
	free (newp->prof);
	free (newp);

	return (0);
}



/* This routine closes the OPROFTAB table. */

static int CloseProfileTab (TblInfo *tabinfo) {

	c_tbtclo (tabinfo->tp);
	if (c_iraferr())
	    return (TABLE_ERROR);

	return (0);
}



