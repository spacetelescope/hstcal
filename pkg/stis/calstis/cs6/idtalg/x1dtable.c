# include <stdio.h>
# include <stdlib.h>

# include "xtables.h"

# include "stis.h"
# include "err.h"
# include "calstis6.h"

/*
   GetX1DTable -- Read X1D table from file.

   Data arrays must be allocated by caller with AllocX1DTable,
   and fred with FreeX1DTable.


   Revision history:
   ----------------
   21 Feb 00  -  Implemented (I.Busko)
*/

int GetX1DTable (TblDesc *tabptr, RowContents **x1d) {

/* arguments
TblDesc *tabptr;	o: _x1d temporary table
RowContents **x1d;	o: row contents
*/
	int irow, i;
	int nwave, ngross, nnet, nextrlocy;

	void FreeX1DTable (RowContents **, int);

	/* Set up pointers to columns. */

	c_tbcfnd1 (tabptr->tp, SPORDER,    &tabptr->sporder);
	c_tbcfnd1 (tabptr->tp, NELEM,      &tabptr->npts);
	c_tbcfnd1 (tabptr->tp, WAVELENGTH, &tabptr->wave);
	c_tbcfnd1 (tabptr->tp, GROSS,      &tabptr->gross);
	c_tbcfnd1 (tabptr->tp, NET,        &tabptr->net);
	c_tbcfnd1 (tabptr->tp, EXTRLOCY,   &tabptr->extrlocy);
	if ((tabptr->sporder    == 0) || (tabptr->npts == 0) ||
            (tabptr->wave == 0) || (tabptr->gross == 0) ||
            (tabptr->net  == 0)) {
	    FreeX1DTable (x1d, tabptr->nrows);
	    c_tbtclo (tabptr->tp);
	    printf ("Column not found in input table.\n");
	    return (TABLE_ERROR);
	}

	/* Alloc memory and read data from each row.

           This is an array of RowContents data structures with
           tabptr->nrows elements. Each element holds data arrays of
           tabptr->npts length.
        */

	for (irow = 0; irow < tabptr->nrows; irow++) {
	    c_tbegts (tabptr->tp, tabptr->sporder, irow+1,
                      &(x1d[irow]->sporder));
	    c_tbegts (tabptr->tp, tabptr->npts, irow+1,
                      &(x1d[irow]->npts));

	    x1d[irow]->wave     = (double *) calloc (x1d[irow]->npts,
                                                        sizeof(double));
	    x1d[irow]->gross    = (float *) calloc (x1d[irow]->npts,
                                                        sizeof(float));
	    x1d[irow]->net      = (float *) calloc (x1d[irow]->npts,
                                                        sizeof(float));
	    x1d[irow]->extrlocy = (float *) calloc (x1d[irow]->npts,
                                                        sizeof(float));
	    if (x1d[irow]->wave == NULL || x1d[irow]->gross    == NULL ||
                x1d[irow]->net  == NULL || x1d[irow]->extrlocy == NULL) {
	        FreeX1DTable (x1d, tabptr->nrows);
	        c_tbtclo (tabptr->tp);
                printf ("Not enough memory to allocate data arrays.\n");
	        return (OUT_OF_MEMORY);
	    }
	    nwave     = c_tbagtd (tabptr->tp, tabptr->wave, irow+1,
                                  x1d[irow]->wave, 1, x1d[irow]->npts);
	    ngross    = c_tbagtr (tabptr->tp, tabptr->gross, irow+1,
                                  x1d[irow]->gross, 1, x1d[irow]->npts);
	    nnet      = c_tbagtr (tabptr->tp, tabptr->net, irow+1,
                                  x1d[irow]->net, 1, x1d[irow]->npts);
	    nextrlocy = c_tbagtr (tabptr->tp, tabptr->extrlocy, irow+1,
                                  x1d[irow]->extrlocy,1,x1d[irow]->npts);
	    if ((nwave     != x1d[irow]->npts) ||
                (ngross    != x1d[irow]->npts) ||
                (nnet      != x1d[irow]->npts) ||
                (nextrlocy != x1d[irow]->npts)) {
	        FreeX1DTable (x1d, tabptr->nrows);
	        c_tbtclo (tabptr->tp);
	        printf ("Insufficent number of elements read from table.\n");
	        return (TABLE_ERROR);
	    }

	    /* Compatibilize extrlocy with 0-indexed image arrays. */

	    for (i = 0; i < x1d[irow]->npts; x1d[irow]->extrlocy[i++] -= 1.0);

	}

	return (0);
}



/*  Opens the X1D table and gets number of rows. */

int OpenX1DTable (char *table, int write, TblDesc *tabptr) {

/* arguments
char *table;		i: input table name
int write;		i: read only or read/write ?
TblDesc *tabptr;	o: _x1d temporary table
*/
	if (write)
	    tabptr->tp = c_tbtopn (table, IRAF_READ_WRITE, 0);
	else
	    tabptr->tp = c_tbtopn (table, IRAF_READ_ONLY, 0);

	if (c_iraferr()) {
	    printf ("_x1d temporary table could not be opened.\n");
	    return (OPEN_FAILED);
	}

	tabptr->nrows = c_tbpsta (tabptr->tp, TBL_NROWS);
	if (c_iraferr()) {
	    c_tbtclo (tabptr->tp);
            printf ("Cannot determine the number of rows in table.\n");
	    return (OPEN_FAILED);
	}
	return (0);
}


/* Alloc memory for a X1D table array. */

RowContents **AllocX1DTable (int nrows) {

/* arguments
int nrows;		i: number of rows
*/
	RowContents **x1d;
	int i;

	x1d = (RowContents **) malloc (nrows * sizeof(RowContents *));
	if (x1d == NULL) {
            printf ("Not enough memory to allocate row pointer.\n");
	    return (NULL);
	}
	for (i = 0; i < nrows; i++) {
	    x1d[i] = (RowContents *) calloc (1, sizeof(RowContents));
	    if (x1d[i] == NULL) {
                printf ("Not enough memory to allocate row pointer.\n");
	        return (NULL);
	    }
	    x1d[i]->wave  = NULL; /* not allocated yet ! */
	    x1d[i]->gross = NULL;
	}
	return (x1d);
}



/* Frees the data areas used by GetX1DTable() */

void FreeX1DTable (RowContents **x1d, int nrows) {

/* arguments
RowContents **x1d;	i: row contents
*/
	int i;

	if (x1d != NULL) {
	    for (i = 0; i < nrows; i++) {
	        if (x1d[i]->gross != NULL)
	            free (x1d[i]->gross);
	        if (x1d[i]->wave != NULL) {
	            free (x1d[i]->extrlocy);
	            free (x1d[i]->net);
	            free (x1d[i]->wave);
	        }
	        free (x1d[i]);
	    }
	    free (x1d);
	    x1d = NULL;
	}
}
