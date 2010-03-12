# include <stdio.h>
# include <stdlib.h>
# include "ctables.h"

# define COL_INCR   10
# define INCREMENT 100
static TblInfoPtr **dataIRAFPointer=NULL;
static int nelem=0;	/* current length of dataIRAFPointer array */

/* public functions (prototypes are in ctables.h):

    tp = init_tp();
    cp = init_cp();
    free_tp (tp);
    columnSpace (tp, newcols);

    tptr = makeTableDescr (t_descr);
    tptr = makeColumnDescr (c_descr);
    tp = saveIrafP (tptr);
    t_descr = getTableDescr (tp);
    c_descr = getColumnDescr (cp);

    tp and cp are integers, indices in the local array dataIRAFPointer.
    tp is returned by c_tbtopn, and it's input to almost every other table
    I/O function.  tp is like a file handle for a table.
    cp is the output variable for c_tbcfnd1 (and the output variable for
    c_tbcfnd is an array of cd).  This is the identifier for a column in
    an open table.
    Internally, t_descr is a pointer to a TableDescr (table descriptor) object,
    and c_descr is a pointer to a ColumnDescr (column descriptor) object.
    tptr is a pointer to a generic structure (TblInfoPtr) that contains both
    a *t_descr and a *c_descr, one of which is NULL; this is used in preference
    to a pointer to void for either a TableDescr or a ColumnDescr.
    Each element of the dataIRAFPointer array is a *tptr, i.e. a pointer to
    a TblInfoPtr object.

    IRAFPointer tp              "pointer" to a table descriptor
    IRAFPointer cp              "pointer" to a column descriptor
    TblInfoPtr *tptr            struct for either a table or column descr
    TableDescr *t_descr         table descriptor
    ColumnDescr *c_descr        column descriptor

    The array of pointers to TblInfoPtr structs is dataIRAFPointer, and
    it is managed by functions in this file:
	findEmptyElement, reallocBuffer
*/

static IRAFPointer findEmptyElement (void);
static int reallocBuffer (int increment);
static int reallocColArray (TableDescr *t_descr, int new_alloc_cols);
static void freeTableP (IRAFPointer tp);
static void freeColumnP (IRAFPointer cp);

/*
saveIrafP        save in list
getTableDescr    get table descriptor from list
getColumnDescr   get column descriptor from list
makeTableDescr   allocate tptr and save a table descriptor
makeColumnDescr  allocate tptr and save a column descriptor
reallocBuffer    reallocate
xxx              deallocate (delete from list)
*/

/* create a tp (table pointer), and initialize to null values */
IRAFPointer init_tp (void) {

	TblInfoPtr *tptr;
	TableDescr *t_descr;
	IRAFPointer tp;

	t_descr = (TableDescr *)calloc (1, sizeof(TableDescr));
	t_descr->tablename = (char *)calloc (SZ_FNAME, sizeof(char));
	t_descr->filename = (char *)calloc (SZ_FNAME, sizeof(char));
	t_descr->brackets = (char *)calloc (SZ_FNAME, sizeof(char));
	t_descr->fptr = NULL;
	t_descr->hdunum = -1;
	t_descr->hdutype = -1;
	t_descr->nrows = 0;
	t_descr->ncols = 0;

	t_descr->alloc_cols = 0;	/* changed by columnSpace */
	t_descr->columns = NULL;
	/* allocate space for at least one column */
	columnSpace (t_descr, 1);

	tptr = makeTableDescr (t_descr);
	tp = saveIrafP (tptr);

	return tp;
}

/* Allocate space for column pointers if necessary. */
void columnSpaceXXX (IRAFPointer tp, int newcols) {

	TableDescr *t_descr;
	t_descr = getTableDescr (tp);
	columnSpace (t_descr, newcols);
}

/* Allocate space for column pointers if necessary. */
void columnSpace (TableDescr *t_descr, int newcols) {

	int total_cols;		/* total number of columns we currently need */
	int new_alloc_cols;

	total_cols = t_descr->ncols + newcols;
	if (total_cols > t_descr->alloc_cols) {		/* need more space */
	    if (newcols == 1)
		new_alloc_cols = COL_INCR;
	    else
		new_alloc_cols = total_cols + COL_INCR;
	    if (reallocColArray (t_descr, new_alloc_cols) > 0)
		setError (ERR_OUT_OF_MEMORY, "out of memory");
	}
}

static int reallocColArray (TableDescr *t_descr, int new_alloc_cols) {

	int alloc_cols;
	int i;

	alloc_cols = t_descr->alloc_cols;
	t_descr->columns = (IRAFPointer *)realloc (t_descr->columns,
			new_alloc_cols * sizeof(IRAFPointer));
	if (t_descr->columns == NULL)
	    return 1;
	t_descr->alloc_cols = new_alloc_cols;

	for (i = alloc_cols;  i < new_alloc_cols;  i++)
	    t_descr->columns[i] = -1;
	return 0;
}

/* create a cp (column pointer), and initialize to null values */
IRAFPointer init_cp (void) {

	TblInfoPtr *tptr;
	ColumnDescr *c_descr;
	IRAFPointer cp;

	c_descr = (ColumnDescr *)calloc (1, sizeof(ColumnDescr));
	c_descr->colnum = -1;
	c_descr->name = (char *)calloc (SZ_FITS_STR, sizeof(char));
	c_descr->tunit = (char *)calloc (SZ_FITS_STR, sizeof(char));
	c_descr->tdisp = (char *)calloc (SZ_FITS_STR, sizeof(char));
	c_descr->dtype = -1;
	c_descr->repeat = 0;

	tptr = makeColumnDescr (c_descr);
	cp = saveIrafP (tptr);
	return cp;
}
/* free a tp (table pointer) */
void free_tp (IRAFPointer tp) {
	freeTableP (tp);
}

/* allocate a TblInfoPtr object, and copy a table descriptor into it */
TblInfoPtr *makeTableDescr (TableDescr *t_descr) {

	TblInfoPtr *tptr;

	tptr = (TblInfoPtr *)calloc (1, sizeof(TblInfoPtr));
	tptr->ptype = TABLE_DESCR;
	tptr->t_descr = t_descr;
	tptr->c_descr = NULL;
	return tptr;
}

/* allocate a TblInfoPtr object, and copy a column descriptor into it */
TblInfoPtr *makeColumnDescr (ColumnDescr *c_descr) {

	TblInfoPtr *tptr;

	tptr = (TblInfoPtr *)calloc (1, sizeof(TblInfoPtr));
	tptr->ptype = COLUMN_DESCR;
	tptr->t_descr = NULL;
	tptr->c_descr = c_descr;
	return tptr;
}

/* save a TblInfoPtr object, which may contain either a table descriptor or
   a column descriptor
*/
IRAFPointer saveIrafP (TblInfoPtr *tptr) {

	IRAFPointer tp;		/* could be either a tp or cp */

	tp = findEmptyElement();
	if (tp >= 0)
	    dataIRAFPointer[tp] = tptr;

	return tp;
}

/* free memory for a table pointer */
static void freeTableP (IRAFPointer tp) {

	TblInfoPtr *tptr;
	TableDescr *t_descr;
	IRAFPointer cp;
	int i;

	tptr = dataIRAFPointer[tp];
	if (tptr == NULL)
	    return;
	if (tptr->ptype != TABLE_DESCR)
	    printf ("debug:  freeTableP was called for a column descriptor");
	t_descr = tptr->t_descr;
	if (t_descr != NULL) {
	    /* free memory for strings */
	    if (t_descr->tablename != NULL)
		free (t_descr->tablename);
	    if (t_descr->filename != NULL)
		free (t_descr->filename);
	    if (t_descr->brackets != NULL)
		free (t_descr->brackets);
	    /* free memory for columns */
	    for (i = 0;  i < t_descr->ncols;  i++) {
		cp = t_descr->columns[i];
		freeColumnP (cp);
	    }
	}
	free (tptr);
	dataIRAFPointer[tp] = NULL;	/* flag this element as empty */
}

/* free memory for a column pointer */
static void freeColumnP (IRAFPointer cp) {

	TblInfoPtr *tptr;
	ColumnDescr *c_descr;

	tptr = dataIRAFPointer[cp];
	if (tptr == NULL)
	    return;
	if (tptr->ptype != COLUMN_DESCR)
	    printf ("debug:  freeColumnP was called for a table descriptor");
	c_descr = tptr->c_descr;
	if (c_descr != NULL) {
	    if (c_descr->name != NULL)
		free (c_descr->name);
	    if (c_descr->tunit != NULL)
		free (c_descr->tunit);
	    if (c_descr->tdisp != NULL)
		free (c_descr->tdisp);
	    free (c_descr);
	}
	free (tptr);
	tptr = NULL;
	dataIRAFPointer[cp] = NULL;	/* flag this element as empty */
}

TableDescr *getTableDescr (IRAFPointer tp) {
	/* should check whether the type is correct */
	return (dataIRAFPointer[tp]->t_descr);
}

ColumnDescr *getColumnDescr (IRAFPointer cp) {
	/* should check whether the type is correct */
	return (dataIRAFPointer[cp]->c_descr);
}

static IRAFPointer findEmptyElement (void) {

	IRAFPointer tp = -1;
	int i;
	int status = 0;

	if (nelem < 1) {
	    if ((status = reallocBuffer (INCREMENT)) > 0)
		return tp;
	}
	for (i = 0;  i < nelem;  i++) {
	    if (dataIRAFPointer[i] == NULL) {
		tp = i;
		break;
	    }
	}
	if (tp < 0) {
	    if ((status = reallocBuffer (INCREMENT)) > 0)
		return tp;
	    for (i = 0;  i < nelem;  i++) {
		if (dataIRAFPointer[i] == NULL) {
		    tp = i;
		    break;
		}
	    }
	}
	return tp;
}

static int reallocBuffer (int increment) {

	int new_nelem;
	int i;

	new_nelem = nelem + increment;
	dataIRAFPointer = (TblInfoPtr **)realloc (dataIRAFPointer,
			new_nelem * sizeof(TblInfoPtr));
	if (dataIRAFPointer == NULL)
	    return 1;
	for (i = nelem;  i < new_nelem;  i++)
	    dataIRAFPointer[i] = NULL;
	nelem = new_nelem;
	return 0;
}
