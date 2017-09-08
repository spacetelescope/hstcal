# include <stdio.h>
# include <stdlib.h>
# include <fitsio.h>
#include "hstcal.h"
# include "ctables.h"

# define COL_INCR   10

/* This is a set of functions used by the tables C interface functions.
   The prototypes for these are in ctables.h.

    TableDescr *tbl_descr       pointer to a table descriptor
    ColumnDescr *col_descr      pointer to a column descriptor
    IRAFPointer tp              (void *)tbl_descr
    IRAFPointer cp              (void *)col_descr

    tp = init_tp();
    cp = init_cp (tp);
    free_tp (tp);
    columnSpace (tp, newcols);

    tbl_descr = (TableDescr *)tp;
    col_descr = (ColumnDescr *)cp;

    When c_tbtopn is called, a tbl_descr is created that points to a struct
    with information about the table.  c_tbtopn returns tp, which is a pointer
    to this object, but cast to void *.  tp is input to almost every other
    table I/O function.

    For each column in an input table and for each new column created, a
    col_descr is created that points to a struct with information about
    that column.  cp is a pointer to a col_descr, but cast to void *.
    Table I/O functions that read or write elements in a column, or that
    get information about the column, use cp to refer to the column.

    The reason that tp and cp are pointers to void is so they can both be
    defined as of type IRAFPointer.
*/

static int reallocColArray (TableDescr *tbl_descr, int new_alloc_cols);
static void freeTableP (IRAFPointer tp);
static void freeColumnP (IRAFPointer cp);

/* create a tp (table pointer), and initialize to default values */
IRAFPointer init_tp (void) {

        TableDescr *tbl_descr;
        IRAFPointer tp;

        tbl_descr = (TableDescr *)calloc (1, sizeof(TableDescr));
        tbl_descr->tablename = (char *)calloc (CHAR_FNAME_LENGTH+1, sizeof(char));
        tbl_descr->fullname = (char *)calloc (CHAR_FNAME_LENGTH+1, sizeof(char));
        tbl_descr->filename = (char *)calloc (CHAR_FNAME_LENGTH+1, sizeof(char));
        tbl_descr->brackets = (char *)calloc (CHAR_FNAME_LENGTH+1, sizeof(char));
        tbl_descr->fptr = NULL;
        tbl_descr->template = NULL;
        tbl_descr->table_exists = 0;
        tbl_descr->iomode = IRAF_READ_ONLY;
        tbl_descr->hdunum = -1;
        tbl_descr->hdutype = -1;
        tbl_descr->nrows = 0;
        tbl_descr->ncols = 0;

        tbl_descr->alloc_cols = 0;      /* changed by columnSpace */
        tbl_descr->columns = NULL;
        /* allocate space for at least one column */
        columnSpace (tbl_descr, 1);

        tp = (void *)tbl_descr;

        return tp;
}

/* Allocate space for column pointers if necessary. */
void columnSpace (TableDescr *tbl_descr, int newcols) {

        /* If the tbl_descr->columns array is not large enough to contain
           newcols additional columns, tbl_descr->columns will be reallocated
           at a larger size.
        */

        int total_cols;         /* total number of columns we currently need */
        int new_alloc_cols;

        total_cols = tbl_descr->ncols + newcols;
        if (total_cols > tbl_descr->alloc_cols) {       /* need more space? */
            new_alloc_cols = total_cols + COL_INCR;
            if (reallocColArray (tbl_descr, new_alloc_cols) > 0)
                setError (ERR_OUT_OF_MEMORY, "out of memory");
        }
}

static int reallocColArray (TableDescr *tbl_descr, int new_alloc_cols) {

        int alloc_cols;
        int i;

        alloc_cols = tbl_descr->alloc_cols;
        tbl_descr->columns = (IRAFPointer *)realloc (tbl_descr->columns,
                        new_alloc_cols * sizeof(IRAFPointer));
        if (tbl_descr->columns == NULL)
            return 1;
        tbl_descr->alloc_cols = new_alloc_cols;

        for (i = alloc_cols;  i < new_alloc_cols;  i++)
            tbl_descr->columns[i] = NULL;
        return 0;
}

/* create a cp (column pointer), and initialize to default values */
IRAFPointer init_cp (IRAFPointer tp) {

        ColumnDescr *col_descr;
        IRAFPointer cp;

        col_descr = (ColumnDescr *)calloc (1, sizeof(ColumnDescr));
        col_descr->tp = tp;
        col_descr->colnum = -1;
        col_descr->name = (char *)calloc (FLEN_VALUE+1, sizeof(char));
        col_descr->tform = (char *)calloc (FLEN_VALUE+1, sizeof(char));
        col_descr->tunit = (char *)calloc (FLEN_VALUE+1, sizeof(char));
        col_descr->tdisp = (char *)calloc (FLEN_VALUE+1, sizeof(char));
        if (col_descr->name == NULL || col_descr->tform == NULL ||
            col_descr->tunit == NULL || col_descr->tdisp == NULL) {
                setError (ERR_OUT_OF_MEMORY, "out of memory");
                return NULL;
        }
        col_descr->typecode = 0;
        col_descr->datatype = 0;
        col_descr->var_length = 0;
        col_descr->repeat = 0;
        col_descr->nelem = 0;
        col_descr->width = 0;

        cp = (void *)col_descr;
        return cp;
}

/* free a tp (table pointer) */
void free_tp (IRAFPointer tp) {
        freeTableP (tp);
}

/* free memory for a table pointer */
static void freeTableP (IRAFPointer tp) {

        TableDescr *tbl_descr;
        IRAFPointer cp;
        int i;

        tbl_descr = (TableDescr *)tp;
        if (tbl_descr != NULL) {
            /* free memory for strings */
            if (tbl_descr->tablename != NULL)
                free (tbl_descr->tablename);
            if (tbl_descr->fullname != NULL)
                free (tbl_descr->fullname);
            if (tbl_descr->filename != NULL)
                free (tbl_descr->filename);
            if (tbl_descr->brackets != NULL)
                free (tbl_descr->brackets);
            /* free memory for columns */
            for (i = 0;  i < tbl_descr->ncols;  i++) {
                cp = tbl_descr->columns[i];
                freeColumnP (cp);
            }
            free (tbl_descr->columns);
        }
        free (tbl_descr);
}

/* free memory for a column pointer */
static void freeColumnP (IRAFPointer cp) {

        ColumnDescr *col_descr;

        col_descr = (ColumnDescr *)cp;
        if (col_descr != NULL) {
            if (col_descr->name != NULL)
                free (col_descr->name);
            if (col_descr->tform != NULL)
                free (col_descr->tform);
            if (col_descr->tunit != NULL)
                free (col_descr->tunit);
            if (col_descr->tdisp != NULL)
                free (col_descr->tdisp);
        }
        free (col_descr);
}
