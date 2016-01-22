# include <fitsio.h>
# include "ctables.h"

void c_tbrudf (IRAFPointer tp, IRAFPointer *cp, int numcols, int row) {

/* Set one or more values in a row to undefined.
arguments:
IRAFPointer tp          i: table descriptor
IRAFPointer *cp         i: array of descriptors for columns to be set to
                           undefined
int numcols             i: the length of the array 'cp'
int row                 i: number of the row to be modified (one indexed)
*/

        TableDescr *tbl_descr;
        ColumnDescr *col_descr;
        long firstelem=1, nelem=1;
        int i;
        int status = 0;

        tbl_descr = (TableDescr *)tp;

        for (i = 0;  i < numcols; i++) {
            col_descr = (ColumnDescr *)cp[i];
            /* fits_write_col_null / ffpclu */
            fits_write_col_null (tbl_descr->fptr, col_descr->colnum,
                        (long)row, firstelem, nelem, &status);
        }
}
