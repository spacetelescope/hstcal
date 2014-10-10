# include <fitsio.h>
# include "ctables.h"

IRAFPointer c_tbcnum (IRAFPointer tp, int colnum) {

/* Return the column pointer for the specified column number.
arguments:
IRAFPointer tp          i: table descriptor
int colnum              i: column number (one indexed)

function value          o: column descriptor, or NULL if colnum is less than
                           or equal to 0 or greater than the number of columns
*/

        TableDescr *tbl_descr;
        ColumnDescr *col_descr;
        IRAFPointer cp;
        int i;
        int foundit;

        tbl_descr = (TableDescr *)tp;

        if (colnum <= 0 || colnum > tbl_descr->ncols)
            return NULL;

        cp = tbl_descr->columns[colnum-1];
        col_descr = (ColumnDescr *)cp;
        if (col_descr->colnum == colnum)
            return cp;

        /* we shouldn't get here */
        foundit = 0;
        for (i = 0;  i < tbl_descr->ncols;  i++) {
            cp = tbl_descr->columns[i];
            col_descr = (ColumnDescr *)cp;
            if (col_descr->colnum == colnum) {
                foundit = 1;
                break;
            }
        }
        if (!foundit)
            cp = NULL;

        return cp;
}
