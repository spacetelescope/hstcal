# include <fitsio.h>
# include "ctables.h"

int c_tbtacc (char *tablename) {

/* Return 1 if the specified table can be opened as a table.  If it exists
   but is actually an image, 0 will be returned.
arguments:
char *tablename         i: name of table; may include extname or HDU number
                           in brackets (HDU = 0 is the primary HDU)

function value          o: 1 if 'tablename' exists and is a table, 0 otherwise
*/

        IRAFPointer tp;
        TableDescr *tbl_descr;
        int value;
        int status;

        tp = c_tbtopn (tablename, IRAF_READ_ONLY, NULL);
        status = checkError();
        if (tp == NULL || status != 0) {
            clearError();
            value = 0;
        } else {

            /* if this is an image, it's not a table */
            tbl_descr = (TableDescr *)tp;
            if (tbl_descr->hdutype == IMAGE_HDU)
                value = 0;
            else
                value = 1;
            c_tbtclo (tp);
        }

        return value;
}
