# include <stdlib.h>
# include <string.h>
# include <fitsio.h>
# include "ctables.h"

int c_tbcigi (IRAFPointer cp, int param) {

/* Get the value of an integer parameter pertaining to a column.
arguments:
IRAFPointer cp          i: column descriptor
int param               i: code for the parameter to be gotten
                           TBL_COL_DATATYPE:  IRAF data type (or -w for a
                              string, where w is the maximum string length)
                           TBL_COL_LENDATA:   length of array (or 1 if scalar)
                           TBL_COL_FMTLEN:    number of char needed to print
function value          o: the value of the parameter
*/

        TableDescr *tbl_descr;
        ColumnDescr *col_descr;
        int value;
        char *errmsg;
        int status = 0;

        col_descr = (ColumnDescr *)cp;

        if (param == TBL_COL_DATATYPE) {
            value = col_descr->datatype;

        } else if (param == TBL_COL_LENDATA) {
            /* number of elements in the array, or 1 if scalar column */
            value = col_descr->nelem;

        } else if (param == TBL_COL_FMTLEN) {
            /* width needed for printing the value */
            tbl_descr = (TableDescr *)col_descr->tp;
            /* fits_get_col_display_width = ffgcdw */
            fits_get_col_display_width (tbl_descr->fptr, col_descr->colnum,
                        &value, &status);

        } else {
            errmsg = (char *)calloc (SZ_ERRMESS+1, sizeof(char));
            sprintf (errmsg,
                "c_tbcigi:  Parameter code %d not recognized.", param);
            setError (ERR_PARAMETER_UNKNOWN, errmsg);
            free (errmsg);
            value = 0;
        }

        return value;
}
