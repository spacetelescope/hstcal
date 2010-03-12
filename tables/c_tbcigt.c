# include <string.h>
# include <fitsio.h>
# include "ctables.h"

void c_tbcigt (IRAFPointer cp, int param, char *outstr, int maxch) {

/* Get the value of a string parameter pertaining to a column.
arguments:
IRAFPointer cp          i: column descriptor
int param               i: code for the parameter to be gotten
                           TBL_COL_NAME:   column name
                           TBL_COL_UNITS:  units for column
                           TBL_COL_FMT:    display format for column
char *outstr            o: the value of the parameter
int maxch               i: maximum length of 'outstr'
*/

        ColumnDescr *col_descr;

        col_descr = (ColumnDescr *)cp;

        if (param == TBL_COL_NAME) {
            copyString (outstr, col_descr->name, maxch);

        } else if (param == TBL_COL_UNITS) {
            copyString (outstr, col_descr->tunit, maxch);

        } else if (param == TBL_COL_FMT) {
            copyString (outstr, col_descr->tdisp, maxch);

        } else {
            char *errmsg;
            errmsg = (char *)calloc (SZ_ERRMESS+1, sizeof(char));
            sprintf (errmsg,
                "c_tbcigt:  Parameter code %d not recognized.", param);
            setError (ERR_PARAMETER_UNKNOWN, errmsg);
            free (errmsg);
            outstr[0] = '\0';
        }
}
