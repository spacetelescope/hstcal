# include <string.h>
# include <fitsio.h>
# include "ctables.h"

int c_tbpsta (IRAFPointer tp, int param) {

/* Get the value of a parameter describing the table.
arguments:
IRAFPointer tp          i: table descriptor
int param               i: code specifying the parameter to be gotten:
                           TBL_NROWS:  number of rows in the table
                           TBL_NCOLS:  number of columns in the table
                           TBL_WHTYPE:  which type of table (TBL_TYPE_FITS)
                           TBL_NPAR:  number of keywords in the header

function value          o: value of parameter
*/

        TableDescr *tbl_descr;
        int value;
        char *errmsg;
        int status = 0;

        tbl_descr = (TableDescr *)tp;

        switch (param) {
        case TBL_NROWS:
            value = tbl_descr->nrows;
            break;

        case TBL_NCOLS:
            value = tbl_descr->ncols;
            break;

        case TBL_WHTYPE:                /* type is FITS */
            value = TBL_TYPE_FITS;
            break;

        case TBL_NPAR:                  /* number of keywords */
            /* fits_get_hdrspace = ffghsp */
            fits_get_hdrspace (tbl_descr->fptr, &value, NULL, &status);
            break;

        default:
            errmsg = (char *)calloc (SZ_ERRMESS+1, sizeof(char));
            sprintf (errmsg,
                "c_tbpsta:  Parameter code %d not recognized.", param);
            setError (ERR_PARAMETER_UNKNOWN, errmsg);
            free (errmsg);
            value = 0;
            break;
        }

        return value;
}
