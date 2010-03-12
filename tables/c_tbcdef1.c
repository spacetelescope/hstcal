# include <string.h>
# include <fitsio.h>
# include "ctables.h"

void c_tbcdef1 (IRAFPointer tp, IRAFPointer *cp,
        char *colname, char *colunits, char *colfmt, int datatype, int nelem) {

/* Create a new column.
arguments:
IRAFPointer tp          i: table descriptor
IRAFPointer cp          i: column descriptor
char *colname           i: column name
char *colunits          i: units for column, or null
char *colfmt            i: display format for column, or null
int datatype            i: IRAF data type for column
int nelem               i: length of array, or 1 for a scalar column
*/

        char ftn_fmt[SZ_FITS_STR+1];    /* colfmt in Fortran style */

        /* It is an error if the column already exists.
           On the other hand, this sets cp for the existing column, so it's
           mostly transparent.
        */
        c_tbcfnd1 (tp, colname, cp);
        if (*cp != NULL) {
            char errmsg[SZ_ERRMESS+1];
            sprintf (errmsg, "column %s already exists", colname);
            setError (ERR_COLUMN_ALREADY_EXISTS, errmsg);
            return;
        }

        /* If the display format colfmt was given in C notation, convert
           to Fortran notation for the FITS file.
        */
        cToFortran (colfmt, ftn_fmt);

        /* create cp, and save info about the column */
        initCol (tp, cp, colname, colunits, ftn_fmt, datatype, nelem);

        /* if the table exists, add the column to the table */
        addCol (tp, *cp, colname, colunits);
}
