# include <string.h>
# include <fitsio.h>
# include "ctables.h"

void c_tbcinf (IRAFPointer cp, int *colnum, char *colname, char *colunits,
        char *colfmt, int *datatype, int *nelem, int *lenfmt) {

/* Get the values of parameters describing a table column.  The buffers
for string values should be at least 71 char in size (including '\0').
arguments:
IRAFPointer cp          i: column descriptor
int *colnum             o: column number (one indexed)
char *colname           o: column name
char *colunits          o: units for column
char *colfmt            o: display format for column
int *datatype           o: IRAF data type of column (or -w for a string,
                           where w is the maximum string length)
int *nelem              o: number of elements in array (or 1)
int *lenfmt             o: width needed for printing one element of column
*/

        ColumnDescr *col_descr;

        col_descr = (ColumnDescr *)cp;

        *colnum = col_descr->colnum;
        strcpy (colname, col_descr->name);      /* TTYPEi */
        strcpy (colunits, col_descr->tunit);    /* TUNITi */
        strcpy (colfmt, col_descr->tdisp);      /* TDISPi */

        /* this is the IRAF data type, not the CFITSIO typecode */
        *datatype = c_tbcigi (cp, TBL_COL_DATATYPE);

        /* the number of elements in the array, or 1 for a scalar column */
        *nelem = c_tbcigi (cp, TBL_COL_LENDATA);

        /* width needed for printing a value (one element) */
        *lenfmt = c_tbcigi (cp, TBL_COL_FMTLEN);
}
