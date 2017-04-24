# include <ctype.h>
# include <stdlib.h>
# include <string.h>
# include <fitsio.h>
# include "ctables.h"

void c_tbcfnd1 (IRAFPointer tp, const char *colname, IRAFPointer *cp) {

/* Find the column pointer for the column with the specified name.
arguments:
IRAFPointer tp          i: table descriptor
char *colname           i: column name
IRAFPointer *cp         o: column descriptor
*/

        *cp = NULL;
        if (!tp || !colname || *colname == '\0')
            return;

        char lc_colname[SZ_FITS_STR+1], lc_name[SZ_FITS_STR+1];
        str_lower (lc_colname, colname);

        TableDescr *tbl_descr = (TableDescr *)tp;

        {unsigned i;
        for (i = 0;  i < tbl_descr->ncols;  ++i) {
            ColumnDescr * col_descr = (ColumnDescr *)tbl_descr->columns[i];
            str_lower (lc_name, col_descr->name);
            if (strcmp (lc_colname, lc_name) == 0) {
                *cp = tbl_descr->columns[i];
                return;
            }
        }}
}

IRAFPointer c_tbcfnd1_retPtr (IRAFPointer tp, const char *colname)
{
    void * ptr = NULL;
    c_tbcfnd1(tp, colname, &ptr);
    return ptr;
}
