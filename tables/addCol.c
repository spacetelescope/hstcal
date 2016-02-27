# include <string.h>
# include <fitsio.h>
# include "ctables.h"

void addCol (IRAFPointer tp, IRAFPointer cp, char *colname, char *colunits) {

/* This is called by c_tbcdef1.  If the table does not exist yet, this
   function will return without doing anything.  If the table does exist,
   the column will be added to the table, and the keywords for units,
   display format, and undefined value (if int or short) will be added to
   the header.
arguments:
IRAFPointer tp          i: table descriptor
IRAFPointer cp          i: column descriptor
char *colname           i: column name
char *colunits          i: units for column, or null
*/
        
        TableDescr *tbl_descr;
        ColumnDescr *col_descr;
        char *keyword;
        int colnum;
        int datatype;           /* IRAF data type */
        int indef_int, indef_short;
        int status = 0;

        tbl_descr = (TableDescr *)tp;
        col_descr = (ColumnDescr *)cp;

        if (tbl_descr->table_exists) {
            keyword = (char *)calloc (SZ_FITS_STR+1, sizeof(char));

            colnum = col_descr->colnum;           /* one indexed */
            datatype = col_descr->datatype;
            indef_int = IRAF_INDEFI;
            indef_short = IRAF_INDEFS;

            /* fits_insert_col = fficol */
            fits_insert_col (tbl_descr->fptr, colnum, colname,
                        col_descr->tform, &status);
            if (status != 0) {
                setError (status, "c_tbcdef1:  couldn't create column");
                return;
            }
            tbl_descr->ncols = colnum;

            if (colunits[0] != '\0') {
                sprintf (keyword, "TUNIT%d", colnum);
                /* fits_update_key = ffuky */
                fits_update_key (tbl_descr->fptr, TSTRING, keyword,
                        colunits, "units for column", &status);
                if (status != 0) {
                    setError (status,
                        "c_tbcdef1:  couldn't add TUNITi keyword");
                    return;
                }
            }

            sprintf (keyword, "TDISP%d", colnum);
            fits_update_key (tbl_descr->fptr, TSTRING, keyword,
                        col_descr->tdisp, "display format for column",
                        &status);
            if (status != 0) {
                setError (status, "c_tbcdef1:  couldn't add TDISPi keyword");
                return;
            }

            if (datatype == IRAF_SHORT || datatype == IRAF_INT) {
                sprintf (keyword, "TNULL%d", colnum);
                if (datatype == IRAF_INT) {
                    int indef_int = IRAF_INDEFI;
                    fits_update_key (tbl_descr->fptr, TINT, keyword,
                        &indef_int, "undefined value for column", &status);
                } else {        /* short */
                    int indef_short = IRAF_INDEFS;
                    fits_update_key (tbl_descr->fptr, TINT, keyword,
                        &indef_short, "undefined value for column", &status);
                }
                if (status != 0) {
                    setError (status,
                        "c_tbcdef1:  couldn't add TNULLi keyword");
                    return;
                }
            }

            free (keyword);
        }
}
