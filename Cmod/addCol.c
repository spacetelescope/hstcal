# include <string.h>
# include <fitsio.h>
# include "ctables.h"

void addCol (IRAFPointer tp, IRAFPointer cp, char *colname,
                char *colunits, char *colfmt, int datatype) {

/* This is called by c_tbcdef1.  If the table does not exist yet, this
   function will return without doing anything.  If the table does exist,
   the column will be added to the table, and the keywords for units and
   display format will be added to the header.  If no display format was
   specified, a default will be assigned based on the data type.
arguments:
IRAFPointer tp          i: table descriptor
IRAFPointer cp          i: column descriptor
char *colname           i: column name
char *colunits          i: units for column, or null
char *colfmt            i: display format for column, or null
int datatype            i: IRAF data type for column
*/
        
        TableDescr *tbl_descr;
        ColumnDescr *c_descr;
        char *keyword, *tdisp;
        int colnum;
        int status = 0;
        char *errmsg;

        tbl_descr = (TableDescr *)tp;
        c_descr = (ColumnDescr *)cp;

        if (tbl_descr->table_exists) {
            keyword = (char *)calloc (SZ_FITS_STR+1, sizeof(char));
            tdisp = (char *)calloc (SZ_FITS_STR+1, sizeof(char));

            if (colfmt[0] == '\0') {
                if (datatype == IRAF_REAL) {
                    strcpy (tdisp, "G15.7");
                } else if (datatype == IRAF_DOUBLE) {
                    strcpy (tdisp, "G25.16");
                } else if (datatype == IRAF_SHORT) {
                    strcpy (tdisp, "I11");
                } else if (datatype == IRAF_INT) {
                    strcpy (tdisp, "I11");
                } else if (datatype == IRAF_BOOL) {
                    strcpy (tdisp, "L6");
                } else if (datatype < 0) {
                    sprintf (tdisp, "A%d", datatype);   /* left justified */
                } else {
                    errmsg = (char *)calloc (SZ_ERRMESS+1, sizeof(char));
                    sprintf (errmsg, "Data type %d not recognized.", datatype);
                    setError (ERR_DATATYPE_UNKNOWN, errmsg);
                    free (errmsg);
                }
            } else {
                strcpy (tdisp, colfmt);
            }

            colnum = c_descr->colnum;           /* one indexed */

            /* fits_insert_col = fficol */
            fits_insert_col (tbl_descr->fptr, colnum, colname, c_descr->tform,
                        &status);
            if (status != 0)
                setError (status, "couldn't create column");
            tbl_descr->ncols = colnum;

            if (colunits[0] != '\0') {
                sprintf (keyword, "TUNIT%d", colnum);
                /* fits_write_key_str = ffpkys */
                fits_write_key_str (tbl_descr->fptr, keyword, colunits,
                        "units for column", &status);
                if (status != 0)
                    setError (status, "couldn't add TUNITi keyword");
            }

            sprintf (keyword, "TDISP%d", colnum);
            /* fits_write_key_str = ffpkys */
            fits_write_key_str (tbl_descr->fptr, keyword, tdisp,
                        "display format for column", &status);
            if (status != 0)
                setError (status, "couldn't add TDISPi keyword");

            free (tdisp);
            free (keyword);
        }
}
