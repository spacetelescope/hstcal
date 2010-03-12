# include <fitsio.h>
# include "ctables.h"

void c_tbeptb (IRAFPointer tp, IRAFPointer cp, int row, Bool buffer) {

/* Write a boolean value to a table column.
arguments:
IRAFPointer tp          i: table descriptor
IRAFPointer cp          i: column descriptor
int row                 i: row number (one indexed)
Bool buffer             i: value (True or False) to write to table
*/

        TableDescr *tbl_descr;
        ColumnDescr *col_descr;
        long firstelem=1, nelem=1;
        char value[11]={'\0'};
        int status = 0;

        tbl_descr = (TableDescr *)tp;
        col_descr = (ColumnDescr *)cp;

        if (buffer == True)
            value[0] = 1;
        else
            value[0] = 0;

        /* fits_write_col_log = ffpcll */
        fits_write_col_log (tbl_descr->fptr, col_descr->colnum,
                        (long)row, firstelem, nelem, value, &status);
        if (status != 0)
            setError (status, "c_tbeptb:  error writing element");
}

void c_tbeptd (IRAFPointer tp, IRAFPointer cp, int row, double buffer) {

/* Write a value of type double to a table column.
arguments:
IRAFPointer tp          i: table descriptor
IRAFPointer cp          i: column descriptor
int row                 i: row number (one indexed)
double buffer           i: value to write to table
*/

        TableDescr *tbl_descr;
        ColumnDescr *col_descr;
        long firstelem=1, nelem=1;
        int status = 0;

        tbl_descr = (TableDescr *)tp;
        col_descr = (ColumnDescr *)cp;

        if (buffer == IRAF_INDEFD) {
            /* fits_write_col_null / ffpclu */
            fits_write_col_null (tbl_descr->fptr, col_descr->colnum,
                        (long)row, firstelem, nelem, &status);
        } else {
            /* fits_write_col_dbl = ffpcld */
            fits_write_col_dbl (tbl_descr->fptr, col_descr->colnum,
                        (long)row, firstelem, nelem, &buffer, &status);
        }
        if (status != 0)
            setError (status, "c_tbeptd:  error writing element");
}

void c_tbeptr (IRAFPointer tp, IRAFPointer cp, int row, float buffer) {

/* Write a value of type float to a table column.
arguments:
IRAFPointer tp          i: table descriptor
IRAFPointer cp          i: column descriptor
int row                 i: row number (one indexed)
float buffer            i: value to write to table
*/

        TableDescr *tbl_descr;
        ColumnDescr *col_descr;
        long firstelem=1, nelem=1;
        int status = 0;

        tbl_descr = (TableDescr *)tp;
        col_descr = (ColumnDescr *)cp;

        if (buffer >= 0.99999 * IRAF_INDEFR &&
            buffer <= 1.00001 * IRAF_INDEFR) {
            fits_write_col_null (tbl_descr->fptr, col_descr->colnum,
                        (long)row, firstelem, nelem, &status);
        } else {
            /* fits_write_col_flt = ffpcle */
            fits_write_col_flt (tbl_descr->fptr, col_descr->colnum,
                        (long)row, firstelem, nelem, &buffer, &status);
        }
        if (status != 0)
            setError (status, "c_tbeptr:  error writing element");
}

void c_tbepti (IRAFPointer tp, IRAFPointer cp, int row, int buffer) {

/* Write an integer value to a table column.
arguments:
IRAFPointer tp          i: table descriptor
IRAFPointer cp          i: column descriptor
int row                 i: row number (one indexed)
int buffer              i: value to write to table
*/

        TableDescr *tbl_descr;
        ColumnDescr *col_descr;
        long firstelem=1, nelem=1;
        int status = 0;

        tbl_descr = (TableDescr *)tp;
        col_descr = (ColumnDescr *)cp;

        if (buffer == IRAF_INDEFI) {
            fits_write_col_null (tbl_descr->fptr, col_descr->colnum,
                        (long)row, firstelem, nelem, &status);
        } else {
            /* fits_write_col_int = ffpclk */
            fits_write_col_int (tbl_descr->fptr, col_descr->colnum,
                        (long)row, firstelem, nelem, &buffer, &status);
        }
        if (status != 0)
            setError (status, "c_tbepti:  error writing element");
}

void c_tbepts (IRAFPointer tp, IRAFPointer cp, int row, short buffer) {

/* Write a short integer value to a table column.
arguments:
IRAFPointer tp          i: table descriptor
IRAFPointer cp          i: column descriptor
int row                 i: row number (one indexed)
short buffer            i: value to write to table
*/

        TableDescr *tbl_descr;
        ColumnDescr *col_descr;
        long firstelem=1, nelem=1;
        int status = 0;

        tbl_descr = (TableDescr *)tp;
        col_descr = (ColumnDescr *)cp;

        if (buffer == IRAF_INDEFS) {
            fits_write_col_null (tbl_descr->fptr, col_descr->colnum,
                        (long)row, firstelem, nelem, &status);
        } else {
            /* fits_write_col_sht = ffpcli */
            fits_write_col_sht (tbl_descr->fptr, col_descr->colnum,
                        (long)row, firstelem, nelem, &buffer, &status);
        }
        if (status != 0)
            setError (status, "c_tbepts:  error writing element");
}

void c_tbeptt (IRAFPointer tp, IRAFPointer cp, int row, char *buffer) {

/* Write a string value to a table column.
arguments:
IRAFPointer tp          i: table descriptor
IRAFPointer cp          i: column descriptor
int row                 i: row number (one indexed)
char *buffer            i: value to write to table
*/

        TableDescr *tbl_descr;
        ColumnDescr *col_descr;
        long firstelem=1, nelem=1;
        int status = 0;

        tbl_descr = (TableDescr *)tp;
        col_descr = (ColumnDescr *)cp;

        /* fits_write_col_str = ffpcls */
        fits_write_col_str (tbl_descr->fptr, col_descr->colnum,
                        (long)row, firstelem, nelem, &buffer, &status);
        if (status != 0)
            setError (status, "c_tbeptt:  error writing element");
}
