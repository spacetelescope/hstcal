# include <stdlib.h>
# include <fitsio.h>
# include "ctables.h"

void c_tbaptb (IRAFPointer tp, IRAFPointer cp, int row, Bool *buffer,
                int first, int nelem) {

/* Write an array of boolean values to a table column.
arguments:
IRAFPointer tp          i: table descriptor
IRAFPointer cp          i: column descriptor
int row                 i: row number (one indexed)
bool *buffer            i: array of values (True or False) to write
int first               i: first element to write (one indexed)
int nelem               i: number of elements to write
*/

        TableDescr *tbl_descr;
        ColumnDescr *col_descr;
        char *value;
        int i;
        int status = 0;

        tbl_descr = (TableDescr *)tp;
        col_descr = (ColumnDescr *)cp;

        value = (char *)calloc (nelem+10, sizeof(char)); /* nelem plus extra */

        for (i = 0;  i < nelem;  i++) {
            if (buffer[i] == True)
                value[i] = 1;
            else
                value[i] = 0;
        }

        /* fits_write_col_log = ffpcll */
        fits_write_col_log (tbl_descr->fptr, col_descr->colnum,
                (long)row, (long)first, (long)nelem, value, &status);
        if (status != 0)
            setError (status, "c_tbaptb:  error writing array to column");

        free (value);
}

void c_tbaptd (IRAFPointer tp, IRAFPointer cp, int row, double *buffer,
                int first, int nelem) {

/* Write an array of double values to a table column.
arguments:
IRAFPointer tp          i: table descriptor
IRAFPointer cp          i: column descriptor
int row                 i: row number (one indexed)
double *buffer          i: array of values to write to the table
int first               i: first element to write (one indexed)
int nelem               i: number of elements to write
*/

        TableDescr *tbl_descr;
        ColumnDescr *col_descr;
        long firstelem, one_elem=1;
        int i;
        int status = 0;

        tbl_descr = (TableDescr *)tp;
        col_descr = (ColumnDescr *)cp;

        /* fits_write_col_dbl = ffpcld */
        fits_write_col_dbl (tbl_descr->fptr, col_descr->colnum,
                (long)row, (long)first, (long)nelem, buffer, &status);
        for (i = 0;  i < nelem;  i++) {
            if (buffer[i] == IRAF_INDEFD) {
                firstelem = i + first;
                /* fits_write_col_null / ffpclu */
                fits_write_col_null (tbl_descr->fptr, col_descr->colnum,
                        (long)row, firstelem, one_elem, &status);
            }
        }
        if (status != 0)
            setError (status, "c_tbaptd:  error writing array to column");
}

void c_tbaptr (IRAFPointer tp, IRAFPointer cp, int row, float *buffer,
                int first, int nelem) {

/* Write an array of float values to a table column.
arguments:
IRAFPointer tp          i: table descriptor
IRAFPointer cp          i: column descriptor
int row                 i: row number (one indexed)
float *buffer           i: array of values to write to the table
int first               i: first element to write (one indexed)
int nelem               i: number of elements to write
*/

        TableDescr *tbl_descr;
        ColumnDescr *col_descr;
        long firstelem, one_elem=1;
        int i;
        int status = 0;

        tbl_descr = (TableDescr *)tp;
        col_descr = (ColumnDescr *)cp;

        /* fits_write_col_flt = ffpcle */
        fits_write_col_flt (tbl_descr->fptr, col_descr->colnum,
                (long)row, (long)first, (long)nelem, buffer, &status);
        for (i = 0;  i < nelem;  i++) {
            if (buffer[i] >= 0.99999 * IRAF_INDEFR &&
                buffer[i] <= 1.00001 * IRAF_INDEFR) {
                firstelem = i + first;
                fits_write_col_null (tbl_descr->fptr, col_descr->colnum,
                        (long)row, firstelem, one_elem, &status);
            }
        }
        if (status != 0)
            setError (status, "c_tbaptr:  error writing array to column");
}

void c_tbapti (IRAFPointer tp, IRAFPointer cp, int row, int *buffer,
                int first, int nelem) {

/* Write an array of integer values to a table column.
arguments:
IRAFPointer tp          i: table descriptor
IRAFPointer cp          i: column descriptor
int row                 i: row number (one indexed)
int *buffer             i: array of values to write to the table
int first               i: first element to write (one indexed)
int nelem               i: number of elements to write
*/

        TableDescr *tbl_descr;
        ColumnDescr *col_descr;
        long firstelem, one_elem=1;
        int i;
        int status = 0;

        tbl_descr = (TableDescr *)tp;
        col_descr = (ColumnDescr *)cp;

        /* fits_write_col_int = ffpclk */
        fits_write_col_int (tbl_descr->fptr, col_descr->colnum,
                (long)row, (long)first, (long)nelem, buffer, &status);
        for (i = 0;  i < nelem;  i++) {
            if (buffer[i] == IRAF_INDEFI) {
                firstelem = i + first;
                fits_write_col_null (tbl_descr->fptr, col_descr->colnum,
                        (long)row, firstelem, one_elem, &status);
            }
        }
        if (status != 0)
            setError (status, "c_tbapti:  error writing array to column");
}

void c_tbapts (IRAFPointer tp, IRAFPointer cp, int row, short *buffer,
                int first, int nelem) {

/* Write an array of short integer values to a table column.
arguments:
IRAFPointer tp          i: table descriptor
IRAFPointer cp          i: column descriptor
int row                 i: row number (one indexed)
short *buffer           i: array of values to write to the table
int first               i: first element to write (one indexed)
int nelem               i: number of elements to write
*/

        TableDescr *tbl_descr;
        ColumnDescr *col_descr;
        long firstelem, one_elem=1;
        int i;
        int status = 0;

        tbl_descr = (TableDescr *)tp;
        col_descr = (ColumnDescr *)cp;

        /* fits_write_col_sht = ffpcli */
        fits_write_col_sht (tbl_descr->fptr, col_descr->colnum,
                (long)row, (long)first, (long)nelem, buffer, &status);
        for (i = 0;  i < nelem;  i++) {
            if (buffer[i] == IRAF_INDEFS) {
                firstelem = i + first;
                fits_write_col_null (tbl_descr->fptr, col_descr->colnum,
                        (long)row, firstelem, one_elem, &status);
            }
        }
        if (status != 0)
            setError (status, "c_tbapts:  error writing array to column");
}

void c_tbaptt (IRAFPointer tp, IRAFPointer cp, int row, char **cbuf,
                int maxch, int first, int nelem) {

/* Write an array of text strings to a table column.
arguments:
IRAFPointer tp          i: table descriptor
IRAFPointer cp          i: column descriptor
int row                 i: row number (one indexed)
char **cbuf             i: array of strings to write to the table
int maxch               i: length of each string (currently ignored)
int first               i: first element to write (one indexed)
int nelem               i: number of elements to write
*/

        TableDescr *tbl_descr;
        ColumnDescr *col_descr;
        int status = 0;

        tbl_descr = (TableDescr *)tp;
        col_descr = (ColumnDescr *)cp;

        /* fits_write_col_str = ffpcls */
        fits_write_col_str (tbl_descr->fptr, col_descr->colnum,
                (long)row, (long)first, (long)nelem, cbuf, &status);
        if (status != 0)
            setError (status, "c_tbaptt:  error writing array to column");
}
