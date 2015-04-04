# include <ctype.h>
# include <string.h>
# include <fitsio.h>
# include "ctables.h"

/* Data type conversion is supported for most cases.  The following
   are NOT supported:

function c_tbeptb:  table columns of type double, float, int or short
function c_tbeptd:  table column of type boolean
function c_tbeptr:  table column of type boolean
*/

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
        char s_value[11]={'\0'};
        int status = 0;

        tbl_descr = (TableDescr *)tp;
        col_descr = (ColumnDescr *)cp;

        if (col_descr->datatype < 0) {
            char cbuf[SZ_FITS_STR+1];
            if (buffer == True)
                strcpy (cbuf, "yes");
            else
                strcpy (cbuf, "no");
            c_tbeptt (tp, cp, row, cbuf);

        } else {

            if (buffer == True)
                s_value[0] = 1;
            else
                s_value[0] = 0;

            /* fits_write_col_log = ffpcll */
            fits_write_col_log (tbl_descr->fptr, col_descr->colnum,
                        (long)row, firstelem, nelem, s_value, &status);
        }
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

        if (col_descr->datatype < 0) {
            if (buffer == IRAF_INDEFD) {
                c_tbeptt (tp, cp, row, "INDEF");
            } else {
                char cbuf[SZ_FITS_STR+1];
                sprintf (cbuf, "%-25.16g", buffer);
                c_tbeptt (tp, cp, row, cbuf);
            }

        } else {

            if (buffer == IRAF_INDEFD) {
                /* fits_write_col_null / ffpclu */
                fits_write_col_null (tbl_descr->fptr, col_descr->colnum,
                        (long)row, firstelem, nelem, &status);
            } else {
                /* fits_write_col_dbl = ffpcld */
                fits_write_col_dbl (tbl_descr->fptr, col_descr->colnum,
                        (long)row, firstelem, nelem, &buffer, &status);
            }
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

        if (col_descr->datatype < 0) {
            if (buffer >= 0.99999 * IRAF_INDEFR &&
                buffer <= 1.00001 * IRAF_INDEFR) {
                c_tbeptt (tp, cp, row, "INDEF");
            } else {
                char cbuf[SZ_FITS_STR+1];
                sprintf (cbuf, "%-15.7g", buffer);
                c_tbeptt (tp, cp, row, cbuf);
            }

        } else {

            if (buffer >= 0.99999 * IRAF_INDEFR &&
                buffer <= 1.00001 * IRAF_INDEFR) {
                fits_write_col_null (tbl_descr->fptr, col_descr->colnum,
                        (long)row, firstelem, nelem, &status);
            } else {
                /* fits_write_col_flt = ffpcle */
                fits_write_col_flt (tbl_descr->fptr, col_descr->colnum,
                        (long)row, firstelem, nelem, &buffer, &status);
            }
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

        if (col_descr->datatype < 0) {
            if (buffer == IRAF_INDEFI) {
                c_tbeptt (tp, cp, row, "INDEF");
            } else {
                char cbuf[SZ_FITS_STR+1];
                sprintf (cbuf, "%d", buffer);
                c_tbeptt (tp, cp, row, cbuf);
            }

        } else if (buffer == IRAF_INDEFI) {
            fits_write_col_null (tbl_descr->fptr, col_descr->colnum,
                        (long)row, firstelem, nelem, &status);
        } else {

            if (col_descr->datatype == IRAF_BOOL) {
                Bool b_value;
                if (buffer)
                    b_value = True;
                else
                    b_value = False;
                c_tbeptb (tp, cp, row, b_value);

            } else if (col_descr->datatype == IRAF_SHORT) {
                c_tbepts (tp, cp, row, (short)buffer);

            } else {
                /* fits_write_col_int = ffpclk */
                fits_write_col_int (tbl_descr->fptr, col_descr->colnum,
                        (long)row, firstelem, nelem, &buffer, &status);
            }
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

        if (col_descr->datatype < 0) {
            if (buffer == IRAF_INDEFS) {
                c_tbeptt (tp, cp, row, "INDEF");
            } else {
                char cbuf[SZ_FITS_STR+1];
                sprintf (cbuf, "%hd", buffer);
                c_tbeptt (tp, cp, row, cbuf);
            }

        } else if (buffer == IRAF_INDEFS) {
            fits_write_col_null (tbl_descr->fptr, col_descr->colnum,
                        (long)row, firstelem, nelem, &status);
        } else {

            if (col_descr->datatype == IRAF_BOOL) {
                Bool b_value;
                if (buffer)
                    b_value = True;
                else
                    b_value = False;
                c_tbeptb (tp, cp, row, b_value);

            } else if (col_descr->datatype == IRAF_INT) {
                c_tbepti (tp, cp, row, (int)buffer);

            } else {
                /* fits_write_col_sht = ffpcli */
                fits_write_col_sht (tbl_descr->fptr, col_descr->colnum,
                        (long)row, firstelem, nelem, &buffer, &status);
            }
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
        char *p;
        double d_value;
        int status = 0;

        tbl_descr = (TableDescr *)tp;
        col_descr = (ColumnDescr *)cp;

        if (col_descr->datatype == IRAF_BOOL) {
            char value[SZ_FITS_STR+1];
            int i, j;
            for (i = 0, j = 0;  buffer[i] != '\0';  i++) {
                if (buffer[i] == ' ')
                    continue;
                if (j >= SZ_FITS_STR)
                    break;
                if (isupper (buffer[i]))
                    value[j] = tolower (buffer[i]);
                else
                    value[j] = buffer[i];
                j++;
            }
            value[j] = '\0';
            if (strcmp (value, "1") == 0 ||
                strcmp (value, "yes") == 0 ||
                strcmp (value, "true") == 0) {
                c_tbeptb (tp, cp, row, True);
            } else if (strcmp (value, "0") == 0 ||
                       strcmp (value, "no") == 0 ||
                       strcmp (value, "false") == 0) {
                c_tbeptb (tp, cp, row, False);
            } else {
                /* neither True nor False, so set the value to undefined */
                c_tbrudf (tp, &cp, 1, row);
            }
        } else if (col_descr->datatype == IRAF_DOUBLE) {
            d_value = strtod (buffer, &p);
            if (*p == '\0')
                c_tbeptd (tp, cp, row, d_value);
            else
                c_tbrudf (tp, &cp, 1, row);
        } else if (col_descr->datatype == IRAF_REAL) {
            d_value = strtod (buffer, &p);
            if (*p == '\0')
                c_tbeptr (tp, cp, row, (float)d_value);
            else
                c_tbrudf (tp, &cp, 1, row);
        } else if (col_descr->datatype == IRAF_INT) {
            int i_value;
            d_value = strtod (buffer, &p);
            i_value = (int)d_value;
            if (*p == '\0' && i_value == d_value)
                c_tbepti (tp, cp, row, i_value);
            else
                c_tbrudf (tp, &cp, 1, row);
        } else if (col_descr->datatype == IRAF_SHORT) {
            short si_value;
            d_value = strtod (buffer, &p);
            si_value = (short)d_value;
            if (*p == '\0' && si_value == d_value)
                c_tbepts (tp, cp, row, si_value);
            else
                c_tbrudf (tp, &cp, 1, row);
        } else {

            /* fits_write_col_str = ffpcls */
            fits_write_col_str (tbl_descr->fptr, col_descr->colnum,
                        (long)row, firstelem, nelem, &buffer, &status);
        }
        if (status != 0)
            setError (status, "c_tbeptt:  error writing element");
}
