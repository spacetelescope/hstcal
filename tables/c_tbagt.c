# include <stdlib.h>
# include <string.h>
# include <fitsio.h>
# include "ctables.h"

int c_tbagtb (IRAFPointer tp, IRAFPointer cp, int row, Bool *buffer,
                int first, int nelem) {

/* Read an array of boolean values from a table column.
arguments:
IRAFPointer tp          i: table descriptor
IRAFPointer cp          i: column descriptor
int row                 i: row number (one indexed)
Bool *buffer            o: array of values (True or False) read from table
int first               i: first element to read (one indexed)
int nelem               i: number of elements to read

function value          o: actual number of elements that were read
*/

        TableDescr *tbl_descr;
        ColumnDescr *col_descr;
        long nret;              /* actual number of elements to get */
        int anynul=0;
        int nulval=0;
        char *value;
        int i;
        int status = 0;

        tbl_descr = (TableDescr *)tp;
        col_descr = (ColumnDescr *)cp;

        if (col_descr->var_length) {
            nret = nelem;
        } else {
            if (nelem <= col_descr->nelem-first+1)
                nret = nelem;
            else
                nret = col_descr->nelem - first + 1;
        }
        value = (char *)calloc (nret+10, sizeof(char)); /* nret plus extra */

        /* fits_read_col_log = ffgcvl */
        fits_read_col_log (tbl_descr->fptr, col_descr->colnum,
                (long)row, (long)first, nret, nulval,
                value, &anynul, &status);
        if (status != 0)
            setError (status, "c_tbagtb:  error reading array from column");
        for (i = 0;  i < nret;  i++) {
            if (value[i])
                buffer[i] = True;
            else
                buffer[i] = False;
        }
        free (value);

        return (int)nret;
}

int c_tbagtd (IRAFPointer tp, IRAFPointer cp, int row, double *buffer,
                int first, int nelem) {

/* Read an array of double values from a table column.
arguments:
IRAFPointer tp          i: table descriptor
IRAFPointer cp          i: column descriptor
int row                 i: row number (one indexed)
double *buffer          o: array of values read from table
int first               i: first element to read (one indexed)
int nelem               i: number of elements to read

function value          o: actual number of elements that were read
*/

        TableDescr *tbl_descr;
        ColumnDescr *col_descr;
        long nret;              /* actual number of elements to get */
        int anynul=0;
        double nulval=IRAF_INDEFD;
        int status = 0;

        tbl_descr = (TableDescr *)tp;
        col_descr = (ColumnDescr *)cp;

        if (col_descr->var_length) {
            nret = nelem;
        } else {
            if (nelem <= col_descr->nelem-first+1)
                nret = nelem;
            else
                nret = col_descr->nelem - first + 1;
        }

        /* fits_read_col_dbl = ffgcvd */
        fits_read_col_dbl (tbl_descr->fptr, col_descr->colnum,
                (long)row, (long)first, nret, nulval,
                buffer, &anynul, &status);
        if (status != 0)
            setError (status, "c_tbagtd:  error reading array from column");

        return (int)nret;
}

int c_tbagtr (IRAFPointer tp, IRAFPointer cp, int row, float *buffer,
                int first, int nelem) {

/* Read an array of float values from a table column.
arguments:
IRAFPointer tp          i: table descriptor
IRAFPointer cp          i: column descriptor
int row                 i: row number (one indexed)
float *buffer           o: array of values read from table
int first               i: first element to read (one indexed)
int nelem               i: number of elements to read

function value          o: actual number of elements that were read
*/

        TableDescr *tbl_descr;
        ColumnDescr *col_descr;
        long nret;              /* actual number of elements to get */
        int anynul=0;
        float nulval=IRAF_INDEFR;
        int status = 0;

        tbl_descr = (TableDescr *)tp;
        col_descr = (ColumnDescr *)cp;

        if (col_descr->var_length) {
            nret = nelem;
        } else {
            if (nelem <= col_descr->nelem-first+1)
                nret = nelem;
            else
                nret = col_descr->nelem - first + 1;
        }

        /* fits_read_col_flt = ffgcve */
        fits_read_col_flt (tbl_descr->fptr, col_descr->colnum,
                (long)row, (long)first, nret, nulval,
                buffer, &anynul, &status);
        if (status != 0)
            setError (status, "c_tbagtr:  error reading array from column");

        return (int)nret;
}

int c_tbagti (IRAFPointer tp, IRAFPointer cp, int row, int *buffer,
                int first, int nelem) {

/* Read an array of integer values from a table column.
arguments:
IRAFPointer tp          i: table descriptor
IRAFPointer cp          i: column descriptor
int row                 i: row number (one indexed)
int *buffer             o: array of values read from table
int first               i: first element to read (one indexed)
int nelem               i: number of elements to read

function value          o: actual number of elements that were read
*/

        TableDescr *tbl_descr;
        ColumnDescr *col_descr;
        long nret;              /* actual number of elements to get */
        int anynul=0;
        int nulval=IRAF_INDEFI;
        int status = 0;

        tbl_descr = (TableDescr *)tp;
        col_descr = (ColumnDescr *)cp;

        if (col_descr->var_length) {
            nret = nelem;
        } else {
            if (nelem <= col_descr->nelem-first+1)
                nret = nelem;
            else
                nret = col_descr->nelem - first + 1;
        }

        /* fits_read_col_int = ffgcvk */
        fits_read_col_int (tbl_descr->fptr, col_descr->colnum,
                (long)row, (long)first, nret, nulval,
                buffer, &anynul, &status);
        if (status != 0)
            setError (status, "c_tbagti:  error reading array from column");

        return (int)nret;
}

int c_tbagts (IRAFPointer tp, IRAFPointer cp, int row, short *buffer,
                int first, int nelem) {

/* Read an array of short integer values from a table column.
arguments:
IRAFPointer tp          i: table descriptor
IRAFPointer cp          i: column descriptor
int row                 i: row number (one indexed)
short *buffer           o: array of values read from table
int first               i: first element to read (one indexed)
int nelem               i: number of elements to read

function value          o: actual number of elements that were read
*/

        TableDescr *tbl_descr;
        ColumnDescr *col_descr;
        long nret;              /* actual number of elements to get */
        int anynul=0;
        short nulval=IRAF_INDEFS;
        int status = 0;

        tbl_descr = (TableDescr *)tp;
        col_descr = (ColumnDescr *)cp;

        if (col_descr->var_length) {
            nret = nelem;
        } else {
            if (nelem <= col_descr->nelem-first+1)
                nret = nelem;
            else
                nret = col_descr->nelem - first + 1;
        }

        /* fits_read_col_sht = ffgcvi */
        fits_read_col_sht (tbl_descr->fptr, col_descr->colnum,
                (long)row, (long)first, nret, nulval,
                buffer, &anynul, &status);
        if (status != 0)
            setError (status, "c_tbagts:  error reading array from column");

        return (int)nret;
}

int c_tbagtt (IRAFPointer tp, IRAFPointer cp, int row, char **cbuf,
                int first, int nelem, int maxch) {

/* Read an array of text strings from a table column.
arguments:
IRAFPointer tp          i: table descriptor
IRAFPointer cp          i: column descriptor
int row                 i: row number (one indexed)
char **cbuf             o: array of strings read from table
int first               i: first element to read (one indexed)
int nelem               i: number of elements to read
int maxch               i: maximum length of each string (not incl NULL)

function value          o: actual number of elements that were read
*/

        TableDescr *tbl_descr;
        ColumnDescr *col_descr;
        long nret;              /* actual number of elements to get */
        long one_element=1;
        long i, j;
        int anynul=0;
        char *value;
        int len;
        int status = 0;

        tbl_descr = (TableDescr *)tp;
        col_descr = (ColumnDescr *)cp;

        if (col_descr->var_length) {
            nret = nelem;
        } else {
            if (nelem <= col_descr->nelem-first+1)
                nret = nelem;
            else
                nret = col_descr->nelem - first + 1;
        }

        if (col_descr->width >= maxch)
            len = col_descr->width;
        else
            len = maxch;
        /* add 5 for extra space */
        value = (char *)calloc (len+5, sizeof(char));

        for (i = first, j = 0;  j < nret;  i++, j++) {
            /* fits_read_col_str = ffgcvs */
            fits_read_col_str (tbl_descr->fptr, col_descr->colnum,
                        (long)row, i, one_element, "INDEF",
                        &value, &anynul, &status);
            copyString (cbuf[j], value, maxch);
        }
        if (status != 0)
            setError (status, "c_tbagtt:  error reading array from column");

        return (int)nret;
}
