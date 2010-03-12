# include <stdlib.h>
# include <string.h>
# include <fitsio.h>
# include "ctables.h"

void c_tbegtb (IRAFPointer tp, IRAFPointer cp, int row, Bool *buffer) {

/* Read a boolean value from a table column.
arguments:
IRAFPointer tp          i: table descriptor
IRAFPointer cp          i: column descriptor
int row                 i: row number (one indexed)
Bool *buffer            o: value (True or False) read from table
*/

        TableDescr *tbl_descr;
        ColumnDescr *col_descr;
        int anynul=0;
        long firstelem=1, nelem=1;
        int nulval=0;
        char value[11]={'\0'};
        int status = 0;

        tbl_descr = (TableDescr *)tp;
        col_descr = (ColumnDescr *)cp;

        /* fits_read_col_log = ffgcvl */
        fits_read_col_log (tbl_descr->fptr, col_descr->colnum,
                (long)row, firstelem, nelem, nulval,
                value, &anynul, &status);
        if (status != 0)
            setError (status, "c_tbegtb:  error reading element");
        if (value[0])
            *buffer = True;
        else
            *buffer = False;
}

void c_tbegtd (IRAFPointer tp, IRAFPointer cp, int row, double *buffer) {

/* Read a value of type double from a table column.
arguments:
IRAFPointer tp          i: table descriptor
IRAFPointer cp          i: column descriptor
int row                 i: row number (one indexed)
double *buffer          o: value read from table
*/

        TableDescr *tbl_descr;
        ColumnDescr *col_descr;
        int anynul=0;
        long firstelem=1, nelem=1;
        double nulval=IRAF_INDEFD;
        int status = 0;

        tbl_descr = (TableDescr *)tp;
        col_descr = (ColumnDescr *)cp;

        /* fits_read_col_dbl = ffgcvd */
        fits_read_col_dbl (tbl_descr->fptr, col_descr->colnum,
                (long)row, firstelem, nelem, nulval,
                buffer, &anynul, &status);
        if (status != 0)
            setError (status, "c_tbegtd:  error reading element");
}

void c_tbegtr (IRAFPointer tp, IRAFPointer cp, int row, float *buffer) {

/* Read a value of type float from a table column.
arguments:
IRAFPointer tp          i: table descriptor
IRAFPointer cp          i: column descriptor
int row                 i: row number (one indexed)
float *buffer           o: value read from table
*/

        TableDescr *tbl_descr;
        ColumnDescr *col_descr;
        int anynul=0;
        long firstelem=1, nelem=1;
        float nulval=IRAF_INDEFR;
        int status = 0;

        tbl_descr = (TableDescr *)tp;
        col_descr = (ColumnDescr *)cp;

        /* fits_read_col_flt = ffgcve */
        fits_read_col_flt (tbl_descr->fptr, col_descr->colnum,
                (long)row, firstelem, nelem, nulval,
                buffer, &anynul, &status);
        if (status != 0)
            setError (status, "c_tbegtr:  error reading element");
}

void c_tbegti (IRAFPointer tp, IRAFPointer cp, int row, int *buffer) {

/* Read an integer value from a table column.
arguments:
IRAFPointer tp          i: table descriptor
IRAFPointer cp          i: column descriptor
int row                 i: row number (one indexed)
int *buffer             o: value read from table
*/

        TableDescr *tbl_descr;
        ColumnDescr *col_descr;
        int anynul=0;
        long firstelem=1, nelem=1;
        int nulval=IRAF_INDEFI;
        int status = 0;

        tbl_descr = (TableDescr *)tp;
        col_descr = (ColumnDescr *)cp;

        /* fits_read_col_int = ffgcvk */
        fits_read_col_int (tbl_descr->fptr, col_descr->colnum,
                (long)row, firstelem, nelem, nulval,
                buffer, &anynul, &status);
        if (status != 0)
            setError (status, "c_tbegti:  error reading element");
}

void c_tbegts (IRAFPointer tp, IRAFPointer cp, int row, short *buffer) {

/* Read a short integer value from a table column.
arguments:
IRAFPointer tp          i: table descriptor
IRAFPointer cp          i: column descriptor
int row                 i: row number (one indexed)
short *buffer           o: value read from table
*/

        TableDescr *tbl_descr;
        ColumnDescr *col_descr;
        int anynul=0;
        long firstelem=1, nelem=1;
        short nulval=IRAF_INDEFS;
        int status = 0;

        tbl_descr = (TableDescr *)tp;
        col_descr = (ColumnDescr *)cp;

        /* fits_read_col_sht = ffgcvi */
        fits_read_col_sht (tbl_descr->fptr, col_descr->colnum,
                (long)row, firstelem, nelem, nulval,
                buffer, &anynul, &status);
        if (status != 0)
            setError (status, "c_tbegts:  error reading element");
}

void c_tbegtt (IRAFPointer tp, IRAFPointer cp, int row, char *buffer,
                int maxch) {

/* Read a string value from a table column.
arguments:
IRAFPointer tp          i: table descriptor
IRAFPointer cp          i: column descriptor
int row                 i: row number (one indexed)
char *buffer            o: value read from table
int maxch               i: maximum length of the string (not incl NULL)
*/

        TableDescr *tbl_descr;
        ColumnDescr *col_descr;
        int anynul=0;
        long firstelem=1, nelem=1;
        char *value;
        int len;
        int status = 0;

        tbl_descr = (TableDescr *)tp;
        col_descr = (ColumnDescr *)cp;

        if (col_descr->width >= maxch)
            len = col_descr->width;
        else
            len = maxch;
        /* add 5 for extra space */
        value = (char *)calloc (len+5, sizeof(char));

        /* fits_read_col_str = ffgcvs */
        fits_read_col_str (tbl_descr->fptr, col_descr->colnum,
                (long)row, firstelem, nelem, "INDEF",
                &value, &anynul, &status);
        if (status != 0)
            setError (status, "c_tbegtt:  error reading element");

        copyString (buffer, value, maxch);
        free (value);
}
