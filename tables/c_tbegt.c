# include <stdlib.h>
# include <string.h>
# include <fitsio.h>
# include "ctables.h"

/* Data type conversion is supported for most cases.  The following
   are NOT supported:

function c_tbegtb:  table columns of type double or float
function c_tbegtd:  table column of type boolean
function c_tbegtr:  table column of type boolean
*/

void c_tbegtb (const IRAFPointer tp, const IRAFPointer cp, int row, Bool *buffer) {

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
        char s_value[11]={'\0'};
        int status = 0;

        tbl_descr = (TableDescr *)tp;
        col_descr = (ColumnDescr *)cp;

        if (col_descr->datatype < 0) {
            char *value;
            int maxch;
            maxch = col_descr->width + 5;
            value = (char *)calloc (maxch, sizeof(char));
            c_tbegtt (tp, cp, row, value, maxch);
            str_lower (value, value);
            if (strcmp (value, "1") == 0) {
                *buffer = True;
            } else if (strcmp (value, "0") == 0) {
                *buffer = False;
            } else if (strncmp (value, "yes", 3) == 0 ||
                       strncmp (value, "true", 4) == 0) {
                *buffer = True;
            } else {
                *buffer = False;
            }
            free (value);

        } else if (col_descr->datatype == IRAF_INT) {
            int i_value;
            c_tbegti (tp, cp, row, &i_value);
            if (i_value == IRAF_INDEFI)
                *buffer = False;
            else if (i_value)
                *buffer = True;
            else
                *buffer = False;

        } else if (col_descr->datatype == IRAF_SHORT) {
            short si_value;
            c_tbegts (tp, cp, row, &si_value);
            if (si_value == IRAF_INDEFS)
                *buffer = False;
            else if (si_value)
                *buffer = True;
            else
                *buffer = False;

        } else {

            /* fits_read_col_log = ffgcvl */
            fits_read_col_log (tbl_descr->fptr, col_descr->colnum,
                (long)row, firstelem, nelem, nulval,
                s_value, &anynul, &status);
            if (status != 0)
                setError (status, "c_tbegtb:  error reading element");
            if (s_value[0])
                *buffer = True;
            else
                *buffer = False;
        }
}

void c_tbegtd (const IRAFPointer tp, const IRAFPointer cp, int row, double *buffer) {

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

        if (col_descr->datatype < 0) {
            char *value;
            int maxch;
            maxch = col_descr->width + 5;
            value = (char *)calloc (maxch, sizeof(char));
            c_tbegtt (tp, cp, row, value, maxch);
            if (strcmp (value, "INDEF") == 0)
                *buffer = IRAF_INDEFD;
            else
                *buffer = atof (value);
            free (value);

        } else if (col_descr->datatype == IRAF_INT) {
            int i_value;
            c_tbegti (tp, cp, row, &i_value);
            if (i_value == IRAF_INDEFI)
                *buffer = IRAF_INDEFD;
            else
                *buffer = i_value;

        } else if (col_descr->datatype == IRAF_SHORT) {
            short si_value;
            c_tbegts (tp, cp, row, &si_value);
            if (si_value == IRAF_INDEFS)
                *buffer = IRAF_INDEFD;
            else
                *buffer = si_value;

        } else {

            /* fits_read_col_dbl = ffgcvd */
            fits_read_col_dbl (tbl_descr->fptr, col_descr->colnum,
                (long)row, firstelem, nelem, nulval,
                buffer, &anynul, &status);
            if (status != 0)
                setError (status, "c_tbegtd:  error reading element");
        }
}

void c_tbegtr (const IRAFPointer tp, const IRAFPointer cp, int row, float *buffer) {

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

        if (col_descr->datatype < 0) {
            char *value;
            int maxch;
            maxch = col_descr->width + 5;
            value = (char *)calloc (maxch, sizeof(char));
            c_tbegtt (tp, cp, row, value, maxch);
            if (strcmp (value, "INDEF") == 0)
                *buffer = IRAF_INDEFR;
            else
                *buffer = atof (value);
            free (value);

        } else if (col_descr->datatype == IRAF_INT) {
            int i_value;
            c_tbegti (tp, cp, row, &i_value);
            if (i_value == IRAF_INDEFI)
                *buffer = IRAF_INDEFR;
            else
                *buffer = i_value;

        } else if (col_descr->datatype == IRAF_SHORT) {
            short si_value;
            c_tbegts (tp, cp, row, &si_value);
            if (si_value == IRAF_INDEFS)
                *buffer = IRAF_INDEFR;
            else
                *buffer = si_value;

        } else {

            /* fits_read_col_flt = ffgcve */
            fits_read_col_flt (tbl_descr->fptr, col_descr->colnum,
                (long)row, firstelem, nelem, nulval,
                buffer, &anynul, &status);
            if (status != 0)
                setError (status, "c_tbegtr:  error reading element");
        }
}

void c_tbegti (const IRAFPointer tp, const IRAFPointer cp, int row, int *buffer) {

/* Read an integer value from a table column.
arguments:
IRAFPointer tp          i: table descriptor
IRAFPointer cp          i: column descriptor
int row                 i: row number (one indexed)
int *buffer             o: value read from table
*/
        int anynul=0;
        long firstelem=1, nelem=1;
        int nulval=IRAF_INDEFI;
        int status = 0;

        TableDescr *tbl_descr = (TableDescr *)tp;
        ColumnDescr *col_descr = (ColumnDescr *)cp;

        if (col_descr->datatype < 0) {
            char *value;
            int maxch;
            maxch = col_descr->width + 5;
            value = (char *)calloc (maxch, sizeof(char));
            c_tbegtt (tp, cp, row, value, maxch);
            if (strcmp (value, "INDEF") == 0)
                *buffer = IRAF_INDEFI;
            else
                *buffer = atoi (value);
            free (value);

        } else if (col_descr->datatype == IRAF_BOOL) {
            Bool b_value;
            c_tbegtb (tp, cp, row, &b_value);
            if (b_value)
                *buffer = 1;
            else
                *buffer = 0;

        } else if (col_descr->datatype == IRAF_SHORT) {
            short si_value;
            c_tbegts (tp, cp, row, &si_value);
            if (si_value == IRAF_INDEFS)
                *buffer = IRAF_INDEFI;
            else
                *buffer = si_value;

        } else {

            /* fits_read_col_int = ffgcvk */
            fits_read_col_int (tbl_descr->fptr, col_descr->colnum,
                (long)row, firstelem, nelem, nulval,
                buffer, &anynul, &status);
            if (status != 0)
                setError (status, "c_tbegti:  error reading element");
        }
}

int c_tbeGetInt (const IRAFPointer tp, const IRAFPointer cp, int row)
{
	int ret = 0;
	c_tbegti(tp, cp, row, &ret);
	return ret;
}
double c_tbeGetDouble (const IRAFPointer tp, const IRAFPointer cp, int row)
{
	double ret = 0;
	c_tbegtd(tp, cp, row, &ret);
	return ret;
}

void c_tbegts (const IRAFPointer tp, const IRAFPointer cp, int row, short *buffer) {

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

        if (col_descr->datatype < 0) {
            char *value;
            int maxch;
            maxch = col_descr->width + 5;
            value = (char *)calloc (maxch, sizeof(char));
            c_tbegtt (tp, cp, row, value, maxch);
            if (strcmp (value, "INDEF") == 0)
                *buffer = IRAF_INDEFS;
            else
                *buffer = atoi (value);
            free (value);

        } else if (col_descr->datatype == IRAF_BOOL) {
            Bool b_value;
            c_tbegtb (tp, cp, row, &b_value);
            if (b_value)
                *buffer = 1;
            else
                *buffer = 0;

        } else if (col_descr->datatype == IRAF_INT) {
            int i_value;
            c_tbegti (tp, cp, row, &i_value);
            if (i_value == IRAF_INDEFI)
                *buffer = IRAF_INDEFS;
            else
                *buffer = (short)i_value;

        } else {

            /* fits_read_col_sht = ffgcvi */
            fits_read_col_sht (tbl_descr->fptr, col_descr->colnum,
                (long)row, firstelem, nelem, nulval,
                buffer, &anynul, &status);
            if (status != 0)
                setError (status, "c_tbegts:  error reading element");
        }
}

void c_tbegtt (const IRAFPointer tp, const IRAFPointer cp, int row, char *buffer,
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

        if (col_descr->datatype < 0) {
            if (col_descr->width >= maxch)
                len = col_descr->width;
            else
                len = maxch;
        } else {
            len = SZ_FITS_STR;
        }
        value = (char *)calloc (len+1, sizeof(char));

        if (col_descr->datatype == IRAF_BOOL) {
            Bool b_value;
            c_tbegtb (tp, cp, row, &b_value);
            if (b_value == True)
                strcpy (value, "yes");
            else
                strcpy (value, "no");

        } else if (col_descr->datatype == IRAF_DOUBLE) {
            double d_value;
            c_tbegtd (tp, cp, row, &d_value);
            if (d_value == IRAF_INDEFD)
                strcpy (value, "INDEF");
            else
                sprintf (value, "%.16g", d_value);

        } else if (col_descr->datatype == IRAF_REAL) {
            float r_value;
            c_tbegtr (tp, cp, row, &r_value);
            if (r_value >= 0.99999 * IRAF_INDEFR &&
                r_value <= 1.00001 * IRAF_INDEFR) {
                strcpy (value, "INDEF");
            } else {
                sprintf (value, "%.7g", r_value);
            }

        } else if (col_descr->datatype == IRAF_INT) {
            int i_value;
            c_tbegti (tp, cp, row, &i_value);
            if (i_value == IRAF_INDEFI)
                strcpy (value, "INDEF");
            else
                sprintf (value, "%d", i_value);

        } else if (col_descr->datatype == IRAF_SHORT) {
            short si_value;
            c_tbegts (tp, cp, row, &si_value);
            if (si_value == IRAF_INDEFS)
                strcpy (value, "INDEF");
            else
                sprintf (value, "%hd", si_value);

        } else {

            /* fits_read_col_str = ffgcvs */
            fits_read_col_str (tbl_descr->fptr, col_descr->colnum,
                (long)row, firstelem, nelem, "INDEF",
                &value, &anynul, &status);
            if (status != 0)
                setError (status, "c_tbegtt:  error reading element");
        }

        copyString (buffer, value, maxch);

        free (value);
}
