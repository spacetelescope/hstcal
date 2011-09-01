# include <stdlib.h>
# include <string.h>
# include <fitsio.h>
# include "ctables.h"

void tbSaveInfo (IRAFPointer tp, int *status) {

/* This is called by c_tbtopn to save information about a table in
   the table descriptor.
arguments:
IRAFPointer tp          i: table descriptor
int *status             o: 0 is OK
*/

        fitsfile *fptr;
        TableDescr *tbl_descr;
        ColumnDescr *col_descr;
        IRAFPointer cp;
        char value[SZ_FITS_STR+1];      /* column name, units, print format */
        char comment[SZ_FITS_STR+1];    /* comment field in a FITS keyword */
        char keyword[SZ_FITS_STR+1];    /* keyword name */
        int hdutype, ncols;
        int typecode, datatype;         /* CFITSIO, IRAF data type codes */
        long nrows;
        long repeat;            /* number of elements in column array */
        long offset;
        long width;
        char *errmsg;
        int i;

        *status = 0;

        tbl_descr = (TableDescr *)tp;
        fptr = tbl_descr->fptr;

        /* fits_get_hdu_type = ffghdt */
        fits_get_hdu_type (fptr, &hdutype, status);
        tbl_descr->hdutype = hdutype;
        if (hdutype == IMAGE_HDU) {
            tbl_descr->nrows = 0;
            tbl_descr->ncols = 0;
            return;
        } else if (hdutype != ASCII_TBL && hdutype != BINARY_TBL) {
            *status = ERR_NOT_A_TABLE;
            return;
        }

        /* fits_get_num_rows = ffgnrw */
        fits_get_num_rows (fptr, &nrows, status);
        tbl_descr->nrows = nrows;

        /* fits_get_num_cols = ffgncl*/
        fits_get_num_cols (fptr, &ncols, status);
        if (*status != 0)
            return;

        /* make sure we have enough space for the column pointers */
        /* (don't update tbl_descr->ncols yet because columnSpace uses it) */
        columnSpace (tbl_descr, ncols);
        if (checkError() != 0)
            return;

        tbl_descr->ncols = ncols;

        /* Get and save info for each column in the table. */
        for (i = 0;  i < ncols;  i++) {

            cp = init_cp (tp);
            col_descr = (ColumnDescr *)cp;

            /* get column name */
            sprintf (keyword, "TTYPE%d", i+1);  /* one indexed column number */
            /* fits_read_key = ffgky */
            fits_read_key (fptr, TSTRING, keyword, value, comment, status);
            if (*status == 0) {
                strcpy (col_descr->name, value);
            } else {
                /* default name is "c" with one indexed column number */
                sprintf (col_descr->name, "c%d", i+1);
                /* fits_clear_errmsg = ffcmsg */
                fits_clear_errmsg();
                *status = 0;
            }

            /* get column data type as a string */
            sprintf (keyword, "TFORM%d", i+1);  /* keyword for data type */
            fits_read_key (fptr, TSTRING, keyword, value, comment, status);
            if (*status == 0) {
                strcpy (col_descr->tform, value);
            } else {
                col_descr->tform[0] = '\0';
                fits_clear_errmsg();
                *status = 0;
            }

            /* get column units */
            sprintf (keyword, "TUNIT%d", i+1);  /* keyword for units */
            fits_read_key (fptr, TSTRING, keyword, value, comment, status);
            if (*status == 0) {
                strcpy (col_descr->tunit, value);
            } else {
                col_descr->tunit[0] = '\0';
                fits_clear_errmsg();
                *status = 0;
            }

            /* get display format */
            sprintf (keyword, "TDISP%d", i+1);  /* keyword for print format */
            fits_read_key (fptr, TSTRING, keyword, value, comment, status);
            if (*status == 0) {
                strcpy (col_descr->tdisp, value);
            } else {
                col_descr->tdisp[0] = '\0';
                fits_clear_errmsg();
                *status = 0;
            }

            /* fits_get_eqcoltype = ffeqty */
            fits_get_eqcoltype (fptr, i+1, &typecode, &repeat, &width, status);
            col_descr->colnum = i + 1;          /* one indexed */
            col_descr->typecode = typecode;
            if (typecode < 0) {
                typecode = -typecode;           /* only change local var. */
                col_descr->var_length = 1;      /* true, variable length */
                /* fits_read_descript / ffgdes */
                fits_read_descript (fptr, i+1, 1, &repeat, &offset, status);
            } else {
                col_descr->var_length = 0;      /* false */
            }
            col_descr->repeat = repeat;
            if (typecode == TSTRING)
                col_descr->nelem = repeat / width;
            else
                col_descr->nelem = repeat;
            col_descr->width = width;

            switch (typecode) {
            case TFLOAT:
                datatype = IRAF_REAL;
                break;
            case TDOUBLE:
                datatype = IRAF_DOUBLE;
                break;
            case TUSHORT:
                datatype = IRAF_USHORT;
                break;
            case TBYTE:
            case TSHORT:
                datatype = IRAF_SHORT;
                break;
            case TINT:
            case TLONG:
            case TLONGLONG:
                datatype = IRAF_INT;
                break;
            case TLOGICAL:
                datatype = IRAF_BOOL;
                break;
            case TCOMPLEX:
            case TDBLCOMPLEX:
                datatype = IRAF_COMPLEX;
                break;
            case TSTRING:
                /* the data type is the negative of the max string length */
                datatype = -col_descr->width;
                break;
            default:
                datatype = 0;
                errmsg = (char *)calloc (SZ_ERRMESS+1, sizeof(char));
                sprintf (errmsg,
                        "c_tbtopn:  data type code %d not supported.",
                        col_descr->typecode);
                setError (ERR_DATATYPE_UNKNOWN, errmsg);
                free (errmsg);
            }
            col_descr->datatype = datatype;

            tbl_descr->columns[i] = cp;
        }

        return;
}
