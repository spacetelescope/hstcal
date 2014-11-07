# include <stdlib.h>
# include <string.h>
# include <fitsio.h>
# include "ctables.h"

void initCol (IRAFPointer tp, IRAFPointer *cp,
        char *colname, char *colunits, char *colfmt, int datatype, int nelem) {

/* This is called by c_tbcdef1.  A column descriptor will be allocated, and
   values will be copied to the descriptor.  The column descriptor will also
   be saved in the table descriptor.
   If the print format colfmt is not specified (i.e. colfmt = ""), a default
   value will be assigned.
arguments:
IRAFPointer tp          i: table descriptor
IRAFPointer *cp         o: column descriptor (created by this function), or
                           NULL in case of error
char *colname           i: column name
char *colunits          i: units for column, or null
char *colfmt           io: display format for column; if null on input,
                                it will be replaced with a default value
int datatype            i: IRAF data type for column
int nelem               i: length of array, or 1 for a scalar column
*/

        TableDescr *tbl_descr;
        ColumnDescr *col_descr;
        int ncols;              /* number of columns in the table */
        int colnum;             /* column number */
        int typecode;
        char *errmsg;

        *cp = init_cp (tp);
        col_descr = (ColumnDescr *)(*cp);

        /* Convert the IRAF data type to TFORM and to the CFITSIO typecode
           for tables.  If colfmt is null, copy a default to colfmt.
        */
        if (datatype == IRAF_REAL) {
            typecode = TFLOAT;
            sprintf (col_descr->tform, "%dE", nelem);
            if (colfmt[0] == '\0')
                strcpy (colfmt, "G15.7");
        } else if (datatype == IRAF_DOUBLE) {
            typecode = TDOUBLE;
            sprintf (col_descr->tform, "%dD", nelem);
            if (colfmt[0] == '\0')
                strcpy (colfmt, "G25.16");
        } else if (datatype == IRAF_SHORT) {
            typecode = TSHORT;
            sprintf (col_descr->tform, "%dI", nelem);
            if (colfmt[0] == '\0')
                strcpy (colfmt, "I11");
        } else if (datatype == IRAF_INT) {
            typecode = TINT;
            sprintf (col_descr->tform, "%dJ", nelem);
            if (colfmt[0] == '\0')
                strcpy (colfmt, "I11");
        } else if (datatype == IRAF_BOOL) {
            typecode = TLOGICAL;
            sprintf (col_descr->tform, "%dL", nelem);
            if (colfmt[0] == '\0')
                strcpy (colfmt, "L6");
        } else if (datatype < 0) {
            typecode = TSTRING;
            if (nelem > 1) {
                sprintf (col_descr->tform, "%dA%d",
                        -datatype * nelem, -datatype);
            } else {
                sprintf (col_descr->tform, "%dA", -datatype);
            }
            if (colfmt[0] == '\0')
                sprintf (colfmt, "A%d", datatype);      /* left justified */
        } else {
            errmsg = (char *)calloc (SZ_ERRMESS+1, sizeof(char));
            sprintf (errmsg, "Data type code %d not recognized.", datatype);
            setError (ERR_DATATYPE_UNKNOWN, errmsg);
            free (errmsg);
            *cp = NULL;
            return;
        }

        strcpy (col_descr->name, colname);
        strcpy (col_descr->tunit, colunits);
        strcpy (col_descr->tdisp, colfmt);

        tbl_descr = (TableDescr *)tp;
        ncols = tbl_descr->ncols;

        colnum = ncols + 1;                     /* one indexed */
        col_descr->colnum = colnum;
        col_descr->typecode = typecode;         /* CFITSIO data type */
        col_descr->datatype = datatype;         /* IRAF data type */
        if (typecode == TSTRING) {
            col_descr->repeat = -datatype * nelem;
            col_descr->nelem = nelem;
            col_descr->width = -datatype;       /* size of one element */
        } else {
            col_descr->repeat = nelem;
            col_descr->nelem = nelem;
            col_descr->width = 1;
        }

        /* make sure we have space for one more column */
        columnSpace (tbl_descr, 1);
        tbl_descr->columns[ncols] = *cp;
        tbl_descr->ncols = colnum;
}
