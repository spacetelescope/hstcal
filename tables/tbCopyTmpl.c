# include <string.h>
# include <fitsio.h>
# include "ctables.h"

void tbCopyTmpl (IRAFPointer tp) {

/* Copy column info from a template table that's already open.
argument:
IRAFPointer tp          i: descriptor for the new table (not the template)
*/

        TableDescr *tbl_descr;          /* table being opened as a NEW_COPY */
        TableDescr *template_descr;     /* table descriptor for template */
        /* for a column in the template and a column in the NEW_COPY table */
        ColumnDescr *tcol_descr, *col_descr;
        IRAFPointer tcp, cp;
        int ncols;              /* number of columns */
        int i;

        tbl_descr = (TableDescr *)tp;
        template_descr = (TableDescr *)tbl_descr->template;

        ncols = template_descr->ncols;

        /* allocate enough space for the column pointers */
        /* (don't update tbl_descr->ncols yet because columnSpace uses it) */
        columnSpace (tbl_descr, ncols);
        if (checkError() != 0)
            return;

        tbl_descr->ncols = ncols;

        /* Copy info for each column in the table. */
        for (i = 0;  i < ncols;  i++) {

            tcp = template_descr->columns[i];
            tcol_descr = (ColumnDescr *)tcp;
            cp = init_cp (tp);
            col_descr = (ColumnDescr *)cp;

            strcpy (col_descr->name, tcol_descr->name);
            strcpy (col_descr->tform, tcol_descr->tform);
            strcpy (col_descr->tunit, tcol_descr->tunit);
            strcpy (col_descr->tdisp, tcol_descr->tdisp);

            col_descr->colnum = i + 1;          /* one indexed */
            col_descr->typecode = tcol_descr->typecode;
            col_descr->datatype = tcol_descr->datatype;
            col_descr->repeat = tcol_descr->repeat;
            col_descr->nelem = tcol_descr->nelem;
            col_descr->width = tcol_descr->width;

            tbl_descr->columns[i] = cp;
        }
}
