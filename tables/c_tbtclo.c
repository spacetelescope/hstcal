# include <fitsio.h>
# include "ctables.h"

void c_tbtclo (IRAFPointer tp) {

/* Close a table.
argument:
IRAFPointer tp          i: table descriptor
*/

        TableDescr *tbl_descr;
        fitsfile *fptr;         /* CFITSIO pointer */
        int status = 0;

        if (tp == NULL)
            return;

        tbl_descr = (TableDescr *)tp;
        fptr = tbl_descr->fptr;

        /* fits_close_file = ffclos */
        fits_close_file (fptr, &status);

        free_tp (tp);
}
