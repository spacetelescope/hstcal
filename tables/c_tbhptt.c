# include <string.h>
# include <fitsio.h>
# include "ctables.h"

void c_tbhptt (IRAFPointer tp, char *keyword, char *text) {

/* Update a keyword, which is supposed to already exist, but for FITS
   tables it doesn't have to.
arguments:
IRAFPointer tp          i: table descriptor
char *keyword           i: keyword name
char *text              i: value for the keyword
*/

        TableDescr *tbl_descr;
        int status = 0;

        tbl_descr = (TableDescr *)tp;

        /* fits_update_key = ffuky */
        fits_update_key (tbl_descr->fptr, TSTRING, keyword, text,
                                NULL, &status);
        if (status != 0)
            setError (status, "c_tbhptt:  error updating keyword");
}
