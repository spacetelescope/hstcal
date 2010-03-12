# include <string.h>
# include <fitsio.h>
# include "ctables.h"

void c_tbhpcm (IRAFPointer tp, char *keyword, char *comment) {

/* Replace the comment for an existing keyword in a header.
arguments:
IRAFPointer tp          i: table descriptor
char *keyword           i: keyword name
char *comment           i: value for the comment
*/

        TableDescr *tbl_descr;
        int status = 0;

        tbl_descr = (TableDescr *)tp;

        /* fits_modify_comment = ffmcom */
        fits_modify_comment (tbl_descr->fptr, keyword, comment, &status);
        if (status != 0)
            setError (status, "c_tbhpcm:  error modifying comment");
}
