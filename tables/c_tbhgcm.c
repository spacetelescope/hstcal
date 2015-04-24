# include <string.h>
# include <fitsio.h>
# include "ctables.h"

void c_tbhgcm (IRAFPointer tp, char *keyword, char *comment, int maxch) {

/* Get the comment for a keyword from a header.
   Note:  The STSDAS tables convention (followed here) is that for HISTORY,
   COMMENT or blank keyword, the following string is the value and there is
   no comment.  The CFITSIO convention is that such keywords have no value
   and the string is the comment.
arguments:
IRAFPointer tp          i: table descriptor
char *keyword           i: keyword name
char *comment           o: value for the comment
int maxch               i: maximum length of the comment (not including '\0')
*/

        TableDescr *tbl_descr;
        char value[SZ_FITS_STR+1], cmt[SZ_FITS_STR+1];
        int status = 0;

        tbl_descr = (TableDescr *)tp;

        /* fits_read_keyword = ffgkey */
        fits_read_keyword (tbl_descr->fptr, keyword, value, cmt, &status);
        if (status != 0)
            setError (status, "c_tbhgcm:  error reading comment");

        if (strcmp (keyword, "HISTORY") == 0 ||
            strcmp (keyword, "COMMENT") == 0 ||
            keyword[0] == ' ') {
            comment[0] = '\0';
        } else {
            copyString (comment, cmt, maxch);
        }
}
