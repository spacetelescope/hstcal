# include <string.h>
# include <fitsio.h>
# include "ctables.h"

void c_tbtnam (IRAFPointer tp, char *tablename, int maxch) {

/* Get the name of an open table.
arguments:
IRAFPointer tp          i: table descriptor
char *tablename         o: name of table; may include extname or HDU number
                           in brackets (HDU = 0 is the primary HDU)
int maxch               i: maximum length of 'tablename' string (not
                           including '\0')
*/

        TableDescr *tbl_descr;
        char hdu_string[SZ_FITS_STR+1];
        int len;

        tbl_descr = (TableDescr *)tp;

        len = strlen (tbl_descr->filename);
        tablename[0] = '\0';
        strncat (tablename, tbl_descr->filename, maxch);

        /* append the extension number (or as specified by the user) */
        if (tbl_descr->brackets[0] == '\0')
            sprintf (hdu_string, "[%d]", tbl_descr->hdunum - 1);
        else
            strcpy (hdu_string, tbl_descr->brackets);

        if (len + strlen (hdu_string) <= maxch) {
            strcat (tablename, hdu_string);
        } else {
            setError (ERR_STRING_TOO_LONG,
                "c_tbtnam:  buffer for table name is too short");
        }
}
