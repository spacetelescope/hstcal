# include <string.h>
# include <fitsio.h>
# include "ctables.h"

void c_tbhgnp (IRAFPointer tp, int parnum,
                char *keyword, int *dtype, char *str) {

/* Get the Nth keyword from a header.
   Note:  This follows the STSDAS tables convention that the string following
   a commentary keyword (HISTORY, COMMENT or blank) is the value, so that
   string will be copied to the output 'str'.  The output buffer 'str' must
   be at least 73 characters long (i.e. 72 plus one for '\0').
   The END record has no value or data type, but something had to be set,
   so 'str' is empty ('\0') and the data type is set to IRAF_CHAR.
arguments:
IRAFPointer tp          i: table descriptor
int parnum              i: number of keyword to be gotten (one indexed)
char *keyword           o: keyword name
int *dtype              o: data type of keyword, one of:
                           IRAF_CHAR, IRAF_BOOL, IRAF_INT, IRAF_DOUBLE
char *str               o: string containing the value of the keyword
*/

        fitsfile *fptr;
        TableDescr *tbl_descr;
        char value[SZ_FITS_STR+1], comment[SZ_FITS_STR+1];
        char *char_implying_float;
        int status = 0;

        tbl_descr = (TableDescr *)tp;
        fptr = tbl_descr->fptr;

        /* fits_read_keyn = ffgkyn */
        fits_read_keyn (fptr, parnum, keyword, value, comment, &status);
        if (status != 0) {
            setError (status, "c_tbhgnp:  couldn't read Nth keyword");
            keyword[0] = '\0';
            *dtype = -1;
            str[0] = '\0';
            return;
        }

        if (strcmp (keyword, "HISTORY") == 0 ||
            strcmp (keyword, "COMMENT") == 0 ||
            keyword[0] == ' ') {
            strcpy (str, comment);
            *dtype = IRAF_CHAR;
        } else if (strcmp (keyword, "END") == 0) {
            str[0] = '\0';              /* no value */
            *dtype = IRAF_CHAR;
        } else {
            if (value[0] == '\'') {
                *dtype = IRAF_CHAR;
            } else if (value[0] == 'T' || value[0] == 'F') {
                *dtype = IRAF_BOOL;
            } else {
                char_implying_float = strpbrk (value, ".EeDd");
                if (char_implying_float == NULL)
                    *dtype = IRAF_INT;
                else
                    *dtype = IRAF_DOUBLE;
            }
            trimString (value);
            strcpy (str, value);
        }
}
