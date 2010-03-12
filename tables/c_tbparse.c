# include <ctype.h>
# include <stdlib.h>
# include <string.h>
# include <fitsio.h>
# include "ctables.h"

int c_tbparse (char *tablename, char *fname, char *extname, int maxch,
                int *hdu) {

/* Extract the file name 'fname' from the table name 'tablename'.  These
   will be the same unless 'tablename' contains an expression in brackets
   following the file name.  The expression in brackets may be just an
   integer value, in which case that value will be assigned to hdu, or the
   expression may be a string (without quotes), in which case the string
   will be copied to extname and hdu will be set to -1.  The expression may
   actually include both extname and extver (integer), separated by a comma,
   but this function will ignore extver.
arguments:
char *tablename         i: name of table; may include extname or HDU number
                           in brackets (HDU = 0 is the primary HDU)
char *fname             o: file name; this may still include an environment
                           variable as part of the directory name
char *extname           o: extension name, if specified
int maxch               i: maximum length of extname (not including '\0')
int *hdu                o: extension number, if specified (0 is primary HDU)

function value          o: length of 'tablename' string, excluding trailing
                           blanks
*/

        /* This implementation is very simple-minded, and it should be
           improved if this function is to be used for anything other
           than getting fname.
        */

        int i, j;
        int nchar;      /* number of non-blank characters in tablename */
        int br;         /* index of '[', or -1 if no brackets */

        /* initial values */
        fname[0] = '\0';
        extname[0] = '\0';
        *hdu = -1;
        nchar = 0;
        br = -1;

        /* Get nchar, and check for blank tablename. */
        for (i = strlen (tablename) - 1;  i >= 0;  i--) {
            if (tablename[i] != ' ') {
                nchar = i + 1;  /* add one to correct for zero indexing */
                break;
            }
        }
        if (nchar <= 0)
            return nchar;

	/* Find the first open bracket, and copy out the file name. */
        for (i = 0;  i < nchar;  i++) {
            if (tablename[i] == '[') {
                br = i;
                break;
            }
            fname[i] = tablename[i];
        }
        /* If an open bracket was found, this truncates fname at that
           point; if no bracket was found, fname will be the same as
           tablename, and this ends the string.
        */
        fname[i] = '\0';

        /* Copy extname or the hdu number, if specified. */
        if (br >= 0) {
            if (isdigit (tablename[br+1])) {
                *hdu = atoi (tablename+br+1);
            } else {
                for (i = br+1, j = 0;  i < nchar;  i++, j++) {
                    if (tablename[i] == ']' ||
                        tablename[i] == ',' || j >= maxch) {
                        extname[j] = '\0';
                        break;
                    } else {
                        extname[j] = tablename[i];
                    }
                }
                extname[j] = '\0';
            }
        }

        return nchar;
}
