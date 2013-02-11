# include <string.h>
# include <fitsio.h>
# include "ctables.h"

int c_tbfres (char *keyword) {

/* Check whether the specified keyword is a reserved keyword for a table.
arguments:
char *keyword           i: header keyword name

function value          o: true (1) if the keyword is reserved for a table
                           header, else 0
*/

        char uc_keyword[SZ_FITS_STR+1]; /* keyword in upper case */
        int value;

        strcpy (uc_keyword, keyword);
        /* fits_uppercase = ffupch */
        fits_uppercase (uc_keyword);

        if (strcmp (uc_keyword, "XTENSION") == 0)
            value = 1;
        else if (strcmp (uc_keyword, "BITPIX") == 0)
            value = 1;
        else if (strncmp (uc_keyword, "NAXIS", 5) == 0)
            value = 1;
        else if (strcmp (uc_keyword, "PCOUNT") == 0)
            value = 1;
        else if (strcmp (uc_keyword, "GCOUNT") == 0)
            value = 1;
        else if (strcmp (uc_keyword, "TFIELDS") == 0)
            value = 1;
        else if (strcmp (uc_keyword, "END") == 0)
            value = 1;
        else if (strncmp (uc_keyword, "TBCOL", 5) == 0)
            value = 1;
        else if (strncmp (uc_keyword, "TFORM", 5) == 0)
            value = 1;
        else if (strncmp (uc_keyword, "TTYPE", 5) == 0)
            value = 1;
        else if (strncmp (uc_keyword, "TUNIT", 5) == 0)
            value = 1;
        else if (strncmp (uc_keyword, "TSCAL", 5) == 0)
            value = 1;
        else if (strncmp (uc_keyword, "TZERO", 5) == 0)
            value = 1;
        else if (strncmp (uc_keyword, "TNULL", 5) == 0)
            value = 1;
        else if (strncmp (uc_keyword, "TDISP", 5) == 0)
            value = 1;
        else if (strncmp (uc_keyword, "TDIM", 4) == 0)
            value = 1;
        else if (strcmp (uc_keyword, "THEAP") == 0)
            value = 1;
        else
            value = 0;

        return value;
}
