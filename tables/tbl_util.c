# include <ctype.h>
# include <stdlib.h>
# include <stdio.h>
# include <string.h>
# include <fitsio.h>
#include "hstio.h"
#include "hstcal.h"
# include "ctables.h"

char *expandfn (char *filename) {
/* Expand an environment variable in a file name.  The variable may be
   either IRAF-style (e.g. abc$def/filename.fits) or
   Unix-style (e.g. $abc/def/filename.fits).

   This function returns a pointer to memory allocated by this function,
   and that memory should be deallocated by the calling routine.
argument:
char *filename          i: name of a file, possibly including an environment
                           variable as part of the directory

function value          o: name of the file with the environment variable
                           expanded to its actual value
*/

        int i, j;
        int dollar, start;
        int done;
        char *var;
        char *value=NULL, *l_value;

        i = 0;
        start = 0;
        dollar = -1;
        done = 0;
        /* look for a '$', indicating an environment variable */
        while (!done) {
            if (filename[i] == '\0') {
                done = 1;
            } else if (filename[i] == '$') {
                done = 1;
                var = (char *)calloc (CHAR_FNAME_LENGTH+1, sizeof(char));
                dollar = i;
                if (dollar == 0) {
                    /* Unix-style file name */
                    for (j = 1;  filename[j] != '\0';  j++) {
                        if (j >= CHAR_FNAME_LENGTH) {
                            var[j] = '\0';
                            start = j + 1;
                            break;
                        }
                        var[j-1] = filename[j];
                        if (filename[j] == '/') {
                            var[j-1] = '\0';
                            start = j + 1;
                            break;
                        }
                    }
                } else {
                    /* IRAF-style file name */
                    start = dollar + 1;
                    if (dollar >= CHAR_FNAME_LENGTH)
                        dollar = CHAR_FNAME_LENGTH;
                    for (j = 0;  j < dollar;  j++)
                        var[j] = filename[j];
                    var[dollar] = '\0';
                }
                value = getenv (var);
                free (var);
            }
            i++;
        }
        l_value = (char *)calloc (CHAR_FNAME_LENGTH+1, sizeof(char));
        if (dollar >= 0) {
            if (value == NULL) {
                strcpy (l_value, filename);
            } else {
                strcpy (l_value, value);
                i = strlen (l_value) - 1;
                if (l_value[i] != '/')
                    strcat (l_value, "/");
                strcat (l_value, filename+start);
            }
        } else {
            strcpy (l_value, filename);
        }
        return l_value;
}

int checkExists (char *filename) {

/* Test whether a file exists.
argument:
char *filename          i: name of a file

function value          o: 1 if the file exists, 0 otherwise
*/

        FILE *fd;
        int file_exists;

        /* check whether the file exists */
        fd = fopen (filename, "r");
        if (fd == NULL) {
            file_exists = 0;
        } else {
            file_exists = 1;
            (void)fcloseWithStatus(&fd);
        }

        return file_exists;
}

void cToFortran (char *c_fmt, char *ftn_fmt) {

/* Convert C-style display (print) format to Fortran-style.
arguments:
char *c_fmt             i: C-style format string (e.g. "%5d")
char *ftn_fmt           o: Fortran-style format string (e.g. "I5")

   If c_fmt does not begin with '%', it will just be copied to ftn_fmt
   without any change.

   Here are some examples of corresponding format strings:

        c_fmt  ftn_fmt      description
        %12.5f  F12.5     floating-point value
        %12.5e  E12.5     floating-point value
        %12.5e  D12.5     floating-point value
        %12.5g  G12.5     general floating-point value
        %12d    I12       integer
        %012d   I12.12    integer padded with '0' on the left
        %12b    L12       boolean
        %17s    A17       character string
*/

        char numpart[SZ_FITS_STR+1];
        int i, j, k;

        if (c_fmt[0] != '%') {
            strcpy (ftn_fmt, c_fmt);
            return;
        }

        k = strlen (c_fmt) - 1;

        /* copy out the numerical portion */
        j = 1;
        if (c_fmt[1] == '0')            /* e.g. "%012d" */
            j = 2;
        for (i = 0;  j < k;  i++, j++)
            numpart[i] = c_fmt[j];
        numpart[i] = '\0';

        if (c_fmt[k] == 'f' || c_fmt[k] == 'e' || c_fmt[k] == 'g') {
            ftn_fmt[0] = toupper (c_fmt[k]);
            ftn_fmt[1] = '\0';
            strcat (ftn_fmt, numpart);
        } else {
            if (c_fmt[k] == 'b') {
                ftn_fmt[0] = 'L';
                ftn_fmt[1] = '\0';
                strcat (ftn_fmt, numpart);
            } else if (c_fmt[k] == 's') {
                ftn_fmt[0] = 'A';
                ftn_fmt[1] = '\0';
                strcat (ftn_fmt, numpart);
            } else if (c_fmt[k] == 'd') {
                ftn_fmt[0] = 'I';
                ftn_fmt[1] = '\0';
                strcat (ftn_fmt, numpart);
                if (c_fmt[1] == '0') {
                    strcat (ftn_fmt, ".");
                    strcat (ftn_fmt, numpart);
                }
            } else {
                /* format code not recognized */
                strcpy (ftn_fmt, c_fmt);
            }
        }
}

void trimString (char *value) {

/* Trim trailing blanks and remove enclosing quotes, in-place.
argument:
char *value             io: string to be modified
*/

        int j;

        j = strlen (value) - 1;
        if (j >= 0 && value[j] == '\'') {
            value[j] = '\0';
            j--;
        }

        while (j >= 0) {
            if (value[j] == ' ') {
                value[j] = '\0';
                j--;
            } else {
                break;
            }
        }

        if (value[0] == '\'') {
            /* value begins with a quote; remove it */
            int done = 0;
            j = 1;
            while (!done) {
                value[j-1] = value[j];
                if (value[j] == '\0')
                    done = 1;
                j++;
            }
        }
}

void copyString (char *output, char *input, int maxch) {

/* copy up to maxch characters from input to output, and then append a null
arguments:
char *output            o: copy of 'input'
char *input             i: string to be copied to 'output'
int maxch               i: maximum length of output string (not including '\0')
*/

        int i;

        for (i = 0;  i < maxch;  i++) {
            output[i] = input[i];
            if (input[i] == '\0')
                break;
        }
        output[i] = '\0';
}

void str_lower (char lc_name[], const char name[]) {
/* Copy name to lc_name, converting to lower case. */

        int i;

        for (i = 0;  name[i] != '\0';  i++) {
            if (isupper (name[i]))
                lc_name[i] = tolower (name[i]);
            else
                lc_name[i] = name[i];
        }
        lc_name[i] = '\0';
}
