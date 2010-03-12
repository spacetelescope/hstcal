# include <stdlib.h>
# include <string.h>
# include <fitsio.h>
# include "ctables.h"

/* public functions:

    int c_iraferr();            return the error code (0 is OK)
    char *c_iraferrmsg();       return the error message
    void clear_cvoserr();       clear the error code and error message

interface functions:

    setError (int status, char *msg);
    clearError();
    error_code = checkError();
*/

static int error_code = 0;
static char error_message[SZ_ERRMESS+1] = "";

void setError (int status, char *msg) {
/* Copy 'status' to error_code, and copy 'msg' to error_message. */

        error_code = status;
        copyString (error_message, msg, SZ_ERRMESS);
}

void clearError (void) {
/* Set the error code to 0 (OK), and erase the error message. */

        error_code = 0;
        error_message[0] = '\0';
}

int checkError (void) {
/* Return the error code (0 is OK). */

        return error_code;
}

int c_iraferr (void) {
/* Return the error code (0 is OK). */

        return error_code;
}

char *c_iraferrmsg (void) {

/* Return a pointer to the error message.  If the error code is greater than
   100, the CFITSIO error stack will first be appended to the current contents
   of the error message.  Note that this clears the CFITSIO error stack.
*/

        if (error_code > 100) {         /* implies cfitsio error code */

            int done=0;
            int len_err;
            int err_num;
            char *err_msg;

            /* concatenate the error stack from CFITSIO */
            err_msg = (char *)calloc (SZ_ERRMESS+1, sizeof(char));
            while (!done) {
                /* fits_read_errmsg = ffgmsg */
                err_num = fits_read_errmsg (err_msg);
                if (err_num == 0)
                    done = 1;
                if (!done) {
                    len_err = strlen (error_message);
                    if (len_err + strlen (err_msg) >= SZ_ERRMESS-1)
                        break;
                    if (error_message[0] != '\0')
                        strcat (error_message, "; ");
                    strcat (error_message, err_msg);
                }
            }
            free (err_msg);
        }

        return error_message;
}

void clear_cvoserr (void) {
/* Set the error code to 0 (OK), and erase the error message. */

        clearError();
}
