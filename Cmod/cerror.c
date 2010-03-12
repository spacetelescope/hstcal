# include <stdlib.h>
# include <string.h>
# include <fitsio.h>
# include "ctables.h"

/* public functions:

    xxx deleted

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
