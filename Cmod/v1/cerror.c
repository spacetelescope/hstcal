# include <stdlib.h>
# include <string.h>
# include "ctables.h"

/* public functions:

    int hstio_err();
    int c_iraferr();
    char *c_iraferrmsg();

interface functions:

    setError (int status, char *msg);
    clearError();
    error_code = checkError();
*/

static int error_code = 0;
static char error_message[SZ_ERRMESS] = "";

void setError (int status, char *msg) {
	error_code = status;
	strncpy (error_message, msg, SZ_ERRMESS-1);
}

void clearError (void) {
	error_code = 0;
	error_message[0] = '\0';
}

int checkError (void) {
	return error_code;
}

int hstio_err (void) {
	return 0;
}

int c_iraferr (void) {
	return error_code;
}

char *c_iraferrmsg (void) {
	return error_message;
}
