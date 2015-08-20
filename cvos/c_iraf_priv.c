/* 
** String defined to allow determination of the CVOS library version 
** from the library file (*.a) or the executable using the library.
*/
const char CVOS_VERSION[] = {"CVOS Version 3.5 (10-April-2001)"};

# include <stdlib.h>
# include <string.h>

/* Error Handling */
static int cvos_errcode = 0;
static char cvos_errmsg[256] = { '\0' };
typedef void (*c_IRAFErrHandler)(void);
c_IRAFErrHandler errhandler[32];
static int max_err_handlers = 32;
static int cvos_errtop = -1;

# include <stdio.h>
void clear_cvoserr(void) { cvos_errcode = 0; cvos_errmsg[0] = '\0'; }
char *c_iraferrmsg(void) { return cvos_errmsg; }
int c_iraferr(void) { return cvos_errcode; }

int c_pusherr(c_IRAFErrHandler x) {
	if (cvos_errtop == (max_err_handlers - 1)) return -1;
	++cvos_errtop;
	errhandler[cvos_errtop] = x;
	return cvos_errtop + 1;	
}

int c_poperr(void) {
	if (cvos_errtop == -1) return -1;
	--cvos_errtop;
	return cvos_errtop + 1;
}

void setError (int status, char *msg) {
/* Copy 'status' to cvos_errcode, and copy 'msg' to cvos_errmsg. */

        cvos_errcode = status;
        strcpy (cvos_errmsg, msg);
}

void clearError (void) {
/* Set the error code to 0 (OK), and erase the error message. */

        cvos_errcode = 0;
        cvos_errmsg[0] = '\0';
}

int checkError (void) {
/* Return the error code (0 is OK). */

        return cvos_errcode;
}
