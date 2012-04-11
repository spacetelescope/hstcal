/* A_ERR: Contains various error handling routines */

# include <stdio.h>
# include <string.h>

# include "hstio.h"	/* defines HST I/O functions */
# include "acs.h"	/* defines ACS data structures */

void asnerror (char *);

extern char MsgText[ACS_LINE+1];

void errchk() {
		
	extern int status;
	
	if (hstio_err()) {
	fprintf (stdout, "\n*** ERROR in HST I/O functions:\n%s\n",
			 hstio_errmsg());
	fflush (stdout);
	status = 1;
	}
}

void asnkwerr (char *keyword, char *file) {
	sprintf (MsgText, "Keyword \"%s\" not found in %s", keyword, file);
	asnerror (MsgText);
}

void asnopenerr (char *name) {
	sprintf (MsgText, "Can't open file %s", name);
	asnerror (MsgText);
}

void asnreaderr (char *name) {
	sprintf (MsgText, "Can't read file %s", name);
	asnerror (MsgText);
}

void asnfilerr (char *name) {
	sprintf (MsgText, "while trying to read file %s", name);
	asnerror (MsgText);
}

void asnmessage (char *message) {
	printf ("%s\n", message);
        fflush(stdout);
}

void asnwarn (char *message) {

	char line[ACS_LINE+1];
	
	/* Use macro for beginning of Warning message */
	sprintf(line,"%s",WARN_PREFIX);
	strcat (line,message);

    asnmessage(line);
}

void asnerror (char *message) {
		
	char line[ACS_LINE+1];
	
	/* Use macro for beginning of Warning message */
	sprintf(line,"%s",ERR_PREFIX);
	strcat (line,message);

    asnmessage(line);
}
