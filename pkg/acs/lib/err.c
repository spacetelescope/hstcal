/* A_ERR: Contains various error handling routines */

# include <stdio.h>
# include <string.h>

#include "hstcal.h"
# include "hstio.h"	/* defines HST I/O functions */
# include "acs.h"	/* defines ACS data structures */

void errchk() {
		
	extern int status;
	
	if (hstio_err()) {
	fprintf (stdout, "\n*** ERROR in HST I/O functions:\n%s\n",
			 hstio_errmsg());
	fflush (stdout);
	status = 1;
	}
}