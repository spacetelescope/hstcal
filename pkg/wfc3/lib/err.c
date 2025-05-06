/* ERR: Contains various error handling routines */

# include <stdio.h>
# include <string.h>

# include "hstio.h"	/* defines HST I/O functions */
# include "msg.h"	/* for MsgText and asnerror */
# include "hstcal.h"
# include "trlbuf.h"	/* for trlerror */

void errchk() {
		
	extern int status;
	
	if (hstio_err()) {
	    status = 1;
	    sprintf (MsgText, " in HST I/O functions:\n%s\n", hstio_errmsg());
	    trlerror (MsgText);
	}
}