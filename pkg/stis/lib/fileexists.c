# include <stdio.h>
# include "hstio.h"	/* for ckNewFile */

# include "stis.h"
# include "hstcalerr.h"

/* This routine takes action if the output file already exists.  If not,
   zero is returned.  If it does, the environment variable imclobber is
   checked.  If imclobber is defined and its value is "yes", the file will
   be deleted.  If imclobber is not defined or has some other value, an
   error message will be printed, and an error code will be returned.

   Phil Hodge, 1998 Oct 5:
	Change status values 1021 and 1023 to GENERIC_ERROR_CODE.
*/

int FileExists (char *fname) {

	int flag;

	flag = ckNewFile (fname);
	if (flag > 0) {
	    if (flag == 1) {
		trlerror("Output file `%s' already exists.", fname);
		return (GENERIC_ERROR_CODE);
	    } else {
		trlerror("Can't clobber `%s'.", fname);
		return (GENERIC_ERROR_CODE);
	    }
	}
	return (0);
}
