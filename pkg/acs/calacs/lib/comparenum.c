# include <string.h>
#include "hstcal.h"
# include "acs.h"
# include "hstcalerr.h"


/* This function checks that the number of input and output files are
   the same, unless no output file was specified.  If the numbers match,
   zero will be returned; if they aren't the same, a value of one will
   be returned. That is, this function returns true for the error condition.
*/

int CompareNumbers (int n_in, int n_out, char *str_out) {
	char buff[ACS_FITS_REC+1];

	buff[0] = '\0';

	if (n_out > 0 && n_out != n_in) {
	    sprintf (MsgText, "You specified %d input file", n_in);
	    if (n_in > 1) {
			strcat (MsgText, "s");
		}
	    sprintf (buff, " but %d %s file", n_out, str_out);
		strcat (MsgText, buff);
	    if (n_out > 1) {
	    	strcat (MsgText, "s");
		}
	    strcat (MsgText, ".\n");

            printf ("%s", MsgText);
		return (1);
	}

	return (0);
}
