/* This file contains SameInt and SameString.

   When comparing a value (e.g. from the input header) with a value read
   from a table row in order to select the appropriate row, the values
   may either match exactly, or the table value may be a wildcard which
   matches anything.  The following routines return one if the values
   are the same (case sensitive for strings) or if the rowvalue is the
   wildcard; otherwise, zero is returned.

   Note that it is rowvalue, not value, that is compared with the
   wildcard.
*/

# include <string.h>

# include "stis.h"
# include "stiswild.h"	/* defines wildcard values */

int SameInt (int rowvalue, int value) {

	if (rowvalue == INT_WILDCARD)
	    return (1);
	else if (rowvalue == value)
	    return (1);
	else
	    return (0);
}

int SameString (char *rowvalue, char *value) {

	if (strcmp (rowvalue, STRING_WILDCARD) == 0)
	    return (1);
	else if (strcmp (rowvalue, value) == 0)
	    return (1);
	else
	    return (0);
}
