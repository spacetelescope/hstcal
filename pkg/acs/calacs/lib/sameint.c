/* This file contains SameInt and SameString.

   When comparing a value (e.g. from the input header) with a value read
   from a table row in order to select the appropriate row, the values
   may either match exactly, or the table value may be a wildcard which
   matches anything.  The following routines return one if the values
   are the same (case sensitive for strings) or if the rowvalue is the
   wildcard; otherwise, zero is returned. 

   Note that it is rowvalue, not value, that is compared with the
   wildcard.
   Additionally, an integer value equivalent to INT_IGNORE will be treated
   as a match, just like a WILDCARD.  As a result, columns with values of 
   INT_IGNORE will be matched as if they were legal values, and let the 
   rest of the selection criteria determine a valid match.  11 Feb 1999 WJH
*/

# include <string.h>
# include "acswild.h"	/* defines wildcard values */

int SameInt (int rowvalue, int value) {

	if (rowvalue == INT_WILDCARD)
	    return (1);
	else if (rowvalue == value)
	    return (1);
	else if (rowvalue == INT_IGNORE)
	    return (1);
	else
	    return (0);
}

int SameFlt (float rowvalue, float value) {

	if (rowvalue == FLT_WILDCARD)
	    return (1);
	else if (rowvalue == value)
	    return (1);
	else if (rowvalue == FLT_IGNORE)
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
