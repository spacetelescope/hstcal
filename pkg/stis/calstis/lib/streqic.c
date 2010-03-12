# include <ctype.h>

# include "stis.h"

/* This function compares two strings without regard to case, returning
   one if the strings are equal.

   Phil Hodge, 1997 Dec 12:
	Function created.
*/

int streq_ic (char *s1, char *s2) {

	int c1, c2;
	int i;

	c1 = 1;
	for (i = 0;  c1 != 0;  i++) {

	    c1 = s1[i];
	    c2 = s2[i];
	    if (isupper(c1))
		c1 = tolower (c1);
	    if (isupper(c2))
		c2 = tolower (c2);
	    if (c1 != c2)
		return (0);
	}
	return (1);
}
