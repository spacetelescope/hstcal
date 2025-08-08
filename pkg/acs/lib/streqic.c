# include <ctype.h>
#include <stddef.h>
#include <string.h>

/* This function compares two strings without regard to case, returning
   one if the strings are equal.

   Phil Hodge, 1997 Dec 12:
	Function created.
*/

int streq_ic (char *s1, char *s2) {
	if (!s1 || !s2)
		return 0;

	const size_t s1_len = strlen(s1);
	const size_t s2_len = strlen(s2);
	if (s1_len != s2_len) {
		return 0;
	}

	for (size_t i = 0; i < s1_len; i++) {
		const char c1 = (char) tolower((unsigned char) s1[i]);
		const char c2 = (char) tolower((unsigned char) s2[i]);
		if (c1 != c2)
			return 0;
	}
	return 1;
}
