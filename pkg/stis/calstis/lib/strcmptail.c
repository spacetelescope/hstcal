# include <string.h>

/* This function checks whether string s1 ends with the string s2.
   If it does, the function value will be zero.  strcmp is used to
   compare the end of s1 with s2, where "end of s1" means the same
   number of characters as the length of s2, or all of s1 if it is
   shorter than s2.  The value returned by strcmp is returned.

   Phil Hodge, 2003 Mar 21:
	Function created.
*/

int strcmptail (char *s1, char *s2) {

	int len_s1, len_s2;
	int start;

	len_s1 = strlen (s1);
	len_s2 = strlen (s2);

	start = (len_s1 < len_s2) ? 0 : (len_s1 - len_s2);

	return strcmp (s1+start, s2);
}
