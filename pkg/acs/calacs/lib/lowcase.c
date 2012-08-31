/* LOWCASE: Contains routines for converting the case of strings */

# include <stdio.h>
# include <string.h>
# include <ctype.h>


/* Convert s2 into lower case and copy out to s1.
**	Based on 'strcpy.c' function.
*/
char *(lowcase) (char *s1, const char *s2) {

	char *s = s1;
	
	for (s = s1; (*s++ = tolower(*s2++)) != '\0'; ) 
		;
	return (s1);
}

