# include <string.h>
# include <ctype.h>

# include "stis.h"		/* for PERFORM and OMIT */
# include "stisomit.h"	/* a list of switches that should be OMIT */
# define SIZE_KEYWORD  9

/* This routine assigns the default value for a calibration switch,
   for the case that no switch was specified on the command line.
*/

int DefSwitch (char *keyword) {

	int i;
	char lckey[SIZE_KEYWORD];	/* keyword converted to lower case */

	for (i = 0;  i < SIZE_KEYWORD;  i++) {
	    if (isupper (keyword[i]))
		lckey[i] = tolower (keyword[i]);
	    else
		lckey[i] = keyword[i];
	    if (keyword[i] == '\0')
		break;
	}

	for (i = 0;  i < STIS_N_OMIT;  i++) {
	    if (strcmp (lckey, omitsw[i]) == 0)
		return (OMIT);
	}

	return (PERFORM);
}
