# include <stdio.h>
# include <stdlib.h>
# include <time.h>
# include <string.h>
# include <ctype.h>

# include "stis.h"

# define SZ_TIMESTRING  100
# define SZ_EIGHTY       80

/* This routine prints the current date and time into a string
   in the format preferred by OPUS, date & time = YYYYDDDHHMMSS.
   Example:

1996212152529-I--------------- Generic Conversion Completed: N2RUSG54B ---------
*/

void TimeStamp (char *message, char *rootname) {

/* arguments:
char *message   i: string to append to date & time info
char *rootname  i: root name to include in printed string
*/

	char *timestring;	/* string to be printed */
	char *uc_rootname;	/* rootname converted to upper case */
	int lenrootname;	/* length of rootname string */
	int i;
	time_t *tp, now;

	if ((timestring = malloc (SZ_TIMESTRING * sizeof (char))) == NULL) {
	    printf ("Warning  (TimeStamp) Out of memory.\n");
	    return;
	}

	lenrootname = strlen (rootname);
	if (lenrootname > 0) {
	    uc_rootname = malloc ((lenrootname+1) * sizeof (char));
	    if (uc_rootname == NULL) {
		printf ("Warning  (TimeStamp) Out of memory.\n");
		return;
	    }
	    strcpy (uc_rootname, rootname);
	    for (i = 0;  i < strlen (rootname);  i++) {
		if (islower (uc_rootname[i]))
		    uc_rootname[i] = toupper (uc_rootname[i]);
	    }
	}

	tp = NULL;
	now = time (tp);
	strftime (timestring, SZ_TIMESTRING, "%Y%j%H%M%S", localtime (&now));
	strcat (timestring, "-I--------------- ");
	strcat (timestring, message);
	if (lenrootname > 0) {
	    strcat (timestring, ": ");
	    strcat (timestring, uc_rootname);
	}
	strcat (timestring, " ");

	/* Fill out to 80 bytes with dashes. */
	for (i = strlen (timestring);  i < SZ_EIGHTY;  i++)
	    timestring[i] = '-';
	timestring[SZ_EIGHTY] = '\0';

	printf ("%s\n", timestring);

	if (lenrootname > 0)
	    free (uc_rootname);
	free (timestring);
}
