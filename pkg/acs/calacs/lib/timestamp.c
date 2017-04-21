# include <stdio.h>
# include <stdlib.h>
# include <time.h>
# include <string.h>
# include <ctype.h>
# include "acs.h"	/* for message output */
# include "hstcalerr.h"

# define SZ_TIMESTRING  100 /* Make sure this is less that ACS_LINE */
# define SZ_EIGHTY       80

/* This routine prints the current date and time into a string
   in the format preferred by OPUS, date & time = YYYYDDDHHMMSS.
   Example:

1996212152529-I--------------- Generic Conversion Completed: N2RUSG54B ---------

	Warren J. Hack, 1998 Nov 17:
		Revised to support trailer file comments.
		
	Warren J. Hack, 1999 Jan 5:
		Revised to work directly with MsgText string, instead of using
			temporary local string...
*/

void TimeStamp (char *message, char *rootname) {

/* arguments:
char *message   i: string to append to date & time info
char *rootname  i: root name to include in printed string
*/

	char *uc_rootname;	/* rootname converted to upper case */
	int lenrootname;	/* length of rootname string */
	int i;
	time_t *tp, now;

	/* Initialize output string */
	MsgText[0] = '\0';
	

	lenrootname = strlen (rootname);
	if (lenrootname > 0) {
	    uc_rootname = malloc ((lenrootname+1) * sizeof (char));
	    if (uc_rootname == NULL) {
			trlerror ("(TimeStamp) Out of memory.");
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
	strftime (MsgText, SZ_TIMESTRING, "%Y%j%H%M%S", localtime (&now));
	strcat (MsgText, "-I--------------- ");
	strcat (MsgText, message);
	if (lenrootname > 0) {
	    strcat (MsgText, ": ");
	    strcat (MsgText, uc_rootname);
	}
	strcat (MsgText, " ");

	/* Fill out to 80 bytes with dashes. */
	for (i = strlen (MsgText);  i < SZ_EIGHTY;  i++)
	    MsgText[i] = '-';
	MsgText[SZ_EIGHTY] = '\0';

	trlmessage (MsgText);
	
	if (lenrootname > 0)
	    free (uc_rootname);

}
