# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <ctype.h>
# include "wf3.h"	/* for trlerror & NOT_APPLICABLE */
# include "err.h"

int GotFileName (char *filename) {

/* arguments:
char *filename  i: name of reference file
*/

	char *fname;		/* name of file converted to lower case */
	int len;
	int i, c;
	int flag;

	len = strlen (filename);

	if ((fname = calloc (len+1, sizeof(char))) == NULL) {
	    trlerror ("Out of memory.\n");
	    return (OUT_OF_MEMORY);
	}

	/* Copy filename to fname, trimming trailing whitespace and
	** converting to lower case. */
	fname[len] = '\0';
	i = len - 1;
	while (i >= 0) {
	    if (isspace (filename[i])) {
		fname[i] = '\0';
		i--;
	    } else {
		break;
	    }
	}
	while (i >= 0) {
	    c = filename[i];
	    if (isupper (c))
		fname[i] = tolower (c);
	    else
		fname[i] = c;
	    i--;
	}

	if (fname[0] == '\0')
	    flag = 0;
	else if (strcmp (fname, NOT_APPLICABLE) == 0)
	    flag = 0;
	else
	    flag = 1;

	free (fname);

	return (flag);
}
