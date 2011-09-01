# include <stdio.h>
# include <ctype.h>
# include <string.h>
# include "c_iraf.h"
# include "hstio.h"

# include "stis.h"
# include "stisdef.h"

/* This routine updates the FILENAME primary header keyword, or adds it
   to the header if it's not already present.  If the input file name
   begins with an environment variable or an explicit directory, that will
   be skipped over before writing the name to the header.
   The section for finding a directory prefix is based on the iraf zfnbrk
   function.
*/

void UFilename (char *filename, Hdr *phdr) {

/* arguments:
char *filename  i: file name, but may begin with directory specification
Hdr *phdr       io: pointer to header; FILENAME will be updated
*/

	int ch;			/* one character in filename */
	int namelen;		/* length of filename */
	int start = 0;		/* start of file name */
	int i;

	namelen = strlen (filename);

	/* If there's a directory prefix, skip over it. */
	for (i = 0;  i < namelen;  i++) {
	    ch = filename[i];
	    if (isalnum(ch) || ch == '_' || ch == '.')
		continue;		/* ordinary character */
	    else if (ch == '\\' && i+1 < namelen)
		i++;			/* skip \ and next character */
	    else
		start = i + 1;	/* not ordinary, so part of directory prefix */
	}

	/* Were we unable to find the filename root? */
	if (start >= namelen - 1)
	    start = 0;

	/* Update the FILENAME keyword, or add it if it doesn't exist. */
	Put_KeyS (phdr, "FILENAME", &filename[start], "name of file");
}
