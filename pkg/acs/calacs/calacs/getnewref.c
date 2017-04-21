# include <stdlib.h>
# include "hstio.h"
# include "acs.h"
# include "err.h"

/* This routine gets the value of a reference file keyword from the
   primary header and appends it to the list of reference file info.

   Phil Hodge, 1997 Dec 10:
	Function created.

   Warren Hack, 26 May 1998:
   	Function edited to work with ACS.
*/

int GetNewRef (Hdr *phdr, char *keyword, RefFileInfo *ref) {

/* arguments:
Hdr *phdr         i: primary header of science file
char *keyword     i: keyword name for reference file
RefFileInfo ref  io: list of keyword,filename pairs
*/

	extern int status;
	char *filename;
	int GetKeyStr (Hdr *, char *, int, char *, char *, int);
	int NewRefFile (RefFileInfo *, char *, char *);

	if ((filename = calloc (ACS_FITS_REC+1, sizeof(char))) == NULL)
	    return (status = OUT_OF_MEMORY);

	if (GetKeyStr (phdr, keyword, USE_DEFAULT, "", filename, ACS_FITS_REC))
	    return (status);

	if (NewRefFile (ref, keyword, filename))
	    return (status);

	free (filename);

	return (status);
}
