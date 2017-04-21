# include <stdlib.h>		/* calloc */
# include "hstio.h"
# include "wf3.h"
# include "hstcalerr.h"

/* This routine gets the value of a header keyword, expected to be the
   name of a reference image or table.

   The RefFileInfo structure ref contains a list of keyword & filename
   pairs.  This list is first searched for a match with the specified
   keyword.  If the keyword is found in the list, the corresponding
   filename is returned as refname, and the header phdr is ignored.
   If the keyword is not found in the list, it is gotten from phdr.
   If it is also not found in phdr, refname is set to ""; this is not
   considered to be an error.

   Warren Hack, 1998 May 26:
	Revised for ACS.  No Major changes from STIS code.

*/

int GetRefName (RefFileInfo *ref, Hdr *phdr, char *keyword, char *refname) {

/* arguments:
RefFileInfo *ref  i: the list of keyword & filename pairs for reference files
Hdr *phdr         i: primary header
char *keyword     i: keyword name
char *refname     o: directory name and reference file name
*/

	extern int status;

	int foundit;			/* true if keyword was found in ref */

	void FindRefFile (RefFileInfo *, char *, char *, int *);
	int GetKeyStr (Hdr *, char *, int, char *, char *, int);

	FindRefFile (ref, keyword, refname, &foundit);

	if (!foundit) {
	    if (GetKeyStr (phdr, keyword, USE_DEFAULT, "", refname, SZ_LINE))
		return (status);
	}

	return (status);
}
