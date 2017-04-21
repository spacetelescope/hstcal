# include <stdlib.h>		/* calloc */
# include "c_iraf.h"
# include "hstio.h"

# include "stis.h"
# include "err.h"
# include "stisdef.h"

/* This routine gets the value of a header keyword, expected to be the
   name of a reference image or table.

   The RefFileInfo structure ref contains a list of keyword & filename
   pairs.  This list is first searched for a match with the specified
   keyword.  If the keyword is found in the list, the corresponding
   filename is returned as refname, and the header phdr is ignored.
   If the keyword is not found in the list, it is gotten from phdr.
   If it is also not found in phdr, refname is set to ""; this is not
   considered to be an error.

   Phil Hodge, 1997 Dec 10:
	Include ref in calling sequence, and search there first;
	missing keyword is not a fatal error, just set the name to null.
*/

int GetRefName (RefFileInfo *ref, Hdr *phdr, char *keyword, char *refname) {

/* arguments:
RefFileInfo *ref  i: the list of keyword & filename pairs for reference files
Hdr *phdr         i: primary header
char *keyword     i: keyword name
char *refname     o: directory name and reference file name
*/

	int status;

	int foundit;			/* true if keyword was found in ref */
	int use_default = 1;		/* missing keyword is not an error */

	FindRefFile (ref, keyword, refname, &foundit);

	if (!foundit) {
	    if ((status = Get_KeyS (phdr, keyword,
                                    use_default, "", refname, STIS_LINE)))
		return (status);
	}

	return (0);
}
