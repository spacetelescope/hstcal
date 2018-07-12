# include <stdlib.h>
# include "c_iraf.h"
# include "hstio.h"

# include "stis.h"
# include "hstcalerr.h"
# include "stisdef.h"

/* This routine gets the value of a reference file keyword from the
   primary header and appends it to the list of reference file info.

   Phil Hodge, 1997 Dec 10:
	Function created.
*/

int GetNewRef (Hdr *phdr, char *keyword, RefFileInfo *ref) {

/* arguments:
Hdr *phdr         i: primary header of science file
char *keyword     i: keyword name for reference file
RefFileInfo ref  io: list of keyword,filename pairs
*/

	int status;
	char *filename;
	int use_default = 1;	/* return default if missing keyword */

	if ((filename = calloc (STIS_FITS_REC+1, sizeof(char))) == NULL)
	    return (status = OUT_OF_MEMORY);

	if ((status = Get_KeyS (phdr, keyword, use_default, "",
                                filename, STIS_FITS_REC)))
	    return (status);

	if ((status = NewRefFile (ref, keyword, filename)))
	    return (status);

	free (filename);

	return (0);
}
