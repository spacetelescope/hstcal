# include <stdio.h>
# include <string.h>

# include "xtables.h"
# include "hstio.h"

# include "stis.h"
# include "stisdef.h"
# include "hstcalerr.h"
# include "calstis6.h"

/*
   Checks status of reference file.

   Revision history:
   ----------------
   07 Aug 00  -  Moved from GetFlags6 (I.Busko)
*/

int GetCheckRef (Hdr *phdr, char *keyword, RefTab *table, int *calswitch,
                 int *missing, int no_fatal) {

	int status;

	/* The default name must be a null string, since this is the
	   way to tell GetPCT that the PC table does not exist.
	*/
	if (Get_KeyS (phdr, keyword, no_fatal, "", table->name, STIS_LINE))
	    return (KEYWORD_MISSING);

	/* TabPedigree opens the table to verify that it exists, and if so,
	   gets pedigree & descrip.
	*/
	if ((status = TabPedigree (table)))
	    return (status);

	if (table->exists != EXISTS_YES) {
	    (*missing)++;
	    if (no_fatal == 0)
	        printf ("ERROR    %s `%s' not found\n", keyword, table->name);
	}
	if (table->goodPedigree != GOOD_PEDIGREE)
	    *calswitch = DUMMY;

	return (0);
}
