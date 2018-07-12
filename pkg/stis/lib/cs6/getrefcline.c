# include <stdio.h>
# include <string.h>

# include "xtables.h"
# include "hstio.h"

# include "stis.h"
# include "stisdef.h"
# include "calstis6.h"
# include "hstcalerr.h"

static int CheckRef (char *, RefTab *);

/*
   Check if names of reference files were provided in the command line,
   and use them to supersede names provided in the header. Also check
   the existence of the reference files and get pedigree and descrip.

   For now on this routine only tests for the two files necessary for
   optimal extraction. It is an error if either, or both, of the files
   cannot be found.




   Revision history:
   ----------------
   16 Sep 98  -  Implemented (IB)
*/

int GetRefCommLine (StisInfo6 *sts) {

	int status;

	if ((status = CheckRef (sts->profilefile, &(sts->pftab))))
	    return (status);

	if ((status = CheckRef (sts->fluxfile, &(sts->pxtab))))
	    return (status);

	return (0);
}


static int CheckRef (char *filename, RefTab *table) {

	int status;

	if (filename[0] != '\0') {

	    /* If a file name was provided, supersede whatever was
               already in the table structure with the new name.
            */
	    strcpy (table->name, filename);

	    /* Check the table consistency. */
	    if ((status = TabPedigree (table)))
	        return (status);

	    /* If table cannot be found, abort. */
	    if (table->exists != EXISTS_YES) {
	        printf ("ERROR    %s not found\n", table->name);
	        return (CAL_FILE_MISSING);
	    }
	}

	return (0);
}
