# include <stdio.h>
# include <string.h>
# include "xtables.h"

# include "stis.h"
# include "calstis6.h"
# include "hstcalerr.h"

/* 
   This routine checks for the existence of both reference files 
   associated with optimal extraction: the profile file and the
   flux file. But the test is performed only if the extraction
   algorithm was explictly defined as "OPTIMAL". It is an error
   if either, or both, files do not exist in this case.

   This routine must be executed after both the input header and
   the command line switches were scanned for the relevant information.




   Revision history:
   ----------------
   16 Sep 98  -  Implemented (IB)
*/

int CheckOptimal (StisInfo6 *sts) {

	if (streq_ic (sts->xtracalg, OPTIMAL)) {

	    if ((sts->pftab.exists == EXISTS_NO) ||
	        (sts->pxtab.exists == EXISTS_NO)) {
	        printf (
"ERROR    Incomplete set of reference files for optimal extraction.\n");
	        return (CAL_FILE_MISSING);
	    }
	}
	return (0);
}


