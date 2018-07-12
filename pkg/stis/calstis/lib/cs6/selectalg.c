# include <stdio.h>
# include <string.h>
# include "xtables.h"

# include "stis.h"
# include "calstis6.h"
# include "hstcalerr.h"

/* 
   Selects the extraction algorithm, based on command-line input and
   the XTRACALG column in XTACTAB. This must be executed for each
   individual spectral order.




   Revision history:
   ----------------
   16 Sep 98  -  Implemented (IB)
   05 Jan 00  -  Scatterd light correction flag (IB)
*/

int SelectAlg (StisInfo6 *sts, XtractInfo *xtractab) {

	int streq_ic (char *, char *);

	/* Default algorithm is UNWEIGHTED, no scattered light correction. */
	sts->optimal = 0;
	sts->scatter = 0;

	/* Get algorithm from extraction table record. */
	if (streq_ic (xtractab->xtracalg, OPTIMAL))
	    sts->optimal = 1;
	else if (streq_ic (xtractab->xtracalg, SCATTER_CORR))
	    sts->scatter = 1;

	/* The table entry can be superseded by the command line choice. */
	if (streq_ic (sts->xtracalg, OPTIMAL)) {
	    sts->optimal = 1;
	    sts->scatter = 0;
	} else if (streq_ic (sts->xtracalg, SCATTER_CORR)) {
	    sts->optimal = 0;
	    sts->scatter = 1;
	}

	/* Test if reference files for optimal extraction are OK. 
           A similar test was already executed at startup time, but
           must be re-executed again here since each spectral order
           has an individual selection of extraction algorithms.

           This behavior was superseded by OPR 37810, but the code
           is left in place just in case (OPRs do come and go...).
        */
	if (sts->optimal) {

	    if (!((sts->pftab.exists == EXISTS_YES) &&
	          (sts->pxtab.exists == EXISTS_YES))) {
	        printf (
"Warning  Incomplete set of reference files for optimal extraction.\n");
	        return (CAL_FILE_MISSING);
	    }
	}

	return (0);
}


