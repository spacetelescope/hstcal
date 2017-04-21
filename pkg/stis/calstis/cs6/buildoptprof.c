# include <stdio.h>
# include <string.h>
# include <math.h>

# include "xtables.h"
# include "hstio.h"

# include "stis.h"
# include "calstis6.h"
# include "err.h"


/*
   This routine is used to build the actual cross-dispersion profiles
   used in optimal extraction. It basically fills up the subsampled
   profile 2-D array by interpolating the subsampled profiles read
   from OPROFTAB.

   This routine assumes that memory for the 2-D subsampled profile that
   lives in the StisInfo6 structure was previously allocated.



   Revision history:
   ----------------
   18 Sep 98  -  Implemented (I.Busko)
   27 Oct 98  -  Warn if any profile value < 0 or extremes > 0.1 flux (IB)
   06 Dec 00  -  Subsampled profiles (IB)
*/

int BuildOptProf (StisInfo6 *sts, ProfileArray **profa) {

/* arguments:
StisInfo6 *sts       io: structure where the 2-D output profile array lives
ProfileArray **profa i:  list with profiles read from OPROFTAB
*/
	int status;

	ProfileArray *profa_x;		/* profile interpolated at x */
	int i, j, negflag, outerflag;
	double sum;

	int InterpProfile (ProfileArray **, int, ProfileArray **, double *);
	void FreeProfileArray (ProfileArray **);

	profa_x = NULL;

	/* Scan entire X length of profile 2-D array. */

	negflag = 0;
	outerflag = 0;
	for (i = 0; i < sts->profile_x ; i++) {

	    /* Get interpolated profile at current X coordinate. */

	    if ((status = InterpProfile (profa, i, &profa_x,
                                         &(sts->profile_offset[i]))))
	        return (status);

	    for (j = 0; j < profa_x->npts; j++)
	        sts->subprof[i][j] += profa_x->prof[j];

	    /* Normalize to unit area and enforce positivity. */

	    sum = 0.0;
	    for (j = 0; j < sts->subprof_size; j++) {
	        if (sts->subprof[i][j] < 0.0) {
	            sts->subprof[i][j] = 0.0;
	            negflag = 1;
	        }
	        sum += sts->subprof[i][j];
	    }
	    if (sum > 0.0)
	        for (j = 0; j < sts->subprof_size; sts->subprof[i][j++] /= sum);
	    /* this test must be revised to account for subsampling. */
	    if ((sts->subprof[i][0] > 0.1) ||
	        (sts->subprof[i][(sts->subprof_size)-1] > 0.1))
	        outerflag = 1;

	}

	FreeProfileArray (&profa_x);

	/* Warn of abnormal conditions. */

	if (negflag)
	    printf ("Warning  Negative value(s) in input profile.\n");
	if (outerflag)
printf ("Warning  Outer pixels in profile have more than 0.1 of the total flux.\n");

	return (0);
}
