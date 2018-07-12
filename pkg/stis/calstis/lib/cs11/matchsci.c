# include <stdio.h>
# include <math.h>	/* fabs */

# include "stis.h"
# include "calstis11.h"

/* This routine finds which imset in the science file is the most
   appropriate to use for subtracting from the wavecal.

   In the current implementation, the imset that is closest in
   time to the wavecal is selected.  That is, option is ignored.

   Note:  It is the imset number that is returned as extver, not
   the array index.

   Phil Hodge, 2000 July 27:
	Add a check on the value of option.
*/

int MatchSci (StisInfo11 *wavecal, StisInfo11 *scidata,
		int option, int *extver) {

/* arguments:
StisInfo11 *wavecal  i: info for wavecal
StisInfo11 *scidata  i: info for science data
int option           i: option for selecting the most appropriate imset
int *extver          o: imset number of best science image to use (one indexed)
*/

	double dt;	/* difference in times */
	double min_dt;	/* minimum value of dt */
	int i;		/* array index */
	int min_i;	/* value of i that gives minimum dt */

	if (option != STIS_NEAREST)
	    printf ("Warning  option = %d will be ignored\n", option);

	if (scidata->nimages == 1) {

	    *extver = 1;

	} else {

	    min_dt = fabs (wavecal->midpt[0] - scidata->midpt[0]);
	    min_i = 0;
	    for (i = 0;  i < scidata->nimages;  i++) {

		dt = fabs (wavecal->midpt[0] - scidata->midpt[i]);
		if (dt < min_dt) {
		    min_dt = dt;
		    min_i = i;
		}
	    }

	    *extver = min_i + 1;	/* extver is one indexed */
	}

	return (0);
}
