# include <stdio.h>
# include <math.h>	/* fabs */
# include "stis.h"
# include "calstis12.h"
# include "cs12.h"		/* for interpolation options */
# include "hstcalerr.h"

/* This routine interpolates the shift in each axis.  The independent
   variable is time, given by midpt for the current science image
   and by the wav_midpt array for the wavecal images.

   If midpt is outside the array of times for the wavecal, the shift
   will be set to the values at the nearer endpoint.

   Phil Hodge, 1998 Oct 5:
	Change status value 1343 to GENERIC_ERROR_CODE.

   Phil Hodge, 2000 Jan 19:
	Get the shift in this routine using interpolation option w_option;
	change calling sequence (only good data, only one shift).
*/

int MatchWav (double wav_midpt[], double wav_shift[], int n,
		double midpt, int w_option, double *shift) {

/* arguments:
double wav_midpt[]   i: array of times for good shifts
double wav_shift[]   i: array of good shifts
int n                i: size of arrays
double midpt         i: midpoint of science data exposure, MJD
int w_option         i: interpolation option
double *shift        o: interpolated shift
*/

	double dt;	/* difference in times */
	double min_dt;	/* minimum value of dt */
	int i;		/* array index */
	int min_i;	/* value of i that gives minimum dt */
	double p, q;	/* for linear interpolation of shift */

	/* If midpt is out of range, return shift at nearest endpoint. */

	if (n < 1) {
	    *shift = UNDEFINED_SHIFT;
	    return (0);
	}

	if (n == 1 || midpt <= wav_midpt[0]) {

	    *shift = wav_shift[0];
	    return (0);
	}

	if (midpt >= wav_midpt[n-1]) {

	    *shift = wav_shift[n-1];
	    return (0);
	}

	/* Return shift at nearest point, or interpolate. */

	if (w_option == STIS_NEAREST) {

	    min_dt = fabs (wav_midpt[0] - midpt);
	    min_i = 0;
	    for (i = 0;  i < n;  i++) {
		dt = fabs (wav_midpt[i] - midpt);
		if (dt < min_dt) {
		    min_dt = dt;
		    min_i = i;
		}
	    }
	    *shift = wav_shift[min_i];

	} else if (w_option == STIS_LINEAR) {

	    for (i = 0;  i < n - 1;  i++) {
		if (midpt >= wav_midpt[i] &&
		    midpt <= wav_midpt[i+1]) {
		    p = (midpt - wav_midpt[i]) /
			(wav_midpt[i+1] - wav_midpt[i]);
		    q = 1. - p;
		    *shift = q * wav_shift[i] +
			     p * wav_shift[i+1];
		    break;
		}
	    }

	} else {

	    printf ("ERROR    interpolation option not supported\n");
	    return (GENERIC_ERROR_CODE);
	}

	return (0);
}
