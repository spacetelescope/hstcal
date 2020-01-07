# include <stdio.h>
# include "stistds.h"

/* This routine generates an array with time-dependent sensitivity 
   correction factors.  The values in the returned array correspond to
   values in the array tds->wl of wavelengths.

   Revision history:
   ----------------
   02 Jan 02  -  Implemented (I. Busko, from code by Phil Hodge)
   27 Dec 04  -  Add temperature dependence; change calling sequence
			(Phil Hodge)
   19 Sep 11  -  Also support COS-like TDS table format, except that the
		 extrapolation is done with the first or last slope, as
		 opposed to the COS extrapolation with zero slope.
			(Phil Hodge)
*/

void TdsCorrection (TdsInfo *tds, double t, double temperature,
		double *factor) {
/* arguments:
TdsInfo *tds;        i: time and temperature dependent sensitivity info
double t;            i: time of exposure (MJD)
double temperature;  i: detector temperature (degrees Celsius)
double *factor;      o: array with correction factors (size = tds->nwl)
*/
	int i, j;
	int nt = tds->nt;	/* number of times in the time array */

	if (tds->format == COS_TDS_FORMAT) {

	    double delta_t;	/* time (days) since reference time */

	    /* Find the time interval (index j) that includes the time of
	       observation.
	    */
	    if (nt == 1 || t >= tds->time[nt-1]) {
		j = nt - 1;
	    } else {
		for (j = 0;  j < nt-1;  j++) {
		    if (t < tds->time[j+1])
			break;
		}
	    }

	    /* tds->time and ref_time are in days.
	       tds->slope is in percent per year.
	    */
	    delta_t = (t - tds->ref_time) / DAYS_PER_YEAR;

	    for (i = 0;  i < tds->nwl;  i++) {
		factor[i] = delta_t * tds->slope[j][i] / 100.0
				    + tds->intercept[j][i];
	    }

	} else {		/* original STIS TDS table format */

	    double slope_j;

	    for (i = 0;  i < tds->nwl;  i++) {

		factor[i] = 1.;			/* initial value */

		for (j = 0;  j < nt;  j++) {

		    /* slopes are expressed in percentage. */

		    slope_j = tds->slope[j][i] / DAYS_PER_YEAR / 100.0;
		    if (j == nt-1 || t <= tds->time[j+1]) {
			factor[i] += (t - tds->time[j]) * slope_j;
			break;
		    } else {
			factor[i] += (tds->time[j+1] - tds->time[j]) * slope_j;
		    }
		}
	    }
	}

	/* Include temperature dependence. */
	for (i = 0;  i < tds->nwl;  i++) {
	    if (tds->ref_temp > 0. && temperature > 0.) {
		factor[i] *= (1. +
			tds->temp_sens[i] * (temperature - tds->ref_temp));
	    }
	}
}
