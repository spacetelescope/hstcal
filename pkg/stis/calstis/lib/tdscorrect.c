# include <stdio.h>
# include "stistds.h"

/* This routine generates an array with time-dependent sensitivity 
   correction factors.

   Revision history:
   ----------------
   02 Jan 02  -  Implemented (I. Busko, from code by Phil Hodge)
   27 Dec 04  -  Add temperature dependence; change calling sequence
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
	double slope_j;

	for (i = 0;  i < tds->nwl;  i++) {

	    factor[i] = 1.;			/* initial value */

	    for (j = 0;  j < tds->nt;  j++) {

	        /* slopes are expressed in percentage. */

		slope_j = tds->slope[j][i] / DAYS_PER_YEAR / 100.0;
		if (j == tds->nt-1 || t <= tds->time[j+1]) {
		    factor[i] += (t - tds->time[j]) * slope_j;
		    break;
		} else {
		    factor[i] += (tds->time[j+1] - tds->time[j]) * slope_j;
		}
	    }

	    /* Include temperature dependence. */
	    if (tds->ref_temp > 0. && temperature > 0.) {
		factor[i] *= (1. +
			tds->temp_sens[i] * (temperature - tds->ref_temp));
	    }
	}
}
