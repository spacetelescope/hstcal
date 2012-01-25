# include <math.h>

# define PI (3.1415926535897932384626433)
# define SEC_PER_DAY (86400.)

/* This function returns the average orbital Doppler shift in low-res
   pixels.  If the exposure time is very short, this is just the shift
   at the middle of the exposure.  Otherwise, the integral of the shift
   over the entire exposure is computed, and the average shift is the
   integral divided by the time interval (exposure time).

   If doppmag is zero, zero will be returned.

   Phil Hodge, 2006 Apr 10:
	Function created.
*/

double orbitalDopp (double t1, double t2,
	double doppzero, double orbitper, double doppmag) {

/* arguments:
double t1, t2           i: MJD at beginning and end of exposure
double doppzero         i: MJD when orbital Doppler shift is zero (increasing)
double orbitper         i: orbital period (seconds)
double doppmag          i: magnitude of Doppler shift (high-res pixels)

The function value is the average orbital Doppler shift in low-res pixels.
*/

	double shift;

	if (doppmag <= 0.)
	    return 0.;

	orbitper /= SEC_PER_DAY;

	if (fabs (t2 - t1) < 1.e-4) {		/* 8.64 seconds */
	    double t = (t1 + t2) / 2.;
	    shift = doppmag * sin ((t - doppzero) * 2*PI / orbitper);
	} else {
	    shift = doppmag * orbitper / (2.*PI) *
		    (cos ((t1 - doppzero) * 2.*PI/orbitper) -
		     cos ((t2 - doppzero) * 2.*PI/orbitper)) / (t2 - t1);
	}

	/* convert to low-res pixels */
	return (shift / 2.);
}
