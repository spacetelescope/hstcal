# include "xtables.h"
# include "stis.h"
# include "calstis6.h"

/* This routine computes the heliocentric correction factor.  An
   observed wavelength should be multiplied by this factor to correct
   for the earth's orbital velocity around the sun:

   heliocentric wavelength = [1 + V/c * cos (theta)] * observed wavelength,

   where V is the speed of the earth in its orbit around the Sun,
   c is the speed of light, and theta is the angle between the earth's
   velocity vector and the direction toward the target.  The expression
   in brackets is returned as the factor.

   Note that the radial velocity is -V * cos (theta).




   Revision history:
   ----------------
   25 Feb 97  -  Borrowed from calstis7 (I.Busko)
   25 Feb 97  -  Rename routine to avoid conflict with cs7 (IB)
   16 Apr 97  -  Changes after code review (IB):
                 - moved define for PI constant to calstis6.h
                 - added some additional comments.
   17 Jan 03  -  Call RadialVel; change the calling sequence (PEH).
*/

double HelioFactor6 (StisInfo6 *sts, double *radvel) {

/* arguments:
StisInfo6 *sts   i: calibration switches and info
double *radvel   o: heliocentric radial velocity

The function value is the heliocentric correction factor.
*/

	/* Take midpoint of exposure as the current date. */

	*radvel = RadialVel (sts->ra_targ, sts->dec_targ,
		(sts->expstart + sts->expend) / 2.);

	return (1. - *radvel / HCFACTOR);
}
