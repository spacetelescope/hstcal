# include "stis.h"
# include "calstis7.h"

# define SPEED_OF_LIGHT 299792.458		/* c in km/sec */

/* This routine computes the heliocentric correction factor.  An
   observed wavelength should be multiplied by this factor to correct
   for the earth's orbital velocity around the sun:

   heliocentric wavelength = [1 + V/c * cos (theta)] * observed wavelength,

   where V is the speed of the earth in its orbit around the Sun,
   c is the speed of light, and theta is the angle between the earth's
   velocity vector and the direction toward the target.  The expression
   in brackets is returned as the factor.

   Phil Hodge, 2000 Aug 18:
	Add comments to EarthVel to explain how the velocity is converted
	to equatorial coordinates.

   Phil Hodge, 2003 Jan 17:
	Call RadialVel to do the work.  Change the calling sequence.
*/

double HelioFactor (StisInfo7 *sts, double *v_helio) {

/* arguments:
StisInfo7 *sts   i: calibration switches and info
double *v_helio  o: heliocentric radial velocity

The function value is the heliocentric correction factor.
*/

	/* Take midpoint of exposure as the current date. */

	*v_helio = RadialVel (sts->ra_targ, sts->dec_targ,
		(sts->expstart + sts->expend) / 2.);

	return (1. - *v_helio / SPEED_OF_LIGHT);
}
