# include <stdio.h>
# include <stdlib.h>		/* for calloc */
# include <math.h>		/* for sin */
# include "hstcalerr.h"

# define PI    3.1415926535897932384626433
# define DT     (1./86400.)	/* time step (one second) in days */

/* nearest integer function */
# define NINT(x)  ((x >= 0.) ? (int) (x + 0.5) : (int) (x - 0.5))

/* This routine computes the Doppler smearing function.
   The calling routine must provide an array ds which may be as large
   as (2 * (doppmag+1) + 1).

Doppler smearing:

          --------- photon would have hit here if no Doppler shift
         |
         |        --------- photon was observed here due to Doppler shift
         V       V
    ___ ___ ___ ___ ___
   |   |   |   |   |   |
   |___|___|___|___|___|
   |   |   |   |   |   |
   |___|___|___|___|___|
   |   | x |   | x |   |
   |___|___|___|___|___|
   |   |   |   |   |   |
   |___|___|___|___|___|
   |   |   |   |   |   |
   |___|___|___|___|___|

The picture above represents a 5 x 5 pixel section of an image, with
dispaxis = 1 (X axis) and dispsign = +1 (wavelength increases to the
right).  The Doppler shift is two pixels to the red; this could
correspond to a time of observation within one half orbital period
following doppzero (see below).

Before applying a flat field (or dark, etc.) to a science file, the
reference file should be shifted two pixels to the left.

Doppzero is a time when the Doppler shift is zero, and HST is near its
closest approach to the target.

   Phil Hodge, 1998 Mar 13:
	Change calling sequence, replacing s0 with d0;
	don't extract the non-zero subset of Doppler smearing array.

   Phil Hodge, 1998 Apr 2:
	Change the sign of the Doppler shift.

   Phil Hodge, 1998 Oct 5:
	Change status value 1001 to GENERIC_ERROR_CODE.
*/

int MakeDopp (double doppzero, double doppmag, double orbitper,
		double expstart, double exptime, int dispsign,
		float *ds, int *nds, int *d0) {

/* arguments:
double doppzero    i: time (MJD) of zero Doppler shift, HST closer to target
double doppmag     i: magnitude of the Doppler shift (hi-res pixels)
double orbitper    i: HST orbital period (seconds)
double expstart    i: start time (MJD) of exposure
double exptime     i: exposure time (seconds)
int dispsign       i: +1 or -1; +1 --> wavelength increases with pixel number
float *ds          o: Doppler smearing array
int *nds           o: size of Doppler smearing array
int *d0            o: index in ds corresponding to zero Doppler shift;
			d0 = (nds - 1) / 2
*/

	int n;			/* amplitude of ds array */
	double t;		/* time at increments of DT during exposure */
	double shift;		/* Doppler shift in pixels at time t */
	float sum;		/* for normalizing smearing function */
	int i;

	if (orbitper <= 0.) {
	    printf ("ERROR    ORBITPER = %.6g.\n", orbitper);
	    return (GENERIC_ERROR_CODE);
	}
	if (doppmag <= 0.) {
	    printf ("ERROR    DOPPMAG = %.6g.\n", doppmag);
	    return (GENERIC_ERROR_CODE);
	}

	/* Array element n is the middle of the relevant portion of ds;
	   it's the element corresponding to zero Doppler shift.
	*/
	n = (int) (doppmag + 1.);		/* amplitude */
	*d0 = n;
	*nds = 2 * n + 1;

	/* Convert units from seconds to days. */
	orbitper /= 86400.;
	exptime /= 86400.;

	/* Loop through the time interval of the exposure at time steps
	   of DT, compute the Doppler shift at each time in units of
	   high-res pixels, and increment the corresponding array element
	   in ds.
	   shift is the number of pixels to add to the observed X
	   coordinate to get the location where the photon would have
	   been observed if the Doppler shift were zero.
	*/
	sum = 0.;
	for (t = expstart-doppzero; t <= expstart-doppzero+exptime; t += DT) {

	    if (dispsign == 1)
		shift = -doppmag * sin (t * 2. * PI / orbitper);
	    else
		shift = doppmag * sin (t * 2. * PI / orbitper);

	    i = NINT(shift);		/* nearest integer */

	    ds[n+i] += 1.;	/* note i = 0 when Doppler shift is zero */
	    sum += 1.;
	}

	/* Normalize ds so the sum of the values is one. */
	for (i = 0;  i < *nds;  i++)
	    ds[i] /= sum;

	return (0);
}
