# include <stdio.h>
# include <math.h>		/* for sqrt */

# include "hstio.h"
# include "stis.h"
# include "calstis1.h"

/* This routine checks whether the error array is all zero, and if so,
   a simple noise model is evaluated and assigned to the error array.

   For the CCD, the noise model is:

	sigma = sqrt ((I - bias) * gain + readnoise^2)  electrons

	      = sqrt ((I - bias) / gain + (readnoise / gain)^2)  dn,

   where I is the science data value in dn, bias is in dn, readnoise is
   the readout noise in electrons, gain is the CCD gain in electrons per dn.
   The value of sigma in dn is what is assigned to the error array.

   We'll subtract sts->err_init_bias as the bias in the above expression.
   This value will either be sts->ccdbias (if blevcorr has not been done)
   or zero (if blevcorr has been done).

   For the MAMAs, the noise model is:

	sigma = sqrt (I).

   Phil Hodge, 1998 June 29:
	Set minimum error to readnoise/gain for CCD, or zero for MAMA.
	Include a separate section for MAMA, rather than relying on
	CCD parameters having been set to appropriate values during
	initialization.

   Phil Hodge, 1999 Nov 2:
	For CCD data, only subtract sts->ccdbias if blevcorr has not
	been done.
*/

int doNoise (StisInfo1 *sts, SingleGroup *x, int *done) {

/* arguments:
StisInfo1 *sts    i: calibration switches and info
SingleGroup *x   io: image to be calibrated; written to in-place
int *done         o: true if we actually did assign error array values
*/

	float bias;		/* bias level to subtract (dn) */
	float gain;		/* electrons per dn */
	float rn2;		/* readnoise^2 (in el^2) */
	float value;		/* dn - bias * gain (i.e. signal in el) */
	int i, j;

	*done = 0;				/* initial value */

	/* First check for a dummy error array.  If it's not dummy,
	   we just return without doing anything.
	*/
	for (j = 0;  j < x->err.data.ny;  j++) {
	    for (i = 0;  i < x->err.data.nx;  i++) {
		if (Pix (x->err.data, i, j) != 0.) {
		    return (0);		/* not a dummy error array */
		}
	    }
	}

	if (sts->ncombine > 1) {
	    printf (
	"Warning  NCOMBINE > 1 before the error array was initialized.\n");
	}

	/* Evaluate the noise model for each pixel. */
	if (sts->detector == CCD_DETECTOR) {

	    bias = sts->err_init_bias;
	    gain = sts->atodgain;
	    rn2 = sts->readnoise * sts->readnoise;

	    /* Subtract the bias and convert to electrons,
	       then include readout noise and convert back to dn.
	    */
	    for (j = 0;  j < x->sci.data.ny;  j++) {
		for (i = 0;  i < x->sci.data.nx;  i++) {
		    value = (Pix (x->sci.data, i, j) - bias) * gain;
		    if (value > 0.)
			Pix (x->err.data, i, j) = sqrt (value + rn2) / gain;
		    else
			Pix (x->err.data, i, j) = sts->readnoise / gain;
		}
	    }

	} else {				/* either MAMA detector */

	    for (j = 0;  j < x->sci.data.ny;  j++) {
		for (i = 0;  i < x->sci.data.nx;  i++) {
		    value = Pix (x->sci.data, i, j);
		    if (value > 0.)
			Pix (x->err.data, i, j) = sqrt (value);
		    else
			Pix (x->err.data, i, j) = 0.;	/* minimum ERR */
		}
	    }
	}

	*done = 1;
	return (0);
}
