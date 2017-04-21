
# include <stdio.h>
# include <string.h>

# include "c_iraf.h"
# include "hstio.h"

# include "stis.h"
# include "calstis7.h"
# include "err.h"
# include "stisdq.h"
# include "stisdef.h"



/* This routine applies the incidence-angle correction, the MAMA
   offset correction, and the MSM shift to the dispersion coefficients.
*/

void AdjustDisp (StisInfo7 *sts, DispRelation *disp_y,
		double delta, InangInfo *iac, int *warn1, int *warn2) {

/* arguments:
StisInfo7 *sts        i: calibration switches and info (MSM offsets)
DispRelation *disp_y  io: dispersion relation interpolated at y
double delta          i: offset of slit from reference
InangInfo *iac        i: incidence-angle correction coeff
*/

	int ncoeff;
	int i;

	/* First apply the incidence-angle correction. */

	if (disp_y->ncoeff < iac->ncoeff1) {
	    ncoeff = disp_y->ncoeff;
	    if (!(*warn1)) {
		printf (
	"Warning  %d dispersion coefficients, but %d incidence-angle coeff.\n",
			disp_y->ncoeff, iac->ncoeff1);
		*warn1 = 1;	/* yes, a warning has been printed */
	    }
	} else {
	    ncoeff = iac->ncoeff1;
	}

	/* first coefficients */
	for (i = 0;  i < ncoeff;  i++)
	    disp_y->coeff[i] += iac->coeff1[i] * delta;

	/* second coefficients */
	if (iac->ncoeff2 > 0)
	    disp_y->coeff[0] += iac->coeff2[0] * delta;

	if (iac->ncoeff2 > 1)
	    disp_y->coeff[0] += iac->coeff2[1] * delta * delta;

	if (iac->ncoeff2 > 2 && !(*warn2)) {
	    printf (
	"Warning  %d incidence-angle second coefficents, limit is 2; \\\n",
			iac->ncoeff2);
	    printf (
	"Warning  the remaining coefficents will not be applied.\n");
	    *warn2 = 1;
	}

	/* Add the MSM offset (which will be zero if we're processing
	   a wavecal).
	*/
	disp_y->coeff[0] += sts->msm_slop[0];
}
