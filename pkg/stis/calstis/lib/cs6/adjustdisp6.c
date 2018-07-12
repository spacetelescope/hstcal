# include <stdio.h>
# include "xtables.h"

# include "stis.h"
# include "calstis6.h"

/* This routine applies the incidence-angle correction and
   the MSM shift to the dispersion coefficients.

   Revision history:
   ----------------
   20 Feb 97  -  Borrowed from calstis7 (I.Busko)
   24 Feb 97  -  Rename routine to avoid conflict with cs7 (IB)
   10 Apr 97  -  Changes after code review (IB):
                 - re-formatted printf statements.
   01 May 97  -  Fixed dither value (IB)
   08 May 97  -  Conform to new _trl standard (IB)
*/

void AdjustDisp6 (StisInfo6 *sts, DispRelation *disp_y,
                  double delta, InangInfo *iac) {

/* arguments:
StisInfo6 *sts        i: calibration switches and info (MSM offsets)
DispRelation *disp_y  io: dispersion relation
double delta          i: offset of slit from reference
InangInfo *iac        i: incidence-angle correction coeff
*/

	int ncoeff;
	int i;

	/* First apply the incidence-angle correction. */

	if (disp_y->ncoeff < iac->ncoeff1) {
	    printf ("Warning  %d dispersion coefficients, ",disp_y->ncoeff);
	    printf ("but %d incidence-angle coeff.\n",iac->ncoeff1);
	    ncoeff = disp_y->ncoeff;
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

	if (iac->ncoeff2 > 2) {
	    printf ("Warning  %d incidence-angle second ", iac->ncoeff2);
	    printf ("coefficents, limit is 2;\n");
	    printf ("Warning  the remaining coefficents ");
	    printf ("will not be applied.\n");
	}

	/* Add the auto-wavecal offset (which will be zero if we're
	   processing a wavecal).
	*/
	disp_y->coeff[0] += sts->msm_offset[0];
}
