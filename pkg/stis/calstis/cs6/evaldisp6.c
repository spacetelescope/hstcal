# include <xtables.h>

# include "../stis.h"
# include "calstis6.h"


/*  Compute dispersion relation and its derivative.




   Revision history:
   ----------------
   20 Feb 97  -  Borrowed from calstis7 (I.Busko)
   20 Feb 97  -  Added derivative computation (IB)
   10 Apr 97  -  Changes after code review (IB):
                 - optimized computations.
*/


void EvalDisp6 (DispRelation *disp, int sporder, double wl,
                      double *ix_r, double *dix_r) {

	double m;		/* spectral order number */
	double w2, m2, mw;

	m = (double) sporder;
	m2 = m * m;
	w2 = wl * wl;
	mw = m * wl;

	*ix_r = disp->coeff[0] +
	        disp->coeff[1] * mw +
	        disp->coeff[2] * m2 * w2 +
	        disp->coeff[3] * m +
	        disp->coeff[4] * wl +
	        disp->coeff[5] * m2 * wl +
	        disp->coeff[6] * m * w2 +
                disp->coeff[7] * m2 * m * w2 * wl;

	(*ix_r)--;		/* make it 0-indexed */

	*dix_r = disp->coeff[1] * m +
	         disp->coeff[2] * 2.0 * m2 * wl +
	         disp->coeff[4] +
	         disp->coeff[5] * m2 +
	         disp->coeff[6] * 2.0 * mw +
                 disp->coeff[7] * 3.0 * m2 * m * w2;
}

