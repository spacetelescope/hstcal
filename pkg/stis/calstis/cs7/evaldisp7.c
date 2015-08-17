# include <stdio.h>

# include "stis.h"
# include "calstis7.h"

# define PIXEL_TOLERANCE  (0.01)	/* for solving disp rel for prism */

static double SolveDisp (DispRelation *, double, double);


/* Functions that handle the dispersion relation.

   Phil Hodge, before 1997 Dec 9:
	Functions created.

   Phil Hodge, 2000 Mar 10:
	Include option to evaluate prism dispersion relation, using
	functions SolveDisp and PrismDisp.

   Phil Hodge, 2000 May 31:
	In EvalDisp, add a seventh term (cubic) in the dispersion relation.

   Phil Hodge, 2001 Feb 22:
	Move prism macros to ../stis.h.

   Phil Hodge, 2001 May 2:
	Move the call to GetDisp from here to Do2Dx, and add disp to
	calling sequence.

   Ivo Busko, 2002 Mar 5:
   	Move EvalDisp to its own file. It is required by the blaze
   	shift algorithm.

   Ivo Busko, 2002 Mar 6:
   	Add GetWavelength function.

   Phil Hodge, 2010 July 30:
	Rename EvalDisp to EvalDisp7, and change the file name to evaldisp7.c.
	Delete PrismDisp; call evalDisp or prismDisp.
*/

void EvalDisp7 (DispRelation *disp_y, int sporder, double wl, int disp_type,
               double *ix_r, double *dix_r) {

	double m;		/* spectral order number */

	if (disp_type == GRATING_DISP) {

	    m = (double) sporder;

	    *ix_r = evalDisp (disp_y->coeff, disp_y->ncoeff, m, wl);

	} else if (disp_type == PRISM_DISP) {

	    *ix_r = SolveDisp (disp_y, wl, PIXEL_TOLERANCE);
	}

	*dix_r = 1.;		/* stub */
}


/* This function evaluates the wavelength, given the pixel number.
   It was put in place to support the blaze shift algorithm. There
   is some code duplication (from SolveDisp); this is unavoidable
   since the approach when introducing new functionality is to keep
   existing (and tested) code untouched. So no refactoring can take
   place here.
*/

double GetWavelength (DispRelation *disp, int sporder, double ix,
                      int disp_type, double low, double high) {

	double wl_high, wl_low, wl_test, x_test;
	double m;

	if (disp_type == GRATING_DISP) {

	    m = (double)sporder;
	    wl_low  = low;
	    wl_high = high;

	    while (wl_high - wl_low > 1.E-5) {

	        wl_test = (wl_low + wl_high) / 2.;
		x_test = evalDisp (disp->coeff, disp->ncoeff, m, wl_test);

	        if (ix < x_test)
		    wl_high = wl_test;
	        else
		    wl_low = wl_test;
	    }

	    wl_test = (wl_low + wl_high) / 2.;

	} else if (disp_type == PRISM_DISP) {

	    wl_test = prismDisp (disp->coeff, ix);

	}

	return (wl_test);
}

/* This function solves the dispersion relation for a prism using a binary
   search.  The function value is the reference pixel number corresponding
   to wl.  If wl is lower or higher than would appear on a full-frame image,
   the pixel number will be 100 pixels out of bounds.

   It is assumed here that the wavelengths increase with pixel number.
*/

static double SolveDisp (DispRelation *disp_y, double wl, double tol) {

	double x_low, x_high;	/* pixel numbers at ends of test range */
	double x_test;		/* pixel number at middle of test range */
	double wl_test;		/* wavelength at x_test */

	/* reference pixel range */
	x_low = 0.;
	x_high = 1023.;

	if (wl < prismDisp (disp_y->coeff, x_low))
	    return (x_low - 100.);		/* out of bounds */

	if (wl > MAX_PRISM_WAVELENGTH)
	    return (x_high + 100.);		/* out of bounds */

	while (x_high - x_low > tol) {

	    x_test = (x_low + x_high) / 2.;
	    wl_test = prismDisp (disp_y->coeff, x_test);

	    if (wl < wl_test)
		x_high = x_test;
	    else
		x_low = x_test;
	}

	x_test = (x_low + x_high) / 2.;

	return (x_test);
}
