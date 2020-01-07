# include <stdio.h>
# include <math.h>

# define MAX_COEFF 9		/* maximum number of coefficients */

/* These three are used by evalInvDisp. */
# define TOO_FAR (1024.)	/* pixels */
# define WL_STEP (100.)		/* Angstroms */
# define MAX_ITER 100

/* Evaluate the dispersion relation at wavelength wl and spectral order m.

   Phil Hodge, 2010 Nov 5:
	Rewrite evalInvDisp to search over the entire range of wavelengths
	for echelle data if the wavelength estimate is too far off.
*/

double evalDisp (double coeff[], int ncoeff, double m, double wl) {

/* arguments:
double coeff[]		i: coefficients of the dispersion relation
int ncoeff		i: number of elements (1-9) in coeff
double m		i: spectral order number (an integer)
double wl		i: wavelength for which we want pixel number

The function value is the zero-indexed pixel number corresponding to wl.
The function will return -1000001. if ncoeff is less than 0 or greater
than 9.
*/

	double ix_r;

	switch (ncoeff) {

	case 1:
	    ix_r = coeff[0];
	    break;
	case 2:
	    ix_r = coeff[0] +
		   coeff[1] * m * wl;
	    break;
	case 3:
	    ix_r = coeff[0] +
		   coeff[1] * m * wl +
		   coeff[2] * m * m * wl * wl;
	    break;
	case 4:
	    ix_r = coeff[0] +
		   coeff[1] * m * wl +
		   coeff[2] * m * m * wl * wl +
		   coeff[3] * m;
	    break;
	case 5:
	    ix_r = coeff[0] +
		   coeff[1] * m * wl +
		   coeff[2] * m * m * wl * wl +
		   coeff[3] * m +
		   coeff[4] * wl;
	    break;
	case 6:
	    ix_r = coeff[0] +
		   coeff[1] * m * wl +
		   coeff[2] * m * m * wl * wl +
		   coeff[3] * m +
		   coeff[4] * wl +
		   coeff[5] * m * m * wl;
	    break;
	case 7:
	    ix_r = coeff[0] +
		   coeff[1] * m * wl +
		   coeff[2] * m * m * wl * wl +
		   coeff[3] * m +
		   coeff[4] * wl +
		   coeff[5] * m * m * wl +
		   coeff[6] * m * wl * wl;
	    break;
	case 8:
	    ix_r = coeff[0] +
		   coeff[1] * m * wl +
		   coeff[2] * m * m * wl * wl +
		   coeff[3] * m +
		   coeff[4] * wl +
		   coeff[5] * m * m * wl +
		   coeff[6] * m * wl * wl +
		   coeff[7] * m * m * m * wl * wl * wl;
	    break;
	case 9:
	    ix_r = coeff[0] +
		   coeff[1] * m * wl +
		   coeff[2] * m * m * wl * wl +
		   coeff[3] * m +
		   coeff[4] * wl +
		   coeff[5] * m * m * wl +
		   coeff[6] * m * wl * wl +
		   coeff[7] * m * m * m * wl * wl * wl +
		   coeff[8] * m * m;
	    break;
	default:
	    ix_r = -1000000.;
	}

	ix_r--;		/* convert to zero indexing */

	return ix_r;
}

/* This function evaluates the dispersion relation for the prism.
   The input is pixel number (reference pixel units, but zero indexed),
   and the function value is the corresponding wavelength in Angstroms.
*/

double prismDisp (double coeff[], double ix_r) {

	double x;	/* pixel number minus first coefficient plus one */
	double wl;	/* wavelength */

	x = ix_r - coeff[0] + 1.;	/* convert to one indexing */

	wl = coeff[5] / x;
	wl = (coeff[4] + wl) / x;
	wl = (coeff[3] + wl) / x;
	wl = (coeff[2] + wl) / x;
	wl += coeff[1];

	return (wl);
}

/* Evaluate the derivative of the dispersion relation at wavelength wl and
   spectral order m.
*/

double evalDerivDisp (double coeff[], int ncoeff, double m, double wl) {

/* arguments:
double coeff[]		i: coefficients of the dispersion relation
int ncoeff		i: number of elements (1-9) in coeff
double m		i: spectral order number (an integer)
double wl		i: wavelength for which we want pixel number

The function value is the slope:  d(pixel) / d(wavelength)
The function will return +1000000. if ncoeff is less than 0 or greater
than 9.
*/

	double slope;

	switch (ncoeff) {

	case 1:
	    slope = 0.;
	    break;
	case 2:
	    slope = coeff[1] * m;
	    break;
	case 3:
	case 4:
	    slope = coeff[1] * m +
		    coeff[2] * 2. * m * m * wl;
	    break;
	case 5:
	    slope = coeff[1] * m +
		    coeff[2] * 2. * m * m * wl +
		    coeff[4];
	    break;
	case 6:
	    slope = coeff[1] * m +
		    coeff[2] * 2. * m * m * wl +
		    coeff[4] +
		    coeff[5] * m * m;
	    break;
	case 7:
	    slope = coeff[1] * m +
		    coeff[2] * 2. * m * m * wl +
		    coeff[4] +
		    coeff[5] * m * m +
		    coeff[6] * 2. * m * wl;
	    break;
	case 8:
	case 9:
	    slope = coeff[1] * m +
		    coeff[2] * 2. * m * m * wl +
		    coeff[4] +
		    coeff[5] * m * m +
		    coeff[6] * 2. * m * wl +
		    coeff[7] * 3. * m * m * m * wl * wl;
	    break;
	default:
	    slope = 1000000.;
	}

	return slope;
}

/* Invert the dispersion relation at pixel number 'pixel' for spectral
   order m.

   If the initial estimate (wl_estimate) of the wavelength is too far off,
   the estimate will be modified in linear steps until the corresponding
   pixel coordinate is no more than 1024 pixels from the specified pixel
   value.  Then Newton's method will be used to find the wavelength to
   within the specified tolerance.
*/

int evalInvDisp (double coeff[], int ncoeff, double m, double pixel,
		double wl_estimate, double tolerance, double *wl) {

/* arguments:
double coeff[]		i: coefficients of the dispersion relation
int ncoeff		i: number of elements (1-9) in coeff
double m		i: spectral order number (an integer)
double pixel		i: zero-indexed pixel coordinate (ref pixel size)
double wl_estimate	i: an estimate of the wavelength; CENWAVE could be
			   given as a starting value
double tolerance	i: the wavelength should be accurate to this value
double *wl		o: the wavelength corresponding to pixel; a value
			   of -1. indicates an error

The function value is status:
	0 is OK
	1 means the number of coefficients (ncoeff) was not valid
	2 means the the maximum number of iterations was exceeded
	3 means the slope of pixel vs wavelength was zero
*/

	double wl_low, wl_high;	/* range of wavelengths for searching */
	double wl_step;		/* wavelength step size for searching */
	double wl_test;		/* current value of wavelength */
	double last_wl_test;	/* last wavelength estimate */
	double x_test;		/* pixel at current value of wavelength */
	double slope;		/* slope of disp. rel. at current wavelength */
	double x_diff;		/* abs value of difference from pixel */
	double min_x_diff;	/* minimum value of x_diff */
	double wl_min_x_diff;	/* wavelength that gives minimum x_diff */
        double save_x_test;     /* x_test at min of x_diff */
	int niter;		/* iteration count */
	int done;		/* boolean flag for terminating a loop */
        int i;

	*wl = -1.;		/* initial value (bad) */

	if (ncoeff < 1 || ncoeff > MAX_COEFF)
	    return 1;

	wl_test = wl_estimate;		/* initial estimate */
	x_test = evalDisp (coeff, ncoeff, m, wl_test);
        save_x_test = x_test;
	slope = evalDerivDisp (coeff, ncoeff, m, wl_test);
	x_diff = fabs (x_test - pixel);
	min_x_diff = x_diff;		/* initial values */
	wl_min_x_diff = wl_test;

	if (slope <= 0. && m < 1.5) {
            /* For first-order data, the wavelength estimate should be good
                enough that the slope will at least be positive.
            */
            printf (
        "Warning:  (evalInvDisp) slope of dispersion relation is negative\n");
        }

	/* If the pixel corresponding to wl_estimate is too far off, or if
	   the slope (pixels per Angstrom) is negative, search for a better
	   estimate of the wavelength.  This section should only be relevant
           for echelle data, hence the test on m.
	*/
	if ((slope <= 0. || x_diff > TOO_FAR) && m > 1.5) {
	    /* The limits wl_low to wl_high were chosen to cover about 20 A
		more than the ranges given in the STIS Instrument Handbook:
		E230M:  cenwave = 1978 - 2707, wl = 1570 - 3111 Angstroms
		E230H:  cenwave = 1763 - 3012, wl = 1620 - 3150
		E140M:  cenwave = 1425,        wl = 1123 - 1710
		E140H:  cenwave = 1234 - 1598, wl = 1133 - 1699
	    */
	    /* search */
	    wl_low = 1123.;
	    wl_high = 3170.;
	    wl_step = WL_STEP;
	    if (slope <= 0.)
		min_x_diff = 1024. * 100;	/* a large number of pixels */
	    else
		min_x_diff = x_diff;
	    wl_min_x_diff = wl_estimate;
	    for (i = 0;  i < 2;  i++) {		/* just two iterations */
		for (wl_test = wl_low;  wl_test <= wl_high;
			wl_test += wl_step) {
		    slope = evalDerivDisp (coeff, ncoeff, m, wl_test);
		    if (slope > 0.) {
			x_test = evalDisp (coeff, ncoeff, m, wl_test);
			x_diff = fabs (x_test - pixel);
			if (x_diff < min_x_diff) {
			    min_x_diff = x_diff;
			    save_x_test = x_test;
			    wl_min_x_diff = wl_test;
			}
		    }
		}
		wl_low = wl_min_x_diff - wl_step;
		wl_high = wl_min_x_diff + wl_step;
		wl_step /= 10.;
	    }
	}

	done = 0;
	niter = 0;
	x_test = save_x_test;
	wl_test = wl_min_x_diff;
	last_wl_test = wl_min_x_diff;
	while (!done) {
	    slope = evalDerivDisp (coeff, ncoeff, m, wl_test);
	    if (slope == 0.) {
		*wl = wl_test;
		return 3;
	    }
	    wl_test = wl_test + (pixel - x_test) / slope;
	    if (fabs (wl_test - last_wl_test) < tolerance)
		done = 1;
	    else if (niter > MAX_ITER) {
		*wl = wl_test;
		return 2;
	    }
	    last_wl_test = wl_test;
	    x_test = evalDisp (coeff, ncoeff, m, wl_test);
	    niter++;
	}

	*wl = wl_test;

	return 0;
}
