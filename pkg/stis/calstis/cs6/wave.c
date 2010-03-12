# include <stdio.h>
# include <string.h>
# include <math.h>

# include <xtables.h>

# include "../stis.h"
# include "../stisdq.h"
# include "calstis6.h"
# include "../stissizes.h"

# define MAXIT 	   100     /* maximum number of Newton-Raphson iterations */
# define PIX_RANGE 10.     /* range in pixels where to look for wavelength */
# define GRATING_DISP  1
# define PRISM_DISP    2

static double CenWave (StisInfo6 *, DispRelation *, int, double *, double *);
static double Newton (DispRelation *, int, double, double, double, double);
static double PrismDisp (DispRelation *, double);


/*
   Compute wavelength array. In the grating case, this is done by solving
   the inverse dispersion relation.


   Revision history:
   ----------------
   25 Feb 97  -  Implemented (I.Busko)
   10 Apr 97  -  Changes after code review (IB):
                 - removed EvalDisp from Newton argument list
   30 Jan 01  -  Prism support (IB)
   26 Feb 01  -  Stop at MAX_PRISM_WAVELENGTH in prism spectra (IB)
   07 Nov 02  -  Move MAX_PRISM_WAVELENGTH to calstis6.h (IB)
*/

void Wave (StisInfo6 *sts, DispRelation *disp, RowContents *rowc) {

/* argument:
StisInfo6 *sts         i:  calibration switches and info
DispRelation *disp     i:  dispersion relation
RowContents rowc       io: row contents
*/

	double rpix;		/* reference pixel */
	double cenwave;		/* central wavelength */
	double wcenter;		/* first guess wavelength */
	double wlow, wupp;	/* wavelength interval for searching root */
	double w0, w1;		/* approximate linear dispersion coeffs. */
        int disp_type;          /* disperser is grating or prism */
	int    i, stopped;

	/* Find disperser type. */
        if (strcmp (sts->opt_elem, "PRISM") == 0)
            disp_type = PRISM_DISP;
        else
            disp_type = GRATING_DISP;

	if (disp_type == GRATING_DISP) {

	    /* Find central wavelength. This must be done here because the
               cenwave value read from the dispersion reference table is
               not correct except for order 1.
            */
	    cenwave = CenWave (sts, disp, rowc->sporder, &rpix, &w1);

	    /* Compute approximate linear dispersion relation. This is
               used to generate bracketing intervals for the root-finding
               algorithm.
            */
	    w1 = 1.0 / w1;
	    w0 = cenwave - w1 * rpix;
	}

	/* Loop over wavelength array. */
        stopped = 0;
	for (i = 0;  i < rowc->npts;  i++) {

	    /* Translate image pixel index into reference pixel index. */
	    rpix = (i - sts->ltv[0]) / sts->ltm[0];

	    if (disp_type == GRATING_DISP) {

	        /* Compute wavelength interval given by approximate
                   wavelength solution using linear dispersion relation.
                   Upper and lower bounds are computed by stepping 20
                   reference pixels in both directions from the center
                   of the interval. For safety, the interval size is chosen
                   to be at least twice as large as the precision used by
                   routine CenWave when finding the central wavelength
                   estimate.
                */
	        wcenter = w0 + w1 * rpix;
	        wlow    = wcenter - w1 * 2.0 * PIX_RANGE;
	        wupp    = wcenter + w1 * 2.0 * PIX_RANGE;

	        /* Solve for wavelength. */
	        rowc->wave[i] = Newton (disp, rowc->sporder, wlow, wupp,
                                        rpix, 0.01);

	    } else if (disp_type == PRISM_DISP) {
                if (stopped) {
	            rowc->wave[i] = MAX_PRISM_WAVELENGTH;
	            rowc->dq[i]  |= CALIBDEFECT ;
                    continue;
	        }
	        rowc->wave[i] = PrismDisp (disp, rpix);
	        if (rowc->wave[i] > MAX_PRISM_WAVELENGTH) {
	            rowc->wave[i] = MAX_PRISM_WAVELENGTH;
	            rowc->dq[i]  |= CALIBDEFECT ;
	            stopped = 1;
	        }
	    }

	    /* Apply heliocentric correction. */
	    if (sts->heliocorr == PERFORM)
		rowc->wave[i] *= sts->hfactor;
	}
}


/*  Find central wavelength. */

static double CenWave (StisInfo6 *sts, DispRelation *disp, int sporder,
                       double *pix, double *w1) {

	double cenwave;		/* central wavelength */
	double cpix;		/* central pixel */
	double step;		/* wavelength search step in Angstrom */
	int    stepsign;	/* step sign */
	int    change;		/* signals change of direction */

	void EvalDisp6 (DispRelation *, int, double, double *, double *);

	/* Initialize. */
	step     = 2.0;
	stepsign = 1;
	change   = 0;
	cenwave  = (double)(sts->cenwave);
	EvalDisp6 (disp, sporder, cenwave, pix, w1);
	if (sts->detector == CCD_DETECTOR) cpix = (double)CCD_NPIX_X  / 2.0;
	else                               cpix = (double)MAMA_NPIX_X / 2.0;

	/* Search. Stop when pixel falls within range. */
	while  ((*pix < (cpix - PIX_RANGE)) || (*pix > (cpix + PIX_RANGE))) {

	    /* Set direction of search for this step. */
	    if (*pix < cpix) {
	        if (stepsign == -1)
	            change = 1;
	        stepsign = +1;
	    } else {
	        if (stepsign == 1)
	            change = 1;
	        stepsign = -1;
	    }

	    /* If direction changed, slow down. */
	    if (change) {
	        step /= 2.0;
	        change = 0;
	    }

	    /* Update cenwave and corresponding pixel number. */
	    cenwave += step * stepsign;
	    EvalDisp6 (disp, sporder, cenwave, pix, w1);
	}

	return (cenwave);
}


/*  Solve dispersion relation using the Newton-Raphson method.
    This is adapted from "Numerical Recipes" to take into account
    the particular functional form, data types and structures used here.
*/

static double Newton (DispRelation *disp, int sporder, double wl1, double wl2,
                      double rpix, double wl_acc) {

	int   j;
	double df, dwl, dwlold, f, fh, fl;
	double temp, wlh, wll, rts;

	void EvalDisp6 (DispRelation *, int, double, double *, double *);

	EvalDisp6 (disp, sporder, wl1, &fl, &df);
	fl -= rpix;
	EvalDisp6 (disp, sporder, wl2, &fh, &df);
	fh -= rpix;
	if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0))
	    return (0.0);
	if (fl == 0.0) return wl1;
	if (fh == 0.0) return wl2;
	if (fl < 0.0) {
	    wll = wl1;
	    wlh = wl2;
	} else {
	    wlh = wl1;
	    wll = wl2;
	}
	rts = 0.5 * (wl1 + wl2);
	dwlold = fabs (wl2 - wl1);
	dwl    = dwlold;
	EvalDisp6 (disp, sporder, rts, &f, &df);
	f -= rpix;
	for (j = 1; j <= MAXIT; j++) {
	    if ((((rts - wlh) * df - f) * ((rts - wll) * df - f) >= 0.0)
	       || (fabs (2.0 * f) > fabs (dwlold * df))) {
	        dwlold = dwl;
	        dwl    = 0.5 * (wlh - wll);
	        rts    = wll + dwl;
	        if (wll == rts) return rts;
	    } else {
	        dwlold = dwl;
	        dwl    = f / df;
	        temp   = rts;
	        rts   -= dwl;
	        if (temp == rts) return rts;
	    }
	    if (fabs (dwl) < wl_acc) return rts;
	    EvalDisp6 (disp, sporder, rts, &f, &df);
	    f -= rpix;
	    if (f < 0.0)
	        wll = rts;
	    else
	        wlh = rts;
	}
	return (0.0);
}


/* This function evaluates the dispersion relation for the prism.
   The input is pixel number (reference pixel units), and the function
   value is the corresponding wavelength in Angstroms.
*/

static double PrismDisp (DispRelation *disp_y, double ix_r) {

        double x;       /* pixel number - first coefficient */
        double wl;      /* wavelength */

        x = ix_r - disp_y->coeff[0] + 1.0;  /* 1-indexed */

        wl = disp_y->coeff[5] / x;
        wl = (disp_y->coeff[4] + wl) / x;
        wl = (disp_y->coeff[3] + wl) / x;
        wl = (disp_y->coeff[2] + wl) / x;
        wl += disp_y->coeff[1];

        return (wl);
}
