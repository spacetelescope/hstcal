# include <stdio.h>
# include <string.h>
# include <math.h>

# include "xtables.h"

# include "stis.h"
# include "stisdq.h"
# include "calstis6.h"
# include "stissizes.h"

# define GRATING_DISP  1
# define PRISM_DISP    2


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
   30 Jul 10  -  Delete function PrismDisp, call prismDisp instead (PEH).
   05 Nov 10  -  For gratings, call evalInvDisp; delete CenWave, Newton (PEH).
*/

void Wave (StisInfo6 *sts, DispRelation *disp, RowContents *rowc) {

/* argument:
StisInfo6 *sts         i:  calibration switches and info
DispRelation *disp     i:  dispersion relation
RowContents rowc       io: row contents
*/

	double rpix;		/* reference pixel */
	/*double cenwave;*/     /* central wavelength */
	double wl;		/* wavelength */
	double wl_estimate;	/* estimate of wavelength */
	double m=1.;		/* spectral order */
	double tolerance=0.001;	/* wavelength error that is acceptable */
        int disp_type;          /* disperser is grating or prism */
	int    i, stopped;
	int status=0;		/* return code from evalInvDisp */

	/* Find disperser type. */
	if (strcmp (sts->opt_elem, "PRISM") == 0) {
	    disp_type = PRISM_DISP;
	} else {
	    disp_type = GRATING_DISP;
	    m = (double)rowc->sporder;
	}

	/* Loop over wavelength array. */
        stopped = 0;
	wl_estimate = (double)(sts->cenwave);
	for (i = 0;  i < rowc->npts;  i++) {

	    /* Translate image pixel index into reference pixel index. */
	    rpix = (i - sts->ltv[0]) / sts->ltm[0];

	    if (disp_type == GRATING_DISP) {
	        /* Solve for wavelength. */
		status = evalInvDisp (disp->coeff, disp->ncoeff,
			m, rpix, wl_estimate, tolerance, &wl);
	        rowc->wave[i] = wl;
		if (status == 0) {
		    wl_estimate = wl;
		} else {
		    printf ("Warning  error code %d from evalInvDisp\n",
			status);
		}

	    } else if (disp_type == PRISM_DISP) {
                if (stopped) {
	            rowc->wave[i] = MAX_PRISM_WAVELENGTH;
	            rowc->dq[i]  |= CALIBDEFECT ;
                    continue;
	        }
	        rowc->wave[i] = prismDisp (disp->coeff, rpix);
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
