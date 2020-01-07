# include <stdio.h>
# include <string.h>
# include <math.h>		/* for sqrt */

# include "hstio.h"
# include "wf3.h"
# include "wf3info.h"
# include "hstcalerr.h"

/* This routine checks whether the error array is all zero, and if so,
   a simple noise model is evaluated and assigned to the error array.
   The noise model is:

	sigma = sqrt ((I - bias) * gain + readnoise^2)  electrons

	      = sqrt ((I - bias) / gain + (readnoise / gain)^2)  dn,

   where I is the science data value in dn, bias is in dn, readnoise is
   the readout noise in electrons, gain is the CCD gain in electrons per dn.
   The value of sigma in dn is what is assigned to the error array.
   
    This has been modified to account for 1/2/4 AMP readout of a CCD
    where each AMP creates regions of the chip with different
    readnoise and gain values.  
    The regions are specified in the CCDTAB table,
    when the READNSE[1-4] values are determined.
    
    For each AMP used, a gain and readnoise value will be non-zero, with
    the AMPX and AMPY parameters specifying the first column affected by
    the second AMP in that dimension.
    
    The boundary for each AMP is modified by the trim values as 
       read in from the LTV1,2 keywords in the header.
       
    For ACS MAMA data, a separate loop initializes the ERR array with a
        value of SQRT(SCI data) or 1, whichever is greater.
        WJH 27 July 1999

    H.Bushouse, 2001 May 8:
	Revisions to keep in sync with calacs: multi-amp bug fixes.

    H.Bushouse, 2001 Nov 16:
	Revisions to keep in sync with calacs: modified to use individual
	bias values for each amp.

    H.Bushouse, 2003 Oct 27:
	Modified to use correct default bias values for CCD chips
	(following CALACS changes).
    H.Bushouse, 2003 Nov 11:
	Modified use of ampx to take into account presence of WFC3
	serial virtual overscan columns.
*/

int doNoise (WF3Info *wf3, SingleGroup *x, int *done) {

/* arguments:
WF3Info *wf3	 i: calibration switches and info
SingleGroup *x	io: image to be calibrated; written to in-place
int *done        o: true if we actually did assign error array values
*/

	extern int status;

	float bias;		/* bias level to subtract (dn) */
	float rn2[NAMPS];       /* square of noise values for observation */
	float gain[NAMPS];	/* gain values for observation */
	float value;		/* dn - bias * gain (i.e. signal in el) */
	int i, j;
	int ampx;		/* border column for 2amp readout regions */
	int ampy;		/* Boundary values corrected for trim regions */
	int dimx, dimy;
	int offsetx, offsety;
	float ccdbias[NAMPS];	/* default ccdbias values for chip */

	void get_nsegn (int, int, int, int, float *, float*, float *, float *);

	*done = 0;				/* initial value */
	
	/* First check for a dummy error array.  If it's not dummy,
	** we just return without doing anything.  */
	dimx = x->err.data.nx;
	dimy = x->err.data.ny;
	for (j = 0;  j < dimy;  j++) {
	    for (i = 0;  i < dimx;  i++) {
		    if (Pix (x->err.data, i, j) != 0.) {
		        return (status);	/* not a dummy error array */
		    }
	    }
	}

    if (wf3->detector != IR_DETECTOR) { 

        /* CCD initialization */
        offsetx = (int)(wf3->offsetx > 0) ? wf3->offsetx : 0;
        offsety = (int)(wf3->offsety > 0) ? wf3->offsety : 0;

	/* Correct the AMP readout boundaries for this offset */
	ampx = ((wf3->ampx == 0) ? 0 : (int)(wf3->ampx + offsetx) );
	ampy = ((wf3->ampy == 0) ? 0 : (int)(wf3->ampy + offsety) );
	ampx += wf3->trimx[2];

	/* Bounds checking to make sure we don't try to apply gain
	** and noise outside the bounds of the image.
	** This would apply if using only 1 AMP on each WFC chip when
	** ampx is given in CCDTAB as something greater than the
	** image size.
	** WJH 8 Sept 2000 (HAB 8 May 2001)
	**
	** We need to make sure that if the ampx value extends into the
	** overscan at the end of the line, ampx gets automatically 
	** moved to cover the whole line. This allows all AMPX and AMPY
	** values to be specified in CCDTAB in trimmed coordinates.
	** WJH 27 Oct 2000 (HAB 8 May 2001)
	*/

	/*if (ampx >= (dimx - wf3->offsetx) || ampx > dimx ) ampx = dimx;*/
	if (ampx >= (dimx - wf3->trimx[1]) || ampx > dimx ) ampx = dimx;
	if (ampy >= (dimy - wf3->offsety) || ampy > dimy ) ampy = dimy;

        /* Set up gain and readnoise arrays for use with chip's data */
        for (i=0; i < NAMPS; i++) {
            gain[i] = 0.;
            rn2[i] = 0.;
        }
        get_nsegn (wf3->detector, wf3->chip, wf3->ampx, wf3->ampy,
		   wf3->atodgain, wf3->readnoise, gain, rn2);

	/* For WFC3 UVIS data, AMPY will always be zero, yet we still need
	** to select the correct ccdbias value to apply. Determine
	** that here. */
	if (wf3->detector == CCD_DETECTOR && wf3->chip == 2) {
	    ccdbias[0] = wf3->ccdbias[2];
	    ccdbias[1] = wf3->ccdbias[3];
	    ccdbias[2] = wf3->ccdbias[2];
	    ccdbias[3] = wf3->ccdbias[3];
	} else {
	    for (i=0; i < NAMPS; i++) ccdbias[i] = wf3->ccdbias[i];
	}

        /* Now square the readnoise */
        for (i = 0; i < NAMPS; i++)
	     rn2[i] = rn2[i] * rn2[i];
        
	if (wf3->ncombine > 1) {
	    trlwarn("NCOMBINE > 1 before the error array was initialized.");
	}

        /* Now apply the initilalization for each AMP used */
        for (j = 0; j < ampy; j++) {
        
            /* This region corresponds to AMP_C, 
                if it is even used for this observation. */
	    /* Let's make sure we actually found a value for the gain
	    **	and readnoise... */

	    if (ampx > 0 && (gain[AMP_C] == 0. || rn2[AMP_C] == 0.)) { 
                trlerror 
		  ("No valid GAIN or READNOISE values to initialize ERR data."); 
		return (status = ERROR_RETURN);
            }
	    bias = ccdbias[2];
	    for (i = 0;  i < ampx;  i++) {
	        /* subtract bias and convert to electrons */
	        value = (Pix (x->sci.data, i, j) - bias) * gain[AMP_C];
	        if (value < 0.)
		    value = 0.;			/* sigma = 0 if signal = 0 */
		/* include readout noise and convert back to dn */
		Pix (x->err.data,i,j) = sqrt(value + rn2[AMP_C]) / gain[AMP_C];

            }

            /* This region corresponds to AMP_D, 
            **  if it is even used for this observation.  */

	    /* Let's make sure we actually found a value for the gain
	    **	and readnoise...  */
	    if (ampx == 0 && (gain[AMP_D] == 0. || rn2[AMP_D] == 0.) ) {
                trlerror 
		  ("No valid GAIN or READNOISE values to initialize ERR data.");
		return (status = ERROR_RETURN);
            }
	    bias = ccdbias[3];
	    for (i = ampx;  i < dimx;  i++) {
		 /* subtract bias and convert to electrons */
		 value = (Pix (x->sci.data, i, j) - bias) * gain[AMP_D];
		 if (value < 0.)
		     value = 0.;		/* sigma = 0 if signal = 0 */
		 /* include readout noise and convert back to dn */
		 Pix (x->err.data,i,j) = sqrt(value + rn2[AMP_D]) / gain[AMP_D];
	    }
        }
        
        
        for (j = ampy; j < dimy; j++) {

            /* This region corresponds to AMP_A, 
            **  if it is even used for this observation.  */

	    /* Let's make sure we actually found a value for the gain
	    **	and readnoise...  */
	    if (ampx > 0 && (gain[AMP_A] == 0. || rn2[AMP_A] == 0.)) {
                trlerror 
		  ("No valid GAIN or READNOISE values to initialize ERR data.");
		return (status = ERROR_RETURN);
            }
	    bias = ccdbias[0];
            /* Only execute this loop for AMP > 0 (multi-amp for line) */
	    for (i = 0;  i < ampx;  i++) {
		 /* subtract bias and convert to electrons */
		 value = (Pix (x->sci.data, i, j) - bias) * gain[AMP_A];
		 if (value < 0.)
		     value = 0.;		/* sigma = 0 if signal = 0 */
		 /* include readout noise and convert back to dn */
		 Pix (x->err.data,i,j) = sqrt(value + rn2[AMP_A]) / gain[AMP_A];
            }

            /* This region corresponds to AMP_B, 
            ** if it is even used for this observation.
            **  
            **  Default 1-AMP loop. AMPX and AMPY are zero.  */

            /* Let's make sure we actually found a value for the gain
            **	and readnoise...  */
            if (ampx == 0 && (gain[AMP_B] == 0. || rn2[AMP_B] == 0.)) {
                trlerror 
		  ("No valid GAIN or READNOISE values to initialize ERR data.");
	        return (status = ERROR_RETURN);
            }
	    bias = ccdbias[1];
	    for (i = ampx;  i < dimx;  i++) {
 
		 /* subtract bias and convert to electrons */
		 value = (Pix (x->sci.data, i, j) - bias) * gain[AMP_B];
		 if (value < 0.)
		     value = 0.;		/* sigma = 0 if signal = 0 */
		 /* include readout noise and convert back to dn */
		 Pix (x->err.data,i,j) = sqrt(value + rn2[AMP_B]) / gain[AMP_B];
	    }
        }

       /* End of CCD initialization */

    } else {

        /* ACS MAMA initialization */
        /* Set error to max of 1 or sqrt(counts) */
        for (j = 0; j < dimy; j++) {
            for (i = 0; i < dimx; i++) {
                value = sqrt(Pix(x->sci.data, i, j));
                Pix (x->err.data, i, j) = (value > 1) ? value : 1;
            }
        }
    }

    *done = 1;
    return (status);

}
