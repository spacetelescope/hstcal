# include <stdio.h>
# include <string.h>
# include <math.h>		/* for sqrt */
# include <time.h>

# include "hstio.h"
# include "acs.h"
# include "acsinfo.h"
# include "hstcalerr.h"

/* This functions converts a a SingleGroup from DN to electrons. */

int to_electrons(ACSInfo *acs, SingleGroup *x) {

  /* arguments:
   ACSInfo *acs    i: calibration switches and info
   SingleGroup *x   io: image to be calibrated; written to in-place
   int *done         o: true if we actually did assign error array values
   */

	extern int status;

	float rn[NAMPS];
  float gain[NAMPS];		    /* gain values for observation */
	int i, j;
	int ampx;		/* border column for 2amp readout regions */
	int ampy;		/* Boundary values corrected for trim regions */
  int dimx, dimy;
  int offsetx, offsety;

  void get_nsegn (int, int, int, int, float *, float*, float *, float *);


  offsetx = (int)(acs->offsetx > 0) ? acs->offsetx : 0;
  offsety = (int)(acs->offsety > 0) ? acs->offsety : 0;
  
  dimx = x->sci.data.nx;
  dimy = x->sci.data.ny;

  /* Correct the AMP readout boundaries for this offset */
  ampx = ((acs->ampx == 0) ? 0 : (int)(acs->ampx + offsetx) );
  ampy = ((acs->ampy == 0) ? 0 : (int)(acs->ampy + offsety) );

  /* Bounds checking to make sure we don't try to apply gain
   and noise outside the bounds of the image.
   This would apply if using only 1 AMP on each WFC chip when
   ampx is given in CCDTAB as something greater image size.
   WJH 8 Sept 2000

   We need to make sure that if the ampx value extends into the
   overscan at the end of the line, ampx gets automatically
   moved to cover the whole line. This allows all AMPX and AMPY
   values to be specified in CCDTAB in trimmed coordinates.
   WJH 27 Oct 2000
   */
  if (ampx >= (dimx - acs->offsetx) || ampx > dimx ) ampx = dimx;
  if (ampy >= (dimy - acs->offsety) || ampy > dimy ) ampy = dimy;

  /* Set up gain and readnoise arrays for use with chip's data */
  for (i=0; i < NAMPS; i++) {
    gain[i] = 0.;
    rn[i] = 0.;
  }
  
  get_nsegn(acs->detector, acs->chip, acs->ampx, acs->ampy, acs->atodgain, 
            acs->readnoise, gain, rn);
  
  /* Now apply the initilalization for each AMP used */
  for (j = 0; j < ampy; j++) {

    /*
     This region corresponds to AMP_C,
     if it is even used for this observation.
     */

    /* Let's make sure we actually found a value for the gain
     */
    if (ampx > 0 && gain[AMP_C] == 0.){
      trlerror ("No valid GAIN value to convert data.");
      return (status = ERROR_RETURN);
    }

    for (i = 0;  i < ampx;  i++) {
      /* convert to electrons */
      Pix(x->sci.data, i, j) = Pix(x->sci.data, i, j) * gain[AMP_C];
      Pix(x->err.data, i, j) = Pix(x->err.data, i, j) * gain[AMP_C];
    }

    /*
     This region corresponds to AMP_D,
     if it is even used for this observation.
     */

    /* Let's make sure we actually found a value for the gain
     **	and readnoise...
     */
    if (ampx == 0 && gain[AMP_D] == 0.){
      trlerror ("No valid GAIN value to convert data.");
      return (status = ERROR_RETURN);
    }

    for (i = ampx;  i < dimx;  i++) {
      /* convert to electrons */
      Pix(x->sci.data, i, j) = Pix(x->sci.data, i, j) * gain[AMP_D];
      Pix(x->err.data, i, j) = Pix(x->err.data, i, j) * gain[AMP_D];
    }
  }

  for (j = ampy; j < dimy; j++) {
    /*
     This region corresponds to AMP_A,
     if it is even used for this observation.
     */
    /* Let's make sure we actually found a value for the gain
     **	and readnoise...
     */
    if (ampx > 0 && gain[AMP_A] == 0.) {
      trlerror ("No valid GAIN value to convert data.");
      return (status = ERROR_RETURN);
    }

    /* Only execute this loop for AMP > 0 (multi-amp for line) */
    for (i = 0;  i < ampx;  i++) {
      /* convert to electrons */
      Pix(x->sci.data, i, j) = Pix(x->sci.data, i, j) * gain[AMP_A];
      Pix(x->err.data, i, j) = Pix(x->err.data, i, j) * gain[AMP_A];
    }

    /*
     This region corresponds to AMP_B,
     if it is even used for this observation.

     Default 1-AMP loop. AMPX and AMPY are zero.
     */
    /* Let's make sure we actually found a value for the gain
     **	and readnoise...
     */
    if (ampx == 0 && gain[AMP_B] == 0.){
      trlerror ("No valid GAIN value to convert data.");
      return (status = ERROR_RETURN);
    }

    for (i = ampx;  i < dimx;  i++) {
      /* convert to electrons */
      Pix(x->sci.data, i, j) = Pix(x->sci.data, i, j) * gain[AMP_B];
      Pix(x->err.data, i, j) = Pix(x->err.data, i, j) * gain[AMP_B];
    }
  }
  
	return (status);
}
