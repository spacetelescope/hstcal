# include <stdio.h>
# include <string.h>
# include <math.h>		/* for sqrt */
# include <time.h>

# include "hstio.h"
# include "acs.h"
# include "acsinfo.h"
# include "err.h"

/* This routine checks whether the error array is all zero, and if so,
 a simple noise model is evaluated and assigned to the error array.
 The noise model is:

 sigma = sqrt (I * gain + readnoise^2)  electrons

 = sqrt (I / gain + (readnoise / gain)^2)  dn,

 where I is the science data value in dn, readnoise is
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

 For MAMA data, a separate loop initializes the ERR array with a
 value of SQRT(SCI data) or 1, whichever is greater.
 WJH 27 July 1999
 */

int doNoise (ACSInfo *acs, SingleGroup *x, int *done) {

  /* arguments:
   ACSInfo *acs    i: calibration switches and info
   SingleGroup *x   io: image to be calibrated; written to in-place
   int *done         o: true if we actually did assign error array values
   */

	extern int status;

	float rn[NAMPS],rn2[NAMPS]; /* square of noise values for observation */
  float gain[NAMPS];          /* gain values for observation */
	float value;                /* signal in e- */
  float err_val;              /* error value in e- */
	int i, j;
	int ampx;		/* border column for 2amp readout regions */
	int ampy;		/* Boundary values corrected for trim regions */
  int dimx, dimy;
  int offsetx, offsety;
  float val;
  
  char targname[ACS_LINE];

	int FindLine(SingleGroup *, SingleGroupLine *, int *, int *,int *, int *, int *);
	int DetCCDChip(char *, int, int, int *);
  int GetKeyStr(Hdr *, char *, int, char *, char *, int);
  int GetKeyDbl(Hdr *, char *, int, double, double *);
  void get_nsegn(int, int, int, int, float *, float*, float *, float *);

	*done = 0;				/* initial value */

  dimx = x->err.data.nx;
  dimy = x->err.data.ny;

  if (acs->detector != MAMA_DETECTOR) {
    /*
     CCD initialization
     */
    offsetx = (int)(acs->offsetx > 0) ? acs->offsetx : 0;
    offsety = (int)(acs->offsety > 0) ? acs->offsety : 0;

    if (GetKeyStr(x->globalhdr, "TARGNAME", USE_DEFAULT, "", targname, ACS_LINE))
      return (status);

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
      rn2[i] = 0.;
    }

    get_nsegn(acs->detector, acs->chip, acs->ampx, acs->ampy, acs->atodgain, acs->readnoise, gain, rn);

    /* Now square the readnoise */
    for (i = 0; i < NAMPS; i++) rn2[i] = rn[i] * rn[i];

    if (acs->ncombine > 1) {
      trlwarn ("NCOMBINE > 1 before the error array was initialized.");
    }

    /* Now apply the initilalization for each AMP used */
    for (j = 0; j < ampy; j++) {

      /*
       This region corresponds to AMP_C,
       if it is even used for this observation.
       */

      /* Let's make sure we actually found a value for the gain
       **	and readnoise...
       */
      if (ampx > 0 && rn[AMP_C] == 0.) {
        trlerror ("No valid READNOISE values to initialize ERR data.");
        return (status = ERROR_RETURN);
      }

      for (i = 0;  i < ampx;  i++) {
        value = Pix(x->sci.data, i, j);
        err_val = Pix(x->err.data, i, j);

        if (value < 0.)
          value = 0.;			/* sigma = 0 if signal = 0 */

        if (strncmp(targname,"BIAS",4) != 0) {
          /* include readout noise and convert back to dn */
          Pix(x->err.data, i, j) = sqrt(value + rn2[AMP_C] + err_val*err_val);
        } else {
          /* BIAS exposure being processed: only set to RN */
          Pix(x->err.data, i, j) = rn[AMP_C];
        }
      } /* End of loop over X for AMP_C */

      /*
       This region corresponds to AMP_D,
       if it is even used for this observation.
       */

      /* Let's make sure we actually found a value for the gain
       **	and readnoise...
       */
      if (ampx == 0 && rn[AMP_D] == 0.) {
        trlerror ("No valid READNOISE values to initialize ERR data.");
        return (status = ERROR_RETURN);
      }

      for (i = ampx;  i < dimx;  i++) {
        value = Pix(x->sci.data, i, j);
        err_val = Pix(x->err.data, i, j);
        
        if (value < 0.)
          value = 0.;			/* sigma = 0 if signal = 0 */
        
        if (strncmp(targname,"BIAS",4) != 0) {
          /* include readout noise and convert back to dn */
          Pix(x->err.data, i, j) = sqrt(value + rn2[AMP_D] + err_val*err_val);
        } else {
          /* BIAS exposure being processed: only set to RN */
          Pix(x->err.data, i, j) = rn[AMP_D];
        }
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
      if (ampx > 0 && rn[AMP_A] == 0.) {
        trlerror ("No valid READNOISE values to initialize ERR data.");
        return (status = ERROR_RETURN);
      }

      /* Only execute this loop for AMP > 0 (multi-amp for line) */
      for (i = 0;  i < ampx;  i++) {
        value = Pix(x->sci.data, i, j);
        err_val = Pix(x->err.data, i, j);
        
        if (value < 0.)
          value = 0.;			/* sigma = 0 if signal = 0 */
        
        if (strncmp(targname,"BIAS",4) != 0) {
          /* include readout noise and convert back to dn */
          Pix(x->err.data, i, j) = sqrt(value + rn2[AMP_A] + err_val*err_val);
        } else {
          /* BIAS exposure being processed: only set to RN */
          Pix(x->err.data, i, j) = rn[AMP_A];
        }
      }

      /*
       This region corresponds to AMP_B,
       if it is even used for this observation.

       Default 1-AMP loop. AMPX and AMPY are zero.
       */
      /* Let's make sure we actually found a value for the gain
       **	and readnoise...
       */
      if (ampx == 0 && rn[AMP_B] == 0.){
        trlerror ("No valid READNOISE values to initialize ERR data.");
        return (status = ERROR_RETURN);
      }

      for (i = ampx;  i < dimx;  i++) {
        value = Pix(x->sci.data, i, j);
        err_val = Pix(x->err.data, i, j);
        
        if (value < 0.)
          value = 0.;			/* sigma = 0 if signal = 0 */

        if (strncmp(targname,"BIAS",4) != 0) {
          /* include readout noise and convert back to dn */
          Pix(x->err.data, i, j) = sqrt(value + rn2[AMP_B] + err_val*err_val);
        } else {
          /* BIAS exposure being processed: only set to RN */
          Pix(x->err.data, i, j) = rn[AMP_B];
        }
      }
    }
    /* End of CCD initialization */
    
  } else {
    /* MAMA initialization */
    /* Set error to max of 1 or sqrt(counts) */
    for (j = 0; j < dimy; j++) {
      for (i = 0; i < dimx; i++) {
        val = Pix(x->sci.data, i, j);
        if (val < 0.0){
          sprintf(MsgText,"Negative value found at (%d,%d) in input MAMA data!",(i+1),(j+1));
          trlwarn (MsgText);
          val = abs(val);
        }
        value = sqrt(val);
        Pix (x->err.data, i, j) = (value > 1) ? value : 1;
      }
    }

  }
	*done = 1;
	return (status);
}


/* check whether the SingleGroup error arrays is all zeros. returns YES
 * if that's the case, NO otherwise. */
int check_zero_noise(SingleGroup *x) {
  int all_zeros = YES;
  
  int i, j, dimx, dimy;
  
  dimx = x->err.data.nx;
  dimy = x->err.data.ny;
  
  for (i = 0; i < dimx; i++) {
    for (j = 0; j < dimy; j++) {
      if (Pix(x->err.data, i, j) != 0.) {
        all_zeros = NO;
        return all_zeros;
        }
      }
    }
  
  return all_zeros;
  }
