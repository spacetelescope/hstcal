/* Contains:
 get_nsegn
 */

# include   <stdio.h>
# include   <string.h>

# include   "hstio.h"
# include   "acs.h"
# include   "acsinfo.h"

/*  This contains the function used to determine the NOISE and GAIN value
 for a given pixel position given a CCDAMP and CHIP configuration.
 
 Description:
 ------------
 These functions determine what gain and noise values are appropriate
 for a given pixel position.  The CCDAMP keyword as well as the OSCNTAB
 table columns, colx and coly, are passed along through the 'multiamp'
 structure for this determination.  In addition, the detector type
 is used to correctly assess which amps are appropriate or possible.
 
 Date            Author			Description
 ----            ------			-----------
 22-Sep-1998     W.J. Hack       initial version, uses multiamp noise,gain
 April 1999      W.J. Hack       simplified to use noise,gain arrays always
 26-Oct-2000 	W.J. Hack		revised the logic for diagonal amp usage
 10-Nov-2000 	W.J. Hack		completed support for all amp configs
 */

/* This function populates the gain and readnoise arrays appropriate
 for the detector used.  In the case of HRC, it provides a copy of both
 from the ACSInfo structure (header).  However, for WFC, it insures
 that the first two elements are the appropriate values for the chip,
 as only 2 amps can be used at most per chip.  The array elements for
 the un-used WFC AMPS for each chip are set to zero.
 
 The logic for what elements are appropriate follow from how AMPX and
 AMPY will be used throughout CALACS.  
 
 IF AMPY > 0., then more than one AMP is used in the Y direction
 If AMPX > 0., then more than one AMP is used in the X direction
 
 So, if line being processed > AMPY then,  
 For all pixels up to AMPX, use AMP_C values 
 If AMPX is ZERO, don't use AMP_C values for any pixel
 For remaining pixels in line (or all), Use AMP_D values 
 
 However, if AMPY is ZERO or the line number > AMPY (and AMPY > 0) then,
 For all pixels up to AMPX, Use AMP_A values
 If AMPX is ZERO, don't use AMP_A values for any pixel
 Then, for remaining pixels in line (or all of them), 
 Use AMP_B values
 
 Therefore, for 1 AMP use, AMPX and AMPY are both ZERO, putting the
 action in the AMP_B slot, regardless of which amp was actually used
 for reading out the chip. This ambiguity is the price to be paid
 for this logic, but the simpler, faster code makes it worthwhile.
 Also, for WFC, AMPY is always zero, resulting in AMP_A and AMP_B
 entries as the only ones needed/used.  All entries which are not
 appropriate will have values of zero.  
 
 For diagonal amp usage,
 this routine will establish symmetry by filling in the opposite
 amp values with the values from the amp on the same side of the chip
 (in the X direction).  With AMPX set to the center and AMPY to 0, 
 we apply the gain values appropriately as if it were read out by
 all 4 amps.  
 
 This routine also assumes that only those amps used will have valid
 values in the ACSInfo structure to begin with, as no further checking
 of amp usage will be required.  
 */  
void get_nsegn (int detector, int chip, int ampx, int ampy, float *hdrgain, float
                *hdrnoise, float *gain, float *rn2) {
  
  /* 
   Parameters:
   int     detector        i: detector ID
   int     chip            i: chip ID
   int     ampx,ampy       i: first columns used by second AMP in x,y axes
   float   *hdrgain          i: atodgn values from image/ccdtab
   float   *hdrnoise         i: readnoise values from image/ccdtab
   float   *gain           o: gain - array of valid values
   float   *rn2            o: readnoise - array of valid values 
   */
  
  int j;
  
  /* Set default values for MAMA */
  for (j = 0; j < NAMPS; j++) {
    rn2[j] = 0.0;
    gain[j] = 1.0;	
  }
  
  if (detector != MAMA_DETECTOR) { 
    /* Start of CCD section */
    if (detector == WFC_CCD_DETECTOR && chip == 2) {
      /* WFC chip 1 */
      gain[0] = hdrgain[AMP_C];
      gain[1] = hdrgain[AMP_D];
      gain[2] = hdrgain[AMP_C];
      gain[3] = hdrgain[AMP_D];
      rn2[0] = hdrnoise[AMP_C];
      rn2[1] = hdrnoise[AMP_D];
      rn2[2] = hdrnoise[AMP_C];
      rn2[3] = hdrnoise[AMP_D];
    } else if (detector == WFC_CCD_DETECTOR && chip == 1) {
      /* WFC chip 2 */
      gain[0] = hdrgain[AMP_A];
      gain[1] = hdrgain[AMP_B];
      gain[2] = hdrgain[AMP_A];
      gain[3] = hdrgain[AMP_B];
      rn2[0] = hdrnoise[AMP_A];
      rn2[1] = hdrnoise[AMP_B];
      rn2[2] = hdrnoise[AMP_A];
      rn2[3] = hdrnoise[AMP_B];
      
    } else if (ampx == 0 && ampy == 0) {
      /* Find the 1 HRC AMP that is used and set it to AMP_B.
       This corresponds to the case for ampx && ampy == 0. 
       */
      for (j = 0; j < NAMPS; j++) {
        if (hdrgain[j] > 0.) {
          rn2[1] = hdrnoise[j];
          gain[1] = hdrgain[j];
          return;	
        }
      }    
    } else {
      /* HRC detector multi-amp case */
			/* For diagonal amp usage (BC or AD), fill 
       in blank amp values with values from other 
       half of chip.  
       For BC, actually build CBCB.
       For AD, actually build ADAD.
       */
			if (hdrgain[0] == 0. && hdrgain[3] == 0.) {
				/* We have CCDAMP == BC... */
				gain[0] = hdrgain[2];
				gain[1] = hdrgain[1];
				gain[2] = hdrgain[2];
				gain[3] = hdrgain[1];
				rn2[0] = hdrnoise[2];
				rn2[1] = hdrnoise[1];
				rn2[2] = hdrnoise[2];
				rn2[3] = hdrnoise[1];
			} else if (hdrgain[1] == 0. && hdrgain[2] == 0.){
				/* We have CCDAMP == AD,... */
				gain[0] = hdrgain[0];
				gain[1] = hdrgain[3];
				gain[2] = hdrgain[0];
				gain[3] = hdrgain[3];
				rn2[0] = hdrnoise[0];
				rn2[1] = hdrnoise[3];
				rn2[2] = hdrnoise[0];
				rn2[3] = hdrnoise[3];
			} else { 
				/* We have an AMP situation which can simply 
         be used as is with AMPX/AMPY... */
        for (j = 0; j < NAMPS; j++) {
          rn2[j] = hdrnoise[j];
          gain[j] = hdrgain[j];	
        }
			}
    }
  } /* End of CCD section */
}
