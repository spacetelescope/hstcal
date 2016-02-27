/* Contains:
    get_nsegn
*/

# include   <stdio.h>
# include   <string.h>

# include   "hstio.h"
# include   "wf3.h"
# include   "wf3info.h"

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
    8 May 2001      H.Bushouse      Revised the logic for diagonal amp usage;
				    completed support for all amp configs
				    (in sync with calacs changes).
    26 Aug 2008     H.Bushouse      Removed unnecessarily confusing ACS/HRC
                                    code that is not needed for WFC3. Also
				    added code necessary to handle WFC3 IR.
*/

/* This function populates the gain and readnoise arrays appropriate
    for the detector used.  In the case of WFC3/IR, it provides a copy of both
    from the WF3Info structure (header).  However, for WFC3/UVIS, it insures
    that the first two elements are the appropriate values for the chip,
    as only 2 amps can be used at most per chip.

    The logic for what elements are appropriate follow from how AMPX and
    AMPY will be used throughout CALWF3.  

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
    Also, for UVIS, AMPY is always zero, resulting in AMP_A and AMP_B
    entries as the only ones needed/used.  All entries which are not
    appropriate will have values of zero.
   
    For diagonal amp usage, this routine will establish symmetry by filling
    in the opposite amp values with the values from the amp on the same
    side of the chip (in the X direction). With AMPX set to the center and
    AMPY to 0, we apply the gain values appropriately as if it were read
    out by all 4 amps.

    This routine also assumes that only those amps used will have valid
    values in the WF3Info structure to begin with, as no further checking
    of amp usage will be required.  
*/  
void get_nsegn (int detector, int chip, int ampx, int ampy, float *hdrgain,
		float *hdrnoise, float *gain, float *rn2) {

/* 
Parameters:
    int     detector        i: detector ID
    int     chip            i: chip ID
    int     ampx,ampy       i: first columns used by second AMP in x,y axes
    float   *hdrgain        i: atodgn values from image/ccdtab
    float   *hdrnoise       i: readnoise values from image/ccdtab
    float   *gain           o: gain - array of valid values
    float   *rn2            o: readnoise - array of valid values 
*/
    
    if (detector == CCD_DETECTOR) { 

        /* Start of CCD section */
        if (chip == 2) {
            /* CCD chip 2: use amp C and D values */
            gain[0] = hdrgain[AMP_C];
            gain[1] = hdrgain[AMP_D];
            gain[2] = hdrgain[AMP_C];
            gain[3] = hdrgain[AMP_D];
            rn2[0] = hdrnoise[AMP_C];
            rn2[1] = hdrnoise[AMP_D];
            rn2[2] = hdrnoise[AMP_C];
            rn2[3] = hdrnoise[AMP_D];

        } else if (chip == 1) {
            /* CCD chip 1: use amp A and B values */
            gain[0] = hdrgain[AMP_A];
            gain[1] = hdrgain[AMP_B];
            gain[2] = hdrgain[AMP_A];
            gain[3] = hdrgain[AMP_B];
            rn2[0] = hdrnoise[AMP_A];
            rn2[1] = hdrnoise[AMP_B];
            rn2[2] = hdrnoise[AMP_A];
            rn2[3] = hdrnoise[AMP_B];
        }
        /* End of CCD section */

    } else {

        /* Start if IR section */
	/* Use all 4 amp values, but need to rearrange the ordering in
	   the arrays because of the different arrangement of the IR
	   amps relative to UVIS: 
		-----          -----
		|A|D|          |A|B|
	    IR: -----    UVIS: -----
		|B|C|          |C|D|
		-----          -----
	*/
        gain[0] = hdrgain[AMP_A];
        gain[1] = hdrgain[AMP_D];
        gain[2] = hdrgain[AMP_B];
        gain[3] = hdrgain[AMP_C];
        rn2[0] = hdrnoise[AMP_A];
        rn2[1] = hdrnoise[AMP_D];
        rn2[2] = hdrnoise[AMP_B];
        rn2[3] = hdrnoise[AMP_C];
    }   /* End of IR section */
}
