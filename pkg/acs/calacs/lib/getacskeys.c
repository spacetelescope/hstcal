# include <stdio.h>
# include <stddef.h>
# include <string.h>
# include <ctype.h>		/* islower, toupper */

# include "hstio.h"

# include "acs.h"
# include "acsinfo.h"
# include "hstcalerr.h"

/* This routine gets keyword values from the primary header.
 
 Warren Hack, 1998 June 8:
 Original ACS version based on Phil Hodge's CALSTIS routine...
 
 24-Sep-1998 WJH - Removed BIAS_REJ as a keyword
 17-Nov-1998 WJH - Revised to support trailer files
 11-Feb-1999 WJH - Read EXPTIME from Primary header instead of SCI hdr
 29-Oct-2001 WJH - Updated to use default value for CCDOFST[A,B,C.D] of 3
 4-Dec-2001 WJH - Read in EXPSTART and EXPEND for computing darktime.
 12-Dec-2012 PLL - Changed FLASHDUR and FLASHSTA defaults.
 */

int GetACSKeys (ACSInfo *acs, Hdr *phdr) {
  
  /* arguments:
   ACSInfo *acs  	io: calibration switches and info
   Hdr *phdr        i: primary header
   */
  
	extern int status;
  
	int nextend;			/* number of FITS extensions */
	int i;
	Bool subarray;
  
	int GetKeyInt (Hdr *, char *, int, int, int *);
	int GetKeyStr (Hdr *, char *, int, char *, char *, int);
	int GetKeyFlt (Hdr *, char *, int, float, float *);
  int GetKeyDbl (Hdr *, char *, int, double, double *);
  int GetKeyBool (Hdr *, char *, int, Bool, Bool *);
  
	/* Get generic parameters. */
  
	if (GetKeyStr (phdr, "ROOTNAME", NO_DEFAULT, "",
                 acs->rootname, ACS_CBUF))
    return (status);
  
	if (GetKeyStr (phdr, "APERTURE", USE_DEFAULT, "", acs->aperture, ACS_CBUF))
    return (status);
	if (GetKeyStr (phdr, "OBSTYPE", USE_DEFAULT, "", acs->obstype, ACS_CBUF))
    return (status);
  if (GetKeyStr (phdr, "JWROTYPE", USE_DEFAULT, "", acs->jwrotype, ACS_CBUF))
    return (status);
  
	if (GetKeyStr (phdr, "DETECTOR", NO_DEFAULT, "", acs->det, ACS_CBUF))
    return (status);
	if (strcmp (acs->det, "SBC") == 0) {
    acs->detector = MAMA_DETECTOR;
	} else if (strcmp (acs->det, "HRC") == 0) {
    acs->detector = HRC_CCD_DETECTOR;
	} else if (strcmp (acs->det, "WFC") == 0) {
    acs->detector = WFC_CCD_DETECTOR;
	} else {
    sprintf (MsgText, "DETECTOR = %s is invalid", acs->det);
    trlerror (MsgText);
		return (status = HEADER_PROBLEM);
	}
  
  /* Grating or mirror name. */
  if (GetKeyStr (phdr, "FILTER1", USE_DEFAULT, "", acs->filter1, ACS_CBUF))
    return (status);
  if (GetKeyStr (phdr, "FILTER2", USE_DEFAULT, "", acs->filter2, ACS_CBUF))
    return (status);
  
	if (GetKeyDbl (phdr, "EXPTIME", NO_DEFAULT, 0., &acs->exptime))
    return (status);
	if (acs->exptime < 0.) {
    sprintf (MsgText,"Exposure time is invalid:  %14.6g.", acs->exptime);
    trlerror (MsgText);
		return (status = INVALID_EXPTIME);
	}
	if (GetKeyDbl (phdr, "EXPSTART", NO_DEFAULT, 0., &acs->expstart))
    return (status);
	if (GetKeyDbl (phdr, "EXPEND", NO_DEFAULT, 0., &acs->expend))
    return (status);
  
	/* Find out how many extensions there are in this file. */
	if (GetKeyInt (phdr, "NEXTEND", USE_DEFAULT, EXT_PER_GROUP, &nextend))
    return (status);
  
	/* Convert number of extensions to number of SingleGroups. */
	acs->nimsets = nextend / EXT_PER_GROUP;
	if (acs->nimsets < 1) {
    sprintf (MsgText, "NEXTEND = %d; must be at least %d.", nextend, EXT_PER_GROUP);
    trlerror (MsgText);
		return (status = INVALID_VALUE);
	}
  
	if (GetKeyBool (phdr, "SUBARRAY", NO_DEFAULT, 0, &subarray))
	  return (status);
	if (subarray)
	  acs->subarray = YES;
	else
	  acs->subarray = NO;
  
	/* Get CCD-specific parameters. */
  
	if (acs->detector != MAMA_DETECTOR) {
    
    if (GetKeyStr (phdr, "CCDAMP", NO_DEFAULT, "",
                   acs->ccdamp, NAMPS))
      return (status);
    
    for (i=0; i < strlen(acs->ccdamp) ; i++) {
			/* Convert each letter in CCDAMP to upper-case. */
      if (islower (acs->ccdamp[i]))
				acs->ccdamp[i] = toupper (acs->ccdamp[i]);
      
			/* Verify that only the letters 'ABCD' are in the string. */
      if (strchr ("ABCD", acs->ccdamp[i]) == NULL) {
				sprintf (MsgText, "CCDAMP = `%s' is invalid.", acs->ccdamp);
				trlerror (MsgText);
				return (status = INVALID_VALUE);
      }
		}
    
    if (GetKeyFlt (phdr, "CCDGAIN", USE_DEFAULT, 1, &acs->ccdgain))
      return (status);
    /*
     ASSUMPTION: if no CCDOFST values are found, assume they were
     taken with the default offset setting of 3.  This will affect
     which row is selected from CCDTAB for setting the default value
     of the bias level in case there is no overscan regions for image.
     */
    if (GetKeyInt (phdr, "CCDOFSTA", USE_DEFAULT, DEFAULT_OFFSET, &acs->ccdoffset[0]))
	    return (status);
    if (GetKeyInt (phdr, "CCDOFSTB", USE_DEFAULT, DEFAULT_OFFSET, &acs->ccdoffset[1]))
	    return (status);
    if (GetKeyInt (phdr, "CCDOFSTC", USE_DEFAULT, DEFAULT_OFFSET, &acs->ccdoffset[2]))
	    return (status);
    if (GetKeyInt (phdr, "CCDOFSTD", USE_DEFAULT, DEFAULT_OFFSET, &acs->ccdoffset[3]))
	    return (status);
    
    if (GetKeyInt (phdr, "BINAXIS1", USE_DEFAULT, 1, &acs->binaxis[0]))
      return (status);
    if (GetKeyInt (phdr, "BINAXIS2", USE_DEFAULT, 1, &acs->binaxis[1]))
      return (status);
    
  if (GetKeyFlt (phdr, "FLASHDUR", USE_DEFAULT, 0.0, &acs->flashdur))
    return (status);
  if (GetKeyStr (phdr, "FLASHSTA", USE_DEFAULT, "", acs->flashstatus, ACS_CBUF))
    return (status);
    
	}
  
	return (status);
}
