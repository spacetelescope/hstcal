# include <stdio.h>
# include <stddef.h>
# include <string.h>
# include <ctype.h>		/* islower, toupper */

# include "hstio.h"

# include "wf3.h"
# include "wf3info.h"
# include "hstcalerr.h"

/* This routine gets keyword values from the primary header.

   Warren Hack, 1998 June 8:
    Original ACS version based on Phil Hodge's CALSTIS routine...

    24-Sep-1998 WJH - Removed BIAS_REJ as a keyword
    17-Nov-1998 WJH - Revised to support trailer files
    11-Feb-1999 WJH - Read EXPTIME from Primary header instead of SCI hdr
    17-Apr-2001 HAB - WF3 IR channel mods: changed exptime from scalar to array;
		      read FILTER keyword instead of FILTER1.
     8-May-2001 HAB - Added support for CCD post-flash keywords;
		      Change UVIS FILTER1 and FILTER2 to just FILTER
		      (WFC3 UVIS will allow only 1 filter at a time).
    16-Nov-2001 HAB - Updated to use default value of 3 for CCDOFST[A,B,C,D]
    21-Jun-2002 HAB - Updated to read in EXPSTART and EXPEND for computing
		      darktime, and also read in OBSTYPE.
    16-Oct-2003 HAB - Updated to use floating-point gain values for WFC3.
    20-Feb-2004 HAB - Eliminated attemp to read BINAXIS keywords from primary
		      header because for WFC3 they're in the sci extension hdr.
    15-Feb-2007 HAB - Updated default gain for IR channel from 2.0 to 2.5.
                      Added 'subtype' to list of IR keywords loaded. Changed
                      default sampzero value to 2.911755 sec, to correspond
		      with new IR timing patterns.
    09-Jan-2009 HAB - Eliminated use of default values for FILTER and CCDGAIN
		      keywords. It will now be an error if they aren't present.
*/

int GetKeys (WF3Info *wf3, Hdr *phdr) {

/* arguments:
WF3Info *wf3  	io: calibration switches and info
Hdr *phdr        i: primary header
*/

	extern int status;

	int nextend;			/* number of FITS extensions */
	int i;
	Bool subarray=0;

	int GetKeyInt (Hdr *, char *, int, int, int *);
	int GetKeyFlt (Hdr *, char *, int, float, float *);
	int GetKeyDbl (Hdr *, char *, int, double, double *);
	int GetKeyStr (Hdr *, char *, int, char *, char *, int);
	int GetKeyBool (Hdr *, char *, int, Bool, Bool *);

	/* Get generic parameters. */

	if (GetKeyStr (phdr, "ROOTNAME", NO_DEFAULT, "", wf3->rootname,SZ_CBUF))
	    return (status);

	if (GetKeyStr (phdr, "APERTURE", USE_DEFAULT, "",wf3->aperture,SZ_CBUF))
	    return (status);

	if (GetKeyStr (phdr, "OBSTYPE", USE_DEFAULT, "",wf3->obstype,SZ_CBUF))
	    return (status);

	if (GetKeyStr (phdr, "DETECTOR", NO_DEFAULT, "", wf3->det, SZ_CBUF))
	    return (status);

	if (GetKeyBool (phdr, "SUBARRAY", NO_DEFAULT, 0, &subarray))
		  return (status);

	if (subarray){
			wf3->subarray=1;
	}

	if (strcmp (wf3->det, "IR") == 0) {
	    wf3->detector = IR_DETECTOR;
	} else if (strcmp (wf3->det, "UVIS") == 0) {
	    wf3->detector = CCD_DETECTOR;
	} else {
	    sprintf (MsgText, "DETECTOR = %s is invalid", wf3->det);
	    trlerror (MsgText);
	    return (status = HEADER_PROBLEM);
	}

	/* Filter or prism/grism name */
	if (GetKeyStr (phdr, "FILTER", NO_DEFAULT, "", wf3->filter, SZ_CBUF))
	    return (status);

	/* Exposure time */
	if (GetKeyDbl (phdr, "EXPTIME", NO_DEFAULT, 0., &(wf3->exptime[0])))
	    return (status);
	if (wf3->exptime[0] < 0.) {
	    sprintf(MsgText,"Exposure time is invalid:  %14.6g.",
		    wf3->exptime[0]);
	    trlerror (MsgText);
	    return (status = INVALID_EXPTIME);
	}
	if (GetKeyDbl (phdr, "EXPSTART", NO_DEFAULT, 0., &wf3->expstart))
	    return (status);
	if (GetKeyDbl (phdr, "EXPEND", NO_DEFAULT, 0., &wf3->expend))
	    return (status);

	/* Find out how many extensions there are in this file. */
	if (GetKeyInt (phdr, "NEXTEND", USE_DEFAULT, EXT_PER_GROUP, &nextend))
	    return (status);

	/* Convert number of extensions to number of SingleGroups. */
	wf3->nimsets = nextend / EXT_PER_GROUP;
	if (wf3->nimsets < 1) {
	    sprintf (MsgText, "NEXTEND = %d; must be at least %d.", nextend,
		     EXT_PER_GROUP);
	    trlerror (MsgText);
	    return (status = INVALID_VALUE);
	}

	/* Get CCD-specific parameters. */

	if (wf3->detector == CCD_DETECTOR) {

	    if (GetKeyStr (phdr, "CCDAMP", NO_DEFAULT, "", wf3->ccdamp, NAMPS))
		return (status);

	    for (i=0; i < strlen(wf3->ccdamp) ; i++) {
			 /* Convert each letter in CCDAMP to upper-case. */
			 if (islower (wf3->ccdamp[i]))
			     wf3->ccdamp[i] = toupper (wf3->ccdamp[i]);

			 /* Verify that only the letters 'ABCD' are in the string. */
			 if (strchr ("ABCD", wf3->ccdamp[i]) == NULL) {
			     sprintf (MsgText, "CCDAMP = `%s' is invalid.",wf3->ccdamp);
			     trlerror (MsgText);
			     return (status = INVALID_VALUE);
			 }
	    }

	    if (GetKeyFlt (phdr, "CCDGAIN", NO_DEFAULT, 0., &wf3->ccdgain))
		return (status);

	    /* ASSUMPTION: if no CCDOFST values are found, assume they were
	    ** taken with the default offset setting of 3. This will affect
	    ** which row is selected from CCDTAB for setting the default value
	    ** of the bias level in case there is no overscan regions for image.
	    */

	    if (GetKeyInt (phdr, "CCDOFSTA", USE_DEFAULT, DEFAULT_OFFSET,
			   &wf3->ccdoffset[0]))
		return (status);
	    if (GetKeyInt (phdr, "CCDOFSTB", USE_DEFAULT, DEFAULT_OFFSET,
			   &wf3->ccdoffset[1]))
		return (status);
	    if (GetKeyInt (phdr, "CCDOFSTC", USE_DEFAULT, DEFAULT_OFFSET,
			   &wf3->ccdoffset[2]))
		return (status);
	    if (GetKeyInt (phdr, "CCDOFSTD", USE_DEFAULT, DEFAULT_OFFSET,
			   &wf3->ccdoffset[3]))
		return (status);
	    if (GetKeyFlt (phdr, "FLASHDUR", USE_DEFAULT, 1.0, &wf3->flashdur))
		return (status);
	    if (GetKeyStr (phdr, "FLASHSTA", NO_DEFAULT, "", wf3->flashstatus,
			   SZ_CBUF))
		return (status);

	} else {

	/* Get IR-specific parameters */

	    sprintf (wf3->ccdamp, "%s", "ABCD");

	    if (GetKeyFlt (phdr, "CCDGAIN", NO_DEFAULT, 0., &wf3->ccdgain))
		return (status);

	    wf3->nsamp = 0;
	    if (GetKeyInt (phdr, "NSAMP", NO_DEFAULT, 0, &wf3->nsamp))
		return (status);
	    wf3->ngroups = wf3->nsamp;

	    wf3->sampseq[0] = '\0';
	    if (GetKeyStr (phdr, "SAMP_SEQ", NO_DEFAULT, "", wf3->sampseq,
			   SZ_CBUF))
		return (status);

	    wf3->subtype[0] = '\0';
	    if (GetKeyStr (phdr, "SUBTYPE", NO_DEFAULT, "", wf3->subtype,
			   SZ_CBUF))
		return (status);

	    wf3->sampzero = 0.0;
	    if (GetKeyDbl (phdr, "SAMPZERO", USE_DEFAULT, 2.911755, &wf3->sampzero))
		return (status);
	}

	return (status);
}
