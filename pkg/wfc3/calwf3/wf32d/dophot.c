# include <stdio.h>
# include <stdlib.h>	/* malloc */
# include <string.h>    /* strchr */
# include <ctype.h>

# include "imphttab.h"
# include "hstio.h"

# include "wf3.h"
# include "wf3info.h"
# include "wf3err.h"

static void  Phot2Obs (char *, char *);

/* This routine calls synphot functions to compute values associated with
   the absolute flux calibration of the science data and populates those
   values into a set of photometry keywords in the science image header.

   Warren Hack, 1998 June 12:
   Initial ACS version.  

   Howard Bushouse, 2000 Aug 29:
   Initial WFC3 version.

   H.Bushouse, 2002 June 17:
   Removed gain normalization of PHOTFLAM keyword value, to be consistent
   with image datatype after flat-field calibration. Also, increased
   ARR_SIZE to match the default value used by SYNPHOT (in accordance
   with CALACS changes).

   H.Bushouse, 2008 Nov 5:
   Added computation of PHOTFNU keyword value, to be consistent with the
   WFC3 IR photcorr process. Also removed some old ACS-specific code that
   is not used for WFC3.

   H.Bushouse, 2011 Sep 7:
   Removed synphot dependencies and replaced with interface to new
   IMPHTTAB reference table, which returns precomputed photometric
   values.

   M. Sosey, 3, July 2013
   Added new FLUXCORR step to the pipeline to make the flux scaling for both chips
   equal. New keywords PHTFLAM1, PHTFLAM2, PHTRATIO added for tracking. PHTRATIO is
   the value chip2 is scaled by and is PHTFLAM1/PHTFLAM2
 */

int doPhot (WF3Info *wf32d, SingleGroup *x) {

	/* arguments:
	   WF3Info *wf32d     i: calibration switches, etc
	   SingleGroup *x    io: image to be calibrated; primary header is modified
	 */

	extern int status;

	PhotPar obs;
	float photfnu;
	float ratio;

	char  photmode[SZ_LINE+1], obsmode[SZ_LINE+1];

	int GetKeyStr (Hdr *, char *, int, char *, char *, int);
	int PutKeyFlt (Hdr *, char *, float, char *);

	if (wf32d->photcorr == DUMMY)
		return (status);

	/* Extract photmode from sci extension header */
	if (GetKeyStr (&x->sci.hdr, "PHOTMODE",USE_DEFAULT,"",photmode,SZ_LINE))
		return (status);

	/* Add commas to PHOTMODE string and change all strings to
	 ** lower case for use in synphot. */
	Phot2Obs (photmode, obsmode);
	if (wf32d->verbose) {
		sprintf (MsgText, "Created PYSYNPHOT obsmode of: %s", obsmode);
		trlmessage (MsgText);
	}

	/* Initialize PhotPar structure */
	InitPhotPar (&obs, wf32d->phot.name, wf32d->phot.pedigree);

	/* Get phot values from IMPHTTAB */
	if (GetPhotTab (&obs, obsmode)) {
		trlerror ("Error return from GetPhotTab.");
		return (status);
	}

	/* Add this information as a HISTORY comment */
	if (wf32d->verbose) {
		sprintf (MsgText, "Retrieved PHOTFLAM value of %g", obs.photflam);
		trlmessage (MsgText);
	}


	/* Update the photometry keyword values in the SCI header. */
	if (PutKeyFlt (&x->sci.hdr, "PHOTFLAM", obs.photflam,
				"inverse sensitivity"))
		return (status);

	/* the value for the phtflam? which is not in the current obsmode
	   should remain 0.0, which is the default. Both are here because
	   the obs structure is updated here, but each chip is processed
	   separately and the imphttab which contains the values we want
	   is read twice, once for each chips obsmode , and time the table
	   is read the obs structure is reset  */

	if (PutKeyFlt (&x->sci.hdr, "PHOTZPT", obs.photzpt, "zero point"))
		return (status);

	if (PutKeyFlt (&x->sci.hdr, "PHOTPLAM", obs.photplam, 
				"pivot wavelength"))
		return (status);

	if (PutKeyFlt (&x->sci.hdr, "PHOTBW", obs.photbw, "RMS bandwidth"))
		return (status);

	photfnu = 3.33564e+4 * obs.photflam * obs.photplam*obs.photplam;

	if (PutKeyFlt (&x->sci.hdr, "PHOTFNU", photfnu, "inverse sensitivity"))
		return (status);

	/*only update the chip that's currently processing,
	  a consequence of the imphttab calling function */    
	if (wf32d->chip == 1){
		if (PutKeyFlt (&x->sci.hdr, "PHTFLAM1", obs.phtflam1,
					"photometry scaling for chip 1")){
			return (status);
		}     
		if (PutKeyFlt (x->globalhdr, "PHTFLAM1", obs.phtflam1,
					"photometry scaling for chip 1")){
			return (status);
		}              
	}

	if (wf32d->chip == 2){
		if (PutKeyFlt (&x->sci.hdr, "PHTFLAM2", obs.phtflam2,
					"photometry scaling for chip 2")) {
			return (status);
		}
	}


	FreePhotPar (&obs);

	return (status);
}

/* This function converts the PHOTMODE string into an OBSMODE string
   suitable for use with synphot functions.
   PHOTMODE - all upper case component names separated by blanks
   OBSMODE - all lower case names separated by commas
NOTE: The operation occurs **in-place** on the photmode string.
 */
static void Phot2Obs (char *photmode, char *obsmode) {

	char blank = 32, comma = 44;
	int i, len, c1;

	len = strlen (photmode);

	for (i=0; i < len; i++) {
		c1 = photmode[i];
		if (c1 == blank) {
			obsmode[i] = comma;
		} else
			obsmode[i] = tolower(c1);
	}
	obsmode[len] = '\0';
}

