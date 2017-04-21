# include <stdio.h>
# include <stdlib.h>	/* calloc */
# include <string.h>

# include "hstio.h"

# include "wf3.h"
# include "wf3info.h"
# include "hstcalerr.h"

/* 
	This function builds the PHOTMODE string for the image header.
	The string will be used as input to 'pysynphot' for determining
	the values of the photometric keywords.
	
   Howard Bushouse, 2001 May 8:
	For WFC3 revised FILTER1,FILTER2 usage to just FILTER, changed
	location of photmode keyword from primary header to SCI extension,
	changed CCD detector name to UVIS, and removed inclusion
	of aperture name in photmode string.
   H.Bushouse, 2001 Nov 16:
	Updated to track CALACS changes and eliminated ACS-specific items
	that aren't needed for WFC3 use (coron, ramp filters, etc.).
   H.Bushouse, 2002 May 21:
	Modified to make universal with WF3IR routines by replacing
	SingleGroup *x input argument with Hdr *hdr.
   H.Bushouse, 2003 Oct 16:
	Modified print format of A2D gain component to accomodate
	non-integer values for WFC3.
   H.Bushouse, 2007 Feb 15:
	Changed UVIS channel detector keyword to always use "UVIS1".
	Changed use of "A2Dx" gain keyword to "DN" and eliminated use of
	it for UVIS images because flatfielding leaves them in units of
	electrons, not counts.
   H.Bushouse, 2009 Jan 14:
	Changed UVIS channel to use UVIS1 and UVIS2 for each chip and to
	use new "cal" keyword for UVIS channel. Removed use of DN keyword
	for IR exposures because they're now in units of electrons.
*/

int PhotMode (WF3Info *wf32d, Hdr *hdr) {

/* arguments:
WF3Info *wf3     i: calibration switches, etc
Hdr *hdr	io: image header to be modified
*/

	extern int status;

	char *photstr;		/* the photmode string */
	char *scratch;		/* scratch space string */

	int PutKeyStr (Hdr *, char *, char *, char *);

	photstr = calloc (SZ_LINE+1, sizeof (char));
	scratch = calloc (SZ_FITS_REC+1, sizeof (char));
	if (photstr == NULL || scratch == NULL)
	    return (status = OUT_OF_MEMORY);

	/* Copy instrument name to PHOTMODE string */
	strcpy (photstr, "WFC3");

	/* Add detector name and, for UVIS, the chip number */
	if (wf32d->detector == IR_DETECTOR)
	    strcat (photstr, " IR");
	else if (wf32d->detector == CCD_DETECTOR) {
	    strcat (photstr, " UVIS");
	    sprintf (scratch, "%d", wf32d->chip);
	    strcat (photstr, scratch);
	}

	/* Add the filter name */
	if (strncmp ("CLEAR", wf32d->filter, 5) != 0) {
	    strcat (photstr, " ");
	    strcat (photstr, wf32d->filter);
	}

	/* Add the CAL keyword for UVIS exposures */
	if (wf32d->detector == CCD_DETECTOR)
	    strcat (photstr, " CAL");

	if (wf32d->verbose) {
	    sprintf (MsgText, "Keyword PHOTMODE built as: %s", photstr);
	    trlmessage (MsgText);
	}

	/* Update PHOTMODE in the header. */
	if (PutKeyStr (hdr, "PHOTMODE", photstr, "list of components"))
	    return (status);

	free (scratch);
	free (photstr);

	return (status);
}
