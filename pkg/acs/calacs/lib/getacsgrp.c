# include <stdio.h>
# include <string.h>
# include <math.h>	/* sqrt */
# include "hstio.h"
# include "acs.h"
# include "acsinfo.h"
# include "hstcalerr.h"

/* This routine gets keyword values from the SCI extension header.

   Warren Hack, 1998 June 2:
   	Revised for ACS data...  Reads in keywords necessary for processing
		from EXT headers.
*/

int GetACSGrp (ACSInfo *acs, Hdr *hdr) {

/* arguments:
ACSInfo *acs   io: calibration switches and info
Hdr *hdr         i: header of current extension
*/

	extern int status;

	int sdqflags;			/* serious data quality flags */
	int GetKeyDbl (Hdr *, char *, int, double, double *);
	int GetKeyInt (Hdr *, char *, int, int, int *);

	/* Get generic parameters. */

	/* Find out which data quality bits are considered serious;
	   default value means all bits are serious.
	*/
	if (GetKeyInt (hdr, "SDQFLAGS", USE_DEFAULT, MAX_DQ, &sdqflags))
	    return (status);
	acs->sdqflags = (short) sdqflags;


	/* Get the pixel size (ignore corner location) from ltm & ltv. */
	/* Set acs->bin here, size of pixel in detector coordinates ... */
	acs->bin[0] = 1;
	acs->bin[1] = 1;
	
	/* Get CCD-specific parameters. */
	if (GetKeyInt (hdr, "CCDCHIP", USE_DEFAULT, 1, &acs->chip))
	    return (status);
    if (GetKeyDbl (hdr, "LTV1", USE_DEFAULT, 0., &acs->offsetx))
        return (status);
    if (GetKeyDbl (hdr, "LTV2", USE_DEFAULT, 0., &acs->offsety))
        return (status);

	/* If images have been combined (e.g. by cosmic-ray rejection),
	   then determine the number of images that were combined together;
	   we need this for bias image subtraction.
	   (This isn't really a CCD-specific keyword, but it does only affect
	   CCD data in the context of acsccd.)
	*/
	if (GetKeyInt (hdr, "NCOMBINE", USE_DEFAULT, 1, &acs->ncombine))
	    return (status);
	if (acs->ncombine < 1) {
	    sprintf (MsgText, "NCOMBINE = %d, reset to one.", acs->ncombine);
	    trlwarn (MsgText);
		acs->ncombine = 1;
	}
	/* Get MAMA-specific parameters. */

	if (acs->detector == MAMA_DETECTOR) {

	    if (GetKeyDbl (hdr, "GLOBRATE", NO_DEFAULT, 0., &acs->globrate))
		return (status);

	}

	return (status);
}
