# include <stdio.h>
# include <string.h>
# include <math.h>	/* sqrt */
# include "hstio.h"
# include "wf3.h"
# include "wf3info.h"
# include "err.h"

/* This routine gets keyword values from the SCI extension header.

   Warren Hack, 1998 June 2:
   	Revised for ACS data...  Reads in keywords necessary for processing
		from EXT headers.

   Howard Bushouse, 2004 Feb 20:
	Eliminated hardwiring wf3->bin to 1 (ACS legacy) and instead read
	the BINAXIS keywords from the sci extension header to populate it.
*/

int GetGrp (WF3Info *wf3, Hdr *hdr) {

/* arguments:
WF3Info *wf3   io: calibration switches and info
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
	wf3->sdqflags = (short) sdqflags;

	/* Get CCD-specific parameters. */
	if (GetKeyInt (hdr, "CCDCHIP", USE_DEFAULT, 1, &wf3->chip))
	    return (status);
	if (GetKeyDbl (hdr, "LTV1", USE_DEFAULT, 0., &wf3->offsetx))
	    return (status);
	if (GetKeyDbl (hdr, "LTV2", USE_DEFAULT, 0., &wf3->offsety))
	    return (status);
	if (GetKeyInt (hdr, "BINAXIS1", USE_DEFAULT, 1, &wf3->bin[0]))
	    return (status);
	if (GetKeyInt (hdr, "BINAXIS2", USE_DEFAULT, 1, &wf3->bin[1]))
	    return (status);

	/* If images have been combined (e.g. by cosmic-ray rejection),
	   then determine the number of images that were combined together;
	   we need this for bias image subtraction.
	   (This isn't really a CCD-specific keyword, but it does only affect
	   CCD data in the context of wf3ccd.)
	*/
	if (GetKeyInt (hdr, "NCOMBINE", USE_DEFAULT, 1, &wf3->ncombine))
	    return (status);
	if (wf3->ncombine < 1) {
	    sprintf (MsgText, "NCOMBINE = %d, reset to one.", wf3->ncombine);
	    trlwarn (MsgText);
	    wf3->ncombine = 1;
	}

	return (status);
}
