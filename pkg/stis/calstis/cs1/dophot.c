# include <ctype.h>
# include <stdio.h>
# include <string.h>

# include "c_iraf.h"
# include "hstio.h"
# include "imphttab.h"
# include "stis.h"
# include "calstis1.h"
# include "stiserr.h"
# include "stisdef.h"

static void Phot2Obs (char *, char *);

/* This routine reads the photmode string from the primary header (which
   must therefore have already been updated) and calls a function that reads
   the IMPHTTAB to compute the inverse sensitivity, reference magnitude
   (actually a constant), pivot wavelength, and RMS bandwidth.  These are
   then written to keywords in the primary header.

   Phil Hodge, 1997 Nov 13:
	Phot table I/O extracted to GetPhot1; also call GetApThr1;
	rename phot photkey and use phot for PhotInfo.

   Phil Hodge, 1998 Oct 5:
	Change status value 1001 to GENERIC_ERROR_CODE.

   Paul Barrett, 2003 Sep 25:
        Add time-dependent sensitivity (TDS) for imaging mode.

   Phil Hodge, 2011 May 5:
	Rewrite to use the functions to read and interpret the imphttab.
	This was based on dophot.c in calacs/acs2d/.

   Phil Hodge, 2011 July 26:
	Print warning messages if phot parameters appear to be garbage.
*/

int doPhot (StisInfo1 *sts, SingleGroup *x) {

/* arguments:
StisInfo1 *sts     i: calibration switches, etc
SingleGroup *x    io: image to be calibrated; primary header is modified
*/

	int status;

	PhotPar obs;

	char photmode[STIS_LINE+1];	/* the photmode string */
	char obsmode[STIS_LINE+1];	/* based on the photmode string */
	int use_default = 1;	/* use default value if keyword is missing */

	/* Get PHOTMODE from the primary header. */
	if ((status = Get_KeyS (x->globalhdr, "PHOTMODE",
		use_default, "", photmode, STIS_LINE)) != 0)
	    return (status);

	Phot2Obs (photmode, obsmode);

	/* Initialize PhotPar struct. */
	InitPhotPar (&obs, sts->phot.name, sts->phot.pedigree);

	/* Get the photometry values. */
	if ((status = GetPhotTab (&obs, obsmode)) != 0) {
	    printf ("Warning  photmode '%s' not found.\n", photmode);
	    FreePhotPar (&obs);
	    return 0;
	}
	if (obs.photbw < 0. || obs.photbw > 1.e5) {
	    printf ("Warning  For photmode '%s', values are strange:\n",
		photmode);
	    printf ("         photflam = %g\n", obs.photflam);
	    printf ("         photplam = %g\n", obs.photplam);
	    printf ("         photbw   = %g\n", obs.photbw);
	    printf ("         photzpt  = %g\n", obs.photzpt);
	    FreePhotPar (&obs);
	    return 0;
	}

	/* Update the photometry keyword values in the primary header. */

	if ((status = Put_KeyF (x->globalhdr, "PHOTFLAM", obs.photflam,
			"inverse sensitivity")) != 0)
	    return (status);

	if ((status = Put_KeyF (x->globalhdr, "PHOTZPT", obs.photzpt,
			"zero point")) != 0)
	    return (status);

	if ((status = Put_KeyF (x->globalhdr, "PHOTPLAM", obs.photplam,
			"pivot wavelength")) != 0)
	    return (status);

	if ((status = Put_KeyF (x->globalhdr, "PHOTBW", obs.photbw,
			"RMS bandwidth")) != 0)
	    return (status);

	FreePhotPar (&obs);

	return (0);
}

/* This function converts the PHOTMODE string into an OBSMODE
   string suitable for use with synphot functions.
   PHOTMODE - all upper case component names separated by blanks
   OBSMODE - all lower case names separated by commas
*/
static void Phot2Obs (char *photmode, char *obsmode) {

	char blank = 32, comma = 44;
	int i, len, c1;

	len = strlen(photmode);

	for (i=0; i < len; i++) {
	    c1 = photmode[i];
	    if (c1 == blank) {
		obsmode[i] = comma;
	    } else {
		obsmode[i] = tolower(c1);
	    }
	}
	obsmode[len] = '\0';
}
