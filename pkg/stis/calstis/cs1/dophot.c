# include <ctype.h>
# include <stdio.h>
# include <string.h>

# include "c_iraf.h"
# include "hstio.h"
# include "imphttab.h"
# include "stis.h"
# include "calstis1.h"
# include "err.h"
# include "stisdef.h"
# include "stistds.h"

static void Phot2Obs (char *, char *);
static int removeA2D (char *, char *);
static double tempFactor (TdsInfo *, double, double);

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

   Phil Hodge, 2011 Aug 31:
	Print warning messages if phot parameters appear to be garbage.
	Set photcorr to IGNORED (which results in PHOTCORR = "SKIPPED" in
	the output header) if photmode not found or parameters are garbage.
	Also set photcorr to IGNORED if the parameters are -9999.

   Phil Hodge, 2011 Nov 8:
	Add function removeA2D to look for "A2D<gain>" in the photmode string
	and, if found, remove it when copying to a temporary photmode.  If
	"A2D" was found in photmode, multiply photflam by atodgain, rather
	than letting the imphttab handle this factor.

   Phil Hodge, 2011 Nov 17:
	Add function tempFactor to interpolate the temperature dependence
	at the pivot wavelength.  Divide photflam by the returned factor.
*/

int doPhot (StisInfo1 *sts, SingleGroup *x) {

/* arguments:
StisInfo1 *sts     i: calibration switches, etc
SingleGroup *x    io: image to be calibrated; primary header is modified
*/

	int status;

	PhotPar obs;

	char photmode[STIS_LINE+1];	/* the photmode string */
	char tempphot[STIS_LINE+1];	/* photmode without A2D */
	char obsmode[STIS_LINE+1];	/* based on the photmode string */
	int use_default = 1;	/* use default value if keyword is missing */
	int found_a2d = 0;	/* true if there's a CCD A2D component */
	TdsInfo tds;		/* time-dep, but use only temperature corr */

	int GetTds1 (char *, char *, TdsInfo *);
	void FreeTds1 (TdsInfo *);

	/* Get PHOTMODE from the primary header. */
	if ((status = Get_KeyS (x->globalhdr, "PHOTMODE",
		use_default, "", photmode, STIS_LINE)) != 0)
	    return (status);

	found_a2d = removeA2D (photmode, tempphot);
	Phot2Obs (tempphot, obsmode);

	/* Initialize PhotPar struct. */
	InitPhotPar (&obs, sts->phot.name, sts->phot.pedigree);

	/* Get the photometry values. */
	if ((status = GetPhotTab (&obs, obsmode)) != 0) {
	    printf ("Warning  photmode '%s' not found.\n", obsmode);
	    FreePhotPar (&obs);
	    sts->photcorr = IGNORED;
	    return 0;
	}

	/* The flag value that indicates that the time of observation is
	   outside the range of times in the imphttab is actually -9999.,
	   but test on -9990. to avoid roundoff problems.
	   Extrapolation is not supported, so the values were not computed.
	   Set photcorr to IGNORED (which results in SKIPPED), but set
	   the keywords to -9999. anyway.
	*/
	if (obs.photflam < -9990. &&
	    obs.photplam < -9990. &&
	    obs.photbw < -9990.) {
	    sts->photcorr = IGNORED;
	}

	/* Make further changes to photflam. */

	if (found_a2d && sts->photcorr == PERFORM)
	    obs.photflam *= sts->atodgain;

	/* Get values from the TDS table */
	if (sts->tdscorr == PERFORM && sts->photcorr == PERFORM) {
	    if ((status = GetTds1(sts->tdstab.name, sts->opt_elem, &tds)))
		return (status);
	    obs.photflam /= tempFactor(&tds, obs.photplam, sts->detector_temp);
	    FreeTds1 (&tds);
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

/* This function looks for "A2D" in the photmode string and, if present,
   removes it.
*/

static int removeA2D (char *photmode, char *tempphot) {

	int i, j, k, n, len;
	int foundit;		/* found "A2D"? */
	int found_space;	/* found a space following "A2D<digit>"? */

	len = strlen(photmode);

	j = len;
	k = len;
	foundit = 0;
	for (n = 0;  n < len-2;  n++) {
	    if (toupper(photmode[n])   == 'A' &&
		photmode[n+1]          == '2' &&
		toupper(photmode[n+2]) == 'D') {
		foundit = 1;
		/* look for a space following "A2D" and digit(s) */
		found_space = 0;
		for (k = n + 3;  k < len;  k++) {
		    if (photmode[k] == ' ') {
			found_space = 1;
			break;
		    }
		}
		if (n == 0) {
		    j = 0;
		    k++;		/* start copying past the space */
		} else if (found_space) {
		    j = n - 1;		/* skip the previous space */
		} else {
		    j = n - 1;
		}
	   }
	}

	i = 0;
	for (n = 0;  n < j;  n++) {
	    tempphot[i] = photmode[n];
	    i++;
	}
	for (n = k;  n < len;  n++) {
	    tempphot[i] = photmode[n];
	    i++;
	}
	tempphot[i] = '\0';

	return foundit;
}

static double tempFactor(TdsInfo *tds, double photplam, double temperature) {

/* arguments:
TdsInfo *tds         i: temperature-dependent sensitivity info
double photplam      i: pivot wavelength (Angstroms)
double temperature   i: detector temperature (degrees Celsius)
   function value:
factor                  divide photflam by this factor
*/

	double sensitivity, factor;
	int starti;
	/*int i;*/

	if (temperature >= 0.) {
	    starti = 0;
	    sensitivity = interp1d (photplam, tds->wl, tds->temp_sens,
				    tds->nwl, &starti);
	    factor = 1. + sensitivity * (temperature - tds->ref_temp);
	} else {
	    factor = 1.;
	}

	return factor;
}

/* This function converts the photmode string into an obsmode string
   suitable for use with the IMPHTTAB.
   photmode - all upper case component names separated by blanks
   obsmode - all lower case names separated by commas
*/
static void Phot2Obs (char *photmode, char *obsmode) {

	char blank = ' ', comma = ',';
	int i, len, c1;

	len = strlen(photmode);

	for (i=0;  i < len;  i++) {
	    c1 = photmode[i];
	    if (c1 == blank) {
		obsmode[i] = comma;
	    } else {
		obsmode[i] = tolower(c1);
	    }
	}
	obsmode[len] = '\0';
}
