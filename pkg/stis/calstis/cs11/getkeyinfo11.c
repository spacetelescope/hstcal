# include <stdio.h>
# include <string.h>

# include "c_iraf.h"
# include "hstio.h"

# include "stis.h"
# include "calstis11.h"
# include "err.h"
# include "stisdef.h"

/* This routine gets keyword values from the primary header.

   Phil Hodge, 1998 Aug 13:
	Get TEXPSTRT.

   Phil Hodge, 1999 Oct 1:
	Get ATODGAIN (or CCDGAIN).
*/

int GetKeyInfo11 (StisInfo11 *sts, Hdr *phdr) {

/* arguments:
StisInfo11 *sts  io: calibration switches and info
Hdr *phdr        i: primary header
*/

	int status;
	int nextend;			/* number of FITS extensions */
	int no_def = 0;			/* missing keyword is fatal error */
	int use_def = 1;		/* use default if keyword is missing */

	if ((status = Get_KeyS (phdr, "ROOTNAME",
                                no_def, "", sts->rootname, STIS_CBUF)))
	    return (status);

	/* Name of calibration lamp (should be HITM1 or HITM2). */
	if ((status = Get_KeyS (phdr, "SCLAMP",
                                no_def, "", sts->sclamp, STIS_CBUF)))
	    return (status);

	if ((status = Get_KeyS (phdr, "OBSMODE",
                                no_def, "", sts->obsmode, STIS_CBUF)))
	    return (status);

	if ((status = Get_KeyS (phdr, "APERTURE",
                                no_def, "", sts->aperture, STIS_CBUF)))
	    return (status);

	if ((status = Get_KeyS (phdr, "OPT_ELEM",
                                no_def, "", sts->opt_elem, STIS_CBUF)))
	    return (status);

	if ((status = Get_KeyS (phdr, "DETECTOR",
                                no_def, "", sts->det, STIS_CBUF)))
	    return (status);

	if (strcmp (sts->det, "NUV-MAMA") == 0)
	    sts->detector = NUV_MAMA_DETECTOR;
	else if (strcmp (sts->det, "FUV-MAMA") == 0)
	    sts->detector = FUV_MAMA_DETECTOR;
	else if (strcmp (sts->det, "CCD") == 0)
	    sts->detector = CCD_DETECTOR;
	else
	    sts->detector = UNKNOWN_DETECTOR;

	/* Central wavelength. */
	if ((status = Get_KeyI (phdr, "CENWAVE", no_def, 0, &sts->cenwave)))
	    return (status);

	/* Find out how many extensions there are in this file. */
	if ((status = Get_KeyI (phdr, "NEXTEND",
                                use_def, EXT_PER_GROUP, &nextend)))
	    return (status);

	/* Convert number of extensions to number of SingleGroups. */
	sts->nimages = nextend / EXT_PER_GROUP;

	/* Get the exposure start time.  We only need to get this from
	   the wavecal header, and only for a CCD HITM exposure, to check
	   whether the external shutter is closed.
	*/
	if ((status = Get_KeyD (phdr, "TEXPSTRT", use_def, 0., &sts->texpstrt)))
	    return (status);

	/* Get the gain.  If the ATODGAIN keyword is not found or is
	   less than or equal to zero, get CCDGAIN instead.
	*/
	if ((status = Get_KeyD (phdr, "ATODGAIN", use_def, 0., &sts->gain)))
	    return (status);
	if (sts->gain <= 0.) {
	    if ((status = Get_KeyD (phdr, "CCDGAIN", use_def, 1., &sts->gain)))
		return (status);
	}

	return (0);
}
