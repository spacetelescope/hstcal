# include <stdio.h>
# include <stdlib.h>	/* calloc */
# include <string.h>

# include "c_iraf.h"
# include "hstio.h"

# include "stis.h"
# include "calstis1.h"
# include "hstcalerr.h"
# include "stisdef.h"

/*
   Phil Hodge, 2000 Nov 1:
	Change the test on whether opt_elem eq "MIRROR"; the test is now
	whether opt_elem begins with "MIR".  The opt_elem value will be
	appended to the photmode string unless opt_elem is any of the
	mirror names.

   Phil Hodge, 2011 May 5:
	Also include MJD#<date> in the photmode string.

   Phil Hodge, 2011 Nov 30:
	Include "MIRROR" in the photmode string.
*/

int PhotMode (StisInfo1 *sts, SingleGroup *x) {

/* arguments:
StisInfo1 *sts    i: calibration switches, etc
SingleGroup *x    io: image to be calibrated; primary header is modified
*/

	int status;

	char *photstr;		/* the photmode string */
	char *scratch;
	int cenwave;		/* central wavelength */
	int use_default = 1;	/* use default if keyword is missing */

	photstr = calloc (STIS_LINE+1, sizeof (char));
	scratch = calloc (STIS_FITS_REC+1, sizeof (char));
	if (photstr == NULL || scratch == NULL)
	    return (OUT_OF_MEMORY);

	strcpy (photstr, "STIS");

	strcat (photstr, " ");
	strcat (photstr, sts->aperture);

	/* If the detector is the CCD, include the gain. */
	if (sts->detector == CCD_DETECTOR) {
	    sprintf (scratch, " A2D%d", sts->ccdgain);
	    strcat (photstr, scratch);
	}

	if (strncmp (sts->opt_elem, "MIR", 3) == 0) {
	    /* Include the mirror name as "MIRROR" in photmode. */
	    strcat (photstr, " MIRROR");
	} else {
	    strcat (photstr, " ");
	    strcat (photstr, sts->opt_elem);
	}

	if (sts->detector == NUV_MAMA_DETECTOR)
	    strcat (photstr, " NUVMAMA");
	else if (sts->detector == FUV_MAMA_DETECTOR)
	    strcat (photstr, " FUVMAMA");
	else if (sts->detector == CCD_DETECTOR)
	    strcat (photstr, " CCD");

	/* Include the central wavelength, if spectroscopic mode. */
	if ((status = Get_KeyS (x->globalhdr, "OBSTYPE",
                use_default, "unknown", scratch, STIS_FITS_REC)))
	    return (status);
	if (strcmp (scratch, "SPECTROSCOPIC") == 0) {
	    if ((status = Get_KeyI (x->globalhdr, "CENWAVE",
                                    use_default, -1, &cenwave)))
		return (status);
	    if (cenwave > 0) {
		sprintf (scratch, " %d", cenwave);
		strcat (photstr, scratch);
	    }
	}

	/* Add 'mjd#' to PHOTMODE string. */
	sprintf (scratch, " MJD#%0.4f", sts->expstart);
	strcat (photstr, scratch);

	/* Update PHOTMODE in the primary header. */
	if ((status = Put_KeyS (x->globalhdr, "PHOTMODE", photstr,
                                "list of components")))
	    return (status);

	free (scratch);
	free (photstr);

	return (0);
}
