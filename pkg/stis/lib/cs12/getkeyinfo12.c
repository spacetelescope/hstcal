# include <stdio.h>
# include <string.h>
# include "c_iraf.h"
# include "hstio.h"

# include "stis.h"
# include "calstis12.h"
# include "hstcalerr.h"
# include "stisdef.h"

/* This routine gets keyword values from the primary header. */

int GetKeyInfo12 (StisInfo12 *sts, Hdr *phdr) {

/* arguments:
StisInfo12 *sts   io: calibration switches and info
Hdr *phdr         i: primary header
*/

	int status;
	int nextend;			/* number of FITS extensions */
	int no_def = 0;			/* missing keyword is fatal error */
	int use_def = 1;		/* use default if keyword is missing */

	if ((status = Get_KeyS (phdr, "ROOTNAME",
                                no_def, "", sts->rootname, STIS_CBUF)))
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
	if ((status = Get_KeyI (phdr, "NEXTEND", use_def, EXT_PER_GROUP,
                                &nextend)))
	    return (status);

	/* Convert number of extensions to number of SingleGroups. */
	sts->nimages = nextend / EXT_PER_GROUP;

	return (0);
}
