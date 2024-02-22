# include <stdio.h>
# include "c_iraf.h"
# include "hstio.h"

# include "stis.h"
# include "calstis12.h"
# include "hstcalerr.h"
# include "stisdef.h"

# define SECONDS_IN_A_DAY  86400.
# define ARCSEC_PER_DEGREE  3600.

/* This routine gets the exposure time and the midpoint of the exposure
   from a science data extension header.

   Phil Hodge, 1998 Mar 18:
	Remove scale from the calling sequence.

   Phil Hodge, 1998 Oct 5:
	Change status value 1001 to GENERIC_ERROR_CODE.
*/

int GetSciGrp (Hdr *hdr, double *exptime, double *midpt) {

/* arguments:
Hdr *hdr         i: header of current imset
double *exptime  o: the exposure time (converted to days)
double *midpt    o: the time (MJD) of the middle of the exposure
*/

	int status;

	double expstart, expend;
	int no_default = 0;	/* missing keyword is fatal error */

	if ((status = Get_KeyD (hdr, "EXPTIME", no_default, 0., exptime)))
	    return (status);
	if (*exptime < 0.) {
	    printf ("ERROR    GetSciGrp:  exposure time = %.6g is invalid.\n",
		*exptime);
	    return (GENERIC_ERROR_CODE);
	}

	*exptime /= SECONDS_IN_A_DAY;

	if ((status = Get_KeyD (hdr, "EXPSTART", no_default, 0., &expstart)))
	    return (status);
	if ((status = Get_KeyD (hdr, "EXPEND", no_default, 0., &expend)))
	    return (status);

	*midpt = (expstart + expend) / 2.;

	return (0);
}
