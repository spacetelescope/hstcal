# include <stdio.h>

# include "c_iraf.h"
# include "hstio.h"

# include "stis.h"
# include "err.h"
# include "stisdef.h"

/* This routine gets the exposure time and the time of the middle of
   the exposure from an extension header.

   Phil Hodge, 1998 Oct 5:
	Change status value 1001 to GENERIC_ERROR_CODE.

   Phil Hodge, 2008 Dec 12:
	Get keyword IMSET_OK (add to calling sequence).
*/

int GetTimes11 (Hdr *hdr, double *exptime, double *midpt, int *imset_ok) {

/* arguments:
Hdr *hdr         i: header of current imset
double *exptime  o: the exposure time (sec)
double *midpt    o: the time (MJD) of the middle of the exposure
int *imset_ok    o: value of header keyword IMSET_OK
*/

	int status;

	double expstart, expend;
	int no_default = 0;	/* missing keyword is fatal error */
	int use_def = 1;	/* use default if keyword is missing */
	Bool value;		/* value of IMSET_OK keyword */

	/* Get IMSET_OK (default is true), which will be false if the current
	   imset was flagged as having zero exptime or constant pixel values.
	*/
	if ((status = Get_KeyB (hdr, "IMSET_OK", use_def, True, &value)) != 0)
	    return status;
	if (value)
	    *imset_ok = 1;
	else
	    *imset_ok = 0;

	if ((status = Get_KeyD (hdr, "EXPTIME", no_default, 0., exptime)))
	    return (status);
	if (*exptime < 0.) {
	    printf ("ERROR    Exposure time = %.6g is invalid.\n",
		*exptime);
	    return (GENERIC_ERROR_CODE);
	}

	if ((status = Get_KeyD (hdr, "EXPSTART", no_default, 0., &expstart)))
	    return (status);
	if ((status = Get_KeyD (hdr, "EXPEND", no_default, 0., &expend)))
	    return (status);

	*midpt = (expstart + expend) / 2.;

	return (0);
}
