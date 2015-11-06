# include <stdio.h>
# include "c_iraf.h"
# include "hstio.h"

# include "stis.h"
# include "stisdef.h"

/* This routine gets the LTV and LTM keyword values.  Note that this
   routine returns the values just as read from the header, which means
   in particular that the LTV values are for one-indexed pixel coordinates.
*/

int GetLT (Hdr *hdr, double *ltm, double *ltv) {

/* arguments:
Hdr *hdr         i: extension header
double ltm[2]    o: diagonal elements of MWCS linear transformation matrix
double ltv[2]    o: MWCS linear transformation vector, one indexed
*/

	int status;
	int use_def = 1;	/* use default value if keyword is missing */

	/* Get linear transformation info for pixel coordinates. */
	if ((status = Get_KeyD (hdr, "LTM1_1", use_def, 1., &ltm[0])))
	    return (status);
	if ((status = Get_KeyD (hdr, "LTM2_2", use_def, 1., &ltm[1])))
	    return (status);
	if ((status = Get_KeyD (hdr, "LTV1", use_def, 0., &ltv[0])))
	    return (status);
	if ((status = Get_KeyD (hdr, "LTV2", use_def, 0., &ltv[1])))
	    return (status);

	if (ltm[0] <= 0. || ltm[1] <= 0.) {
	    printf ("ERROR:  (LTM1_1, LTM2_2) = (%.8g, %.8g) is invalid\n",
                        ltm[0], ltm[1]);
	    return (2030);
	}

	return (0);
}
