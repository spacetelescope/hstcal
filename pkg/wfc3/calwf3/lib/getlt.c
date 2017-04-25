# include <stdio.h>
# include "hstio.h"
# include "wf3.h"
# include "hstcalerr.h"

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

	extern int status;
	int GetKeyDbl (Hdr *, char *, int, double, double *);

	/* Get linear transformation info for pixel coordinates. */
	if (GetKeyDbl (hdr, "LTM1_1", USE_DEFAULT, 1., &ltm[0]))
	    return (status);
	if (GetKeyDbl (hdr, "LTM2_2", USE_DEFAULT, 1., &ltm[1]))
	    return (status);
	if (GetKeyDbl (hdr, "LTV1", USE_DEFAULT, 0., &ltv[0]))
	    return (status);
	if (GetKeyDbl (hdr, "LTV2", USE_DEFAULT, 0., &ltv[1]))
	    return (status);

	if (ltm[0] <= 0. || ltm[1] <= 0.) {
	    sprintf (MsgText, "(LTM1_1, LTM2_2) = (%.8g, %.8g) is invalid",
		     ltm[0], ltm[1]);
	    trlerror (MsgText);
        status = INVALID_VALUE;
	    return (status);
	}

	return (status);
}

/* This routine gets the LTV and LTM keyword values and converts LTV
   from one to zero indexing.
*/

int GetLT0 (Hdr *hdr, double *ltm, double *ltv) {

/* arguments:
Hdr *hdr         i: extension header
double ltm[2]    o: diagonal elements of MWCS linear transformation matrix
double ltv[2]    o: MWCS linear transformation vector, zero indexed
*/

	extern int status;

	int GetLT (Hdr *, double *, double *);

	/* Get the ltm and ltv values from the header. */
	if (GetLT (hdr, ltm, ltv))
	    return (status);

	/* The LTV values in the header are for one-indexed pixel
	   coords.  Since we use zero indexing, adjust LTV accordingly.
	*/
	ltv[0] += (ltm[0] - 1.);
	ltv[1] += (ltm[1] - 1.);

	return (status);
}
