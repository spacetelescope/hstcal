# include <stdio.h>
# include "c_iraf.h"
# include "hstio.h"

# include "stis.h"
# include "stisdef.h"

/* This routine gets the LTV and LTM keyword values and converts LTV
   from one to zero indexing.
*/

int GetLT0 (Hdr *hdr, double *ltm, double *ltv) {

/* arguments:
Hdr *hdr         i: extension header
double ltm[2]    o: diagonal elements of MWCS linear transformation matrix
double ltv[2]    o: MWCS linear transformation vector, zero indexed
*/

	int status;

	/* Get the ltm and ltv values from the header. */
	if ((status = GetLT (hdr, ltm, ltv)))
	    return (status);

	/* The LTV values in the header are for one-indexed pixel
	   coords.  Since we use zero indexing, adjust LTV accordingly.
	*/
	ltv[0] += (ltm[0] - 1.);
	ltv[1] += (ltm[1] - 1.);

	return (0);
}
