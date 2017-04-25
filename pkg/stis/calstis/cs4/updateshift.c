# include <stdio.h>
# include <string.h>

# include "c_iraf.h"
# include "hstio.h"

# include "stis.h"
# include "calstis4.h"
# include "hstcalerr.h"
# include "stisdef.h"

/* This routine updates the values of keywords SHIFTA1 and SHIFTA2 in the
   extension header.

   Phil Hodge, 2000 Jan 7
*/

int UpdateShift (StisInfo4 *sts, int extver, double w_shift, double s_shift) {

/* arguments:
StisInfo4 *sts           i: info
int extver               i: imset number
double w_shift, s_shift  i: the shifts that were determined
*/

	IODescPtr im;		/* image */
	Hdr hdr;		/* extension header */
	int status;

	/* Open the extension header. */
	initHdr (&hdr);
	im = openUpdateImage (sts->input, "SCI", extver, &hdr);
	if (hstio_err())
	    return (HEADER_PROBLEM);

	if (sts->dispaxis == 1) {
	    if ((status = Put_KeyD (&hdr, "SHIFTA1", w_shift,
                                    "MSM shift in dispersion direction")))
		return (status);
	    if ((status = Put_KeyD (&hdr, "SHIFTA2", s_shift,
                                    "MSM shift in spatial direction")))
		return (status);
	} else if (sts->dispaxis == 2) {
	    if ((status = Put_KeyD (&hdr, "SHIFTA1", s_shift,
                                    "MSM shift in spatial direction")))
		return (status);
	    if ((status = Put_KeyD (&hdr, "SHIFTA2", w_shift,
                                    "MSM shift in dispersion direction")))
		return (status);
	}

	putHeader (im);
	if (hstio_err())
	    return (OPEN_FAILED);
	closeImage (im);
	freeHdr (&hdr);

	return (0);
}
