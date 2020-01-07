# include <stdio.h>
# include <math.h>	/* fabs */

# include "c_iraf.h"
# include "hstio.h"

# include "stis.h"
# include "calstis12.h"
# include "hstcalerr.h"
# include "stisdef.h"

static int UpdateShift (Hdr *, double, double);

/* This should be smaller in absolute value than UNDEFINED_SHIFT. */
# define MUCH_TOO_BIG  2048.	/* much too large to be a good shift */

/* This routine records the shifts in the keywords SHIFTA1 and SHIFTA2
   for the SCI, ERR, and DQ science data extension headers.

   If either or both shifts could not be determined, they will have values
   UNDEFINED_SHIFT.  (Values in the science file header, on the other hand,
   should never have the UNDEFINED_SHIFT value.)  In this case, the SHIFTAi
   value(s) will be set to zero.

   The name of this routine goes back to the time when the TARGPOSi keywords
   were updated.

   Phil Hodge, 1998 Mar 18:
	Also print values of shifta1 and shifta2.

   Phil Hodge, 2000 June 16:
	Don't modify CRPIX1 or CRPIX2, and don't add HISTORY records.
*/

int TargPos (StisInfo12 *sci, int extver, double shift1, double shift2) {

/* arguments:
StisInfo12 *sci        i: info about science data file
int extver             i: EXTVER number of extensions to update
double shift1, shift2  i: shift to be assigned to SHIFTAi keywords
*/

	int status;

	IODescPtr im;		/* descriptor for an image */
	Hdr hdr;		/* header for an extension */

	/* Update SCI extension. */

	initHdr (&hdr);
	im = openUpdateImage (sci->input, "SCI", extver, &hdr);
	if (hstio_err())
	    return (OPEN_FAILED);

	/* Update SHIFTAi. */
	if ((status = UpdateShift (&hdr, shift1, shift2)))
	    return (status);

	putHeader (im);
	if (hstio_err())
	    return (OPEN_FAILED);
	closeImage (im);
	freeHdr (&hdr);

	/* Update ERR extension. */

	im = openUpdateImage (sci->input, "ERR", extver, &hdr);
	if (hstio_err())
	    return (OPEN_FAILED);

	if ((status = UpdateShift (&hdr, shift1, shift2)))
	    return (status);

	putHeader (im);
	if (hstio_err())
	    return (OPEN_FAILED);
	closeImage (im);
	freeHdr (&hdr);

	/* Update DQ extension. */

	im = openUpdateImage (sci->input, "DQ", extver, &hdr);
	if (hstio_err())
	    return (OPEN_FAILED);

	if ((status = UpdateShift (&hdr, shift1, shift2)))
	    return (status);

	putHeader (im);
	if (hstio_err())
	    return (OPEN_FAILED);
	closeImage (im);
	freeHdr (&hdr);

	/* write to trailer */
	if (fabs (shift1) < MUCH_TOO_BIG)	/* shift is OK */
	    printf ("         SHIFTA1 set to %.6g\n", shift1);
	if (fabs (shift2) < MUCH_TOO_BIG)
	    printf ("         SHIFTA2 set to %.6g\n", shift2);

	return (0);
}

/* This routine updates SHIFTAi in the current output header. */

static int UpdateShift (Hdr *hdr, double shift1, double shift2) {

/* arguments:
Hdr *hdr                   i: pointer to header struct
double shift1, shift2      i: shift to be assigned to SHIFTAi keywords
*/

	int status;

	if (fabs (shift1) < MUCH_TOO_BIG)	/* shift is OK */
	    status = Put_KeyD (hdr, "SHIFTA1", shift1, "");
	else
	    status = Put_KeyD (hdr, "SHIFTA1", 0., "");

	if (status)
	    return (status);

	if (fabs (shift2) < MUCH_TOO_BIG)	/* shift is OK */
	    status = Put_KeyD (hdr, "SHIFTA2", shift2, "");
	else
	    status = Put_KeyD (hdr, "SHIFTA2", 0., "");

	if (status)
	    return (status);

	return (0);
}
