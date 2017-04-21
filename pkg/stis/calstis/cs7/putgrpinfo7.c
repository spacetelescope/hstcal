# include <stdio.h>

# include "c_iraf.h"
# include "hstio.h"

# include "stis.h"
# include "calstis7.h"
# include "err.h"
# include "stisdef.h"

# define ARCSEC_TO_DEGREES    (1./3600.)

static int PutGrpHdr7 (Hdr *, CoordInfo *, int);

/* This routine copies input headers to output and updates or adds
   keywords to the extension headers.  If extver is one, the primary
   header is copied as well as the extension headers.  The keywords
   written to the extension header are the coordinate parameters and
   order number.  The CD matrix values are obtained from cdelt1 and
   cdelt2, with zero for the cross terms.  CTYPE is assumed to be set
   correctly (i.e. LAMBDA and ANGLE for spectroscopic mode) and will
   not be changed.
   The LTV keywords will be reset to their default values of zero, since
   the output can have an arbitrary border.  The LTM keywords will retain
   their input values, however, because the output scale (i.e. binning)
   should be approximately the same as the input scale.
   The spectral order number will be added to the extension headers
   with keyword SPORDER, if obstype is spectroscopic.

   Phil Hodge, 1997 Nov 13:
	Replace dispaxis with obstype in calling sequence.

   Phil Hodge, 1998 Jan 26:
	Write CUNIT1 = angstrom to the extension headers.

   Phil Hodge, 1998 July 9:
	Write SPORDER only to SCI extension.

   Phil Hodge, 1998 Oct 5:
	Change status value 1071 to HEADER_PROBLEM.

   Phil Hodge, 2000 June 30:
	Only write SPORDER for spectroscopic data.
*/

int PutGrpInfo7 (SingleGroup *in, SingleGroup *out,
	CoordInfo *coords, int extver, int obstype) {

/* arguments:
SingleGroup *in    i: input data
SingleGroup *out   o: output data
CoordInfo *coords  i: coordinate parameters
int extver         i: extension number
int obstype        i: spectroscopic or imaging
*/

	int status;

	/* Copy the primary header without changes. */
	if (extver == 1) {
	    copyHdr (out->globalhdr, in->globalhdr);
	    if (hstio_err())
		return (HEADER_PROBLEM);
	}

	/* Copy extension headers. */
	copyHdr (&out->sci.hdr, &in->sci.hdr);
	if (hstio_err())
	    return (HEADER_PROBLEM);
	copyHdr (&out->err.hdr, &in->err.hdr);
	if (hstio_err())
	    return (HEADER_PROBLEM);
	copyHdr (&out->dq.hdr, &in->dq.hdr);
	if (hstio_err())
	    return (HEADER_PROBLEM);

	/* Update keyword values. */

	if (obstype == SPECTROSCOPIC_TYPE) {
	    if ((status = Put_KeyI (&out->sci.hdr, "SPORDER", coords->sporder,
                                    "spectral order number")))
	    return (status);
	}

	if ((status = PutGrpHdr7 (&out->sci.hdr, coords, obstype)))
	    return (status);
	if ((status = PutGrpHdr7 (&out->err.hdr, coords, obstype)))
	    return (status);
	if ((status = PutGrpHdr7 (&out->dq.hdr, coords, obstype)))
	    return (status);

	return (0);
}

/* This does the work for an individual header. */

static int PutGrpHdr7 (Hdr *hdr, CoordInfo *coords, int obstype) {

/* arguments:
Hdr *hdr           o: output extension header
CoordInfo *coords  i: coordinate parameters
int obstype        i: spectroscopic or imaging
*/

	int status;

	/* Set the offset LTV to zero, but leave the scale LTM unchanged. */
	if ((status = Put_KeyD (hdr, "LTV1", 0., "")))
	    return (status);
	if ((status = Put_KeyD (hdr, "LTV2", 0., "")))
	    return (status);

	if (obstype == SPECTROSCOPIC_TYPE) {

	    /* Note:  add one to crpix to convert to one indexing. */
	    if ((status = Put_KeyD (hdr, "CRPIX1", coords->crpix[0]+1., "")))
		return (status);
	    if ((status = Put_KeyD (hdr, "CRPIX2", coords->crpix[1]+1., "")))
		return (status);
	    if ((status = Put_KeyD (hdr, "CRVAL1", coords->crval[0], "")))
		return (status);
	    if ((status = Put_KeyD (hdr, "CRVAL2", coords->crval[1], "")))
		return (status);
	    if ((status = Put_KeyD (hdr, "CD1_1", coords->cdelt[0], "")))
		return (status);
	    if ((status = Put_KeyD (hdr, "CD1_2", 0., "")))
		return (status);
	    if ((status = Put_KeyD (hdr, "CD2_1", 0., "")))
		return (status);
	    if ((status = Put_KeyD (hdr, "CD2_2",
                                    coords->cdelt[1] * ARCSEC_TO_DEGREES, "")))
		return (status);
	    if ((status = Put_KeyS (hdr, "CUNIT1", "angstrom",
                                    "units for first axis coordinates")))
		return (status);
	}

	return (0);
}
