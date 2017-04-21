# include <stdio.h>
# include <string.h>
# include <ctype.h>

# include "c_iraf.h"
# include "hstio.h"

# include "stis.h"
# include "calstis7.h"
# include "err.h"
# include "stisdef.h"

/* This routine gets keyword values from the primary header.

   Phil Hodge, 2000 Jan 13:
	Set sts->mama_offset to zero, rather than getting them from MOFFSETi.
	Add one to buffer size (buf).

   Phil Hodge, 2000 Aug 9:
	Remove references to mama_offset.

   Phil Hodge, 2000 Dec 18:
	If SCLAMP is NONE (the observation is not a wavecal), check whether
	PROPAPER ends in "E1", and if so, append "E1" to the value of APERTURE.

   Phil Hodge, 2001 Aug 8:
	Get FILTER, and convert to upper case (because the header values
	can be mixed case, e.g. Clear, Crystal_Quartz, Strontium_Fluoride).

   Phil Hodge, 2003 Mar 21:
	Call pseudoap to check whether PROPAPER ends in one of the pseudo-
	aperture suffixes, and if so, append that string to APERTURE.
	(See also the change made on 2000 Dec 18).
*/

int GetKeyInfo7 (StisInfo7 *sts, Hdr *phdr) {

/* arguments:
StisInfo7 *sts  io: calibration switches and info
Hdr *phdr       i: primary header
*/

	int status;
	char buf[STIS_FITS_REC+1];	/* for string parameters */
	char sclamp[STIS_CBUF+1];
	char propaper[STIS_CBUF+1];
	int nextend;			/* number of FITS extensions */
	int ccdgain;
	int i, ch;
	int no_def = 0;			/* missing keyword is fatal error */
	int use_default = 1;		/* use default if keyword is missing */

	if ((status = Get_KeyS (phdr, "ROOTNAME",
                                no_def, "", sts->rootname, STIS_CBUF)))
	    return (status);

	if ((status = Get_KeyS (phdr, "OBSMODE",
                                no_def, "", sts->obsmode, STIS_CBUF)))
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

	if ((status = Get_KeyS (phdr, "OBSTYPE", no_def, "", buf, STIS_FITS_REC)))
	    return (status);

	if (strcmp (buf, "SPECTROSCOPIC") == 0)
	    sts->obstype = SPECTROSCOPIC_TYPE;
	else if (strcmp (buf, "IMAGING") == 0)
	    sts->obstype = IMAGING_TYPE;
	else
	    sts->obstype = UNKNOWN_TYPE;

	/* Slit or aperture name, and aperture name from the proposal. */
	if ((status = Get_KeyS (phdr, "APERTURE",
                                no_def, "", sts->aperture, STIS_CBUF)))
	    return (status);

	/* Name of filter. */
	if ((status = Get_KeyS (phdr, "FILTER",
                                use_default, "CLEAR", sts->filter, STIS_CBUF)))
	    return (status);
	for (i = 0, ch = 1;  ch != 0;  i++) {
	    ch = sts->filter[i];
	    if (islower(ch))
		sts->filter[i] = toupper (ch);
	}

	/* Unless a lamp was on, check whether PROPAPER ends in one of the
	   "pseudo-aperture" suffixes, and if it does, append that string
	   to APERTURE.  (The correct way to check for a wavecal is with
	   ASN_MTYP, but that's in the extension header.)
	*/
	if ((status = Get_KeyS (phdr, "SCLAMP",
                                use_default, "NONE", sclamp, STIS_CBUF)))
	    return (status);
	if (strcmp (sclamp, "NONE") == 0) {	/* not a wavecal */
	    if ((status = Get_KeyS (phdr, "PROPAPER",
                    use_default, sts->aperture, propaper, STIS_CBUF)))
		return (status);
	    /* Append a suffix to sts->aperture, if it's found on propaper. */
	    pseudoap (propaper, sts->aperture, sts->verbose);
	}

	/* Grating name. */
	if ((status = Get_KeyS (phdr, "OPT_ELEM",
                                no_def, "", sts->opt_elem, STIS_CBUF)))
	    return (status);

	if (sts->obstype == SPECTROSCOPIC_TYPE) {

	    /* Central wavelength. */
	    if ((status = Get_KeyI (phdr, "CENWAVE", no_def, 0, &sts->cenwave)))
		return (status);

	    /* Target location. */
	    if ((status = Get_KeyD (phdr, "RA_TARG", no_def, 0., &sts->ra_targ)))
		return (status);
	    if ((status = Get_KeyD (phdr, "DEC_TARG", no_def, 0., &sts->dec_targ)))
		return (status);
	}

	/* Find out how many extensions there are in this file. */
	if ((status = Get_KeyI (phdr, "NEXTEND",
                                use_default, EXT_PER_GROUP, &nextend)))
	    return (status);

	/* Convert number of extensions to number of SingleGroups. */
	sts->nimages = nextend / EXT_PER_GROUP;

	if (sts->detector == CCD_DETECTOR) {

	    /* Get CCD-specific info. */
	    if ((status = Get_KeyD (phdr, "ATODGAIN",
                                    use_default, -1., &sts->atodgain)))
		return (status);
	    if (sts->atodgain < 0.) {	/* keyword wasn't found in header */
		if ((status = Get_KeyI (phdr, "CCDGAIN",
                                        use_default, 1, &ccdgain)))
		    return (status);
		sts->atodgain = (double)ccdgain;
	    }
	}

	return (0);
}
