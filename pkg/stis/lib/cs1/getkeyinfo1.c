# include <stdio.h>
# include <stddef.h>
# include <string.h>
# include <ctype.h>		/* islower, toupper */

# include "c_iraf.h"
# include "hstio.h"

# include "stis.h"
# include "calstis1.h"
# include "hstcalerr.h"
# include "stisdef.h"

/* This routine gets keyword values from the primary header.

   Phil Hodge, 1997 Oct 20:
	Get NCOMBINE instead of CRSPLIT * NRPTEXP, as the number of images
	that have been combined (e.g. cr-reject) to make the input image.

   Phil Hodge, 1998 Mar 10:
	Error if detector is not recognized.

   Phil Hodge, 1998 May 4:
	Don't get NCOMBINE from the primary header, since it's now in
	the SCI extension header.

   Phil Hodge, 1998 June 8:
	Get OBSTYPE and APER_FOV.  Print a warning if OBSTYPE is neither
	IMAGING nor SPECTROSCOPIC.  Don't get SUBARRAY.

   Phil Hodge, 1998 July 23:
	Get TARGNAME, and set the bias_or_dark flag accordingly.

   Phil Hodge, 1998 July 30:
	Size of ROOTNAME buffer is only STIS_CBUF, not STIS_FNAME.

   Phil Hodge, 1998 Sept 24:
	Remove bias_rej.

   Phil Hodge, 1998 Oct 7:
	Change status values 1111 and 1003 to GENERIC_ERROR_CODE.
	Include GetTemperature, and call it for NUV MAMA data.

   Phil Hodge, 1998 Oct 30:
	In GetTemperature, change the suffix for the wavecal support
	file name from _spt to _wsp.

   Phil Hodge, 1999 Mar 29:
	Move GetTemperature from this file to getgrpinfo1.c.

   Phil Hodge, 2003 Jan 20:
	Get RA_TARG & DEC_TARG.
	Set sts->wavecal if targname is either waveline or wavehitm.

   Phil Hodge, 2007 May 9:
	Set sts->bias_exposure to 1 (true) if targname is 'BIAS'.
*/

int GetKeyInfo1 (StisInfo1 *sts, Hdr *phdr) {

/* arguments:
StisInfo1 *sts  io: calibration switches and info
Hdr *phdr        i: primary header
*/

	int status;

	int nextend;			/* number of FITS extensions */
	char targname[STIS_CBUF+1];	/* target name, check if BIAS or DARK */
	char crcorr[STIS_CBUF+1];       /* cosmic-ray rejection keyword */

	int use_def = 1;		/* use default if missing keyword */
	int no_default = 0;		/* missing keyword is fatal error */

	/* Get generic parameters. */

	if ((status = Get_KeyS (phdr, "ROOTNAME", no_default, "",
                                sts->rootname, STIS_CBUF)))
	    return (status);

	if ((status = Get_KeyS (phdr, "OBSMODE", use_def, "unknown",
                                sts->obsmode, STIS_CBUF)))
	    return (status);

	if ((status = Get_KeyS (phdr, "APERTURE", use_def, "",
                                sts->aperture, STIS_CBUF)))
	    return (status);

	if ((status = Get_KeyS (phdr, "DETECTOR", no_default, "",
                                sts->det, STIS_CBUF)))
	    return (status);
	if (strcmp (sts->det, "NUV-MAMA") == 0) {
	    sts->detector = NUV_MAMA_DETECTOR;
	} else if (strcmp (sts->det, "FUV-MAMA") == 0) {
	    sts->detector = FUV_MAMA_DETECTOR;
	} else if (strcmp (sts->det, "CCD") == 0) {
	    sts->detector = CCD_DETECTOR;
	} else {
	    printf ("ERROR    DETECTOR = %s is invalid\n", sts->det);
	    return (HEADER_PROBLEM);
	}

	/* Grating or mirror name. */
	if ((status = Get_KeyS (phdr, "OPT_ELEM",
                                use_def, "", sts->opt_elem, STIS_CBUF)))
	    return (status);

	/* If TARGNAME is BIAS or DARK, set the flag to indicate this. */
	if ((status = Get_KeyS (phdr, "TARGNAME",
                                use_def, "", targname, STIS_CBUF)))
	    return (status);
	if (strcmp (targname, "BIAS") == 0)
	    sts->bias_exposure = 1;
	else
	    sts->bias_exposure = 0;
	if (strcmp (targname, "BIAS") == 0 || strcmp (targname, "DARK") == 0)
	    sts->bias_or_dark = 1;
	else
	    sts->bias_or_dark = 0;
	/* If TARGNAME indicates that this is a wavecal observation, set
	   the wavecal flag.
	*/
	if (strcmp (targname, "WAVELINE") == 0 ||
	    strcmp (targname, "WAVEHITM") == 0) {
	    sts->wavecal = 1;
	}

	if ((status = Get_KeyD (phdr, "RA_TARG", no_default, 0.,
				&sts->ra_targ)))
	    return (status);
	if ((status = Get_KeyD (phdr, "DEC_TARG", no_default, 0.,
				&sts->dec_targ)))
	    return (status);

	/* Check if CRCORR is complete */
	if ((status = Get_KeyS (phdr, "CRCORR", use_def, "",
                                crcorr, STIS_CBUF)))
	    return (status);
	if (strcmp(crcorr, "COMPLETE") == 0) {
	    sts->crcorr = COMPLETE;
	} else {
	    sts->crcorr = PERFORM;
	}

	/* Find out how many extensions there are in this file. */
	if ((status = Get_KeyI (phdr, "NEXTEND",
                                use_def, EXT_PER_GROUP, &nextend)))
	    return (status);

	/* Convert number of extensions to number of SingleGroups. */
	sts->nimages = nextend / EXT_PER_GROUP;
	if (sts->nimages < 1) {
	    printf ("ERROR    NEXTEND = %d; must be at least three.\n",
		sts->nimages);
	    return (GENERIC_ERROR_CODE);
	}

	/* If obstype is spectroscopic and photcorr is set to perform,
	   we will reset photcorr.
	*/
	if ((status = Get_KeyS (phdr, "OBSTYPE",
                                use_def, "", sts->obstype, STIS_CBUF)))
	    return (status);
	if (strcmp (sts->obstype, "IMAGING") != 0 &&
	    strcmp (sts->obstype, "SPECTROSCOPIC") != 0) {
	    printf (
	"Warning  OBSTYPE = %s; should be either IMAGING or SPECTROSCOPIC\n",
			sts->obstype);
	}

	/* For CCD observations, we need the APER_FOV in order to flag
	   regions beyond the aperture for such apertures as F28X50LP.
	*/
	if ((status = Get_KeyS (phdr, "APER_FOV",
                                use_def, "", sts->aper_fov, STIS_CBUF)))
	    return (status);

	/* Get MAMA-specific parameters. */
		/* there are none, at present */

	/* Get CCD-specific parameters. */

	if (sts->detector == CCD_DETECTOR) {

	    if ((status = Get_KeyS (phdr, "CCDAMP", no_default, "",
                                    sts->ccdamp, STIS_CBUF-1)))
		return (status);
	    if (sts->ccdamp[0] != '\0' && sts->ccdamp[1] != '\0') {
		printf (
		"Warning  Multiple amp readout `%s' is not supported.\n",
		    sts->ccdamp);
	    }
	    if (islower (sts->ccdamp[0]))
		sts->ccdamp[0] = toupper (sts->ccdamp[0]);

	    if (strchr ("ABCD", sts->ccdamp[0]) == NULL) {
		printf ("ERROR    CCDAMP = `%s' is invalid.\n", sts->ccdamp);
		return (GENERIC_ERROR_CODE);
	    }

	    if ((status = Get_KeyI (phdr, "CCDGAIN",
                                    use_def, 1, &sts->ccdgain)))
		return (status);

	    if ((status = Get_KeyI (phdr, "CCDOFFST",
                                    use_def, 1, &sts->ccdoffset)))
		return (status);

	    if ((status = Get_KeyI (phdr, "BINAXIS1",
                                    use_def, 1, &sts->binaxis[0])))
		return (status);
	    if ((status = Get_KeyI (phdr, "BINAXIS2",
                                    use_def, 1, &sts->binaxis[1])))
		return (status);
	}

	return (0);
}
