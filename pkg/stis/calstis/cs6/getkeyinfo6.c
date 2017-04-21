# include <stdio.h>
# include <string.h>

# include "xtables.h"
# include "hstio.h"

# include "stis.h"
# include "stisdef.h"
# include "calstis6.h"
# include "err.h"

/*
   Get keyword values from the primary header.




   Revision history:
   ----------------
   20 Feb 97  -  Adapted from similar routine in calstis7 (I.Busko)
   08 May 97  -  Conform to new _trl standard (IB)
   09 May 97  -  Added OBSMODE keyword input (IB)
   03 Sep 97  -  Add ATODGAIN correction (IB)
   14 Nov 97  -  Add POSTARG (IB)
   08 Dec 97  -  Remove POSTARG (IB)
   18 Sep 98  -  Read CCD-related keywords (IB)
   18 Dec 00  -  Append "E1" to aperture name, if appropriate (PEH)
   21 Mar 03  -  Append any pseudo-aperture to aperture name (PEH)
   17 Jun 03  -  Adde CCDOFFST, BINAXIS1, BINAXIS2 keywords (PB)

*/

int GetKeyInfo6 (StisInfo6 *sts, Hdr *phdr) {

/* arguments:
StisInfo6 *sts  io: calibration switches and info
Hdr *phdr       i: primary header
*/

	int status;
	int nextend;			/* number of FITS extensions */
	char sclamp[STIS_CBUF+1];
	char propaper[STIS_CBUF+1];
	int no_def = 0;			/* missing keyword is fatal error */
	int use_default = 1;		/* use default if keyword is missing */

	/* Get generic parameters. */

	if ((status = Get_KeyS (phdr, "ROOTNAME", no_def, "", sts->rootname, STIS_CBUF)))
	    return (status);

	if ((status = Get_KeyS (phdr, "DETECTOR", no_def, "", sts->det, STIS_CBUF)))
	    return (status);
	if (strcmp (sts->det, "NUV-MAMA") == 0)
	    sts->detector = NUV_MAMA_DETECTOR;
	else if (strcmp (sts->det, "FUV-MAMA") == 0)
	    sts->detector = FUV_MAMA_DETECTOR;
	else if (strcmp (sts->det, "CCD") == 0)
	    sts->detector = CCD_DETECTOR;
	else
	    sts->detector = UNKNOWN_DETECTOR;

	/* Only spectra are meaningful for calstis6. */
	if ((status = Get_KeyS (phdr, "OBSTYPE",
                                no_def, "", sts->obsmode, STIS_FITS_REC)))
	    return (status);
	if (strcmp (sts->obsmode, "SPECTROSCOPIC") != 0)
	    return (NOTHING_TO_DO);

	/* This is required for _trl reporting only. */
	if ((status = Get_KeyS (phdr, "OBSMODE",
                                no_def, "", sts->obsmode, STIS_FITS_REC)))
	    return (status);

	/* Find out how many extensions there are in this file. */
	if ((status = Get_KeyI (phdr, "NEXTEND",
                                use_default, EXT_PER_GROUP, &nextend)))
	    return (status);

	/* Convert number of extensions to number of SingleGroups. */
	sts->nimages = nextend / EXT_PER_GROUP;
	if (nextend != sts->nimages * EXT_PER_GROUP) {
	    printf (
            "ERROR    NEXTEND must be a multiple of %d\n", EXT_PER_GROUP);
	    return (HEADER_PROBLEM);
	}

	/* Target location. */
	if ((status = Get_KeyD (phdr, "RA_TARG", no_def, 0., &sts->ra_targ)))
	    return (status);
	if ((status = Get_KeyD (phdr, "DEC_TARG", no_def, 0., &sts->dec_targ)))
	    return (status);

	/* POSTARG */
/*	if ((status = Get_KeyD (phdr, "POSTARG1", no_def, 0., &sts->pos_targ[0])))
	    return (status);
	if ((status = Get_KeyD (phdr, "POSTARG2", no_def, 0., &sts->pos_targ[1])))
	    return (status);
*/
	/* Get CCD-specific info. */
	if ((status = Get_KeyS (phdr, "CCDAMP",
                                use_default, "", sts->ccdamp, STIS_CBUF)))
	    return (status);
	if ((status = Get_KeyD (phdr, "READNSE",
                                use_default, 0., &sts->readnse)))
	    return (status);
	if ((status = Get_KeyD (phdr, "ATODGAIN",
                                use_default, 1., &sts->atodgain)))
	    return (status);
 	if ((status = Get_KeyI (phdr, "CCDGAIN", use_default, 1, &sts->ccdgain)))
	    return (status);
 	if ((status = Get_KeyI (phdr, "CCDOFFST", use_default, 1,
                                &sts->ccdoffset)))
	    return (status);
 	if ((status = Get_KeyI (phdr, "BINAXIS1", use_default, 1,
                                &(sts->binaxis[0]))))
	    return (status);
 	if ((status = Get_KeyI (phdr, "BINAXIS2", use_default, 1,
                                &(sts->binaxis[1]))))
	    return (status);
 	if ((status = Get_KeyI (phdr, "CRSPLIT", use_default, 1, &sts->crsplit)))
	    return (status);

	/* Get MAMA "dither", i.e. the deliberate offset to prevent
	   exposures from always being at the same place.
	*/
	/*
	if ((status = Get_KeyD (phdr, "MOFFSET1", use_default, 0.,
                     &sts->dither_offset[0])))
	    return (status);
	if ((status = Get_KeyD (phdr, "MOFFSET2", use_default, 0.,
                    &sts->dither_offset[1])))
	    return (status);
	*/

	/* Zeroed for good... */
        sts->dither_offset[0] = 0.;
        sts->dither_offset[1] = 0.;

	/* Get calstis6-specific parameters. */

	/* Slit name. */
	if ((status = Get_KeyS (phdr, "APERTURE",
                                no_def, "", sts->aperture, STIS_CBUF)))
	    return (status);

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

	/* Central wavelength. */
	if ((status = Get_KeyI (phdr, "CENWAVE", no_def, 0., &sts->cenwave)))
	    return (status);

	return (0);
}
