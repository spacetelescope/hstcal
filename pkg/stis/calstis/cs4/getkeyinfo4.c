# include <stdio.h>
# include <string.h>

# include "c_iraf.h"
# include "hstio.h"

# include "stis.h"
# include "calstis4.h"
# include "hstcalerr.h"
# include "stisdef.h"

/* This routine gets keyword values from the primary header.

   Phil Hodge, 2000 Jan 5:
	Get keywords for global echelle offset.

   Phil Hodge, 2001 Feb 27:
	Set disp_type, based on opt_elem.  Remove section where cenwave
	was set to zero for non-echelle data.

   Phil Hodge, 2001 Apr 2:
	If opt_elem is PRISM, set sts->sclamp to PRISM, rather than reading
	the sclamp keyword from the header.

   Phil Hodge, 2001 May 3:
	Get cenwave as an integer, rather than as a double.
	Get moffset2 for echelle data.

   Phil Hodge, 2011 Jan 5:
	Delete moffset2.
*/

int GetKeyInfo4 (StisInfo4 *sts, Hdr *phdr) {

/* arguments:
StisInfo4 *sts  io: calibration switches and info
Hdr *phdr       i: primary header
*/

	int status;
	int nextend;			/* number of FITS extensions */
	int no_def = 0;			/* missing keyword is fatal error */
	int use_def = 1;		/* use default if keyword is missing */

	/* Get generic parameters. */

	if ((status = Get_KeyS (phdr, "ROOTNAME",
                                no_def, "", sts->rootname, STIS_CBUF)))
	    return (status);

	if ((status = Get_KeyS (phdr, "OBSMODE",
                                no_def, "", sts->obsmode, STIS_CBUF)))
	    return (status);

	/* Get the aperture name, and interpret it as to slit type. */
	if ((status = Get_KeyS (phdr, "APERTURE",
                                no_def, "", sts->aperture, STIS_CBUF)))
	    return (status);

	/* aperture field of view */
	if ((status = Get_KeyS (phdr, "APER_FOV",
                                use_def, "", sts->aper_fov, STIS_CBUF)))
	    return (status);

	/* Get the grating name, and check whether it's an echelle or prism. */
	if ((status = Get_KeyS (phdr, "OPT_ELEM",
                                no_def, "", sts->opt_elem, STIS_CBUF)))
	    return (status);
	if (strcmp (sts->opt_elem, "PRISM") == 0) {
	    sts->disp_type = PRISM_DISP;
	} else if (sts->opt_elem[0] == 'E' || sts->opt_elem[0] == 'e') {
	    sts->disp_type = ECHELLE_DISP;
	} else {
	    sts->disp_type = RECTIFIED;	/* assume first-order is rectified */
	}

	if ((status = Get_KeyS (phdr, "DETECTOR",
                                no_def, "", sts->det, STIS_CBUF)))
	    return (status);

	if (strncmp (sts->det, "NUV-MAMA", 8) == 0)
	    sts->detector = NUV_MAMA_DETECTOR;
	else if (strncmp (sts->det, "FUV-MAMA", 8) == 0)
	    sts->detector = FUV_MAMA_DETECTOR;
	else if (strncmp (sts->det, "CCD", 3) == 0)
	    sts->detector = CCD_DETECTOR;
	else
	    sts->detector = UNKNOWN_DETECTOR;

	/* Find out how many extensions there are in this file. */
	if ((status = Get_KeyI (phdr, "NEXTEND",
                                use_def, EXT_PER_GROUP, &nextend)))
	    return (status);

	/* Convert number of extensions to number of image sets. */
	sts->nimages = nextend / EXT_PER_GROUP;

	/* Central wavelength. */
	if ((status = Get_KeyI (phdr, "CENWAVE", no_def, 0, &sts->cenwave)))
	    return (status);

	/* Spectral calibration lamp name.  Setting sclamp to PRISM for
	   prism observations is a way to select the correct row in the
	   template lamp table, since the row would otherwise not be
	   unique.
	*/
	if (sts->disp_type == PRISM_DISP) {
	    strcpy (sts->sclamp, "PRISM");
	} else {
	    if ((status = Get_KeyS (phdr, "SCLAMP",
                                    no_def, "", sts->sclamp, STIS_CBUF)))
		return (status);
	}

	/* Current used for calibration lamp. */
	if ((status = Get_KeyS (phdr, "LAMPSET", no_def, "",
                                sts->lampset, STIS_CBUF)))
	    return (status);

	return (0);
}
