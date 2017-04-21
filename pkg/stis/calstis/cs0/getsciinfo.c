# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include "c_iraf.h"
# include "hstio.h"

# include "stis.h"
# include "calstis0.h"
# include "err.h"
# include "stisdef.h"

/* This routine gets info from the primary header of the science file,
   and it calls routines to get calibration switches and reference file
   names.

   Phil Hodge, 1997 Dec 10:
	Replace missing with sciref in calling sequences of this routine
	and SciFlags; get binning and gain info.

   Phil Hodge, 1998 Mar 10:
	Error if detector is not recognized.

   Phil Hodge, 1998 May 1:
	Call CheckVolt.

   Phil Hodge, 1998 July 21:
	Get wavecal name before calling SciFlags; if wavecal name is "N/A",
	reset wavecorr to SKIPPED.

   Phil Hodge, 1999 Sept 23:
	If TEXPTIME is less than or equal to zero, and we're not processing
	a BIAS image, return NOTHING_TO_DO.

   Phil Hodge, 2001 Apr 30:
	Move the tests on opt_elem to set the echelle and prism flags
	from GetWavInfo to this function.

   Phil Hodge, 2001 Aug 27:
	Remove the call to CheckVolt; move this to cs1/Do2D.
*/

int GetSciInfo (StisInfo *sts, CalSwitch *sci_sw, RefFileInfo *sciref) {

/* arguments:
StisInfo *sts         i: calibration flags and other info
CalSwitch *sci_sw     o: all calibration switches (0 or 1) for science file
RefFileInfo *sciref  io: list of keyword,filename pairs for science file
*/

	int status;

	IODescPtr im;		/* descriptor for an image */
	Hdr phdr;		/* primary header */
	char *buf;		/* for string parameters */
	double texptime;	/* total exposure time */
	int nextend;		/* number of FITS extensions in rawfile */
	int no_def = 0;		/* missing keyword is fatal error */
	int use_default = 1;	/* use default if keyword is missing */

	int GetFlags (CalSwitch *, Hdr *);
	int SciFlags (StisInfo *, CalSwitch *, Hdr *, RefFileInfo *);

	/* Read primary header of rawfile into phdr. */
	initHdr (&phdr);
	im = openInputImage (sts->rawfile, "", 0);
	if (hstio_err())
	    return (OPEN_FAILED);
	getHeader (im, &phdr);		/* get primary header */
	if (hstio_err())
	    return (OPEN_FAILED);
	closeImage (im);

	if ((buf = calloc (STIS_FITS_REC+1, sizeof(char))) == NULL)
	    return (OUT_OF_MEMORY);

	/* Get generic parameters. */

	if ((status = Get_KeyS (&phdr, "ROOTNAME",
				no_def, "", sts->rootname, STIS_CBUF)))
	    return (status);

	if ((status = Get_KeyS (&phdr, "DETECTOR",
				no_def, "", buf, STIS_FITS_REC)))
	    return (status);

	if (strcmp (buf, "NUV-MAMA") == 0) {
	    sts->detector = NUV_MAMA_DETECTOR;
	} else if (strcmp (buf, "FUV-MAMA") == 0) {
	    sts->detector = FUV_MAMA_DETECTOR;
	} else if (strcmp (buf, "CCD") == 0) {
	    sts->detector = CCD_DETECTOR;
	} else {
	    printf ("ERROR    DETECTOR = %s is invalid\n", buf);
	    return (HEADER_PROBLEM);
	}

	/* Quit now if the total exposure time is zero, unless
	   this is a bias observation.
	*/
	if ((status = Get_KeyD (&phdr, "TEXPTIME", use_default, 0., &texptime)))
	    return (status);
	if (texptime <= 0.) {
	    strcpy (buf, "not a bias obs");
	    if (sts->detector == CCD_DETECTOR) {
		if ((status = Get_KeyS (&phdr, "TARGNAME",
                                        use_default, "", buf, STIS_FITS_REC))) {
		    return (status);
		}
	    }
	    if (strcmp (buf, "BIAS") != 0) {
		printf ("Warning  Total exposure time = %.6g\n", texptime);
		free (buf);
		freeHdr (&phdr);
		return (NOTHING_TO_DO);
	    }
	}

	if ((status = Get_KeyS (&phdr, "OBSTYPE",
				no_def, "", buf, STIS_FITS_REC)))
	    return (status);

	if (strcmp (buf, "SPECTROSCOPIC") == 0)
	    sts->obstype = SPECTROSCOPIC_TYPE;
	else if (strcmp (buf, "IMAGING") == 0)
	    sts->obstype = IMAGING_TYPE;
	else
	    printf ("Warning  Unknown OBSTYPE = '%s'\n", buf);

	/* Get opt_elem in order to check for echelle or prism. */
	if ((status = Get_KeyS (&phdr, "OPT_ELEM",
				use_default, "", buf, STIS_FITS_REC)))
	    return (status);
	sts->echelle = (buf[0] == 'E' || buf[0] == 'e');
	if (strcmp (buf, "PRISM") == 0)
	    sts->prism = 1;
	else
	    sts->prism = 0;

	/* Find out how many extensions there are in this file. */
	if ((status = Get_KeyI (&phdr, "NEXTEND",
				use_default, EXT_PER_GROUP, &nextend)))
	    return (status);
	sts->nimages = nextend / EXT_PER_GROUP;

	/* Get binning and gain info.  We really only need this for the CCD. */
	if ((status = Get_KeyI (&phdr, "BINAXIS1",
				use_default, 1, &sts->scibin[0])))
	    return (status);
	if ((status = Get_KeyI (&phdr, "BINAXIS2",
				use_default, 1, &sts->scibin[1])))
	    return (status);
	if ((status = Get_KeyI (&phdr, "CCDGAIN", use_default, 1, &sts->scigain)))
	    return (status);
	sts->samebin = 1;	/* default */

	/* Get calibration switches. */
	if ((status = GetFlags (sci_sw, &phdr)))
	    return (status);

	if (sci_sw->wavecorr == PERFORM) {
	    /* Was wavecal name not specified on command line? */
	    if (sts->wavfile[0] == '\0') {
		if ((status = Get_KeyS (&phdr, "WAVECAL", no_def, "",
                                        sts->wavfile, STIS_LINE)))
		    return (status);
	    }
	    /* Reset wavecorr if wavecal name is "N/A". */
	    if (strcmp (sts->wavfile, "N/A") == 0) {
		printf (
		"Warning  WAVECORR will be skipped because WAVECAL is N/A.\n");
		sci_sw->wavecorr = SKIPPED;
	    }
	}

	/* Check that reference files exist. */
	if ((status = SciFlags (sts, sci_sw, &phdr, sciref)))
	    return (status);

	freeHdr (&phdr);
	free (buf);
	return (0);
}
