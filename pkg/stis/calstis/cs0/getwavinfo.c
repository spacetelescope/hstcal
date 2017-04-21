# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include "c_iraf.h"
# include "hstio.h"

# include "stis.h"
# include "calstis0.h"
# include "err.h"
# include "stisshutter.h"	/* for EXT_SHUTTER_CLOSED */
# include "stisdef.h"

/* This routine gets switches and reference file names from the primary
   header of the wavecal file.

   Phil Hodge, 1997 Dec 10:
	Replace missing with wavref in calling sequences of this routine
	and WavFlags; get binning and gain info, and set samebin.

   Phil Hodge, 1998 Mar 10:
	Error if detector is not recognized.

   Phil Hodge, 1998 May 1:
	Call CheckVolt.

   Phil Hodge, 1998 July 30:
	Get TEXPSTRT and compare with EXT_SHUTTER_CLOSED, and only set
	wav_subscicorr to perform for earlier dates.

   Phil Hodge, 2000 Jan 5:
	Get opt_elem in order to set echelle flag.

   Phil Hodge, 2001 Feb 22:
	Set prism flag, based on opt_elem.

   Phil Hodge, 2001 Apr 30:
	Move the tests on opt_elem for echelle and prism from this
	function to GetSciInfo.

   Phil Hodge, 2001 Aug 27:
	Remove the call to CheckVolt; move this to cs1/Do2D.
*/

int GetWavInfo (StisInfo *sts, CalSwitch *wav_sw, RefFileInfo *wavref) {

/* arguments:
StisInfo *sts         i: calibration flags and other info
CalSwitch *wav_sw     o: all calibration switches (0 or 1) for wavecal
RefFileInfo *sciref  io: list of keyword,filename pairs for science file
*/

	int status;

	IODescPtr im;		/* descriptor for an image */
	Hdr phdr;		/* primary header */
	char *buf;		/* for string parameters */
	double texpstrt;	/* start time (MJD) of first exposure */
	int nextend;		/* number of FITS extensions in wavecal */
	int mismatch;		/* flag for comparing science and wavecal */
	int no_def = 0;		/* missing keyword is fatal error */
	int use_default = 1;	/* use default if keyword is missing */

	int GetFlags (CalSwitch *, Hdr *);
	int WavFlags (StisInfo *, CalSwitch *, Hdr *, RefFileInfo *);

	/* Read primary header of wavfile into phdr. */
	initHdr (&phdr);
	im = openInputImage (sts->wavfile, "", 0);
	if (hstio_err())
	    return (OPEN_FAILED);
	getHeader (im, &phdr);		/* get primary header */
	if (hstio_err())
	    return (OPEN_FAILED);
	closeImage (im);

	if ((buf = calloc (STIS_FITS_REC+1, sizeof(char))) == NULL)
	    return (OUT_OF_MEMORY);

	/* Get generic parameters and compare with science file. */

	/* Get detector name and compare with the science file.  It is
	   an error if they do not match.
	*/
	if ((status = Get_KeyS (&phdr, "DETECTOR",
				no_def, "", buf, STIS_FITS_REC)))
	    return (status);
	mismatch = 0;
	if (strcmp (buf, "NUV-MAMA") == 0) {
	    if (sts->detector != NUV_MAMA_DETECTOR)
		mismatch = 1;
	} else if (strcmp (buf, "FUV-MAMA") == 0) {
	    if (sts->detector != FUV_MAMA_DETECTOR)
		mismatch = 1;
	} else if (strcmp (buf, "CCD") == 0) {
	    if (sts->detector != CCD_DETECTOR)
		mismatch = 1;
	} else {
	    printf ("ERROR    (in wavecal) DETECTOR = %s is invalid\n", buf);
	    return (HEADER_PROBLEM);
	}
	if (mismatch) {
	    printf (
	"ERROR    DETECTOR mismatch between science and wavecal.\n");
	    return (HEADER_PROBLEM);
	}

	/* Compare OBSTYPE in science and wavecal. */
	if ((status = Get_KeyS (&phdr, "OBSTYPE",
				no_def, "", buf, STIS_FITS_REC)))
	    return (status);
	mismatch = 0;
	if (strcmp (buf, "SPECTROSCOPIC") == 0) {
	    if (sts->obstype != SPECTROSCOPIC_TYPE)
		mismatch = 1;
	} else if (strcmp (buf, "IMAGING") == 0) {
	    if (sts->obstype != IMAGING_TYPE)
		mismatch = 1;
	} else {
	    printf ("Warning  Unknown OBSTYPE = '%s'\n", buf);
	}
	if (mismatch)
	    printf (
	"Warning  OBSTYPE mismatch between science and wavecal.\n");

	/* Was HITM used for wavecal?  Note that we're setting a switch
	   here (to run calstis11 or not), but it's based on SCLAMP and
	   DETECTOR, not gotten directly from a header switch.
	*/
	if (sts->detector == CCD_DETECTOR) {
	    if ((status = Get_KeyS (&phdr, "SCLAMP",
                                    no_def, "", buf, STIS_FITS_REC)))
		return (status);
	    if (strcmp (buf, "HITM1") == 0 || strcmp (buf, "HITM2") == 0) {
		if ((status = Get_KeyD (&phdr, "TEXPSTRT",
                                        use_default, 0., &texpstrt)))
		    return (status);
		/* Later than this date the shutter was closed, so we don't
		   need to subtract the science image from the wavecal.
		*/
		if (texpstrt < EXT_SHUTTER_CLOSED)
		    sts->wav_subscicorr = PERFORM;
	    }
	}

	/* Check that the number of extensions is a multiple of three. */
	if ((status = Get_KeyI (&phdr, "NEXTEND",
				use_default, EXT_PER_GROUP, &nextend)))
	    return (status);

	/* Get binning and gain info.  We really only need this for the CCD. */
	if ((status = Get_KeyI (&phdr, "BINAXIS1",
				use_default, 1, &sts->wavbin[0])))
	    return (status);
	if ((status = Get_KeyI (&phdr, "BINAXIS2",
				use_default, 1, &sts->wavbin[1])))
	    return (status);
	if ((status = Get_KeyI (&phdr, "CCDGAIN", use_default, 1, &sts->wavgain)))
	    return (status);
	sts->samebin = (sts->scibin[0] == sts->wavbin[0] &&
			sts->scibin[1] == sts->wavbin[1] &&
			sts->scigain == sts->wavgain);

	/* Get calibration switches. */
	if ((status = GetFlags (wav_sw, &phdr)))
	    return (status);

	/* Get reference file names. */
	if ((status = WavFlags (sts, wav_sw, &phdr, wavref)))
	    return (status);

	freeHdr (&phdr);
	free (buf);
	return (0);
}
