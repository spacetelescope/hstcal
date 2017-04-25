/* calstis11 -- subtract science from wavecal

   This file contains:
	CalStis11
	StisInit11
	PrintInfo11
*/

# include <stdio.h>
# include <string.h>

# include "hstio.h"

# include "stis.h"
# include "calstis11.h"
# include "hstcalerr.h"
# include "stisshutter.h"	/* for EXT_SHUTTER_CLOSED */

static void StisInit11 (StisInfo11 *, StisInfo11 *);

/* This routine subtracts the science image from the wavecal, if the
   external shutter was open during the wavecal and the HITM system
   was used.

   Phil Hodge, 1998 Aug 13:
	Add a check on the date (TEXPSTRT) of the observation.  Prior to
	a certain date, the shutter was open for CCD HITM wavecals, but
	afterwards the shutter was closed for all wavecals.  If the science
	image should not be subtracted (either for this reason, or because
	the detector is not CCD, or the lamp was not a HITM), then this
	function will return immediately with status = NOTHING_TO_DO.
*/

int CalStis11 (char *inwav, char *insci, char *output,
		int printtime, int verbose) {

	int status;

	StisInfo11 wavecal, scidata;	/* calibration switches, etc. */

	IODescPtr imWav;	/* descriptor for input wavecal */
	IODescPtr imSci;	/* descriptor for input science file */
	Hdr phdrWav;		/* primary header for input wavecal */
	Hdr phdrSci;		/* primary header for input science file */
	int subscicorr;		/* PERFORM if CCD and sclamp is HITM1 or 2 */

	int GetKeyInfo11 (StisInfo11 *, Hdr *);
	int SubSci (StisInfo11 *, StisInfo11 *);

	PrBegin (11);

	if (printtime)
	    TimeStamp ("CALSTIS-11 started", "");

	/* Initialize structure containing calstis information. */
	StisInit11 (&wavecal, &scidata);

	/* Copy command-line arguments into wavecal & scidata. */
	strcpy (wavecal.input, inwav);
	strcpy (scidata.input, insci);
	strcpy (wavecal.output, output);
	wavecal.printtime = printtime;
	scidata.printtime = printtime;
	wavecal.verbose = verbose;
	scidata.verbose = verbose;

	PrFileName ("wavecal", wavecal.input);
	PrFileName ("science", scidata.input);
	PrFileName ("output", wavecal.output);

	initHdr (&phdrWav);
	initHdr (&phdrSci);

	/* Check whether the output file already exists. */
	if ((status = FileExists (wavecal.output)))
	    return (status);

	/* Read primary header of input wavecal. */
	imWav = openInputImage (wavecal.input, "", 0);
	if (hstio_err())
	    return (OPEN_FAILED);
	getHeader (imWav, &phdrWav);
	if (hstio_err())
	    return (OPEN_FAILED);
	closeImage (imWav);

	/* Get keyword values from wavecal primary header. */
	if ((status = GetKeyInfo11 (&wavecal, &phdrWav)))
	    return (status);

	freeHdr (&phdrWav);

	/* Print information about the input wavecal. */
	PrHdrInfo (wavecal.obsmode, wavecal.aperture,
		wavecal.opt_elem, wavecal.det);

	/* Do we need to subtract the science image from the wavecal? */
	subscicorr = PERFORM;			/* initial value */
	if (wavecal.detector != CCD_DETECTOR) {
	    subscicorr = OMIT;
	    printf ("Warning  Detector is %s\n", wavecal.det);
	}
	if (strcmp (wavecal.sclamp, "HITM1") != 0 &&
		   strcmp (wavecal.sclamp, "HITM2") != 0) {
	    subscicorr = OMIT;
	    printf ("Warning  Wavecal SCLAMP is `%s', not HITM1 or HITM2\n",
			wavecal.sclamp);
	}
	if (wavecal.texpstrt >= EXT_SHUTTER_CLOSED) {
	    subscicorr = OMIT;
	    printf (
	"Warning  TEXPSTRT=%.2f implies external shutter is closed.\n",
		wavecal.texpstrt);
	}

	if (subscicorr != PERFORM) {
	    printf (
	"Warning  Science data will not be subtracted from wavecal.\n");
	    return (NOTHING_TO_DO);
	}

	/* Read primary header of input science file. */
	imSci = openInputImage (scidata.input, "", 0);
	if (hstio_err())
	    return (OPEN_FAILED);
	getHeader (imSci, &phdrSci);
	if (hstio_err())
	    return (OPEN_FAILED);
	closeImage (imSci);

	if (wavecal.printtime)
	    TimeStamp ("Begin processing", wavecal.rootname);

	/* Get keyword values from science file primary header. */
	if ((status = GetKeyInfo11 (&scidata, &phdrSci)))
	    return (status);

	freeHdr (&phdrSci);

	/* Detector, central wavelength, grating, and aperture must be
	   the same in the wavecal and science file.
	*/
	if (wavecal.detector != scidata.detector ||
	    wavecal.cenwave != scidata.cenwave ||
	    strcmp (wavecal.opt_elem, scidata.opt_elem) != 0 ||
	    strcmp (wavecal.aperture, scidata.aperture) != 0) {

	    printf ("Warning  Wavecal and science file do not match; \\\n");
	    printf ("Warning  the science file will not be subtracted.\n");
	    return (NOTHING_TO_DO);
	}

	/* Subtract the science image from the wavecal. */
	if ((status = SubSci (&wavecal, &scidata)))
	    return (status);

	printf ("\n");
	PrEnd (11);

	if (wavecal.printtime)
	    TimeStamp ("CALSTIS-11 completed", wavecal.rootname);

	return (0);
}

/* Initialize the calstis11 structure. */

static void StisInit11 (StisInfo11 *wavecal, StisInfo11 *scidata) {

	wavecal->input[0] = '\0';
	wavecal->output[0] = '\0';
	wavecal->rootname[0] = '\0';

	wavecal->sclamp[0] = '\0';
	wavecal->obsmode[0] = '\0';
	wavecal->aperture[0] = '\0';
	wavecal->opt_elem[0] = '\0';
	wavecal->det[0] = '\0';
	wavecal->detector = UNKNOWN_DETECTOR;
	wavecal->nimages = 1;
	wavecal->cenwave = 0;

	wavecal->exptime = NULL;
	wavecal->midpt = NULL;

	/* some of the following will not be used at all */

	scidata->input[0] = '\0';
	scidata->output[0] = '\0';
	scidata->rootname[0] = '\0';

	scidata->sclamp[0] = '\0';
	scidata->obsmode[0] = '\0';
	scidata->aperture[0] = '\0';
	scidata->opt_elem[0] = '\0';
	scidata->det[0] = '\0';
	scidata->detector = UNKNOWN_DETECTOR;
	scidata->nimages = 1;
	scidata->cenwave = 0;

	scidata->exptime = NULL;
	scidata->midpt = NULL;
}
