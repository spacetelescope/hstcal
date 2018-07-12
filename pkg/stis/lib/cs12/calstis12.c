/* calstis12 -- update SHIFTAi keywords in science file extensions

   This file contains:
	CalStis12
	StisInit12
	WavOption
*/

# include <stdio.h>
# include <string.h>	/* strcpy, strcmp */

# include "hstio.h"

# include "stis.h"
# include "calstis12.h"
# include "hstcalerr.h"

static void StisInit12 (StisInfo12 *, StisInfo12 *);

/* This routine updates SHIFTAi keywords in the science file.

   Phil Hodge, 1997 Dec 11:
	Include refnames in calling sequence, and pass to GetFlags12.

   Phil Hodge, 1998 Mar 18:
	Remove RefTab apdes, don't call GetFlags12 or GetApDes12, since
	the aperture description table is not actually needed.

   Phil Hodge, 1998 Oct 5:
	Change status value 1201 to GENERIC_ERROR_CODE.

   Phil Hodge, 2000 June 16:
	Update a comment.
*/

int CalStis12 (char *inwav, char *insci,
		int w_option, int printtime, int verbose) {

/* arguments:
char *inwav            i: name of input wavecal
char *insci            i: name of science file (header modified in-place)
int w_option           i: option for selecting wavecal ("nearest" is supported)
int printtime          i: print timing info?
int verbose            i: print additional info?
*/

	int status;

	StisInfo12 wavecal, scidata;	/* calibration switches, etc. */

	IODescPtr imSci;	/* descriptor for input science file */
	IODescPtr imWav;	/* descriptor for input wavecal */
	Hdr phdrSci;		/* primary header for input science file */
	Hdr phdrWav;		/* primary header for input wavecal */

	int AddShifts (StisInfo12 *, StisInfo12 *, int);
	int GetKeyInfo12 (StisInfo12 *, Hdr *);

	PrBegin (12);

	if (printtime)
	    TimeStamp ("CALSTIS-12 started", "");

	/* Initialize structure containing calstis information. */
	StisInit12 (&wavecal, &scidata);

	/* Copy command-line arguments into wavecal & scidata. */
	strcpy (wavecal.input, inwav);
	strcpy (scidata.input, insci);
	wavecal.printtime = printtime;
	scidata.printtime = printtime;
	wavecal.verbose = verbose;
	scidata.verbose = verbose;

	PrFileName ("wavecal", wavecal.input);
	PrFileName ("science", scidata.input);

	/* Read primary header of input science file. */
	initHdr (&phdrSci);
	imSci = openInputImage (scidata.input, "", 0);
	if (hstio_err())
	    return (OPEN_FAILED);
	getHeader (imSci, &phdrSci);
	if (hstio_err())
	    return (OPEN_FAILED);
	closeImage (imSci);

	/* Get keyword values from science data primary header. */
	if ((status = GetKeyInfo12 (&scidata, &phdrSci)))
	    return (status);

	/* Print information about the science image. */
	PrHdrInfo (scidata.obsmode, scidata.aperture,
		scidata.opt_elem, scidata.det);

	freeHdr (&phdrSci);

	/* Read primary header of input wavecal. */
	initHdr (&phdrWav);
	imWav = openInputImage (wavecal.input, "", 0);
	if (hstio_err())
	    return (OPEN_FAILED);
	getHeader (imWav, &phdrWav);
	if (hstio_err())
	    return (OPEN_FAILED);
	closeImage (imWav);

	/* Get keyword values from wavecal primary header. */
	if ((status = GetKeyInfo12 (&wavecal, &phdrWav)))
	    return (status);

	freeHdr (&phdrWav);

	/* Detector, central wavelength, and grating, must be the same
	    in the science file and wavecal.
	*/
	if (wavecal.detector != scidata.detector ||
	    wavecal.cenwave != scidata.cenwave ||
	    strcmp (wavecal.opt_elem, scidata.opt_elem) != 0) {

	    printf ("ERROR    Science file and wavecal do not match.\n");
	    return (GENERIC_ERROR_CODE);
	}

	if (wavecal.printtime)
	    TimeStamp ("Begin processing", scidata.rootname);

	/* Update the target position keywords in the science file header. */
	if ((status = AddShifts (&scidata, &wavecal, w_option)))
	    return (status);

	if (scidata.printtime)
	    TimeStamp ("CALSTIS-12 completed", scidata.rootname);

	printf ("\n");
	PrEnd (12);

	return (0);
}

/* Initialize the calstis12 structure. */

static void StisInit12 (StisInfo12 *wavecal, StisInfo12 *scidata) {

	wavecal->input[0] = '\0';
	wavecal->rootname[0] = '\0';

	wavecal->det[0] = '\0';
	wavecal->obsmode[0] = '\0';
	wavecal->aperture[0] = '\0';
	wavecal->opt_elem[0] = '\0';
	wavecal->nimages = 1;
	wavecal->detector = UNKNOWN_DETECTOR;
	wavecal->cenwave = 0;

	wavecal->midpt = NULL;
	wavecal->shift1 = NULL;
	wavecal->shift2 = NULL;

	scidata->input[0] = '\0';
	scidata->rootname[0] = '\0';

	scidata->det[0] = '\0';
	scidata->obsmode[0] = '\0';
	scidata->aperture[0] = '\0';
	scidata->opt_elem[0] = '\0';
	scidata->nimages = 1;
	scidata->detector = UNKNOWN_DETECTOR;
	scidata->cenwave = 0;

	/* These four are actually not used for the science data. */
	scidata->midpt = NULL;
	scidata->shift1 = NULL;
	scidata->shift2 = NULL;
}
