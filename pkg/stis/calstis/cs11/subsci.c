# include <stdio.h>
# include <stdlib.h>	/* malloc */

# include "c_iraf.h"
# include "hstio.h"

# include "stis.h"
# include "calstis11.h"
# include "err.h"
# include "stisdef.h"

/* This routine opens the input wavecal and science files, selects the
   best (closest in time) science image for each wavecal, and calls a
   routine to scale and subtract the science image from the wavecal.
   The difference is written to the output file.

   Phil Hodge, 1998 Oct 5:
	Change status value 1001 to GENERIC_ERROR_CODE.

   Phil Hodge, 1999 Oct 1:
	Include the ratio of the gains in ratio; previously it was only
	the ratio of exposure times.

   Phil Hodge, 2000 Jan 13:
	Add one to history buffer size.

   Phil Hodge, 2008 Nov 3:
	Skip science image sets for which IMSET_OK is false.
*/

int SubSci (StisInfo11 *wavecal, StisInfo11 *scidata) {

/* arguments:
StisInfo11 *wavecal  i: info for wavecal
StisInfo11 *scidata  i: info for science data
*/

	int status;

	IODescPtr im;		/* descriptor for an image */
	Hdr hdr;		/* header for SCI extension */
	SingleGroup wav, sci;
	char history[STIS_LINE+1];
	double ratio;		/* ratio of exposure times and gains */
	double exptime, midpt;	/* temporary variables */
	int imset_ok;		/* value of header keyword IMSET_OK */
	int n;
	int k;			/* loop index for number of good sci images */
	int extverWav;		/* imset number for wavecal */
	int extverSci;		/* imset number for science data */
	int index;		/* extverSci-1 (index into arrays) */
	int option;		/* selection option for best science imset */

	int BinSubtract (SingleGroup *, SingleGroup *, double, int);
	int GetTimes11 (Hdr *, double *, double *, int *);
	int MatchSci (StisInfo11 *, StisInfo11 *, int, int *);

	/* Allocate space for imset-specific info. */

	/* We only need one element for wavecals. */
	if ((wavecal->exptime = malloc (sizeof (double))) == NULL)
	    return (OUT_OF_MEMORY);
	if ((wavecal->midpt = malloc (sizeof (double))) == NULL)
	    return (OUT_OF_MEMORY);

	n = scidata->nimages;
	if ((scidata->exptime = malloc (n * sizeof (double))) == NULL)
	    return (OUT_OF_MEMORY);
	if ((scidata->midpt = malloc (n * sizeof (double))) == NULL)
	    return (OUT_OF_MEMORY);

	/* Get info from each extension in the science file; save
	   exptime and midpt info in scidata.  But skip an imset if
	   keyword IMSET_OK is false.
	*/
	initHdr (&hdr);
	k = 0;
	for (extverSci = 1;  extverSci <= scidata->nimages;  extverSci++) {
	    im = openInputImage (scidata->input, "SCI", extverSci);
	    if (hstio_err())
		return (OPEN_FAILED);
	    getHeader (im, &hdr);
	    if (hstio_err())
		return (OPEN_FAILED);
	    closeImage (im);
	    if ((status = GetTimes11 (&hdr, &exptime, &midpt, &imset_ok)) != 0)
		return (status);
	    if (imset_ok) {
		scidata->exptime[k] = exptime;
		scidata->midpt[k] = midpt;
		k++;
	    }
	    freeHdr (&hdr);
	}
	scidata->nimages = k;	/* update to be the number with imset_ok */

	initSingleGroup (&wav);
	initSingleGroup (&sci);

	/* Process each wavecal. */

	for (extverWav = 1;  extverWav <= wavecal->nimages;  extverWav++) {

	    PrGrpBegin ("imset", extverWav);

	    getSingleGroup (wavecal->input, extverWav, &wav);
	    if (hstio_err())
		return (OPEN_FAILED);

	    /* Get keyword values from wavecal header.  Note that we put
		exptime and midpt in the first array element, since there
		is only one element allocated.
	    */
	    if ((status = GetTimes11 (&wav.sci.hdr,
			wavecal->exptime, wavecal->midpt, &imset_ok)) != 0)
		return (status);
	    if (!imset_ok) {
		freeSingleGroup (&wav);
		continue;
	    }

	    /* Find the appropriate imset in science file. */
	    option = STIS_NEAREST;
	    if ((status = MatchSci (wavecal, scidata, option, &extverSci)))
		return (status);
	    index = extverSci - 1;		/* array index for exptime */

	    getSingleGroup (scidata->input, extverSci, &sci);
	    if (hstio_err())
		return (OPEN_FAILED);

	    /* There's an additional factor of the binning, but that will be
		taken into account in BinSubtract.
	    */
	    ratio = wavecal->exptime[0] / scidata->exptime[index] *
			scidata->gain / wavecal->gain;

	    /* Do the subtraction (in-place in wav). */
	    if ((status = BinSubtract (&wav, &sci, ratio, scidata->verbose)))
		return (status);

	    freeSingleGroup (&sci);

	    if (extverWav == 1) {
		sprintf (history, "Science file %s was subtracted.\n",
			scidata->input);
		addHistoryKw (wav.globalhdr, history);
		if (hstio_err())
		    return (HEADER_PROBLEM);
	    }

	    printf ("         %s[EXTVER=%d] was subtracted.\n",
			scidata->input, extverSci);

	    sprintf (history, "%s[EXTVER=%d] was subtracted\n",
			scidata->input, extverSci);
	    addHistoryKw (&wav.sci.hdr, history);
	    if (hstio_err())
		return (HEADER_PROBLEM);

	    UCalVer (wav.globalhdr);			/* update CAL_VER */
	    UFilename (wavecal->output, wav.globalhdr);	/* update FILENAME */
	    putSingleGroup (wavecal->output, extverWav, &wav, 0);
	    if (hstio_err()) {
		printf ("ERROR    Couldn't write imset %d of %s.\n",
			extverWav, wavecal->output);
		return (GENERIC_ERROR_CODE);
	    }

	    freeSingleGroup (&wav);

	    PrGrpEnd ("imset", extverWav);
	}

	free (wavecal->exptime);
	free (wavecal->midpt);
	free (scidata->exptime);
	free (scidata->midpt);

	return (0);
}
