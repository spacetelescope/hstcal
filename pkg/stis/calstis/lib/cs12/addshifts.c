# include <stdio.h>
# include <stdlib.h>	/* malloc */
# include <math.h>	/* fabs */
# include "hstio.h"
# include "stis.h"
# include "calstis12.h"
# include "cs12.h"		/* for interpolation options */
# include "hstcalerr.h"

/* This routine gets the shift for each axis, by interpolation
   within the array of wavecals, and copies those values to the
   SHIFTA1 & SHIFTA2 keywords for each extension header.

   Phil Hodge, 1998 Mar 18:
	Remove the section that added the offsets between the wavecal
	and science apertures to the shifts.

   Phil Hodge, 2000 Jan 19:
	Don't call AvgShift, since echelle data no longer need to be
	averaged; check that wavecals are time ordered; get the shifts
	in MatchWav instead of AvgShift.

   Phil Hodge, 2000 June 16:
	Change comments to remove references to CRPIX and to HISTORY
	in the extension headers.
*/

int AddShifts (StisInfo12 *scidata, StisInfo12 *wavecal, int w_option) {

/* arguments:
StisInfo12 *scidata  io: info for science data
StisInfo12 *wavecal  i: info for wavecal
int w_option         i: option for selecting wavecal
*/

	int status;

	IODescPtr im;		/* descriptor for an image */
	Hdr phdr;		/* primary header */
	Hdr hdr;		/* header for an extension */
	double exptime, midpt;	/* values from science data SCI header */
	double shift1, shift2;	/* average shift in each axis */
	int n;
	int extverSci;		/* imset number for science data */
	int extverWav;		/* imset number for wavecal */
	/* these are the good shifts and their times */
	double *wav_midpt1, *wav_shift1, *wav_midpt2, *wav_shift2;
	int n1, n2;		/* number of good shift1 and shift2 */

	int GetWavGrp (StisInfo12 *, Hdr *, int);
	int GetSciGrp (Hdr *, double *, double *);
	int History12 (Hdr *, char *);
	int MatchWav (double *, double *, int, double, int, double *);
	int TargPos (StisInfo12 *, int, double, double);

	/* Allocate space for imset-specific info. */
	n = wavecal->nimages;
	if ((wavecal->midpt = malloc (n * sizeof (double))) == NULL)
	    return (OUT_OF_MEMORY);
	if ((wavecal->shift1 = malloc (n * sizeof (double))) == NULL)
	    return (OUT_OF_MEMORY);
	if ((wavecal->shift2 = malloc (n * sizeof (double))) == NULL)
	    return (OUT_OF_MEMORY);

	/* Get and save info from each SCI extension in the wavecal file. */
	initHdr (&hdr);
	for (extverWav = 1;  extverWav <= wavecal->nimages;  extverWav++) {
	    im = openInputImage (wavecal->input, "SCI", extverWav);
	    if (hstio_err())
		return (OPEN_FAILED);
	    getHeader (im, &hdr);
	    if (hstio_err())
		return (OPEN_FAILED);
	    closeImage (im);
	    if ((status = GetWavGrp (wavecal, &hdr, extverWav)))
		return (status);
	    freeHdr (&hdr);
	}

	/* Check that the wavecal images are sorted. */
	if (w_option != STIS_NEAREST) {
	    midpt = wavecal->midpt[0];
	    for (n = 1;  n < wavecal->nimages;  n++) {
		if (midpt >= wavecal->midpt[n]) {
		    printf (
"Warning  Times of wavecal exposures are not strictly increasing; \\\n");
		    printf (
"Warning  interpolation option will be reset to nearest neighbor.\n");
		    w_option = STIS_NEAREST;
		    break;
		}
	    }
	}

	/* Extract the good shifts to temporary arrays. */
	wav_midpt1 = malloc (wavecal->nimages * sizeof (double));
	wav_midpt2 = malloc (wavecal->nimages * sizeof (double));
	wav_shift1 = malloc (wavecal->nimages * sizeof (double));
	wav_shift2 = malloc (wavecal->nimages * sizeof (double));
	if (wav_midpt1 == NULL || wav_midpt2 == NULL ||
	    wav_shift1 == NULL || wav_shift2 == NULL)
	    return (OUT_OF_MEMORY);
	n1 = 0;
	n2 = 0;
	for (n = 0;  n < wavecal->nimages;  n++) {
	    if (fabs (wavecal->shift1[n]) < UNREASONABLE_SHIFT) {
		wav_midpt1[n1] = wavecal->midpt[n];
		wav_shift1[n1] = wavecal->shift1[n];
		n1++;
	    }
	    if (fabs (wavecal->shift2[n]) < UNREASONABLE_SHIFT) {
		wav_midpt2[n2] = wavecal->midpt[n];
		wav_shift2[n2] = wavecal->shift2[n];
		n2++;
	    }
	}

	if (wavecal->verbose && (n1 > 1 || n2 > 1)) {
	    if (w_option == STIS_NEAREST)
		printf (
"         The wavecal that is nearest in time will be selected.\n");
	    else if (w_option == STIS_LINEAR)
		printf (
"         Linear interpolation will be used to get the output shifts.\n");
	}

	/* Process each imset in the science data file. */

	for (extverSci = 1;  extverSci <= scidata->nimages;  extverSci++) {

	    PrGrpBegin ("imset", extverSci);

	    /* Open the science file SCI extension header read-only,
		and get exposure time and midpoint of the exposure.
	    */
	    im = openInputImage (scidata->input, "SCI", extverSci);
	    if (hstio_err())
		return (OPEN_FAILED);
	    getHeader (im, &hdr);
	    if (hstio_err())
		return (OPEN_FAILED);
	    if ((status = GetSciGrp (&hdr, &exptime, &midpt)))
		return (status);
	    closeImage (im);
	    freeHdr (&hdr);

	    /* Get the shift in each axis for the current science file imset. */
	    if ((status = MatchWav (wav_midpt1, wav_shift1, n1,
                                    midpt, w_option, &shift1)))
		return (status);
	    if ((status = MatchWav (wav_midpt2, wav_shift2, n2,
                                    midpt, w_option, &shift2)))
		return (status);

	    /* Record the shifts in SHIFTA1 & SHIFTA2. */
	    if ((status = TargPos (scidata, extverSci, shift1, shift2)))
		return (status);

	    if (extverSci == 1) {		/* update primary header */
		initHdr (&phdr);
		im = openUpdateImage (scidata->input, "", 0, &phdr);
		if (hstio_err())
		    return (OPEN_FAILED);
		if ((status = History12 (&phdr, wavecal->input)))
		    return (status);
		putHeader (im);
		if (hstio_err())
		    return (OPEN_FAILED);
		closeImage (im);
		freeHdr (&phdr);
	    }

	    PrGrpEnd ("imset", extverSci);
	}

	free (wav_midpt1);
	free (wav_midpt2);
	free (wav_shift1);
	free (wav_shift2);

	free (wavecal->midpt);
	free (wavecal->shift1);
	free (wavecal->shift2);

	return (0);
}
