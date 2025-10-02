/* This file contains:
	WaveCal
	FreeLampSpec
*/

# include <stdio.h>
# include <stdlib.h>
# include <string.h>

# include "c_iraf.h"
# include "hstcal.h"
# include "hstio.h"
# include "stis.h"
# include "calstis4.h"
# include "hstcalerr.h"
# include "stisdq.h"
# include "stisdef.h"

static void ScaleTrim (StisInfo4 *);
static void ScaleOne (int, double, int, int []);
static void SaveDispCoeff (DispRelation *);
static void PrintWCP (StisInfo4 *);
static void PrintSection (StisInfo4 *);
static void PrintRefMsg (StisInfo4 *);
static void FreeLampSpec (LampInfo *);

/* This routine finds the shift in the wavecal spectrum in both the
   wavelength and spatial directions and assigns header keyword values
   with those shifts.  The keywords are in the SCI extension headers.
   The primary header will be updated with history information, and
   the calibration switch WAVECAL in the header will be reset from
   "PERFORM" to "COMPLETE".

   Phil Hodge, 1997 Sept 11:
	Print pixel size, if verbose.

   Phil Hodge, 1998 Dec 11:
	Get info from WCPTAB table; print parameters if verbose.
	Add extver to the calling sequence of FlagCR.

   Phil Hodge, 1999 Sept 23:
	Don't try to process an imset if EXPTIME <= 0.

   Phil Hodge, 2000 Jan 14:
	Do global offset for echelle data.

   Phil Hodge, 2000 Mar 10:
	Change the format for printing the shifts from %.5g to %.3f.

   Phil Hodge, 2000 July 21:
	cmplx.h is now in ../

   Phil Hodge, 2001 Mar 7:
	Also get dispersion coeff and trace info for first-order data;
	pass disp to WaveShift and trace to SpatialShift.  Print reference
	table info in PrintRefMsg.

   Phil Hodge, 2001 May 3:
	For echelle data, get MAMA offset info, and apply correction to the
	dispersion coefficients.

   Phil Hodge, 2004 July 23:
	Print slit_angle, if it's non-zero, for echelle data.

   Phil Hodge, 2007 March 9:
	Felix Stoehr of the ST-ECF pointed out that the exposure start
	time was being used by GetTrace4 before it was read by GetGrpInfo4.
	To fix this bug, move the call to GetTrace4 to inside the loop over
	extver, to a point after the call to GetGrpInfo4.

   Phil Hodge, 2007 May 2:
	Remove extver from the calling sequence for EchShift.  Initialize
	ref_names_printed to 0, and test on this rather than imset == 1
	to determine whether to call PrintRefMsg.  Add function getMinMax
	to check whether data are constant.

   Phil Hodge, 2008 Nov 3:
	Delete function getMinMax, and use sts->imset_ok.

   Phil Hodge, 2011 Jan 5:
	Add function SaveDispCoeff.  Move clamp from this function to EchShift.

   Phil Hodge, 2011 Feb 2:
        Initialize disp.ncoeff to 0.
*/

int WaveCal (StisInfo4 *sts) {

/* arguments:
StisInfo4 *sts    i: calibration switches and info
*/

	int status;

	IODescPtr iim;		/* descriptor for input image */
	Hdr phdr;		/* primary header for input image */
	Hdr hdr;		/* a SCI extension header */
	SingleGroup in;		/* input data */
	int extver;		/* loop index for extension version number */
	LampInfo lamp;		/* reference spectrum of cal lamp */
	ApInfo slit;		/* description of slit */
	DispRelation disp;	/* dispersion relation */
	char ref_aper[STIS_CBUF+1];	/* name of reference aperture */
	double angle;		/* incidence angle */
	double w_shift, s_shift; /* shifts in wavelength and slit direction */
	/* The 1-D spectrum computed by WaveShift by adding along the
	    spatial direction is used as weights by SpatialShift.
	*/
	/* set to true after the ref file names have been printed */
	int ref_names_printed = 0;

	int FlagCR (StisInfo4 *, SingleGroup *, int);
	int GetAngle4 (StisInfo4 *, char *, double *);
	int GetApDes4 (StisInfo4 *, ApInfo *);
	int GetDisp4 (StisInfo4 *, DispRelation *, char *);
	int GetGrpInfo4 (StisInfo4 *, Hdr *);
	int GetInang4 (StisInfo4 *, DispRelation *, double);
	int GetLamp (StisInfo4 *, LampInfo *);
	int GetTrace4 (StisInfo4 *, SpTrace **);
	int GetWCP (StisInfo4 *);
	int History4 (StisInfo4 *, Hdr *);
	int EchShift (StisInfo4 *, SingleGroup *,
		LampInfo *, DispRelation *, SpTrace *,
		double *, double *);
	void ScaleRef (StisInfo4 *, double *, double *);
	int SpatialShift (StisInfo4 *, ApInfo *, SpTrace *,
		SingleGroup *, double *);
	int UpdateShift (StisInfo4 *, int, double, double);
	int WaveShift (StisInfo4 *, ApInfo *, DispRelation *, LampInfo *,
		SingleGroup *, double *);

	SpTrace *trace;		/* list of spectrum traces */
	void FreeTrace4 (SpTrace **);

	/* memory not allocated yet */
	lamp.allocated = 0;
	trace = NULL;

        /* the dispersion relation is not used for first-order data */
        disp.ncoeff = 0;

	/* Get parameters that control wavecal processing. */
	if (sts->wcptab.exists == EXISTS_YES) {
	    if ((status = GetWCP (sts)))
		return (status);
	}
	if (sts->verbose)
	    PrintWCP (sts);

	/* Get the reference spectrum for the calibration lamp used. */
	if ((status = GetLamp (sts, &lamp)))
	    return (status);

	if (sts->disp_type == ECHELLE_DISP || sts->disp_type == PRISM_DISP) {

	    /* get dispersion relation, and get reference aperture name */
	    if ((status = GetDisp4 (sts, &disp, ref_aper)))
		return (status);

	    /* get incidence angle from apdestab */
	    if ((status = GetAngle4 (sts, ref_aper, &angle)))
		return (status);

	    /* Get incidence angle coefficients, and update the dispersion coeff
	       for the current angle.  The disp coeff are modified in-place.
	    */
	    if ((status = GetInang4 (sts, &disp, angle)))
		return (status);
	}

	SaveDispCoeff (&disp);		/* copy coeff to coeff_save */

	if (sts->disp_type != ECHELLE_DISP) {

	    /* Get aperture description from apdestab. */
	    if ((status = GetApDes4 (sts, &slit)))
		return (status);
	}

	initHdr (&hdr);
	for (extver = 1;  extver <= sts->nimages;  extver++) {

	    trlmessage("");
	    PrGrpBegin ("imset", extver);

	    if (sts->dbg != NULL)
		fprintf (sts->dbg, "# Begin imset %d ###\n", extver);

	    initSingleGroup (&in);

	    getSingleGroup (sts->input, extver, &in);
	    if (hstio_err())
		return (OPEN_FAILED);

	    sts->nx = in.sci.data.nx;
	    sts->ny = in.sci.data.ny;

	    /* Get keyword values from extension header. */
	    if ((status = GetGrpInfo4 (sts, &in.sci.hdr)))
		return (status);

	    if (sts->disp_type == ECHELLE_DISP ||
		sts->disp_type == PRISM_DISP) {
		/* Get the list of spectral traces, one for each order. */
		if ((status = GetTrace4 (sts, &trace)))
		    return (status);
		if (sts->verbose && sts->trace_rotation != 0.) {
		  trlmessage("         trace was rotated by = %.6g degree.",
                       sts->trace_rotation);
		}
	    }

	    if (sts->imset_ok) {

		/* Get the image sections to use when finding the shifts. */
		ScaleTrim (sts);
		if (sts->verbose)
		    PrintSection (sts);

		trlmessage("");
		PrSwitch ("wavecorr", PERFORM);

		if (!ref_names_printed) {
		    PrintRefMsg (sts);
		    ref_names_printed = 1;
		}

		/* Flag cosmic rays in DQ (if that bit is set). */
		if (sts->detector == CCD_DETECTOR &&
			sts->sdqflags & DATAREJECT) {
		    trlmessage("");
		    PrSwitch ("flagcr", PERFORM);
		    if ((status = FlagCR (sts, &in, extver)))
			return (status);
		    PrSwitch ("flagcr", COMPLETE);
		}

		if (sts->disp_type == ECHELLE_DISP) {

		    if ((status = EchShift (sts, &in,
				&lamp, &disp, trace,
                                &w_shift, &s_shift)))
			return (status);

		} else {	/* first order data */

		    /* Determine the shifts in each direction separately. */
		    if ((status = WaveShift (sts, &slit, &disp, &lamp, &in,
                                             &w_shift))) {
			if (status != NO_GOOD_DATA)
			    return (status);
			status = 0;
			w_shift = UNDEFINED_SHIFT;
		    }
		    if ((status = SpatialShift (sts, &slit, trace, &in,
                                                &s_shift))) {
			if ((status != NO_GOOD_DATA))
			    return (status);
			status = 0;			/* not fatal */
			s_shift = UNDEFINED_SHIFT;
		    }
		}

		if (w_shift == UNDEFINED_SHIFT) {
		    trlwarn("Shift in dispersion direction could not be determined.");
		} else {
		    trlwarn("         Shift in dispersion direction is %.3f pixels.", w_shift);
		}
		if (s_shift == UNDEFINED_SHIFT) {
		    trlwarn("Shift in spatial direction could not be determined.");
		} else {
		    trlwarn("         Shift in spatial direction is %.3f pixels.", s_shift);
		}

	    } else {

		char msg1[81], msg2[81];
		sprintf (msg2,
		    "wavecal imset %d skipped (IMSET_OK = F)\n", extver);
		strcpy (msg1, "Warning  ");
		strcat (msg1, msg2);
		trlmessage("%s", msg1);
		if (sts->dbg != NULL) {
		    strcpy (msg1, "# Warning:  ");
		    strcat (msg1, msg2);
		    fprintf (sts->dbg, "%s", msg1);
		}
		w_shift = UNDEFINED_SHIFT;
		s_shift = UNDEFINED_SHIFT;
	    }

	    /* Scale the shifts to reference pixel size. */
	    ScaleRef (sts, &w_shift, &s_shift);

	    PrSwitch ("wavecorr", COMPLETE);

	    /* Update current SCI extension header with the shift values. */
	    if ((status = UpdateShift (sts, extver, w_shift, s_shift)))
		return (status);

	    /* Write history records to primary header. */
	    if (extver == 1) {
		initHdr (&phdr);
		iim = openUpdateImage (sts->input, "", 0, &phdr);
		if (hstio_err())
		    return (OPEN_FAILED);
		if ((status = History4 (sts, &phdr)))
		    return (status);
		putHeader (iim);
		if (hstio_err())
		    return (OPEN_FAILED);
		closeImage (iim);
		freeHdr (&phdr);
	    }

	    freeSingleGroup (&in);
	    if (trace != NULL) {
		FreeTrace4 (&trace);
		trace = NULL;
	    }
	    PrGrpEnd ("imset", extver);
	    if (sts->printtime)
		TimeStamp ("Ending current imset", sts->rootname);
	}

	FreeLampSpec (&lamp);			/* free memory */

	return (0);
}

/* This routine adjusts the scales the size of the trim regions depending
   on binning, and computes the image sections to use for determining
   the shift in the wavelength and spatial directions.
*/
static void ScaleTrim (StisInfo4 *sts) {

/* argument:
StisInfo4 *sts      i: calibration switches and info

the following (zero indexed) values are assigned in sts:
    for determining the shift in the wavelength direction:
int wl_sect1[0], wl_sect1[1]     first and last pixels to use in first axis
int wl_sect2[0], wl_sect2[1]     first and last pixels to use in second axis
    for determining the shift in the spatial direction:
int sp_sect1[0], sp_sect1[1]     first and last pixels to use in first axis
int sp_sect2[0], sp_sect2[1]     first and last pixels to use in second axis
*/

	/* shift in dispersion direction, first and second axes */
	ScaleOne (sts->wl_trim1, sts->scale[0], sts->nx, sts->wl_sect1);
	ScaleOne (sts->wl_trim2, sts->scale[1], sts->ny, sts->wl_sect2);

	/* shift in cross dispersion direction, first and second axes */
	ScaleOne (sts->sp_trim1, sts->scale[0], sts->nx, sts->sp_sect1);
	ScaleOne (sts->sp_trim2, sts->scale[1], sts->ny, sts->sp_sect2);
}

static void ScaleOne (int trim, double scale, int npix, int section[]) {

	double btrim;		/* trim size scaled by binning */
	int itrim;		/* nearest integer to btrim */

	/* In zero-indexed reference pixels, the first and last pixels
	   to use are trim and (npix - 1 - trim) respectively.
	   We assume the wavecal is not a subarray, but it could be
	   binned, so adjust for the bin size.
	*/
	btrim = (double)trim / scale;
	itrim = NINT (btrim);

	section[0] = itrim;
	section[1] = npix - 1 - itrim;

	/* Make sure the endpoints are within the image. */
	if (section[0] < 0)
	    section[0] = 0;
	if (section[1] > npix - 1)
	    section[1] = npix - 1;

	/* Make sure the trim region wasn't too large. */
	if (section[1] <= section[0]) {
	    section[0] = 0;
	    section[1] = npix - 1;
	}
}

/* This copies the dispersion coefficients from coeff to coeff_save. */

static void SaveDispCoeff (DispRelation *disp) {

	int i;

	for (i = 0;  i < disp->ncoeff;  i++)
	    disp->coeff_save[i] = disp->coeff[i];

	for (i = disp->ncoeff;  i < MAX_DISP_COEFF;  i++)
	    disp->coeff_save[i] = 0.;
}

/* This routine prints values read from the wavecal parameters table. */

static void PrintWCP (StisInfo4 *sts) {

	if (sts->disp_type == ECHELLE_DISP) {
	    if (sts->slit_angle != 0.) {
		trlmessage("         Slit angle = %.5g degrees",
			sts->slit_angle / DEGREES_TO_RADIANS);
	    }
	    trlmessage("         Wavecal parameters are:");
	    trlmessage("         WL_TRIM1 = %d", sts->wl_trim1);
	    trlmessage("         WL_TRIM2 = %d", sts->wl_trim2);
	} else {
	    if (sts->slit_angle != 0.) {
		trlwarn("Slit angle was specified for non-echelle data, will be ignored.");
	    }
	    trlmessage("         Wavecal parameters are:");
	    trlmessage("         WL_TRIM1 = %d", sts->wl_trim1);
	    trlmessage("         WL_TRIM2 = %d", sts->wl_trim2);
	    trlmessage("         SP_TRIM1 = %d", sts->sp_trim1);
	    trlmessage("         SP_TRIM2 = %d", sts->sp_trim2);
	    trlmessage("         WL_RANGE = %d", sts->wl_range);
	    trlmessage("         SP_RANGE = %d", sts->sp_range);
	}
	if (sts->detector == CCD_DETECTOR) {
	    trlmessage("");
	    trlmessage("         NSIGMA_CR = %.6g", sts->nsigma_cr);
	    trlmessage("         NSIGMA_ILLUM = %.6g",
				sts->nsigma_illum);
	    trlmessage("         MAD_REJECT = %.6g", sts->mad_reject);
	    trlmessage("         MIN_MAD = %.6g", sts->min_mad);
	} else {
	    trlmessage("");
	}
}

/* Only a subset of the full image is used for finding the shifts.
   The region to exclude is given by the WL_TRIMi and SP_TRIMi values
   read from the WCP table, but these can be ignored if they're too
   large, and they can get scaled depending on binning.  More meaningful
   values from the user's perspective are the image sections of the
   data that will actually be used; note that for echelle data one
   section is used, since both shifts will be determined simultaneously.
   This routine prints the image sections that will be used when
   finding the shifts in dispersion and cross dispersion directions.
   It also prints info about the pixel size.
*/

static void PrintSection (StisInfo4 *sts) {

	if (sts->disp_type == ECHELLE_DISP) {

	    trlmessage("         Section to use for finding shift:");
	    trlmessage("             %d:%d, %d:%d",
		sts->wl_sect1[0]+1, sts->wl_sect1[1]+1,
		sts->wl_sect2[0]+1, sts->wl_sect2[1]+1);

	} else {

	    trlmessage("         Section to use for finding shift in dispersion direction: ");
	    trlmessage("             %d:%d, %d:%d",
		sts->wl_sect1[0]+1, sts->wl_sect1[1]+1,
		sts->wl_sect2[0]+1, sts->wl_sect2[1]+1);

	    trlmessage("         Section to use for finding shift in spatial direction: ");
	    trlmessage("             %d:%d, %d:%d",
		sts->sp_sect1[0]+1, sts->sp_sect1[1]+1,
		sts->sp_sect2[0]+1, sts->sp_sect2[1]+1);
	}

	if (sts->scale[0] == 1. && sts->scale[1] == 1.) {
	    trlmessage("         Image pixels are reference pixel size.");
	} else {
	    trlmessage("         Image pixel size is %.2g by %.2g reference pixels.",
		sts->scale[0], sts->scale[1]);
	}

	fflush (stdout);
}

/* This routine prints the names, etc., of the reference tables. */

static void PrintRefMsg (StisInfo4 *sts) {

	if (sts->wcptab.exists == EXISTS_YES) {
		PrRefInfo ("wcptab", sts->wcptab.name,
		sts->wcptab.pedigree,
		sts->wcptab.descrip, sts->wcptab.descrip2);
	}
	PrRefInfo ("lamptab",
		sts->lamptab.name, sts->lamptab.pedigree,
		sts->lamptab.descrip, sts->lamptab.descrip2);

	PrRefInfo ("apdestab", sts->apdestab.name,
		sts->apdestab.pedigree,
		sts->apdestab.descrip, sts->apdestab.descrip2);

	if (sts->disp_type == ECHELLE_DISP || sts->disp_type == PRISM_DISP) {
	    PrRefInfo ("disptab",
		sts->disptab.name, sts->disptab.pedigree,
		sts->disptab.descrip, sts->disptab.descrip2);
	    PrRefInfo ("inangtab",
		sts->inangtab.name, sts->inangtab.pedigree,
		sts->inangtab.descrip, sts->inangtab.descrip2);
	    PrRefInfo ("sptrctab",
		sts->sptrctab.name, sts->sptrctab.pedigree,
		sts->sptrctab.descrip, sts->sptrctab.descrip2);
	}
	if (sts->disp_type == PRISM_DISP) {
	    PrRefInfo ("sdctab",
		sts->sdctab.name, sts->sdctab.pedigree,
		sts->sdctab.descrip, sts->sdctab.descrip2);
	}

	fflush (stdout);
}

/* This routine frees memory for the calibration lamp spectrum,
   if it has been allocated, and resets the flag to indicate that
   memory is no longer allocated.
*/

static void FreeLampSpec (LampInfo *lamp) {

	if (lamp->allocated) {
	    free (lamp->wl);
	    free (lamp->flux);
	    lamp->allocated = 0;
	}
}
