/* This file contains:
	doDark
	MedSciVal
*/

# include <stdio.h>
# include <stdlib.h>		/* calloc */
# include <math.h>		/* fabs */

# include "c_iraf.h"
# include "hstio.h"
# include "stis.h"
# include "calstis1.h"
# include "hstcalerr.h"
# include "stistemperature.h"	/* defines DARKRATE */
# include "stisdef.h"

static int getDarkParam (Hdr *, double *, double *);
static double CCDFactor (double, double, double);
static int MedSciVal (SingleGroup *, float *);
static int meanNUVDark(SingleGroup *, double *);
static void OverrideFactor (StisInfo1 *, int, double *, int *);


/* This routine subtracts the dark image from x (in-place).
   For CCD data, the dark image is multiplied by the exposure time and
   divided by the gain before subtracting.  The dark time is just the
   exposure time; it DOES NOT INCLUDE the idle time since the last
   flushing of the chip or the readout time.

   For MAMA data, the dark image is just multiplied by the exposure time
   before subtracting.

   Phil Hodge, 1997 Oct 1:
	Change avg from 1 to 0, so bin2d will sum within bins instead of
	averaging.

   Phil Hodge, 1998 Feb 6:
	Remove the include for "../stissizes.h".

   Phil Hodge, 1998 Mar 13:
	Change calling sequence of DoppConv.

   Phil Hodge, 1998 Aug 6:
	Add border=0 to calling sequence of DoppConv.

   Phil Hodge, 1998 Oct 7:
	Include NUVFactor, and scale the dark by that factor.

   Phil Hodge, 1998 Nov 13:
	In NUVFactor, scale the expression from the Instrument Handbook
	by a constant, rather than by the mean of the dark image.  The
	calling sequence was changed, also.

   Phil Hodge, 1999 Mar 10:
	DARKRATE (used by NUVFactor) is now defined in ../stistemperature.h
	instead of in this file.

   Ivo Busko, 2001 Oct 26:
	Include CCDFactor to scale side 2 darks by CCD housing temperarure.

   Ivo Busko, 2002 Mar 20:
        Dark scale factor from command line.

   Ivo Busko, 2002 Mar 29:
        NUV dark time-dependency.

   Phil Hodge, 2003 Aug 26:
	Use median instead of mean for meandark.

   Paul Barrett, 2003 Sep 18:
        Added CCD_DETECTOR check before CCDFactor().

   Phil Hodge, 2004 Dec 27:
	Use detector_temp instead of temperature or ccd_temperature.

   Phil Hodge, 2010 May 6:
	Add function getDarkParam to read reference temperature and slope
	from the dark file primary header.  Add these arguments to the
	calling sequence of CCDFactor.

   Phil Hodge, 2012 July 6:
	Add function medianDark; include median_dark in the call to
	GetTdcCorr.

   Phil Hodge, 2013 Oct 3:
	Replace function medianDark with meanNUVDark; change median_dark to
	mean_nuv_dark in the call to GetTdcCorr.  meanNUVDark actually does
	compute the mean rather than the median of the NUV-MAMA dark
	reference image.
*/

int doDark (StisInfo1 *sts, SingleGroup *x, float *meandark, int sci_extver) {

/* arguments:
StisInfo1 *sts     i: calibration switches, etc
SingleGroup *x    io: image to be calibrated; written to in-place
float *meandark    o: median of dark image values subtracted
int sci_extver     i: IMSET number in input (science) file
*/

	int status;

	SingleGroup y, z;	/* y and z are scratch space */
	float *ds;		/* Doppler smearing array */
	double ref_temp;	/* reference temperature for CCD dark */
	double drk_vs_t;	/* slope of dark vs temperature */
	double factor;		/* scale factor, depending on temperature */
	int nds, d0;		/* size of ds, index of zero point */
	int extver = 1;		/* get this imset from dark image */
	int rx, ry;		/* for binning dark down to size of x */
	int x0, y0;		/* offsets of sci image */
	int same_size;		/* true if no binning of ref image required */
	int high_res;		/* true if high-res pixels in dispersion dir */
	int avg = 0;		/* bin2d should sum values within each bin */
	int override;		/* override flag for darkscale warning msg */

	int FindBin (StisInfo1 *, SingleGroup *, SingleGroup *,
		int *, int *, int *, int *, int *, int *);
	int GetTdcCorr (StisInfo1 *, double, double *);
	int MakeDopp (double, double, double, double, double, int,
		float *, int *, int *);
	int DoppConv (SingleGroup *, int, float *, int, int);

	if (sts->darkcorr != PERFORM)
	    return (0);

	initSingleGroup (&y);

	/* Get the dark image data. */
	getSingleGroup (sts->dark.name, extver, &y);
	if (hstio_err())
	    return (OPEN_FAILED);

	/* Compare binning of science image and reference image;
	   get same_size and high_res flags, and get info about
	   binning and offset for use by bin2d.
	*/
	if ((status = FindBin (sts, x, &y,
                               &same_size, &high_res, &rx, &ry, &x0, &y0)))
	    return (status);

	/* Do we need to do Doppler convolution? */
	if (sts->doppcorr == PERFORM) {

	    if (!high_res) {
		printf (
		"ERROR    Doppler convolution (DOPPCORR) was specified, \\\n");
		printf ("ERROR    but %s is binned to low-res pixels.\n",
			sts->dark.name);
		return (SIZE_MISMATCH);
	    }

	    /* Allocate space for the Doppler smearing array, making it
		larger than we will need.  The actual size nds will be
		updated by MakeDopp.
	    */
	    nds = 2 * (sts->doppmag + 1) + 21;
	    ds = (float *) calloc (nds, sizeof (float));

	    if ((status = MakeDopp (sts->doppzero, sts->doppmag, sts->orbitper,
                                    sts->expstart, sts->exptime, sts->dispsign,
                                    ds, &nds, &d0)))
		return (status);

	    /* Convolve y with the Doppler smearing function. */
	    if ((status = DoppConv (&y, 0, ds, nds, d0)))
		return (status);

	    free (ds);
	}

	/* Multiply the dark image by the exposure time and divide by the
	   atodgain (or just by exposure time for the MAMAs), and
	   subtract it from x.  For NUV MAMA data, there's an additional
	   scaling factor that depends on the temperature.
	*/

	/* Get temperature scaling factor. */

	factor = 1;

	if (sts->detector == NUV_MAMA_DETECTOR) {
	    double mean_nuv_dark;
	    if ((status = meanNUVDark(&y, &mean_nuv_dark)))
	        return (status);
	    if ((status = GetTdcCorr(sts, mean_nuv_dark, &factor)))
	        return (status);
	}
        else if (sts->detector == CCD_DETECTOR) {
	    if ((status = getDarkParam (y.globalhdr,
				&ref_temp, &drk_vs_t)) != 0)
		return (status);
            factor *= CCDFactor (sts->detector_temp, ref_temp, drk_vs_t);
        }

	/* Override it by optional command line argument. */

	OverrideFactor (sts, sci_extver, &factor, &override);

	if (same_size) {

	    /* No binning required. */

	    if ((status = multk2d (&y, factor * sts->exptime / sts->atodgain)))
		return (status);
	    if ((status = MedSciVal (&y, meandark)))
		return (status);
	    if ((status = sub2d (x, &y))) {
		printf ("ERROR    (darkcorr) size mismatch.\n");
		return (status);
	    }
	    freeSingleGroup (&y);

	} else {

	    /* Bin the dark image down to the actual size of x. */

	    initSingleGroup (&z);
	    allocSingleGroup (&z, x->sci.data.nx, x->sci.data.ny, True);
	    if ((status = bin2d (&y, x0, y0, rx, ry, avg, &z))) {
		printf ("ERROR    (darkcorr) size mismatch.\n");
		return (status);
	    }
	    freeSingleGroup (&y);			/* done with y */

	    if ((status = multk2d (&z, factor * sts->exptime / sts->atodgain)))
		return (status);
	    if ((status = MedSciVal (&z, meandark)))
		return (status);
	    if ((status = sub2d (x, &z)))
		return (status);
	    freeSingleGroup (&z);			/* done with z */
	}

	if (sts->verbose && (factor != 1. || override)) {
	    printf (
"         Dark reference image was scaled by the factor %.6g, \\\n", factor);
	    printf ("         in addition to the exposure time.\n");
	    if (override)
	        printf (
"Warning  Default dark scaling overriden by DARKSCALE parameter.\n");
	}

	return (0);
}

/* This function reads the primary header of the dark reference file to
   get the reference temperature (default 18 C) and the slope of the dark
   vs temperature (default 0.07).
*/

static int getDarkParam (Hdr *phdr, double *ref_temp, double *drk_vs_t) {

/* arguments:
Hdr phdr            i: primary header of dark reference file
double ref_temp     o: reference temperature
double drk_vs_t     o: slope of dark vs temperature
The function value is the status (0 is OK).
*/

	int use_def = 1;		/* use default if missing keyword */
	int status = 0;

	if ((status = Get_KeyD (phdr, "REF_TEMP", use_def, CCD_REF_TEMP,
				ref_temp)) != 0)
	    return (status);

	if ((status = Get_KeyD (phdr, "DRK_VS_T", use_def, 0.07,
				drk_vs_t)) != 0)
	    return (status);

	return (0);
}


/* The dark count rate for the side 2 CCD varies with temperature. For this
   detector, sts->detector_temp was gotten earlier from a header keyword.  In
   this routine we will determine the scale factor by which the dark image
   should be scaled.  The factor will be one for other detectors or if
   the temperature could not be obtained from the header; those cases are
   flagged by a value of temperature less than or equal to zero.
*/

static double CCDFactor (double temperature,
			double ref_temp, double drk_vs_t) {

/* arguments:
double temperature  i: from header keyword, in degrees Celsius
double ref_temp     i: reference temperature, read from primary header of
                        dark file
double drk_vs_t     i: slope of dark vs temperature, read from primary header
                        of dark file
The function value is the factor by which the dark image should be multiplied.
*/

	double factor;

	if (temperature <= 0.)
	    factor = 1.;
	else
	    factor = 1.0 + drk_vs_t * (temperature - ref_temp);

	return (factor);
}


/*  Overrides darkscale factor with command line value. */

static void OverrideFactor (StisInfo1 *sts, int extver, double *factor,
                            int *override){
	int i;

	*override = 0;

	if (sts->ndarkscale > 0) {
	    i = extver - 1;
	    if (extver > sts->ndarkscale)
	        i = sts->ndarkscale - 1;

	    *factor = sts->darkscale[i];

	    *override = 1;
	}
}


/* This routine computes the median of the science data values
   that are not flagged as bad by the data quality array.
*/

static int MedSciVal (SingleGroup *y, float *meandark) {

	float *dark;		/* a copy of the unflagged values in the dark */
	int numgood;		/* number of good pixels */
	int i, j;
	int inplace = 1;	/* flag:  compute the median in-place */

	if ((dark = (float *) calloc (y->sci.data.nx * y->sci.data.ny,
		sizeof (float))) == NULL) {
	    printf ("ERROR    (darkcorr) out of memory in MedSciVal.\n");
	    return (OUT_OF_MEMORY);
	}

	numgood = 0;
	for (j = 0;  j < y->sci.data.ny;  j++) {
	    for (i = 0;  i < y->sci.data.nx;  i++) {
		if (DQPix (y->dq.data, i, j) == 0) {
		    dark[numgood] = Pix (y->sci.data, i, j);
		    numgood++;
		}
	    }
	}

	if (numgood > 0)
	    *meandark = MedianFloat (dark, numgood, inplace);
	else
	    *meandark = 0.;

	free (dark);

	return (0);
}

/* This function computes the average value of the dark reference image,
   ignoring bad pixels.  The result is then multiplied by a factor
   (should be 4) to account for binning of the dark to low-res pixels.

   This function is only used for NUV-MAMA data.
*/

static int meanNUVDark(SingleGroup *y, double *mean_nuv_dark) {

    /*int status;*/
	int numgood;		/* number of good pixels */
	int i, j;
	double sum, mean_value;

	numgood = 0;
	sum = 0.;
	for (j = 0;  j < y->sci.data.ny;  j++) {
	    for (i = 0;  i < y->sci.data.nx;  i++) {
		if (DQPix(y->dq.data, i, j) == 0) {
		    sum += Pix(y->sci.data, i, j);
		    numgood++;
		}
	    }
	}

	if (numgood > 0)
	    mean_value = sum / (double)numgood;
	else
	    mean_value = 0.;

	*mean_nuv_dark = mean_value * (double)(y->sci.data.nx) / 1024. *
	                              (double)(y->sci.data.ny) / 1024.;

	return 0;
}
