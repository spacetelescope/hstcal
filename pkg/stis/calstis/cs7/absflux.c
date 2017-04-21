# include <stdio.h>

# include "c_iraf.h"
# include "hstio.h"

# include "stis.h"
# include "calstis7.h"
# include "stisdq.h"
# include "stispht.h"
# include "stisdef.h"
# include "stistds.h"
# include "hstcalerr.h"

# define  HST_AREA  45238.93416		/* cm**2, including obscured areas */
# define  H_PLANCK  (6.6260755e-27)	/* Planck's constant, erg * sec */
# define  C_LIGHT   29979245800.	/* cm / sec */
# define  CM_PER_ANGSTROM  (1.e-8)	/* cm / angstrom */

/* This routine corrects the output SCI and ERR data (in-place) for
   telescope throughput and detector sensitivity.  The data before
   correction are expected to be in counts; after correction the BUNIT
   keyword in the science and error extension headers will be set to
   "erg /s /cm**2 /angstrom /arcsec**2".

   Counts are converted to specific intensity by multiplying by:

				h * c * atodgain
	--------------------------------------------------------------
	exptime * R(lambda) * lambda * hst_area * disp * scale * width

   where:

	h = Planck's constant (erg * sec); c = speed of light (cm / sec)
	atodgain = analog to digital gain (electrons / dn)
	exptime = exposure time (sec)
	lambda = wavelength (cm)
	R(lambda) = instrumental response (e.g. output from synphot.calcband)
	hst_area = area of HST (cm**2), _including_ obscured regions
	disp = dispersion (angstrom / pixel)
	scale = image scale along slit (arcsec / pixel)
	width = slit width (arcsec)

   The R(lambda) that we got from the phottab is for an "infinite"
   extraction box height, which is what we need for a diffuse source.

   The DIFF2PT keyword is assigned the value (scale * width * pcorr) / Ta,
   where Ta is the aperture throughput for a point source, scale and width
   are as above, and pcorr is the factor read from the PCTAB, the ratio of
   R(lambda) for an infinite extraction box height divided by R(lambda)
   for the default extraction box.  DIFF2PT will be written to the SCI
   extension header.  Since Ta and pcorr are functions of wavelength and
   we only write one value to the header, we average the throughput over
   the wavelength range of the current extension when computing DIFF2PT.
### (Actually, this range is too large because it covers the full output
   image in the dispersion direction rather than just the pixels mapped
   from points within the input image.)
   For a point source, to convert the specific intensity to flux density
   Flambda, sum the SCI image values along the length of the slit and
   multiply the sum by DIFF2PT.

   Phil Hodge, 1997 Sept 5:
	Include atodgain as a factor.
	plate_scale is now an array of two elements, to allow for unequal
	binning in the two axes.

   Phil Hodge, 1997 Nov 13:
	Assume dispaxis = 1, and remove it from the calling sequence.
	Correct R(lambda) for default extraction box height using pcorr
	in the phot struct, and take account of this in the DIFF2PT factor.

   Phil Hodge, 1998 Jan 26:
	Remove the pcorr factor from the calibration for diffuse source;
	that is, assume R(lambda) is for an infinite extraction box height.
	(We've been assuming all along that R(lambda) has been corrected
	for the aperture throughput.)

   Phil Hodge, 2001 Feb 20:
	Check for zero or negative wavelength (needed for prism).

   Ivo busko, 2002 Jan 10:
	Time-dependent sensitivity correction.

   Ivo busko, 2002 Mar 6:
	MSM/blaze shift correction.

   Ivo Busko, 2002 Apr 24:
        Output blaze shift value.

   Phil Hodge, 2004 Dec 27:
	Change the calling sequence of TdsCorrection, to include
	temperature dependence in addition to time.

   Phil Hodge, 2010 Mar 24:
	Mike Droettboom discovered by using valgrind that tds_starti was not
	initialized.  Instead, thr_starti was initialized twice, so the
	second statement has been changed to initialize tds_starti.
*/

int AbsFlux (StisInfo7 *sts, SingleGroup *out, CoordInfo *coord_o,
		PhotInfo *phot, ApInfo *slit, TdsInfo *tds,
		double *plate_scale, double atodgain,
		double exptime, double hfactor, int sporder,
		double *blazeshift) {

/* arguments:
StisInfo7 *sts         i: calibration switches and info
SingleGroup *out       io: output data
CoordInfo *coord_o     i: coordinate parameters
PhotInfo *phot         i: photometry info
ApInfo *slit           i: slit info (for throughput)
TdsInfo *tds           i: time-dependent sensitivity info
double plate_scale[2]  i: arcsec / pixel for each axis
double atodgain        i: electrons / adu for CCD
double exptime         i: exposure time in seconds
double hfactor         i: divide to convert back to observed wavelength
int sporder            i: spectral order number
double *blazeshift;    o: blaze shift value actually used
*/

	int status;

	double wl;		/* heliocentric wavelength at a pixel */
	double response;	/* instrumental response at a wavelength */
	double pct_factor;	/* correction to infinite extr box height */
	double tds_factor;	/* correction for time-dependent sens. */
	double throughput;	/* factor to correct for slit throughput */
	double sum;		/* for getting average throughput */
	double photfactor;	/* to get inverse sensitivity */
	float correction;	/* combined correction factor */
	float diff_pt;		/* conversion from diffuse source to point */
	float cont_eml;		/* conversion from continuum to emission line */
	int i, j;
	int abs_starti;		/* index to begin search in interp1d */
	int pct_starti;
	int thr_starti;
	int tds_starti;
	double dispersion;
	double *tds_factors;	/* array with time-dependent sensitivity
	                           correction factors */
	double *wlt;		/* temp array for blaze shift */

	void TdsCorrection (TdsInfo *, double, double, double *);
	void BlazeCorr (PhotInfo *, double, double, int, double, double *,
	                double);

	abs_starti = 1;				/* initial values */
	pct_starti = 1;
	thr_starti = 1;
	tds_starti = 1;

	sum = 0.;

	/* Generate time-dependent sensitivity correction factors. */
	if (sts->tdscorr == PERFORM) {
	    tds_factors = (double *) malloc (tds->nwl * sizeof (double));
	    if (tds_factors == NULL)
	        return (OUT_OF_MEMORY);

	    TdsCorrection (tds, sts->expstart, sts->detector_temp,
				tds_factors);
	}

	/* incomplete; will also divide by average throughput */
	diff_pt = coord_o->cdelt[1] * slit->width[0];

	cont_eml = coord_o->cdelt[0] * slit->width[0] /
			plate_scale[0];

	/* The photometry table gives the instrumental response for
	   a point source.  This factor is part of the conversion
	   from count rate to specific intensity.
	*/
	photfactor = H_PLANCK * C_LIGHT * atodgain / (HST_AREA * exptime *
		coord_o->cdelt[0] * coord_o->cdelt[1] * slit->width[0]);

	/* Apply MSM/blaze correction to throughput array. Save the
	   original wavelength array for later use.
	 */
	wlt = (double *) malloc (phot->nelem * sizeof (double));
	if (wlt == NULL)
	    return (OUT_OF_MEMORY);
	for (i = 0;  i < phot->nelem; i++)
	    wlt[i] = phot->wl[i];

	/* Get dispersion from data. This is required in the case the
	   reference values from the _pht table are invalid, AND we are
	   passing a blazeshift value thru the command line. In this case,
	   we use the dispersion at the mid point in the current spectral
	   order.
	*/
	i = out->sci.data.nx / 2;
	wl         = (((double)i     - coord_o->crpix[0]) * coord_o->cdelt[0] +
			coord_o->crval[0]);
	dispersion = (((double)(i+1) - coord_o->crpix[0]) * coord_o->cdelt[0] +
			coord_o->crval[0]) - wl;
	dispersion /= hfactor;

	if (phot->wpos != 0.0 || sts->blazeshift != NO_VALUE) {
	    BlazeCorr (phot, sts->expstart, sts->expend, sporder,
                       sts->blazeshift, blazeshift, dispersion);
	}

	/* for each pixel in the dispersion direction ... */
	for (i = 0;  i < out->sci.data.nx; i++) {

	    /* convert back to observed wavelength */
	    wl = (((double)i - coord_o->crpix[0]) * coord_o->cdelt[0] +
			coord_o->crval[0]) / hfactor;	/* Angstroms */
	    if (wl <= 0.) {
		for (j = 0;  j < out->dq.data.ny;  j++) {
		    DQSetPix (out->dq.data, i, j,
				DQPix(out->dq.data,i,j) | CALIBDEFECT);
		}
		continue;
	    }

	    /* Interpolate in throughput array with shifted blaze. */
	    response = interp1d (wl, phot->wl, phot->thru,
			phot->nelem, &abs_starti);

	    /* Interpolate remaining factors with non-shifted blaze. */
	    pct_factor = interp1d (wl, wlt, phot->pcorr,
			 phot->nelem, &pct_starti);
	    throughput = interp1d (wl, slit->wl, slit->thr,
			 slit->nelem, &thr_starti);
	    if (sts->tdscorr == PERFORM)
	        tds_factor  = interp1d (wl, tds->wl, tds_factors,
	                                tds->nwl, &tds_starti);
	    else
	        tds_factor = 1;

	    /* We need the average throughput for diff2pt. */
	    sum += (throughput / pct_factor);
	    if (response <= 0.)
		continue;

	    wl *= CM_PER_ANGSTROM;	/* convert from Angstroms to cm */
	    correction = (float) (photfactor / (response * wl));

	    correction /= tds_factor;

	    /* for each pixel along the slit ... */
	    if (correction <= 0.) {
		for (j = 0;  j < out->dq.data.ny;  j++) {
		    DQSetPix (out->dq.data, i, j,
				DQPix(out->dq.data,i,j) | CALIBDEFECT);
		}
	    } else {
		for (j = 0;  j < out->sci.data.ny;  j++) {
		    Pix (out->sci.data, i, j) *= correction;
		    Pix (out->err.data, i, j) *= correction;
		}
	    }
	}

	free (wlt);

	/* average throughput */
	throughput = sum / (double)(out->sci.data.nx);

	/* Update BUNIT in the science and error extension headers. */
	if ((status = Put_KeyS (&out->sci.hdr, "BUNIT", phot->bunit,
                                "units for flux-calibrated data")))
	    return (status);
	if ((status = Put_KeyS (&out->err.hdr, "BUNIT", phot->bunit,
                                "units for flux-calibrated data")))
	    return (status);

	if (throughput > 0.)
	    diff_pt /= throughput;
	else
	    diff_pt = -9999.;
	if ((status = Put_KeyF (&out->sci.hdr, "DIFF2PT", diff_pt,
                                "diffuse to point source conversion factor")))
	    return (status);

	if ((status = Put_KeyF (&out->sci.hdr, "CONT2EML", cont_eml,
                                "continuum to emission line conversion factor")))
	    return (status);
	if ((status = Put_KeyF (&out->sci.hdr, "SCALE_A1", plate_scale[0],
                                "arcsec / pixel in dispersion direction")))
	    return (status);
	if ((status = Put_KeyD (&out->sci.hdr, "OMEGAPIX",
                plate_scale[0] * plate_scale[1], "pixel area, arcsec^2")))
	    return (status);

	if (sts->tdscorr == PERFORM)
	    free (tds_factors);

	return (0);
}
