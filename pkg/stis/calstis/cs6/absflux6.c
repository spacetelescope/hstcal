# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <math.h>

# include "xtables.h"
# include "stis.h"
# include "calstis6.h"
# include "stispht.h"
# include "stisdq.h"
# include "stissizes.h"
# include "stistds.h"
# include "hstcalerr.h"

static int getHalo (StisInfo6 *, int, double [], float [], CTICorrInfo *);
static float CtiCorr(StisInfo6 *, CTICorrInfo *,
        double, float, float, float, float);
static double getDispersion (RowContents *, int);
static float QFactor(StisInfo6 *, double, float, float);
double interp1d (double, double *, double *, int, int *);

/* Generates the FLUX array from the NET array, and also modifies the ERROR
   array in-place for telescope throughput and detector sensitivity. The
   data before correction are expected to be in c/s; after correction its
   units will be "erg /s /cm**2 /angstrom".

   Counts are converted to specific intensity by multiplying by:

                          h * c * H
            -----------------------------------------
            R(lambda) * Ta * lambda * hst_area * disp

   where:

        h = Planck's constant (erg * sec); c = speed of light (cm / sec)
        lambda = wavelength (cm)
        R(lambda) = instrumental response (e.g. output from synphot.calcband)
        hst_area = area of HST (cm**2), including obscured regions
        disp = dispersion (angstrom / pixel)
        Ta = aperture throughput for a point source
        H = correction to infinite aperture from PCTAB



   Revision history:
   ----------------
   26 Feb 97  -  Implemented, with some code adapted from a similar
                 routine in calstis7 (I.Busko)
   10 Apr 97  -  Changes after code review (IB):
                 - C_LIGHT in scientific notation
   03 Sep 97  -  Add ATODGAIN correction (IB)
   28 Jan 98  -  Add extraction box height correction (IB)
   09 Apr 99  -  Check for zeroed throughput - OPR 38749 (IB)
   17 Dec 99  -  Flux normalization factor in opt. extr. rejection (IB)
   17 Aug 00  -  Fix dispersion computation in hires data (IB)
   26 Feb 01  -  Got rid of EvalDisp (IB)
   06 Mar 01  -  Correct wavelengths back to observed values (IB)
   09 Jan 02  -  Time-dependent sensitivity (IB)
   04 Feb 02  -  Blaze shift correction (IB)
   17 Apr 02  -  Output blaze shift value (IB)
   24 Jun 03  -  Correct error array for effect of Q factor (PB).
   19 Aug 03  -  Add CTI correction to flux calibration (PB).
   15 Dec 03  -  Correct Q factor algorithm for extraction width (PB).
   27 Dec 04  -  Change calling sequence of TdsCorrection, to include
                 temperature dependence as well as time (PEH).
   08 Apr 05  -  Include GAC, the grating-aperture correction (PEH).
   28 Apr 05  -  In CtiCorr, include ltv offset for y (PEH).
   20 Feb 06  -  Include getHalo(); in CtiCorr, use the halpar term (PEH).
   20 Apr 06  -  For pct_factor, use wavelengths that have not been
                 modified for the blaze shift.
                 In getHalo, return if ctecorr is not PERFORM (PEH).
   10 Feb 12  -  Apply the Q factor to the NET_ERROR column (PEH).
*/

int AbsFlux6 (StisInfo6 *sts, RowContents *out,
              PhotInfo *phot, PhotInfo *photc, ApInfo *slit, TdsInfo *tds,
              CTICorrInfo *cti, double hfactor, double atodgain, int optimal,
              double *pfactor, double *blazeshift)
{

    /* arguments:
       StisInfo6 *sts         i: calibration switches and info
       RowContents *out       io: output data
       PhotInfo *phot         i: photometry info
       PhotInfo *photc        i: for extraction box height correction
       ApInfo *slit           i: slit info (for throughput)
       TdsInfo *tds           i: time-dependent sensitivity info
       double hfactor         i: divide to convert back to observed wavelength
       double atodgain        i: CCD ATODGAIN value
       int optimal            i: is optimal extraction ?
       double *pfactor;       i: flux normalization factor from opt. extr.
                                 profile
       double *blazeshift;    o: blaze shift value actually used
    */

        double wl;              /* observed wavelength at a pixel */
        double response;        /* instrumental response at a wavelength */
        double pct_factor;      /* correction to infinite extr box height */
        double pct_factor2;     /* correction for current extr box height */
        double tds_factor;      /* correction for time-dependent sens. */
        double gac_factor;      /* grating-aperture correction */
        double throughput;      /* factor to correct for slit throughput */
        double photfactor;      /* to get inverse sensitivity */
        double dispersion;      /* dispersion in A/pix */
        float correction;       /* combined correction factor */
        float ecorrection;      /* combined correction factor for the error */
        int i;
        int status;
        int abs_starti;         /* index to begin search in interp1d */
        int pct_starti;
        int thr_starti;
        int gac_starti;
        int tds_starti;
        double hf;
        double *tds_factors;    /* array with time-dependent sensitivity
                                   correction factors */
        double *wl_no_blaze;    /* phot->wl before blaze shift correction */

        void BlazeCorr (PhotInfo *, double, double, int, double, double *,
                        double);
        void TdsCorrection (TdsInfo *, double, double, double *);

        abs_starti = 1;                         /* initial values */
        pct_starti = 1;
        thr_starti = 1;
        gac_starti = 1;
        tds_starti = 1;

        /* Get the CCD halo factor (haloterm), for opt_elem G750L or G750M. */
        cti->allocated = 0;
        if ((status = getHalo (sts,
                out->npts, out->wave, out->net, cti)) != 0) {
            return (status);
        }

        photfactor = H_PLANCK * C_LIGHT / HST_AREA;

        if (sts->heliocorr == PERFORM)
            hf = sts->hfactor;
        else
            hf = 1.0;

        /* Generate time-dependent sensitivity correction factors. */
        if (sts->tdscorr == PERFORM) {
            tds_factors = (double *) malloc (tds->nwl * sizeof (double));
            if (tds_factors == NULL)
                return (OUT_OF_MEMORY);

            TdsCorrection (tds, sts->expstart, sts->detector_temp,
                           tds_factors);
        }

        /* Apply MSM/blaze correction to throughput array. Dispersion
           is required in the case the reference values from the _pht
           table are invalid, AND we are passing a blazeshift value
           thru the command line. In this case, we use the dispersion
           at the mid point in the current spectral order.

           First we have to save the phot->wl array, because when we
           interpolate to get pct_factor we must use wavelengths that
           do not include the blaze shift correction.
        */
        wl_no_blaze = malloc (phot->nelem * sizeof (double));
        if (wl_no_blaze == NULL)
            return (OUT_OF_MEMORY);
        for (i = 0;  i < phot->nelem;  i++)
            wl_no_blaze[i] = phot->wl[i];

        if (phot->wpos != 0.0 || sts->blazeshift != NO_VALUE) {

            i = out->npts / 2;
            dispersion = getDispersion (out, i);
            dispersion /= hf;

            BlazeCorr (phot, sts->expstart, sts->expend, out->sporder,
                       sts->blazeshift, blazeshift, dispersion);
        }

        /* Loop over flux array. */
        for (i = 0; i < out->npts; i++) {

            /* Convert back to observed wavelength. */
            wl = out->wave[i] / hfactor;

            /* Interpolate in the reference tables. */
            response    = extrap1d (wl, phot->wl, phot->thru, phot->nelem,
                                    &abs_starti);
            pct_factor  = interp1d (wl, wl_no_blaze, phot->pcorr, phot->nelem,
                                    &pct_starti);
            pct_factor2 = interp1d (wl, photc->wl, photc->pcorr, photc->nelem,
                                    &pct_starti);
            throughput  = interp1d (wl, slit->wl, slit->thr, slit->nelem,
                                    &thr_starti);
            if (sts->gaccorr == PERFORM) {
                gac_factor = interp1d (wl, slit->gac_wl, slit->gac_thr,
                                slit->gac_nelem, &gac_starti);
            } else {
                gac_factor = 1.;
            }
            if (sts->tdscorr == PERFORM)
                tds_factor  = interp1d (wl, tds->wl, tds_factors,
                                        tds->nwl, &tds_starti);
            else
                tds_factor = 1.;

            /* Check for zeroed throughput (OPR 38749) */
            if (throughput <= 0.0) {
                out->flux[i]  = 0.0F;
                out->error[i] = 0.0F;
                out->dq[i]   |= CALIBDEFECT;
                continue;
            }

            /* Compute dispersion at current pixel. */

            dispersion = getDispersion (out, i);

            if (dispersion <= 0.0) {
                out->flux[i]  = 0.0F;
                out->error[i] = 0.0F;
                out->dq[i]   |= CALIBDEFECT;
                continue;
            }
            dispersion /= hf;

            /* Compute flux. */

            if (response <= 0.0) {
                out->flux[i]  = 0.0F;
                out->error[i] = 0.0F;
                out->dq[i]   |= CALIBDEFECT;
            } else {
                float qfactor = 1.0;
                float net   = out->net[i];

                correction = (float) (photfactor / (response * throughput *
                                      wl * dispersion * CM_PER_ANGSTROM));
                correction *= atodgain;  /* added 9/3/97 (IB) */
                correction *= pct_factor / pct_factor2;
                correction /= gac_factor;
                correction /= tds_factor;
                ecorrection = correction;

                if (sts->detector == CCD_DETECTOR) {
                    qfactor = QFactor(sts, wl, out->extrsize, out->gross[i]);
                    if (sts->ctecorr == PERFORM) {
                        if (out->gross[i] > 0.0) {
                            net = out->gross[i] *
                                CtiCorr(sts, cti, cti->haloterm[i],
                                        out->extrsize, out->extrlocy[i],
                                        out->gross[i], out->back[i])
                                - out->back[i];
                        } else {
                            out->dq[i] |= NOT_CTI_CORR;
                        }
                    }
                }
                if (optimal)
                    correction *= pfactor[i];
                if (ecorrection <= 0.) {
                    out->flux[i]  = 0.0F;
                    out->error[i] = 0.0F;
                    out->dq[i]   |= CALIBDEFECT;
                } else {
                    out->flux[i] = net * correction;
                    out->error[i] *= ecorrection * qfactor;
                    out->net_error[i] *= qfactor;
                }
            }
        }

        if (sts->tdscorr == PERFORM)
            free (tds_factors);

        if (cti->allocated)
            free (cti->haloterm);

        free (wl_no_blaze);

        return (0);
}


static double getDispersion (RowContents *out, int index) {

        double dispersion;

        if (index > 0 && index < out->npts - 1)
            dispersion = (out->wave[index+1] - out->wave[index-1]) / 2.0;
        else if (index == 0)
            dispersion = out->wave[1] - out->wave[0];
        else
            dispersion = out->wave[out->npts-1] - out->wave[out->npts-2];

        return (dispersion);
}

static int getHalo (StisInfo6 *sts,
                int nelem, double wavelength[], float net[],
                CTICorrInfo *cti) {

/* arguments:
StisInfo6 *sts          i: calibration switches and info
int nelem               i: size of arrays
double wavelength[]     i: wavelength at each pixel
float net[]             i: net count rate
CTICorrInfo *cti        io: nelem and haloterm will be assigned
*/

        PhotInfo phot_200;      /* PCT info for extrheight = 200 */
        double extrheight = 200.;       /* 200 pixel extraction height */
        double net_to_electrons;
        int i;
        int warn = WARN;        /* yes, print a warning if appropriate */
        int status;
        void FreePhot6 (PhotInfo *);
        int GetPCT6 (StisInfo6 *, PhotInfo *, double, int);

        cti->nelem = nelem;
        cti->haloterm = calloc (nelem, sizeof (double));
        if (cti->haloterm == NULL)
            return (OUT_OF_MEMORY);
        cti->allocated = 1;

        if (sts->ctecorr != PERFORM)
            return (0);

        if (strcmp (sts->opt_elem, "G750L") != 0 &&
            strcmp (sts->opt_elem, "G750M") != 0) {
            return (0);
        }

        /* wl, thru, and error must be allocated, and the wl array must
           have the wavelengths we want.
        */
        phot_200.nelem = nelem;
        phot_200.wl    = malloc (nelem * sizeof (double));
        phot_200.thru  = malloc (nelem * sizeof (double));
        phot_200.error = malloc (nelem * sizeof (double));
        phot_200.pcorr = NULL;          /* allocated in GetPCT6 */
        phot_200.allocated = 1;
        if (phot_200.wl == NULL || phot_200.thru == NULL ||
            phot_200.error == NULL) {
            return (OUT_OF_MEMORY);
        }
        /* copy the wavelengths, converting back to observed wavelength */
        for (i = 0;  i < nelem;  i++)
            phot_200.wl[i] = wavelength[i] / sts->hfactor;

        /* This function will get the PC factors and interpolate them onto
           the wavelength scale we just assigned in phot_200.
        */
        if ((status = GetPCT6 (sts, &phot_200, extrheight, warn)) != 0) {
            FreePhot6 (&phot_200);
            return (status);
        }

        net_to_electrons = sts->atodgain * sts->exptime / sts->ncombine;
        for (i = 0;  i < nelem;  i++) {
            cti->haloterm[i] =
                ((phot_200.pcorr[i] - 1.) / 2. - cti->halominfrac) *
                                net[i] * net_to_electrons;
            if (cti->haloterm[i] < 0.)
                cti->haloterm[i] = 0.;
            else
                cti->haloterm[i] *= cti->halofac;
        }

        FreePhot6 (&phot_200);

        return (0);
}


static float CtiCorr(StisInfo6 *sts, CTICorrInfo *cti, double haloterm,
                float width, float y, float gross, float back)
{
    /*
       Arguments:

       StisInfo6 *sts         i:  calibration switches and info
       CtiCorrInfo *cti       i:  CTI correction info
       double haloterm        i:  halo correction term for current wavelength
       float width            i:  extraction width in cross dispersion
                                  direction
       float y                i:  extraction location in cross dispersion
                                  direction
       float gross            i:  input uncorrected gross count rate, but
                                  will be converted to electrons
       float back             i:  input background count rate, but will be
                                  converted to electrons
    */

    double to_elec;                     /* electron conversion factor */
    double gross_p;                     /* gross counts/row */
    double back_p;                      /* background counts/pixel */
    double cti_corr;                    /* charge transfer inefficiency */

    to_elec = sts->atodgain/sts->ncombine;
    /*
       Calculate gross counts per extraction row (in electrons).
       Include dark counts and spurious charge, which are per pixel;
       so multiply by extraction width.
    */
    gross_p = (gross*sts->exptime + sts->meandark*width)*to_elec +
        cti->spurcharge*width;
    /*
       Calculate background counts per pixel (in electrons)!
       Include dark counts and spurious charge.
    */
    back_p  = (back*sts->exptime/width + sts->meandark)*to_elec +
        cti->spurcharge;

    cti_corr = cti->ctinorm * pow(gross_p, -cti->ctigpower) *
        (1. + cti->ctitimfac*(sts->expstart - cti->ctirefmjd)/DAYS_PER_YEAR);

    if (back_p > 0.0)
        cti_corr *= exp(cti->ctibgfac *
                pow((back_p + haloterm) / gross_p, cti->ctibgpower));

    return (float) (1./pow(1. - cti_corr,
                CCD_NPIX_Y - (y - sts->ltv[1]) / sts->ltm[1]));
}


static float QFactor(StisInfo6 *sts, double wl, float width, float gross) {

    /*

      The Q factor is a correction to the error array when more than
      one electron is created per photon.  This effect occurs in the
      CCD at wavelengths shorter than ~3400 Angstroms.  The flux array
      already takes this effect into account via the sensitivity
      array.

                           (S*t + D*w)*g/n + R^2*w
      Q_correction = sqrt(-------------------------)
                          (S*t/Q + D*w)*g/n + R^2*w

      where S are the gross counts, D the mean dark background counts
      per pixel, R the read noise counts per pixel, G the gain, N the
      number of combined exposures, t the exposure time, and w the
      width of the extraction region.

    */

    double Q_wave[]  = {1127.09, 2066.33, 2213.93, 2339.25, 3424.86, 11000.00};
    double Q_value[] = {   3.00,    1.50,    1.35,    1.40,    1.00,     1.00};
    int    Q_num   = 6;
    int    wl_start = 0;

    double S = gross > 0? gross: 0.0;   /* gross counts/row */
    double D = sts->meandark;           /* mean dark counts per pixel */
    double R = sts->readnse;            /* read noise counts per pixel */
    double g = sts->atodgain;           /* electron conversion factor */
    double n = sts->ncombine;           /* number of combined exposures */
    double t = sts->exptime;            /* exposure time */
    double w = width;                   /* width of extraction region */
    double Q = interp1d(wl, Q_wave, Q_value, Q_num, &wl_start);

    return (float) sqrt(((S*t + D*w)*g/n + R*R*w)/((S*t/Q + D*w)*g/n + R*R*w));
}
