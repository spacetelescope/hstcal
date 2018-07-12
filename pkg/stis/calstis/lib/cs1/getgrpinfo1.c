/* This file contains:
	GetGrpInfo1
	UpdateCCDTemperature
*/

# include <stdio.h>
# include <string.h>
# include <math.h>	/* sqrt */

# include "c_iraf.h"
# include "hstio.h"

# include "stis.h"
# include "calstis1.h"
# include "hstcalerr.h"
# include "stistemperature.h"
# include "stisdef.h"

static int UpdateCCDTemperature(StisInfo1 *, Hdr *);

/* This routine gets keyword values from the SCI extension header.

   Phil Hodge, 1998 May 4:
	Get NCOMBINE from the SCI extension header.  Previously we were
	getting it from the primary header.

   Phil Hodge, 1998 June 8:
	Get CD matrix and compute the pixel size in each axis.

   Phil Hodge, 1998 July 30:
	Check value of ASN_MTYP to see if the input file is a wavecal,
	and set sts->wavecal instead of sts->assoc_typ.

   Phil Hodge, 1998 Oct 5:
	Change status value 1001 to GENERIC_ERROR_CODE.

   Phil Hodge, 1999 Mar 29:
	Move GetTemperature from getkeyinfo1.c to this file.
	In GetTemperature, read the [UDL,1] extension (instead of the
	primary header) of the support file to get OM2CAT; use MIN_OM2CAT
	as a minimum valid value for om2cat; if there's a problem, print a
	warning regardless of the verbose switch.

   Phil Hodge, 1999 July 26:
	Initialize the temperature to -1.

   Ivo Busko, 2001 Oct 26:
	Get CCD housing temperature from OCCDHTAV..

   Ivo Busko, 2002 Apr 03:
        Get EXPSTART if NUV-MAMA and darkcor = perform.

   Phil Hodge, 2003 Jan 20:
	Get EXPSTART and EXPEND.
	Don't reset sts->wavecal if asn_mtyp doesn't contain the string
	"WAVECAL".  This is because the current observation might not be
	part of an association, and the wavecal flag may have already been
	set by GetKeyInfo1 based on targname.

   Paul Barrett, 2003 Sep 18:
        Get CCD housing temperature from Engineering Parameter Calibration
        (EPC) file and calculate a time-averaged temperature.  Update
        header keyword with new value. If EPC data is not available, use
        value from OCCDHTAV keyword.

   Phil Hodge, 2004 Dec 27:
	Use detector_temp instead of temperature and ccd_temperature.
	Delete GetTemperature (replaced by GetDetTemp).
	Rename GetCCDTemperature to UpdateCCDTemperature, and modify it
	so that it only updates OCCDHTAV (from the EPC table), then get
	the temperature afterwards by calling GetDetTemp.
*/

int GetGrpInfo1 (StisInfo1 *sts, Hdr *hdr) {

/* arguments:
StisInfo1 *sts   io: calibration switches and info
Hdr *hdr         i: header of current extension
*/

	int status;

	char *buf;			/* scratch for keyword value */
	int sdqflags;			/* serious data quality flags */
	int rsize;			/* 1 for CCD, 2 for MAMA */
	int corner[2];		/* subarray start points, detector coords */
	double cd11, cd12, cd21, cd22;	/* CD matrix */
	int doppon;			/* Doppler correction done on-board? */
	int use_def = 1;		/* use default if missing keyword */
	int no_default = 0;		/* missing keyword is fatal error */
	int GetDetTemp (Hdr *, int, double *);
	int GetEPCTab (StisInfo1 *, float);

	/* Get generic parameters. */

	/* Check whether we're processing a science file or wavecal,
	   based on ASN_MTYP.  Note that we won't reset the wavecal flag
	   to false (it may have been set in GetKeyInfo1) if ASN_MTYP
	   doesn't indicate that the observation is a wavecal, because
	   this file might not be part of an association.
	*/
	if ((buf = calloc (STIS_FNAME+1, sizeof(char))) == NULL)
	    return (OUT_OF_MEMORY);
	if ((status = Get_KeyS (hdr, "ASN_MTYP",
                                use_def, "unknown", buf, STIS_FNAME)))
	    return (status);
	/* Possible values for wavecals are "AUTO-WAVECAL" and "GO-WAVECAL" */
	if (strstr (buf, "WAVECAL") != NULL)
	    sts->wavecal = 1;
	free (buf);

	if ((status = Get_KeyD (hdr, "EXPTIME", no_default, 0., &sts->exptime)))
	    return (status);
	if (sts->exptime < 0.) {
	    printf ("ERROR    Exposure time is invalid:  %14.6g.\n",
		sts->exptime);
	    return (GENERIC_ERROR_CODE);
	}
	if ((status = Get_KeyD (hdr, "EXPSTART", no_default, 0., &sts->expstart)))
	    return (status);
	if ((status = Get_KeyD (hdr, "EXPEND", no_default, 0., &sts->expend)))
	    return (status);

	/* Find out which data quality bits are considered serious;
	   default value means all bits are serious.
	*/
	if ((status = Get_KeyI (hdr, "SDQFLAGS", use_def, 32767, &sdqflags)))
	    return (status);
	sts->sdqflags = (short) sdqflags;

	/* Get the pixel size (ignore corner location) from ltm & ltv. */
	rsize = (sts->detector == CCD_DETECTOR) ? 1 : 2;
	if ((status = GetCorner (hdr, rsize, sts->bin, corner)))
	    return (status);

	/* For spectroscopic data, we want the dispersion axis and the
	   sign of the dispersion.  We'll get the latter from one element
	   of the CD matrix.
	   We also want the pixel size, which we compute from the CD matrix.
	*/
	sts->dispaxis = 1;		/* initial values */
	sts->dispsign = 1;

	if (strcmp (sts->obstype, "IMAGING") == 0) {

	    if ((status = Get_KeyD (hdr, "CD1_1", use_def, 1., &cd11)))
		return (status);
	    if ((status = Get_KeyD (hdr, "CD1_2", use_def, 0., &cd12)))
		return (status);
	    if ((status = Get_KeyD (hdr, "CD2_1", use_def, 0., &cd21)))
		return (status);
	    if ((status = Get_KeyD (hdr, "CD2_2", use_def, 1., &cd22)))
		return (status);
	    sts->cdelt[0] = sqrt (cd11 * cd11 + cd21 * cd21);
	    sts->cdelt[1] = sqrt (cd12 * cd12 + cd22 * cd22);

	} else if (strcmp (sts->obstype, "SPECTROSCOPIC") == 0) {

	    if ((status = Get_KeyI (hdr, "DISPAXIS", use_def, 1, &sts->dispaxis)))
		return (status);
	    if ((status = Get_KeyD (hdr, "CD1_1", use_def, 1., &cd11)))
		return (status);
	    if ((status = Get_KeyD (hdr, "CD2_2", use_def, 1., &cd22)))
		return (status);
	    if (sts->dispaxis == 1) {
		if (cd11 >= 0.)
		    sts->dispsign = 1;
		else
		    sts->dispsign = -1;
		sts->cdelt[0] = cd22;		/* assume square pixels */
		sts->cdelt[1] = cd22;
	    } else if (sts->dispaxis == 2) {
		if (cd22 >= 0.)
		    sts->dispsign = 1;
		else
		    sts->dispsign = -1;
		sts->cdelt[0] = cd11;
		sts->cdelt[1] = cd11;
	    }
	}

	sts->detector_temp = -1.;	/* initial value (not defined) */

	/* Get MAMA-specific parameters. */

	if (sts->detector == NUV_MAMA_DETECTOR ||
	    sts->detector == FUV_MAMA_DETECTOR) {

	    if ((status = Get_KeyD (hdr, "GLOBRATE",
                                    no_default, 0., &sts->globrate)))
		return (status);

	    /* Get info if we need to do Doppler convolution of ref files. */

	    if (sts->doppcorr == PERFORM) {

		/* Was Doppler correction done on-board? */
		if ((status = Get_KeyI (hdr, "DOPPON", use_def, 0, &doppon)))
		    return (status);

		/* doppon could be False in timetag mode. */
		if (!doppon) {
		    if (strcmp (sts->obsmode, "TIME-TAG") == 0) {
			if ((status = Get_KeyD (hdr, "DOPPMAG", use_def, -1.,
                                                &sts->doppmag)))
			    return (status);
			doppon = (sts->doppmag > 0.);
		    }
		}

		if (doppon) {
		    if ((status = Get_KeyD (hdr, "DOPPZERO", no_default, 0.,
                                            &sts->doppzero)))
			return (status);
		    if ((status = Get_KeyD (hdr, "DOPPMAG", no_default, 0.,
                                            &sts->doppmag)))
			return (status);
		    if ((status = Get_KeyD (hdr, "ORBITPER", no_default, 0.,
                                            &sts->orbitper)))
			return (status);
		} else {
		    /* Silently reset the switch. */
		    sts->doppcorr = OMIT;
		}
	    }
	}

	/* Get CCD-specific parameters. */

        if (sts->detector == CCD_DETECTOR) {
	    float occdhtav = -1.;

	    /*  if OCCDHTAV is present and > 0, then data is from Side B,
	     *  so look for EPC file.  Otherwise, data is from Side A and
	     *  EPC file and temperature can be ignored.
	     */

	    if (((status = Get_KeyF(hdr, "OCCDHTAV", use_def, -1.,
				    &occdhtav)) == 0) && (occdhtav >= 0.)) {
		char expname[81], epcname[81];
		expname[0] = epcname[0] = '\0';
		if ((status = Get_KeyS(hdr, "EXPNAME", no_default, "",
                                       expname, 80)))
		    return (status);
		if (sts->crcorr != COMPLETE) {
		    strncpy(epcname, expname, 8); epcname[8] = '\0';
		    sprintf(sts->epctab.name, "%sj_epc.fits", epcname);
		    PrFileName("epcfile", sts->epctab.name);
		    status = GetEPCTab(sts, 0.40);
		    if (status && status != OPEN_FAILED)
			return (status);
		}
		if ((status = UpdateCCDTemperature(sts, hdr)))
		    return (status);
	    }
        }

	/* Get the detector temperature (or housing temperature, if CCD).
	   Note that for side-2 CCD data this function must be called after
	   calling UpdateCCDTemperature.
	*/
	if ((status = GetDetTemp (hdr, sts->detector, &sts->detector_temp))) {
	    return (status);
	}

	/* If images have been combined (e.g. by cosmic-ray rejection),
	   then determine the number of images that were combined together;
	   we need this for bias image subtraction.
	   (This isn't really a CCD-specific keyword, but it does only affect
	   CCD data in the context of calstis1.)
	*/
	if ((status = Get_KeyI (hdr, "NCOMBINE", use_def, 1, &sts->ncombine)))
	    return (status);

	if (sts->ncombine < 1) {
	    printf ("Warning  NCOMBINE = %d, reset to one.\n", sts->ncombine);
	    sts->ncombine = 1;
	}

	return (0);
}

/*
  This function updates the OCCDHTAV keyword, which gives the CCD housing
  temperature that is used for accurate calibration of the side 2 dark
  images.  A time-averaged temperature is calculated using data from the
  Engineering Parameter Calibration file, if this file is available.
  If cosmic-ray rejection has already been done, OCCDHTAV is assumed
  to have been updated previously (during the first call to CalStis1),
  and in that case this function returns without doing anything.
*/

int UpdateCCDTemperature(StisInfo1 *sts, Hdr *hdr) {

    /* arguments:
       StisInfo1 *sts  io: calibration switches and info
       Hdr *hdr i: header of current extension */

    int status = 0;
    int row;
    int row_start = -1;            /*  beginning row of EPC table data  */
    int row_end = -1;              /*  ending row of EPC table data  */
    double ccd_temperature = -1.;  /*  average from EPC table  */

    /*  If CRCORR is complete we don't need to update OCCDHTAV */
    if (sts->crcorr == COMPLETE)
	return (0);

    /*  Finding starting and ending row of EPC table data for this
        image extension ONLY if the EPC data bounds the exposure data.
    */
    for (row = 0; row < sts->epc_rows; row++) {
        if (sts->epc_mjd[row] < sts->expstart)
            row_start = row;
        if (sts->epc_mjd[row] > sts->expend) {
            row_end = row;
            break;
        }
    }

    /*  Determine CCD temperature  */

    if (row_start >= 0 && sts->exptime == 0.) {
        ccd_temperature = sts->epc_temp[row_start];
        status = Put_KeyF(hdr, "OCCDHTAV", (float) ccd_temperature, "");
    }
    else if (row_start >= 0 && row_end >= 0) {
        /*  From time-averaged CCD temperatures in EPC table.  */
        double avtemp = 0.;

        for (row = row_start; row < row_end; row++) {
            double fmjd, ftemp;

            if (row_start == row_end-1) {
                fmjd  = sts->expend - sts->expstart;
                ftemp = sts->epc_temp[row];
            }
            else if (row == row_end-1) {
                fmjd  = sts->expend - sts->epc_mjd[row];
                ftemp = 0.5*(sts->epc_temp[row+1] - sts->epc_temp[row]) *
                    (sts->expend - sts->epc_mjd[row]) /
                    (sts->epc_mjd[row+1] - sts->epc_mjd[row]) +
                    sts->epc_temp[row];
            }
            else if (row == row_start) {
                fmjd  = sts->epc_mjd[row+1] - sts->expstart;
                ftemp = 0.5*(sts->epc_temp[row+1] - sts->epc_temp[row]) *
                    (sts->expstart - sts->epc_mjd[row]) /
                    (sts->epc_mjd[row+1] - sts->epc_mjd[row]) +
                    sts->epc_temp[row];
            }
            else {
                fmjd  = sts->epc_mjd[row+1] - sts->epc_mjd[row];
                ftemp = 0.5*(sts->epc_temp[row+1] + sts->epc_temp[row]);
            }
            avtemp += fmjd*ftemp;
        }
        ccd_temperature = avtemp/(sts->expend-sts->expstart);

        /*  Update the header keyword with the average temperature  */

        status = Put_KeyF(hdr, "OCCDHTAV", (float) ccd_temperature, "");
    }
    return (status);
}
