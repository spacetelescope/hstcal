# include <stdio.h>
# include <stdlib.h>
# include <string.h>		/* strcmp */

# include "hstio.h"
# include "xtables.h"

# include "stis.h"
# include "stisdef.h"
# include "calstis6.h"
# include "hstcalerr.h"
# include "stisdq.h"

/*
   Get keyword values from an extension header.




   Revision history:
   ----------------
   20 Feb 97  -  Adapted from similar routine in calstis7 (I.Busko)
   09 Apr 97  -  Changes after code review (IB)
   10 Apr 97  -  Changes after code review (IB):
                 - replaced literals by INVALID constant.
                 - explicit cast for calloc-returned pointer.
   08 May 97  -  Conform to new _trl standard (IB)
   06 Apr 98  -  Removed MOFFSET input (IB)
   29 Sep 98  -  Set sdqflags = 0 (IB)
   17 Nov 98  -  Store sdqflags in separate StisInfo6 structure element (IB)
   17 Dec 98  -  Get CRVAL1, CRPIX1 for geocoronal Lya detection (IB)
   17 Dec 98  -  Separate sdqflags value for crosscorr function (IB)
   17 Nov 00  -  Ignore smallblemish and hotpix sdqflags (IB)
   16 Jun 03  -  Get MEANDARK for CTI correction (PB)
   27 Dec 04  -  Get detector_temp (PEH)
   10 Apr 06  -  Get Doppler shift info, needed for blaze shift corr. (PEH)
*/

int GetGrpInfo6 (StisInfo6 *sts, Hdr *hdr) {

/* arguments:
StisInfo6 *sts  io: calibration switches and info
Hdr *hdr        i: header of current IMSET
*/

	int status;

	int sdqflags;			/* serious data quality flags */
	char *buf;			/* scratch for keyword value */
	int use_def = 1;		/* use default if missing keyword */
	int no_default = 0;		/* missing keyword is fatal error */
	int GetDetTemp (Hdr *, int, double *);

	/* Get the dispersion axis. */
	if ((status = Get_KeyI (hdr, "DISPAXIS", use_def, 1, &sts->dispaxis)))
	    return (status);
	if (sts->dispaxis < 1 || sts->dispaxis > 2) {
	    printf ("ERROR    Dispaxis = %d is invalid\n", sts->dispaxis);
	    return (INVALID);
	}

	/* This should be SCIENCE or WAVECAL. */
	if ((buf = (char *) calloc (STIS_FNAME+1, sizeof(char))) == NULL)
	    return (OUT_OF_MEMORY);
	if ((status = Get_KeyS (hdr, "ASN_MTYP",
                                use_def, "unknown", buf, STIS_FNAME)))
	    return (status);
	sts->wavecal = (strcmp (buf, "WAVECAL") == 0);
	free (buf);

	/* Exposure info. */
	if ((status = Get_KeyD (hdr, "EXPTIME", no_default, 0., &sts->exptime)))
	    return (status);
	if (sts->exptime < 0.) {
	    printf ("ERROR    Exposure time is invalid:  %14.6g\n",
		sts->exptime);
	    return (INVALID);
	}
	if ((status = Get_KeyD (hdr, "EXPSTART", no_default, 0., &sts->expstart)))
	    return (status);
	if ((status = Get_KeyD (hdr, "EXPEND", no_default, 0., &sts->expend)))
	    return (status);
	if ((status = Get_KeyD (hdr, "DOPPMAG", use_def, 0., &sts->doppmag)))
	    return (status);
	if (sts->doppmag > 0. && strcmp (sts->obsmode, "TIME-TAG") != 0) {
	    /* doppon could be False in timetag mode; otherwise, reset
		doppzero to zero if on-board Doppler correction was not done.
	    */
	    int doppon;
	    if ((status = Get_KeyI (hdr, "DOPPON", use_def, 0, &doppon)))
		return (status);
	    if (!doppon)
		sts->doppmag = 0.;
	}
	if ((status = Get_KeyD (hdr, "DOPPZERO", use_def, 0., &sts->doppzero)))
	    return (status);
	if ((status = Get_KeyD (hdr, "ORBITPER", use_def, 1., &sts->orbitper)))
	    return (status);

	/* Find out which data quality bits are considered serious;
	   default value means all bits are serious.
	*/
	if ((status = Get_KeyI (hdr, "SDQFLAGS", use_def, ALL_BITS, &sdqflags)))
	    return (status);
	sdqflags &= ~HOTPIX;		/* hot (but probably just warm) pixel */
	sdqflags &= ~SMALLBLEM;		/* too small to be concerned about */
	sts->sdqflags_orig = (short) sdqflags;
	sts->cc_sdqflags   = (short) sdqflags;

	/* This turns off the data quality flag checking everywhere
           except  in crosscorrelation module.
        */
	sts->sdqflags = 0;

	/* Get the "plate scale" coefficients. */
	if ((status = Get_KeyD (hdr, "CD1_1", no_default, 0., &sts->cd[0])))
	    return (status);
	if ((status = Get_KeyD (hdr, "CD2_2", no_default, 0., &sts->cd[1])))
	    return (status);
	if ((status = Get_KeyD (hdr, "CRPIX1", no_default, 1., &sts->crpix[0])))
	    return (status);
	if ((status = Get_KeyD (hdr, "CRVAL1", no_default, 0., &sts->crval[0])))
	    return (status);

	/* Get the linear transformation (zero indexed) from the reference
	   coordinate system to current image pixel coordinates.
	*/
	if ((status = GetLT0 (hdr, sts->ltm, sts->ltv)))
	    return (status);

	/* Get the detector temperature (or housing temperature, if CCD). */
	if (sts->fluxcorr == PERFORM) {
	    if ((status = GetDetTemp (hdr, sts->detector,
			&sts->detector_temp))) {
		return (status);
	    }
	} else {
	    sts->detector_temp = -1.;
	}

	/* Get the MSM slop. */
	if (sts->wavecal) {
	    sts->msm_offset[0] = 0.0;
	    sts->msm_offset[1] = 0.0;
	} else {
	    if ((status = Get_KeyD (hdr, "SHIFTA1",
                                    use_def, 0., &sts->msm_offset[0])))
		return (status);
	    if ((status = Get_KeyD (hdr, "SHIFTA2",
                                    use_def, 0., &sts->msm_offset[1])))
		return (status);
	}

        /* Get the mean dark count rate and the number of
           combined images for CTI correction. */
        if ((status = Get_KeyD(hdr, "MEANDARK", no_default, 0.,
                               &sts->meandark)))
            return status;
        if ((status = Get_KeyI(hdr, "NCOMBINE", no_default, 0.,
                               &sts->ncombine)))
            return status;

	return (0);
}
