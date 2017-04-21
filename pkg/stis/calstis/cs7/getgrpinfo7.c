# include <stdio.h>
# include <stdlib.h>
# include <string.h>		/* strstr */

# include "c_iraf.h"
# include "hstio.h"

# include "stis.h"
# include "calstis7.h"
# include "err.h"
# include "stisdef.h"

/* This routine gets keyword values from an extension header.

   Phil Hodge, 1998 Oct 5:
	Change status value 1001 to GENERIC_ERROR_CODE.

   Phil Hodge, 2000 Aug 4:
	Get CRPIXi for spectroscopic type as well as imaging.

   Phil Hodge, 2004 Dec 27:
	Get detector_temp.
*/

int GetGrpInfo7 (StisInfo7 *sts, Hdr *hdr) {

/* arguments:
StisInfo7 *sts  io: calibration switches and info
Hdr *hdr        i: header of current imset
*/

	int status;

	int sdqflags;			/* serious data quality flags */
	char *buf;			/* scratch for keyword value */
	int use_def = 1;		/* use default if missing keyword */
	int no_default = 0;		/* missing keyword is fatal error */
	int GetDetTemp (Hdr *, int, double *);

	/* Get generic parameters. */

	/* This should be SCIENCE or WAVECAL. */
	if ((buf = calloc (STIS_FNAME+1, sizeof(char))) == NULL)
	    return (OUT_OF_MEMORY);
	if ((status = Get_KeyS (hdr, "ASN_MTYP",
                                use_def, "unknown", buf, STIS_FNAME)))
	    return (status);
	/* Possible values for wavecals are "AUTO-WAVECAL" and "GO-WAVECAL" */
	sts->wavecal = (strstr (buf, "WAVECAL") != NULL);
	free (buf);

	if ((status = Get_KeyD (hdr, "EXPTIME", no_default, 0., &sts->exptime)))
	    return (status);
	if (sts->exptime < 0.) {
	    printf ("ERROR    Exposure time %.6g is invalid.\n",
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

	/* Get the linear transformation (zero indexed) from the reference
	   coordinate system to current image pixel coordinates.
	*/
	if ((status = GetLT0 (hdr, sts->ltm, sts->ltv)))
	    return (status);

	if ((status = Get_KeyD (hdr, "CRPIX1", use_def, 0., &sts->crpix[0])))
	    return (status);
	if ((status = Get_KeyD (hdr, "CRPIX2", use_def, 0., &sts->crpix[1])))
	    return (status);
	sts->crpix[0]--;	/* convert to zero-indexed */
	sts->crpix[1]--;

	/* Get the detector temperature (or housing temperature, if CCD). */
	if (sts->fluxcorr == PERFORM) {
	    if ((status = GetDetTemp (hdr, sts->detector,
			&sts->detector_temp))) {
		return (status);
	    }
	} else {
	    sts->detector_temp = -1.;
	}

	/* Get the MSM offset found during wavecal processing. */
	if (sts->wavecal) {

	    sts->msm_slop[0] = 0.;
	    sts->msm_slop[1] = 0.;

	} else {

	    if ((status = Get_KeyD (hdr, "SHIFTA1",
                                    use_def, 0., &sts->msm_slop[0])))
		return (status);

	    if (sts->wx2dcorr == COMPLETE) {
		/* shifta2 has already been accounted for by wx2d */
		sts->msm_slop[1] = 0.;
	    } else {
		if ((status = Get_KeyD (hdr, "SHIFTA2",
                                        use_def, 0., &sts->msm_slop[1])))
		    return (status);
	    }
	}

	if (sts->obstype == SPECTROSCOPIC_TYPE) {

	    /* Get the dispersion axis. */
	    if ((status = Get_KeyI (hdr, "DISPAXIS",
                                    use_def, 1, &sts->dispaxis)))
		return (status);
	    if (sts->dispaxis < 1 || sts->dispaxis > 2) {
		printf ("ERROR    Dispaxis = %d is invalid.\n", sts->dispaxis);
		return (GENERIC_ERROR_CODE);
	    }

	} else {

	    sts->dispaxis = 0;
	}

	return (0);
}
