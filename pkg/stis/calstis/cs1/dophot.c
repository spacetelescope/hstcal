# include <stdio.h>
# include <stdlib.h>	/* malloc */

# include <c_iraf.h>
# include <hstio.h>
# include <xsynphot.h>		/* for c_phopar */
# include "../stis.h"
# include "calstis1.h"
# include "../stiserr.h"
# include "../stisdef.h"
# include "../stistds.h"

# define NPHOT      4		/* size of phot returned by c_phopar */

static void MultFilter (PhotInfo *);

/* This routine gets the absolute flux conversion from PHOTTAB and the
   aperture throughput (filter throughput, really) from APERTAB.  These
   are multiplied together, and then the synphot routine phopar is called
   to determine the inverse sensitivity, reference magnitude (actually a
   constant), pivot wavelength, and RMS bandwidth.  These are written to
   keywords in the primary header.

   Phil Hodge, 1997 Nov 13:
	Phot table I/O extracted to GetPhot1; also call GetApThr1;
	rename phot photkey and use phot for PhotInfo.

   Phil Hodge, 1998 Oct 5:
	Change status value 1001 to GENERIC_ERROR_CODE.

   Paul Barrett, 2003 Sep 25:
        Add time-dependent sensitivity (TDS) for imaging mode.
*/

int doPhot (StisInfo1 *sts, SingleGroup *x) {

/* arguments:
StisInfo1 *sts     i: calibration switches, etc
SingleGroup *x    io: image to be calibrated; primary header is modified
*/

	int status;

	PhotInfo phot;		/* QE and filter throughput, vs wavelength */

	float photkey[NPHOT];	/* values of photometry keywords */

	int GetPhot1 (StisInfo1 *, PhotInfo *);
	int GetApThr1 (StisInfo1 *, PhotInfo *);
        int GetTds1 (StisInfo1 *, PhotInfo *);

	/* Get values from the photometry table. */
	if (status = GetPhot1 (sts, &phot))
	    return (status);
	if (sts->photcorr == DUMMY)
	    return (0);

	/* Get values from the aperture throughput table. */
	if (sts->filtcorr == PERFORM) {		/* hasn't been reset? */
	    if (status = GetApThr1 (sts, &phot))
		return (status);
	}

	/* Combine aperture throughput with photometry table info. */
	if (sts->filtcorr == PERFORM)		/* still hasn't been reset? */
	    MultFilter (&phot);

        /* Get values from the TDS table */
        if (sts->tdscorr == PERFORM) {
            if (status = GetTds1(sts, &phot))
                return (status);
        }

        /* Combine TDS value with photometry table info. */
        if (sts->tdscorr == PERFORM) {
            MultFilter(&phot);
        }

	/* Compute photometry values. */
	c_phopar (phot.p_nelem, phot.p_wl, phot.p_thru, photkey);
	if (c_iraferr()) {
	    printf ("ERROR    Error return from c_phopar\n");
	    return (GENERIC_ERROR_CODE);
	}

	/* If the detector is the CCD, multiply PHOTFLAM by the gain. */
	if (sts->detector == CCD_DETECTOR)
	    photkey[0] *= sts->atodgain;

	/* Update the photometry keyword values in the primary header. */

	if (status = Put_KeyF (x->globalhdr, "PHOTFLAM", photkey[0],
			"inverse sensitivity"))
	    return (status);

	if (status = Put_KeyF (x->globalhdr, "PHOTZPT", photkey[1],
			"zero point"))
	    return (status);

	if (status = Put_KeyF (x->globalhdr, "PHOTPLAM", photkey[2],
			"pivot wavelength"))
	    return (status);

	if (status = Put_KeyF (x->globalhdr, "PHOTBW", photkey[3],
			"RMS bandwidth"))
	    return (status);

	free (phot.p_wl);
	free (phot.p_thru);
	if (sts->filtcorr == PERFORM || sts->tdscorr == PERFORM) {
	    free (phot.f_wl);
	    free (phot.f_thru);
	}

	return (0);
}

/* This routine interpolates the aperture throughput (APERTAB) to the
   same wavelengths as the QE (PHOTTAB), and then multiplies the
   QE throughput in-place by the aperture throughput.
*/

static void MultFilter (PhotInfo *phot) {

	double wl;			/* a wavelength in QE array */
	double filt_throughput;		/* interpolated filter throughput */
	int filt_starti;		/* index to begin search in interp1d */
	int i;

	filt_starti = 1;			/* initial value */

	for (i = 0;  i < phot->p_nelem;  i++) {

	    wl = phot->p_wl[i];

	    filt_throughput = interp1d (wl, phot->f_wl, phot->f_thru,
			phot->f_nelem, &filt_starti);

	    phot->p_thru[i] *= filt_throughput;
	}
}
