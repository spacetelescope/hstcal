/* Do basic 2-D image reduction for the current set of extensions.
   If we're processing the first set (extver = 1), then the primary
   header will be updated with history information, and the calibration
   switches in the header will be reset from "PERFORM" to "COMPLETE".

   Phil Hodge, 1997 Oct 28:
	Change calling sequences of doBlev & blevHistory to include driftcorr.

   Phil Hodge, 1997 Nov 13:
	Move the call to PhotMsg;
	in PhotMsg, remove the call to PrSwitch and include apertab.

   Phil Hodge, 1998 Jan 16:
	Print "uncertainty" instead of "error" for doNoise message;
	in PhotMsg, check filtcorr before calling PrRefInfo for apertab.

   Phil Hodge, 1998 July 16:
	For MAMA detectors, don't print readnoise, etc, for error init.

   Phil Hodge, 1998 July 30:
	After calling doBias, check that biascorr hasn't been reset to SKIPPED.
	Reset CRCORR to SKIPPED rather than OMIT, if only one imset.

   Phil Hodge, 1998 Sept 24:
	Include UpdatePlateSc, and call it after doLoRes.

   Phil Hodge, 1998 Oct 5:
	Change status value 1001 to IO_ERROR.

   Phil Hodge, 1999 Nov 2:
	Do blevcorr before error array initialization (move the calls to
	doNoise and doBlev); assign a value to err_init_bias.
	Use SingleGroup *x and x = calloc ..., rather than SingleGroup x;
	pass &x to doBlev and doLoRes instead of a second SingleGroup.

   Phil Hodge, 2001, Aug 27:
	Move the call to CheckVolt from cs0/GetSciInfo and GetWavInfo
	to here.

   Ivo Busko, 2002 Feb 01:
        Disable call to CheckVolt as per OPR # 45138

   Ivo Busko, 2002 Mar 20:
        IMSET number passed to doDark to allow scale factors from command line.

   Phil Hodge, 2003 Jan 17:
	Compute heliocentric radial velocity, and update v_helio keyword.

   Phil Hodge, 2003 Feb 24:
	Include units for v_helio, in keyword comment and in trailer output.

   Paul Barrett, 2003 Sep 18:
        Print EPC table information.

   Paul Barrett, 2003 Sep 25:
        PrRefInfo to PhotMsg for TDS table.

   Paul Barrett, 2004 Jun 25:
        Removed EPC table information, moved to GetGrpInfo1.

   Phil Hodge, 2007 May 9:
	For a wavecal exposure or science exposure that is not a bias
	(targname != 'BIAS'), check for zero exposure time or constant data
	values, and in this case do not write the current imset to the output
	file.  Note that if all imsets are bad, the output file will not be
	written at all, so this case needs to be caught at an earlier stage
	of processing.  Add output_extver to the calling sequence.
	(The getMinMax function in this file is currently the same as in
	cs4/wavecal.c.)

   Phil Hodge, 2007 Oct 4:
	Add function dummyKw to add a few DUMMY keywords to a DQ header to
	work around a FITS kernel bug.

   Phil Hodge, 2008 Nov 3:
	Rename variable output_extver to ngood_extver.
	Add keyword IMSET_OK = F if an imset has zero EXPTIME (except for
	bias images) or the data values are constant, and write the imset
	but without doing any calibration.  Use extver instead of
	*output_extver in most places.

   Phil Hodge, 2011 May 9:
	Move the call to PhotMode to a point just before the call to doPhot.
	In PhotMsg, change keyword phottab to imphttab, and don't print info
	about apertab or tdstab.

   Phil Hodge, 2011 July 27:
	After calling doPhot, check sts->photcorr (may be SKIPPED).

   Phil Hodge, 2012 Oct 15:
	Add a call to PrRefInfo for the TDCTAB, for detector = NUV-MAMA.
	Delete function dummyKw, since it should no longer be necessary.
*/

# include <math.h>	/* for fabs and sqrt */
# include <string.h>
# include <stdio.h>

# include "c_iraf.h"
# include "hstio.h"
# include "stis.h"
# include "calstis1.h"
# include "err.h"
# include "stisdef.h"

static void getMinMax (SingleGroup *, float *, float *);
static void AtoDMsg (StisInfo1 *, int);
static void BiasMsg (StisInfo1 *, int);
static void BlevMsg (StisInfo1 *);
static void DarkMsg (StisInfo1 *, int);
static void DoppMsg (char *, int *, char *);
static void dqiMsg (StisInfo1 *, int);
static void FlatMsg (StisInfo1 *, int);
static void LoResMsg (StisInfo1 *);
static void NonLinMsg (StisInfo1 *, int);
static void PhotMsg (StisInfo1 *);
static void ShadMsg (StisInfo1 *, int);
static int UpdatePlateSc (StisInfo1 *, SingleGroup *);

int Do2D (StisInfo1 *sts, int extver, int *ngood_extver) {

/* arguments:
StisInfo1 *sts      i: calibration switches and info
int extver          i: "imset" number, the current set of extensions
int ngood_extver   io: incremented unless the current imset has zero
                       exposure time or constant pixel values
*/

	int status;

	SingleGroup *x;		/* used for both input and output */
	int option = 0;
	float meanblev;		/* mean value of overscan bias (for history) */
	float meandark;		/* mean value of dark (for history) */
	int crcorr;		/* for checking if CRCORR is inappropriate */
	int gsat, lsat;		/* > 0 if saturated pixels found */
	int driftcorr;	/* true means bias level was corrected for drift */
	int done;	/* true means the input SingleGroup has been freed */
	char *doppstr;		/* message regarding Doppler correction */
	int doppcount = 0;	/* Doppler convolution not applied yet */
	double v_helio;		/* radial velocity in direction of target */

	int doAtoD (StisInfo1 *, SingleGroup *);
	int atodHistory (StisInfo1 *, Hdr *);
	int doBias (StisInfo1 *, SingleGroup *);
	int biasHistory (StisInfo1 *, Hdr *);
	int doBlev (StisInfo1 *, SingleGroup **, int,
		float *, int *, int *);
	int blevHistory (StisInfo1 *, Hdr *, int, int);
	int CCDHistory (StisInfo1 *, Hdr *);
	int doDark (StisInfo1 *, SingleGroup *, float *, int);
	int darkHistory (StisInfo1 *, Hdr *);
	int doDQI (StisInfo1 *, SingleGroup *);
	int dqiHistory (StisInfo1 *, Hdr *);
	int doFlat (StisInfo1 *, SingleGroup *);
	int flatHistory (StisInfo1 *, Hdr *);
	int doLoRes (StisInfo1 *, SingleGroup **, int *);
	int loresHistory (StisInfo1 *, Hdr *, int);
	int doNonLin (StisInfo1 *, SingleGroup *, int *, int *);
	int nonlinHistory (StisInfo1 *, Hdr *);
	int doNoise (StisInfo1 *, SingleGroup *, int *);
	int noiseHistory (Hdr *);
	int doPhot (StisInfo1 *, SingleGroup *);
	int PhotMode (StisInfo1 *, SingleGroup *);
	int photHistory (StisInfo1 *, Hdr *);
	int doShad (StisInfo1 *, SingleGroup *);
	int shadHistory (StisInfo1 *, Hdr *);
	int GetGrpInfo1 (StisInfo1 *, Hdr *);

	if ((x = calloc (1, sizeof (SingleGroup))) == NULL)
	    return (OUT_OF_MEMORY);
	initSingleGroup (x);

	/* Open the input image. */
	getSingleGroup (sts->input, extver, x);
	if (hstio_err())
	    return (OPEN_FAILED);
	if (sts->printtime)
	    TimeStamp ("Input read into memory", sts->rootname);

	/* Get header info that varies from imset to imset. */
	if ((status = GetGrpInfo1 (sts, &x->sci.hdr)))
	    return (status);

	if (strcmp (sts->obstype, "SPECTROSCOPIC") == 0 && !sts->wavecal) {
	    v_helio = RadialVel (sts->ra_targ, sts->dec_targ,
			(sts->expstart + sts->expend) / 2.);
	    if ((status = Put_KeyD (&x->sci.hdr, "V_HELIO", v_helio,
                                    "heliocentric radial velocity (km/s)")))
		return (status);
	    if (sts->verbose) {
		printf ("         Heliocentric radial velocity = %.3f (km/s)\n",
				v_helio);
	    }
	}

	/* For anything other than a bias image, check for zero exposure time
	   or constant data values.
	*/
	if (sts->wavecal || !sts->bias_exposure) {
	    float minval=0., maxval=1.;
	    if (sts->exptime > 0.)
		getMinMax (x, &minval, &maxval);
	    if (sts->exptime <= 0. || maxval <= minval) {
		char msg1[81], wavecal_str[81];
		if (sts->wavecal)
		    strcpy (wavecal_str, "wavecal ");
		else
		    wavecal_str[0] = '\0';
		if (sts->exptime <= 0.) {
		    sprintf (msg1,
	    "Warning  %simset %d flagged as bad because exptime = %.6g\n",
				wavecal_str, extver, sts->exptime);
		} else {
		    sprintf (msg1,
	    "Warning  %simset %d flagged as bad because all values = %.6g\n",
				wavecal_str, extver, maxval);
		}
		printf ("%s", msg1);
		if ((status = Put_KeyB (&x->sci.hdr, "IMSET_OK", 0,
				"is the current imset good?")) != 0)
		    return (status);
		if (extver == 1) {
		    UCalVer (x->globalhdr);
		    UFilename (sts->output, x->globalhdr);
		}
		putSingleGroup (sts->output, extver, x, option);
		if (hstio_err()) {
		    printf ("ERROR    Couldn't write imset %d.\n", extver);
		    return (IO_ERROR);
		}
		freeSingleGroup (x);
		free (x);
		return (0);
	    }
	}
	/* After incrementing, this will normally be the same as extver. */
	(*ngood_extver)++;

	/* Print a warning if voltages are too low. */
	/* Removed 2/1/02 as per OPR # 45138 */

	/* For the CCD, update two primary header keywords.
	   Also reset CRCORR if there's only one image set.
	*/
	if (sts->detector == CCD_DETECTOR && *ngood_extver == 1) {

	    if ((status = Put_KeyF (x->globalhdr, "ATODGAIN",
                                    sts->atodgain, "")))
		return (status);
	    if ((status = Put_KeyF (x->globalhdr, "READNSE",
                                    sts->readnoise, "")))
		return (status);
	    printf ("\n");
	    PrRefInfo ("ccdtab", sts->ccdpar.name, sts->ccdpar.pedigree,
		sts->ccdpar.descrip, sts->ccdpar.descrip2);
	    if ((status = CCDHistory (sts, x->globalhdr)))
		return (status);

	    if (sts->nimages == 1) {		/* only one imset? */
		if ((status = GetSwitch (x->globalhdr, "CRCORR", &crcorr)))
		    return (status);
		if (crcorr == PERFORM) {
		    printf (
	"Warning  Only one imset, so CRCORR will be reset to SKIPPED.\n");
		    if ((status = Put_KeyS (x->globalhdr,
                                            "CRCORR", "SKIPPED", "")))
			return (status);
		}
	    }
	}

	if (sts->detector != CCD_DETECTOR) {
	    /* Allocate string for Doppler message. */
	    if ((doppstr = calloc (STIS_LINE+1, sizeof (char))) == NULL)
		return (OUT_OF_MEMORY);
	}

	/* Data quality initialization and (for the CCD) check saturation. */
	dqiMsg (sts, *ngood_extver);
	if ((sts->dqicorr == PERFORM) ||
	    ((sts->dqicorr == DUMMY) && (sts->detector == CCD_DETECTOR))) {
	    if ((status = doDQI (sts, x)))
		return (status);
	    PrSwitch ("dqicorr", COMPLETE);
	    if (sts->printtime)
		TimeStamp ("DQICORR complete", sts->rootname);
	}
	if (*ngood_extver == 1 && !OmitStep (sts->dqicorr))
	    if ((status = dqiHistory (sts, x->globalhdr)))
		return (status);
	if (sts->dqicorr == PERFORM && sts->doppcorr == PERFORM)
	    DoppMsg (doppstr, &doppcount, "DQICORR");

	/* Analog to digital correction. */
	AtoDMsg (sts, *ngood_extver);
	if (sts->atodcorr == PERFORM) {
	    if ((status = doAtoD (sts, x)))
		return (status);
	    PrSwitch ("atodcorr", COMPLETE);
	    if (sts->printtime)
		TimeStamp ("ATODCORR complete", sts->rootname);
	}
	if (*ngood_extver == 1 && !OmitStep (sts->atodcorr))
	    if ((status = atodHistory (sts, x->globalhdr)))
		return (status);

	/* Subtract bias level, as determined from overscan regions.
	   Note that x is freed and reallocated by doBlev.
	*/
	BlevMsg (sts);
	if (sts->blevcorr == PERFORM) {
	    if ((status = doBlev (sts, &x, extver,
                                  &meanblev, &done, &driftcorr)))
		return (status);
	    if ((status = Put_KeyF (&x->sci.hdr, "MEANBLEV", meanblev,
                                    "mean of bias levels subtracted")))
		return (status);
	    if (done) {
		printf (
"         Bias level from overscan has been subtracted; \\\n");
		printf ("         mean of bias levels subtracted was %.6g.\n",
			meanblev);
	    } else {
		printf ("         Default bias level %.6g was subtracted.\n",
			meanblev);
	    }
	    PrSwitch ("blevcorr", COMPLETE);
	    if (sts->printtime)
		TimeStamp ("BLEVCORR complete", sts->rootname);

	    /* We've just subtracted the bias, so nothing further needs to be
		subtracted during error array initialization.
	    */
	    sts->err_init_bias = 0.F;

	} else {

	    /* If blevcorr has not been done, we must subtract a bias
		during error array initialization.  If blevcorr was done
		on a previous run of calstis1, the error array will already
		have been initialized, and this value will not be used.
	    */
	    sts->err_init_bias = sts->ccdbias;
	}
	if (*ngood_extver == 1 && !OmitStep (sts->blevcorr))
	    if ((status = blevHistory (sts, x->globalhdr, done, driftcorr)))
		return (status);

	/* Fill in the error array, if it currently contains all zeros. */
	if (sts->noisecorr == PERFORM) {
	    if ((status = doNoise (sts, x, &done)))
		return (status);
	    if (done) {
		if (*ngood_extver == 1) {
		    if ((status = noiseHistory (x->globalhdr)))
			return (status);
		}
		printf ("         Uncertainty array initialized");
		if (sts->detector == CCD_DETECTOR) {
		    if (sts->err_init_bias > 0.) {
			printf (", readnoise=%.5g, gain=%.5g, bias=%.5g\n",
			sts->readnoise, sts->atodgain, sts->err_init_bias);
		    } else {
			printf (", readnoise=%.5g, gain=%.5g\n",
			sts->readnoise, sts->atodgain);
		    }
		} else {
		    printf (".\n");
		}
		if (sts->printtime)
		    TimeStamp ("Uncertainty array initialized", sts->rootname);
	    }
	}

	/* Convert MAMA data from high-res to low-res.
	   Note that x can be freed and reallocated by doLoRes.
	*/
	LoResMsg (sts);
	if (sts->lorscorr == PERFORM) {
	    if ((status = doLoRes (sts, &x, &done)))
		return (status);
	    if (done) {
		if (*ngood_extver == 1) {	/* update PLATESC keyword */
		    if ((status = UpdatePlateSc (sts, x)))
			return (status);
		}
		PrSwitch ("lorscorr", COMPLETE);
		if (sts->printtime)
		    TimeStamp ("LORSCORR complete", sts->rootname);
	    } else {
		PrSwitch ("lorscorr", SKIPPED);
	    }
	    if (*ngood_extver == 1)
		if ((status = loresHistory (sts, x->globalhdr, done)))
		    return (status);
	}

	/* Check (possibly correct) for nonlinearity. */
	NonLinMsg (sts, *ngood_extver);
	if (sts->glincorr == PERFORM || sts->lflgcorr == PERFORM) {
	    if ((status = doNonLin (sts, x, &gsat, &lsat)))
		return (status);
	    if (sts->glincorr == PERFORM)
		PrSwitch ("glincorr", COMPLETE);
	    if (sts->lflgcorr == PERFORM)
		PrSwitch ("lflgcorr", COMPLETE);
	    if (sts->printtime)
		TimeStamp ("Nonlinearity corr. complete", sts->rootname);
	}
	if (!OmitStep (sts->glincorr) || !OmitStep (sts->lflgcorr)) {
	    if (*ngood_extver == 1)
		if ((status = nonlinHistory (sts, x->globalhdr)))
		    return (status);
	}

	/* Subtract bias image. */
	BiasMsg (sts, *ngood_extver);
	if (sts->biascorr == PERFORM) {
	    if ((status = doBias (sts, x)))
		return (status);
	    if (sts->biascorr == PERFORM) {	/* not reset to SKIPPED? */
		PrSwitch ("biascorr", COMPLETE);
		if (sts->printtime)
		    TimeStamp ("BIASCORR complete", sts->rootname);
	    } else {
		PrSwitch ("biascorr", sts->biascorr);
	    }
	}
	if (*ngood_extver == 1 && !OmitStep (sts->biascorr))
	    if ((status = biasHistory (sts, x->globalhdr)))
		return (status);

	/* Subtract dark image. */
	DarkMsg (sts, *ngood_extver);
	if (sts->darkcorr == PERFORM) {
	    if ((status = doDark (sts, x, &meandark, extver)))
		return (status);
	    if ((status = Put_KeyF (&x->sci.hdr, "MEANDARK", meandark,
                                    "mean of dark values subtracted")))
		return (status);
	    PrSwitch ("darkcorr", COMPLETE);
	    if (sts->printtime)
		TimeStamp ("DARKCORR complete", sts->rootname);
	}
	if (*ngood_extver == 1 && !OmitStep (sts->darkcorr))
	    if ((status = darkHistory (sts, x->globalhdr)))
		return (status);
	if (sts->darkcorr == PERFORM && sts->doppcorr == PERFORM)
	    DoppMsg (doppstr, &doppcount, "DARKCORR");

	/* Multiply by flat field(s). */
	FlatMsg (sts, *ngood_extver);
	if (sts->flatcorr == PERFORM) {
	    if ((status = doFlat (sts, x)))
		return (status);
	    PrSwitch ("flatcorr", COMPLETE);
	    if (sts->printtime)
		TimeStamp ("FLATCORR complete", sts->rootname);
	}
	if (*ngood_extver == 1 && !OmitStep (sts->flatcorr))
	    if ((status = flatHistory (sts, x->globalhdr)))
		return (status);
	if (sts->flatcorr == PERFORM && sts->doppcorr == PERFORM)
	    DoppMsg (doppstr, &doppcount, "FLATCORR");

	/* Apply shutter shading correction. */
	ShadMsg (sts, *ngood_extver);
	if (sts->shadcorr == PERFORM) {
	    if ((status = doShad (sts, x)))
		return (status);
	    PrSwitch ("shadcorr", COMPLETE);
	    if (sts->printtime)
		TimeStamp ("SHADCORR complete", sts->rootname);
	}
	if (*ngood_extver == 1 && !OmitStep (sts->shadcorr))
	    if ((status = shadHistory (sts, x->globalhdr)))
		return (status);

	/* Assign values to photometry keywords.
	   Note that this is only done for the first imset.  This is
	   because we're updating keywords in the primary header.
	*/
	if (*ngood_extver == 1) {
	    /* Update the PHOTMODE keyword regardless of PHOTCORR. */
	    if ((status = PhotMode (sts, x)))
		return (status);

	    printf ("\n");
	    PrSwitch ("photcorr", sts->photcorr);
	    if (sts->photcorr == PERFORM) {
		if ((status = doPhot (sts, x)))
		    return (status);
		PhotMsg (sts);
		if (sts->photcorr == PERFORM || sts->photcorr == COMPLETE) {
		    PrSwitch ("photcorr", COMPLETE);
		    if (sts->printtime)
			TimeStamp ("PHOTCORR complete", sts->rootname);
		} else {
		    PrSwitch ("photcorr", SKIPPED);
		    if (sts->printtime)
			TimeStamp ("PHOTCORR skipped", sts->rootname);
		}
	    }
	    if (!OmitStep (sts->photcorr))
		if ((status = photHistory (sts, x->globalhdr)))
		    return (status);
	}

	if (sts->detector != CCD_DETECTOR) {
	    /* Doppler message and history. */
	    printf ("\n");
	    if (doppcount > 0) {
		printf ("%s\n", doppstr);
		addHistoryKw (x->globalhdr, doppstr);
	    } else if (sts->doppcorr != PERFORM) {
		PrSwitch ("doppcorr", sts->doppcorr);	/* OMIT */
	    }
	    free (doppstr);
	}

	/* Compute min, max, mean, etc. of good science data. */
	if (sts->statcorr == PERFORM) {
	    printf ("\n");
	    PrSwitch ("statflag", PERFORM);
	    if ((status = doStat (x, sts->sdqflags)))
		return (status);
	    PrSwitch ("statflag", COMPLETE);
	    if (sts->printtime)
		TimeStamp ("STATFLAG complete", sts->rootname);
	}

	/* Write this imset to the output file.  If this is the first imset,
	   CAL_VER and FILENAME keywords will be updated, and the primary
	   header will be written.
	*/
	if (extver == 1) {
	    UCalVer (x->globalhdr);
	    UFilename (sts->output, x->globalhdr);
	}
	putSingleGroup (sts->output, extver, x, option);
	if (hstio_err()) {
	    printf ("ERROR    Couldn't write imset %d.\n", extver);
	    return (IO_ERROR);
	}
	if (sts->printtime)
	    TimeStamp ("Output written to disk", sts->rootname);

	freeSingleGroup (x);
	free (x);

	return (0);
}

static void getMinMax (SingleGroup *x, float *minval, float *maxval) {

	int nx, ny;
	int i, j;
	float value;

	nx = x->sci.data.nx;
	ny = x->sci.data.ny;

	*minval = Pix (x->sci.data, 0, 0);
	*maxval = *minval;
	for (i = 0;  i < nx;  i++) {
	    for (j = 0;  j < ny;  j++) {
		value = Pix (x->sci.data, i, j);
		if (value < *minval)
		    *minval = value;
		if (value > *maxval)
		    *maxval = value;
	    }
	}
}

static void AtoDMsg (StisInfo1 *sts, int extver) {

	if (sts->detector == CCD_DETECTOR) {
	    printf ("\n");
	    PrSwitch ("atodcorr", sts->atodcorr);
	}

	if (extver == 1 && !OmitStep (sts->atodcorr)) {

	    PrRefInfo ("atodtab", sts->atod.name, sts->atod.pedigree,
			sts->atod.descrip, sts->atod.descrip2);
	}
}

static void BiasMsg (StisInfo1 *sts, int extver) {

	if (sts->detector == CCD_DETECTOR) {
	    printf ("\n");
	    PrSwitch ("biascorr", sts->biascorr);
	}

	if (extver == 1 && !OmitStep (sts->biascorr)) {

	    PrRefInfo ("biasfile", sts->bias.name, sts->bias.pedigree,
			sts->bias.descrip, "");
	}
}

static void BlevMsg (StisInfo1 *sts) {

	if (sts->detector == CCD_DETECTOR) {
	    printf ("\n");
	    PrSwitch ("blevcorr", sts->blevcorr);
	}
}

static void DarkMsg (StisInfo1 *sts, int extver) {

	printf ("\n");
	PrSwitch ("darkcorr", sts->darkcorr);

	if (extver == 1 && !OmitStep (sts->darkcorr)) {

	    PrRefInfo ("darkfile", sts->dark.name, sts->dark.pedigree,
			sts->dark.descrip, "");
	    if (sts->detector == NUV_MAMA_DETECTOR) {
		PrRefInfo ("tdctab", sts->tdctab.name, sts->tdctab.pedigree,
			    sts->tdctab.descrip, "");
	    }
	}
}

/* This routine builds a message string for DOPPCORR, one step at a time.
   doppcount should be initialized to zero at the beginning of processing for
   a given imset.  This function should then be called after each calibration
   step for which Doppler convolution was applied to the reference data.
   doppcount will be incremented with each call.  If doppcount is greater
   than zero at the end of processing, the doppstr message should be printed
   (note that a newline is not included) and written as a history record.
*/

static void DoppMsg (char *doppstr, int *doppcount, char *calswitch) {

/* arguments:
char *doppstr    io: message regarding DOPPCORR
int *doppcount   io: the number of cal steps to which DOPPCORR has been applied
char *calswitch  i: name of current calibration step (e.g. FLATCORR)
*/

	if (*doppcount <= 0)
	    strcpy (doppstr, "DOPPCORR applied to ");
	else
	    strcat (doppstr, ", ");

	strcat (doppstr, calswitch);

	(*doppcount)++;
}

static void dqiMsg (StisInfo1 *sts, int extver) {

	printf ("\n");
	PrSwitch ("dqicorr", sts->dqicorr);

	if (extver == 1 && !OmitStep (sts->dqicorr)) {

	    PrRefInfo ("dqitab", sts->bpix.name, sts->bpix.pedigree,
			sts->bpix.descrip, sts->bpix.descrip2);
	}
}

static void FlatMsg (StisInfo1 *sts, int extver) {

	printf ("\n");
	PrSwitch ("flatcorr", sts->flatcorr);

	if (extver == 1 && !OmitStep (sts->flatcorr)) {

	    if (GotFileName (sts->pflt.name)) {		/* pixel-to-pixel */
		PrRefInfo ("pfltfile", sts->pflt.name,
			sts->pflt.pedigree, sts->pflt.descrip, "");
	    }
	    if (GotFileName (sts->dflt.name)) {		/* delta flat */
		PrRefInfo ("dfltfile", sts->dflt.name,
			sts->dflt.pedigree, sts->dflt.descrip, "");
	    }
	    if (GotFileName (sts->lflt.name)) {		/* low-order flat */
		PrRefInfo ("lfltfile", sts->lflt.name,
			sts->lflt.pedigree, sts->lflt.descrip, "");
	    }
	}
}

static void LoResMsg (StisInfo1 *sts) {

	if (sts->detector != CCD_DETECTOR) {
	    printf ("\n");
	    PrSwitch ("lorscorr", sts->lorscorr);
	}
}

static void NonLinMsg (StisInfo1 *sts, int extver) {

	if (sts->detector != CCD_DETECTOR) {
	    printf ("\n");
	    PrSwitch ("glincorr", sts->glincorr);
	    PrSwitch ("lflgcorr", sts->lflgcorr);
	}

	if (!OmitStep (sts->glincorr) || !OmitStep (sts->lflgcorr)) {

	    if (extver == 1) {

		PrRefInfo ("mlintab", sts->mlin.name, sts->mlin.pedigree,
			sts->mlin.descrip, sts->mlin.descrip2);
	    }
	}
}

static void PhotMsg (StisInfo1 *sts) {

	if (!OmitStep (sts->photcorr)) {

	    PrRefInfo ("imphttab", sts->phot.name, sts->phot.pedigree,
			sts->phot.descrip, sts->phot.descrip2);
	}
}

static void ShadMsg (StisInfo1 *sts, int extver) {

	if (sts->detector == CCD_DETECTOR) {
	    printf ("\n");
	    PrSwitch ("shadcorr", sts->shadcorr);
	}

	if (extver == 1 && !OmitStep (sts->shadcorr)) {

	    PrRefInfo ("shadfile", sts->shad.name, sts->shad.pedigree,
			sts->shad.descrip, "");
	}
}

# define ARCSEC_PER_DEGREE  3600.

/* This routine recomputes and updates the PLATESC keyword in the primary
   header, based on the values in the CD matrix in the SCI extension header.
*/

static int UpdatePlateSc (StisInfo1 *sts, SingleGroup *x) {

/* arguments:
StisInfo1 *sts   i: calibration switches and info
SingleGroup *x   i: image that was binned
*/

	int status;

	double cd[4];		/* cd1_1, cd1_2, cd2_1, cd2_2 */
	float scale;		/* plate scale */
	float scale_x, scale_y;	/* scales in X and Y, used for MIRROR */
	int use_def = 1;	/* use default if missing keyword */

	/* Get the CD matrix from the SCI header. */
	if ((status = Get_KeyD (&x->sci.hdr, "CD1_1", use_def, 1., &cd[0])))
	    return (status);
	if ((status = Get_KeyD (&x->sci.hdr, "CD1_2", use_def, 0., &cd[1])))
	    return (status);
	if ((status = Get_KeyD (&x->sci.hdr, "CD2_1", use_def, 0., &cd[2])))
	    return (status);
	if ((status = Get_KeyD (&x->sci.hdr, "CD2_2", use_def, 1., &cd[3])))
	    return (status);

	if (sts->opt_elem[0] == 'G' || sts->opt_elem[0] == 'E' ||
	    strcmp (sts->opt_elem, "PRISM") == 0) {

	    scale = fabs (cd[3]) * ARCSEC_PER_DEGREE;		/* CD2_2 */

	} else if (sts->opt_elem[0] == 'X') {

	    scale = fabs (cd[0]) * ARCSEC_PER_DEGREE;		/* CD1_1 */

	} else if (strncmp (sts->opt_elem, "MIR", 3) == 0) {

	    /* scale in X direction, based on CD1_1 & CD2_1 */
	    scale_x = sqrt (cd[0] * cd[0] + cd[2] * cd[2]) * ARCSEC_PER_DEGREE;

	    /* scale in Y direction, based on CD1_2 & CD2_2 */
	    scale_y = sqrt (cd[1] * cd[1] + cd[3] * cd[3]) * ARCSEC_PER_DEGREE;

	    /* take whichever is larger */
	    scale = (scale_x >= scale_y) ? scale_x : scale_y;

	} else {			/* no change */
	    return (0);
	}

	/* Update the value in the primary header. */
	if ((status = Put_KeyF (x->globalhdr, "PLATESC", scale,
                                "plate scale (arcsec/pixel)")))
	    return (status);

	return (0);
}
