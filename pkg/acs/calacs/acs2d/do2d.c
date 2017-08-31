/* Do basic 2-D image reduction for the current set of extensions.
 If we're processing the first set (extver = 1), then the primary
 header will be updated with history information, and the calibration
 switches in the header will be reset from "PERFORM" to "COMPLETE".

 Warren Hack, 1998 June 10:
 Initial ACS version.
 Warren Hack, 1999 Jan 7:
 Revised to not proceed unless it detects a overscan-trimmed image.
 The keywords LTV1, LTV2 should be zero for an overscan-trimmed
 image, otherwise it will return without processing.
 Warren Hack, 2002 Jan 30:
 Added update to BUNIT keyword after flat-field correction to reflect
 change in datatype for science data after calibration.
 Warren Hack, 2002 Apr 9:
 Changed BUNIT update to work on SCI and ERR headers, instead of Primary.
 Warren Hack, 2002 Apr 17:
 Removed all references to 'statcorr' and STATFLAG. Always perform 'doStat'.
 Warren Hack, 2008 Jan 18:
 Added logic to trap negative values in MAMA data.
 Pey Lian Lim, 2012 Dec 11:
 Moved FLASHCORR stuff from ACSCCD.
 */

# include <string.h>
# include <stdio.h>

#include "hstcal.h"
# include "hstio.h"
# include "acs.h"
# include "acsinfo.h"
# include "hstcalerr.h"

static void DarkMsg (ACSInfo *, int, int);
static void FlashMsg (ACSInfo *, int, int);
static void dqiMsg (ACSInfo *, int);
static void NonLinMsg (ACSInfo *, int);
static void PhotMsg (ACSInfo *);
static void ShadMsg (ACSInfo *, int);
static void FlatMsg (ACSInfo *, int);

static int OscnTrimmed (Hdr*, Hdr *);

int Do2D (ACSInfo *acs2d, int extver) {

  /* arguments:
   ACSInfo *acs2d   i: calibration switches and info
   int extver       i: "imset" number, the current set of extensions
   */

	extern int status;

	SingleGroup x;	/* used for both input and output */
	int option = 0;
	float meandark;		/* mean value of dark (for history) */
        float meanflash; /* mean value of post-flash image (for history) */
	int gsat, lsat;		/* > 0 if saturated pixels found */
	int done;	/* true means error array was initialized in donoise */
	int i;

  int logit;
	char buff[ACS_FITS_REC+1];
  Bool subarray;

  char units[CHAR_LINE_LENGTH];

  int to_electrons(ACSInfo *, SingleGroup *);
  int check_zero_noise(SingleGroup *);
	int CCDHistory (ACSInfo *, Hdr *);
	int doDark (ACSInfo *, SingleGroup *, float *);
	int darkHistory (ACSInfo *, Hdr *);
        int doFlash (ACSInfo *, SingleGroup *, float *);
        int flashHistory (ACSInfo *, Hdr *);
	int doDQI (ACSInfo *, SingleGroup *);
	int dqiHistory (ACSInfo *, Hdr *);
	int doFlat (ACSInfo *, int, SingleGroup *);
	int flatHistory (ACSInfo *, Hdr *);
	int doNonLin (ACSInfo *, SingleGroup *, int *, int *);
	int nonlinHistory (ACSInfo *, Hdr *);
	int doNoise (ACSInfo *, SingleGroup *, int *);
	int noiseHistory (Hdr *);
	int doPhot (ACSInfo *, SingleGroup *);
	int PhotMode (ACSInfo *, SingleGroup *);
	int photHistory (ACSInfo *, Hdr *);
	int doShad (ACSInfo *, int, SingleGroup *);
	int shadHistory (ACSInfo *, Hdr *);
	int doStat (SingleGroup *, short);
	int GetACSGrp (ACSInfo *, Hdr *);
	int OmitStep (int);
	int PutKeyFlt (Hdr *, char *, float, char *);
  int PutKeyStr(Hdr *, char *, char *, char *);
  int GetKeyStr (Hdr *, char *, int, char *, char *, int);
	void PrSwitch (char *, int);
	void PrRefInfo (char *, char *, char *, char *, char *);
	void TimeStamp (char *, char *);
	void UCalVer (Hdr *);
	void UFilename (char *, Hdr *);
	int GetCCDTab (ACSInfo *, int, int);
        int GetKeyBool (Hdr *, char *, int, Bool, Bool *);
        int UpdateSwitch (char *, int, Hdr *, int *);

	initSingleGroup (&x);

	/* Open the input image. */
	getSingleGroup (acs2d->input, extver, &x);
	if (hstio_err())
    return (status = OPEN_FAILED);
	if (acs2d->printtime)
    TimeStamp ("Input read into memory", acs2d->rootname);

	/* Get header info that varies from imset to imset. */
	if (GetACSGrp (acs2d, &x.sci.hdr)) {
		freeSingleGroup (&x);
    return (status);
	}

	/* For the CCD, update key primary header keywords.	*/
	if (acs2d->detector != MAMA_DETECTOR) {

    /* Has the CCD image been overscan trimmed? */
    if (OscnTrimmed (x.globalhdr, &x.sci.hdr)) {
      freeSingleGroup (&x);
      return (status);
    }

    /* Read in keywords from primary header... */
    if (GetKeyBool (x.globalhdr, "SUBARRAY", NO_DEFAULT, 0, &subarray))
      return (status);
    if (subarray)
      acs2d->subarray = YES;
    else
      acs2d->subarray = NO;

		/* Get values from tables, using same function used in ACSCCD. */
    if (GetCCDTab (acs2d, x.sci.data.nx, x.sci.data.ny)) {
      freeSingleGroup (&x);
			return (status);
    }

    if (PutKeyFlt (x.globalhdr, "ATODGNA", acs2d->atodgain[0], ""))
      return (status);
    if (PutKeyFlt (x.globalhdr, "ATODGNB", acs2d->atodgain[1], ""))
      return (status);
		if (PutKeyFlt (x.globalhdr, "ATODGNC", acs2d->atodgain[2], ""))
      return (status);
    if (PutKeyFlt (x.globalhdr, "ATODGND", acs2d->atodgain[3], ""))
      return (status);
    if (PutKeyFlt (x.globalhdr, "READNSEA", acs2d->readnoise[0], ""))
      return (status);
    if (PutKeyFlt (x.globalhdr, "READNSEB", acs2d->readnoise[1], ""))
      return (status);
    if (PutKeyFlt (x.globalhdr, "READNSEC", acs2d->readnoise[2], ""))
      return (status);
    if (PutKeyFlt (x.globalhdr, "READNSED", acs2d->readnoise[3], ""))
      return (status);
    trlmessage ("\n");

    PrRefInfo ("ccdtab", acs2d->ccdpar.name, acs2d->ccdpar.pedigree,
               acs2d->ccdpar.descrip, acs2d->ccdpar.descrip2);

		if (extver == 1 ){
      if (CCDHistory (acs2d, x.globalhdr))
				return (status);
		}

	}

	/* Data quality initialization and (for the CCD) check saturation. */
	dqiMsg (acs2d, extver);

	if (acs2d->dqicorr == PERFORM ||
	    (acs2d->dqicorr == DUMMY && acs2d->detector != MAMA_DETECTOR)) {
    if (doDQI (acs2d, &x))
      return (status);
    PrSwitch ("dqicorr", COMPLETE);
    if (acs2d->printtime)
      TimeStamp ("DQICORR complete", acs2d->rootname);
	}

  if (extver == 1 && !OmitStep (acs2d->dqicorr))
    if (dqiHistory (acs2d, x.globalhdr))
      return (status);

	/* Check (possibly correct) for nonlinearity. */
	NonLinMsg (acs2d, extver);
	if (acs2d->glincorr == PERFORM || acs2d->lflgcorr == PERFORM) {
    if (doNonLin (acs2d, &x, &gsat, &lsat))
      return (status);
    if (acs2d->glincorr == PERFORM)
      PrSwitch ("glincorr", COMPLETE);
    if (acs2d->lflgcorr == PERFORM)
      PrSwitch ("lflgcorr", COMPLETE);
    if (acs2d->printtime)
      TimeStamp ("Nonlinearity corr. complete", acs2d->rootname);
	}

  if (!OmitStep (acs2d->glincorr) || !OmitStep (acs2d->lflgcorr)) {
    if (extver == 1)
      if (nonlinHistory (acs2d, x.globalhdr))
		    return (status);
	}

  /* if the data aren't already in electrons, do it here. */
  if (GetKeyStr(&x.sci.hdr, "BUNIT", USE_DEFAULT, "", units, CHAR_LINE_LENGTH)) {
    freeSingleGroup(&x);
    return status;
  }

  if (strncmp(units, "ELECTRONS", 6) != 0) {
    if (to_electrons(acs2d, &x)) {
      freeSingleGroup(&x);
      return (status);
    }

    if (PutKeyStr (&x.sci.hdr, "BUNIT", "ELECTRONS", "")) {
      freeSingleGroup(&x);
      return (status);
    }
    if (PutKeyStr (&x.err.hdr, "BUNIT", "ELECTRONS", "")) {
      freeSingleGroup(&x);
      return (status);
    }
  }

  /* Fill in the error array, if it initially contains all zeros. */
	if (acs2d->noisecorr == PERFORM && check_zero_noise(&x) == YES) {
    if (doNoise (acs2d, &x, &done))
      return (status);

    if (done) {
      if (extver == 1) {
		    if (noiseHistory (x.globalhdr))
          return (status);
      }

      trlmessage("         Uncertainty array initialized,");
      buff[0] = '\0';

	    if (acs2d->detector != MAMA_DETECTOR) {
        sprintf(MsgText, "    readnoise =");

        for (i=0; i < NAMPS-1; i++) {
          if (acs2d->readnoise[i] > 0) {
            sprintf (buff, "%.5g,",acs2d->readnoise[i]);
            strcat (MsgText, buff);
          }
        }

        if (acs2d->readnoise[NAMPS-1] > 0) {
          sprintf(buff, "%.5g",acs2d->readnoise[NAMPS-1]);
          strcat (MsgText, buff);
        }

        trlmessage(MsgText);

        sprintf(MsgText, "    gain =");
        for (i=0; i < NAMPS-1; i++) {
          if (acs2d->atodgain[i] > 0) {
            sprintf(buff, "%.5g,",acs2d->atodgain[i]);
            strcat (MsgText, buff);
          }
        }

        if (acs2d->atodgain[NAMPS-1] > 0) {
          sprintf(buff, "%.5g",acs2d->atodgain[NAMPS-1]);
          strcat (MsgText, buff);
        }

        trlmessage(MsgText);

      }

      if (acs2d->printtime)
		    TimeStamp ("Uncertainty array initialized", acs2d->rootname);
    }
	}

	/* Subtract dark image. */
	DarkMsg (acs2d, extver, acs2d->pctecorr);
	if (acs2d->darkcorr == PERFORM) {
    if (doDark (acs2d, &x, &meandark))
			return (status);

		sprintf(MsgText,"Mean of dark image (MEANDARK) = %g",meandark);
		trlmessage(MsgText);

		if (PutKeyFlt (&x.sci.hdr, "MEANDARK", meandark,
                   "mean of dark values subtracted"))
			return (status);
		PrSwitch ("darkcorr", COMPLETE);

		if (acs2d->printtime)
			TimeStamp ("DARKCORR complete", acs2d->rootname);
	}

  if (extver == 1 && !OmitStep (acs2d->darkcorr))
    if (darkHistory (acs2d, x.globalhdr))
      return (status);

    /**************************************************************************/
    /* Subtract post-flash image. */
    FlashMsg(acs2d, extver, acs2d->pctecorr);

    if (acs2d->flashcorr == PERFORM) {
        if (doFlash(acs2d, &x, &meanflash))
            return (status);

        sprintf(MsgText, "Mean of post-flash image (MEANFLSH) = %g", meanflash);
        trlmessage(MsgText);

        if (PutKeyFlt(&x.sci.hdr, "MEANFLSH", meanflash,
                      "mean of post-flash values subtracted"))
            return (status);

        PrSwitch ("flshcorr", COMPLETE);

        if (acs2d->printtime)
            TimeStamp("FLSHCORR complete", acs2d->rootname);
    }

    if (extver == 1 && !OmitStep (acs2d->flashcorr))
        if (flashHistory(acs2d, x.globalhdr))
            return (status);
    /**************************************************************************/

	/*  These corrections are not yet ready to be run...  They require
   the use of Sections...

   Multiply by flat field(s).
   */
	FlatMsg (acs2d, extver);
	if (acs2d->flatcorr == PERFORM) {
    if (doFlat (acs2d, extver, &x))
      return (status);
    PrSwitch ("flatcorr", COMPLETE);
    if (acs2d->printtime)
      TimeStamp ("FLATCORR complete", acs2d->rootname);
    if (PutKeyStr (&x.sci.hdr, "BUNIT", "ELECTRONS", ""))
      return (status);
    if (PutKeyStr (&x.err.hdr, "BUNIT", "ELECTRONS", ""))
      return (status);

	}

  if (extver == 1 && !OmitStep (acs2d->flatcorr))
    if (flatHistory (acs2d, x.globalhdr))
      return (status);

	/*
   Apply shutter shading correction.
   */
	ShadMsg (acs2d, extver);
	if (acs2d->shadcorr == PERFORM) {
    if (doShad (acs2d, extver, &x))
      return (status);
    PrSwitch ("shadcorr", COMPLETE);
    if (acs2d->printtime)
      TimeStamp ("SHADCORR complete", acs2d->rootname);
	}
	if (extver == 1 && !OmitStep (acs2d->shadcorr))
    if (shadHistory (acs2d, x.globalhdr))
      return (status);

	/* Assign values to photometry keywords.	*/
	trlmessage ("\n");
	PrSwitch ("photcorr", acs2d->photcorr);
	/* Update the PHOTMODE keyword regardless of PHOTCORR. */
	if (PhotMode (acs2d, &x))
    return (status);

  /* Calculate the photometry keyword values if set to PERFORM... */
	if (acs2d->photcorr == PERFORM) {
    status = doPhot (acs2d, &x);
    if (status == ACS_OK) {
      PhotMsg (acs2d);
      PrSwitch ("photcorr", COMPLETE);
      if (acs2d->printtime)
        TimeStamp ("PHOTCORR complete", acs2d->rootname);
    } else if (status == CAL_STEP_NOT_DONE) {
      acs2d->photcorr = IGNORED;
      PrSwitch("photcorr", SKIPPED);
    } else {
      return status;
    }
	}

	if (extver == 1 && !OmitStep (acs2d->photcorr))
    if (photHistory (acs2d, x.globalhdr))
      return (status);

	/* Compute min, max, mean, etc. of good science data. */
	if (doStat (&x, acs2d->sdqflags))
		return (status);
	if (acs2d->printtime)
		TimeStamp ("Image statistics computed", acs2d->rootname);

	/* Write this imset to the output file.  If extver is one, the
   CAL_VER and FILENAME keywords will be updated, and the primary
   header will be written.
   */
	if (extver == 1) {
    UCalVer (x.globalhdr);
    UFilename (acs2d->output, x.globalhdr);
	}

  /* Update EXPSCORR switch upon completion */
  /* Added update to EXPSCORR switch to record when it
   has been completed... 4-June-2003 WJH */
  logit = 0;
  if (UpdateSwitch ("EXPSCORR", acs2d->expscorr, x.globalhdr, &logit))
    return(status);

	putSingleGroup (acs2d->output, extver, &x, option);
	if (hstio_err()) {
    sprintf (MsgText, "Couldn't write imset %d.", extver);
    trlerror (MsgText);
		return (status = 1001);
	}

	if (acs2d->printtime)
    TimeStamp ("Output written to disk", acs2d->rootname);

	freeSingleGroup (&x);

	return (status);
}

static void DarkMsg (ACSInfo *acs2d, int extver, int pctecorr) {

	int OmitStep (int);
	void PrSwitch (char *, int);
	void PrRefInfo (char *, char *, char *, char *, char *);

	trlmessage ("\n");
	PrSwitch ("darkcorr", acs2d->darkcorr);

	if (extver == 1 && !OmitStep (acs2d->darkcorr)) {

    if (pctecorr == PERFORM) {
      PrRefInfo ("darkfile", acs2d->darkcte.name, acs2d->darkcte.pedigree,
                 acs2d->darkcte.descrip, "");
    } else {
      PrRefInfo ("darkfile", acs2d->dark.name, acs2d->dark.pedigree,
                 acs2d->dark.descrip, "");
    }
	}
}

static void FlashMsg (ACSInfo *acs2d, int extver, int pctecorr) {

    int OmitStep (int);
    void PrSwitch (char *, int);
    void PrRefInfo (char *, char *, char *, char *, char *);

    trlmessage ("\n");
    PrSwitch ("flshcorr", acs2d->flashcorr);

    if (extver == 1 && !OmitStep (acs2d->flashcorr)) {

        if (pctecorr == PERFORM) {
            PrRefInfo ("flshfile", acs2d->flashcte.name,
                       acs2d->flashcte.pedigree, acs2d->flashcte.descrip, "");
        } else {
            PrRefInfo ("flshfile", acs2d->flash.name,
                       acs2d->flash.pedigree, acs2d->flash.descrip, "");
        }
    }
}

static void dqiMsg (ACSInfo *acs2d, int extver) {

	int OmitStep (int);
	void PrSwitch (char *, int);
	void PrRefInfo (char *, char *, char *, char *, char *);

	trlmessage ("\n");
	PrSwitch ("dqicorr", acs2d->dqicorr);

	if (extver == 1 && !OmitStep (acs2d->dqicorr)) {

    PrRefInfo ("dqitab", acs2d->bpix.name, acs2d->bpix.pedigree,
               acs2d->bpix.descrip, acs2d->bpix.descrip2);
	}
}

static void FlatMsg (ACSInfo *acs2d, int extver) {

	int GotFileName (char *);
	int OmitStep (int);
	void PrSwitch (char *, int);
	void PrRefInfo (char *, char *, char *, char *, char *);

	trlmessage ("\n");
	PrSwitch ("flatcorr", acs2d->flatcorr);

	if (extver == 1 && !OmitStep (acs2d->flatcorr)) {

    if (GotFileName (acs2d->pflt.name)) {		/* pixel-to-pixel */
      PrRefInfo ("pfltfile", acs2d->pflt.name,
                 acs2d->pflt.pedigree, acs2d->pflt.descrip, "");
    }
    if (GotFileName (acs2d->dflt.name)) {		/* delta flat */
      PrRefInfo ("dfltfile", acs2d->dflt.name,
                 acs2d->dflt.pedigree, acs2d->dflt.descrip, "");
    }
    if (GotFileName (acs2d->lflt.name)) {		/* low-order flat */
      PrRefInfo ("lfltfile", acs2d->lflt.name,
                 acs2d->lflt.pedigree, acs2d->lflt.descrip, "");
    }
    if (GotFileName (acs2d->cflt.name)) {		/* coronagraphic flat */
      PrRefInfo ("cfltfile", acs2d->cflt.name,
                 acs2d->cflt.pedigree, acs2d->cflt.descrip, "");
      /* spot position table */
	    PrRefInfo ("spottab", acs2d->spot.name, acs2d->spot.pedigree,
                 acs2d->spot.descrip, acs2d->spot.descrip2);
    }
	}
}

static void NonLinMsg (ACSInfo *acs2d, int extver) {

	int OmitStep (int);
	void PrSwitch (char *, int);
	void PrRefInfo (char *, char *, char *, char *, char *);

	if (acs2d->detector == MAMA_DETECTOR) {
    trlmessage ("\n");
    PrSwitch ("glincorr", acs2d->glincorr);
    PrSwitch ("lflgcorr", acs2d->lflgcorr);
	}

	if (!OmitStep (acs2d->glincorr) || !OmitStep (acs2d->lflgcorr)) {

    if (extver == 1) {

      PrRefInfo ("mlintab", acs2d->mlin.name, acs2d->mlin.pedigree,
                 acs2d->mlin.descrip, acs2d->mlin.descrip2);
    }
	}
}

static void PhotMsg (ACSInfo *acs2d) {

	int OmitStep (int);
	void PrRefInfo (char *, char *, char *, char *, char *);

	if (!OmitStep (acs2d->photcorr)) {
    PrRefInfo ("imphttab", acs2d->phot.name, acs2d->phot.pedigree,
               acs2d->phot.descrip, acs2d->phot.descrip2);
	}
}

static void ShadMsg (ACSInfo *acs2d, int extver) {

	int OmitStep (int);
	void PrSwitch (char *, int);
	void PrRefInfo (char *, char *, char *, char *, char *);

	if (acs2d->detector != MAMA_DETECTOR) {
    trlmessage ("\n");
    PrSwitch ("shadcorr", acs2d->shadcorr);
	}

	if (extver == 1 && !OmitStep (acs2d->shadcorr)) {

    PrRefInfo ("shadfile", acs2d->shad.name, acs2d->shad.pedigree,
               acs2d->shad.descrip, "");
	}
}


/* This function verifies whether the overscan region has been
 trimmed from the image during BLEVCORR.  If not, it tells the user,
 and shuts down gracefully.
 *****************
 This test needs to be made robust against sub-arrays!!!!!
 *****************
 */

static int OscnTrimmed (Hdr *phdr, Hdr *hdr) {

	extern int status;

	double ltv1, ltv2;
  int blevcorr;

	int GetKeyDbl (Hdr *, char *, int, double, double *);
	int GetSwitch (Hdr *, char *, int *);

  if (GetSwitch (phdr, "BLEVCORR", &blevcorr))
    return(status);

	if (GetKeyDbl (hdr, "LTV1", USE_DEFAULT, 0., &ltv1))
    return (status);
	if (GetKeyDbl (hdr, "LTV2", USE_DEFAULT, 0., &ltv2))
    return (status);

	/* If there is any overscan region still specified in the
   image header, then we can not process this image yet...
   */
	if ( blevcorr != COMPLETE || (ltv1 > 0. || ltv2 > 0.)) {
		trlerror ("Overscan region not trimmed.  Please perform BLEVCORR.");
		status = CAL_STEP_NOT_DONE;
	}

	return (status);
}
