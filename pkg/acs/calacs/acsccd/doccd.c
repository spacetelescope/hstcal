/* Do basic CCD image reduction for the current set of extensions.
   If we're processing the first set (extver = 1), then the primary
   header will be updated with history information, and the calibration
   switches in the header will be reset from "PERFORM" to "COMPLETE".

   Warren Hack, 1998 June 1:
   	Revised to only perform CCD operations on ACS observations. 
   Warren Hack, 1999 Jan 7:
   	Revised to only output overscan-trimmed image when BLEVCORR is 
	performed successfully.  In other cases, full image is written out.
   Warren Hack, 2000 Sept 12:
   	Added support for new processing step: post-flash. 
	
   29-Oct-2001 WJH: Updated to use CCDBIAS[A,B,C,D] for default bias levels
        and to use CCDOFST to select CCDBIAS from revised CCDTAB file. Also,
        reports bias levels for each AMP. 

   10-Dec-2001 WJH: Added update of SIZAXIS keywords to reflect trimmed size.

** This code is a trimmed down version of CALSTIS1 do2d.c.  
*/

# include <string.h>
# include <stdio.h>
# include "hstio.h"

# include "acs.h"
# include "acsinfo.h"
# include "acserr.h"

static void AtoDMsg (ACSInfo *, int);
static void BiasMsg (ACSInfo *, int);
static void FlashMsg (ACSInfo *, int);
static void BlevMsg (ACSInfo *, int);
static void dqiMsg (ACSInfo *, int);

int DoCCD (ACSInfo *acs, int extver) {

/* arguments:
ACSInfo *acs   i: calibration switches and info
int extver       i: "imset" number, the current set of extensions
*/

    extern int status;

    SingleGroup x;	/* used for both input and output */
    int option = 0;
    float meanblev;		/* mean value of overscan bias (for history) */
    float meanflash;	/* mean value of post-flash image (for history) */
    int driftcorr;	/* true means bias level was corrected for drift */
    int done;	/* true means the input SingleGroup has been freed */
    int sizex, sizey; 	/* size of output image */

    int i;		/* loop index */
    int overscan;
    int blevcorr;
    char buff[ACS_FITS_REC+1];
    Bool subarray;

    int doAtoD (ACSInfo *, SingleGroup *);
    int atodHistory (ACSInfo *, Hdr *);
    int doBias (ACSInfo *, SingleGroup *);
    int biasHistory (ACSInfo *, Hdr *);
    int doFlash (ACSInfo *, SingleGroup *, float *);
    int flashHistory (ACSInfo *, Hdr *);
    int doBlev (ACSInfo *, SingleGroup *, int, float *, int *, int *);
    int blevHistory (ACSInfo *, Hdr *, int, int);
    int CCDHistory (ACSInfo *, Hdr *);
    int doDQI (ACSInfo *, SingleGroup *);
    int dqiHistory (ACSInfo *, Hdr *);
    int doNoise (ACSInfo *, SingleGroup *, int *);
    int noiseHistory (Hdr *);
    int GetACSGrp (ACSInfo *, Hdr *);
    int OmitStep (int);
    int PutKeyFlt (Hdr *, char *, float, char *);
    int PutKeyInt (Hdr *, char *, int, char *);	
    void PrSwitch (char *, int);
    void PrRefInfo (char *, char *, char *, char *, char *);
    void TimeStamp (char *, char *);
    void UCalVer (Hdr *);
    void UFilename (char *, Hdr *);
    int FindOverscan (ACSInfo *, int, int, int *);
    int GetCCDTab (ACSInfo *, int, int);
    int GetKeyBool (Hdr *, char *, int, Bool, Bool *);
    
/*========================Start Code here =======================*/	
	initSingleGroup (&x);
	if (acs->printtime)
	    TimeStamp ("Open SingleGroup now...", "");

	/* Open the input image. 
	** FOR ACS: Keep this in memory throughout processing.
	** Read in reference images line-by-line within the individual
	** processing step functions and pass along modified input image
	** from one step to the next...
	*/
	getSingleGroup (acs->input, extver, &x);
	if (hstio_err())
	    return (status = OPEN_FAILED);
	if (acs->printtime)
	    TimeStamp ("Input read into memory", acs->rootname);

	/* Get header info that varies from imset to imset necessary
	**	for reading CCDTAB. 
	*/
	if (GetACSGrp (acs, &x.sci.hdr)) {
		freeSingleGroup (&x);
		return (status);
	}
    
    /* Read in keywords from primary header... */

	if (GetKeyBool (x.globalhdr, "SUBARRAY", NO_DEFAULT, 0, &subarray))
        return (status);
    if (subarray)
        acs->subarray = YES;
    else
        acs->subarray = NO;
    
	/* For the CCD, update primary header keywords.
	   Also reset CRCORR if there's only one image set.
	*/
	/* Get values from tables, using same function used in ACSCCD. */
	if (GetCCDTab (acs, x.sci.data.nx, x.sci.data.ny)) {
        freeSingleGroup (&x);
		return (status);
    }

	if (PutKeyFlt (x.globalhdr, "ATODGNA", acs->atodgain[0], ""))
		return (status);
	if (PutKeyFlt (x.globalhdr, "ATODGNB", acs->atodgain[1], ""))
		return (status);
	if (PutKeyFlt (x.globalhdr, "ATODGNC", acs->atodgain[2], ""))
		return (status);
	if (PutKeyFlt (x.globalhdr, "ATODGND", acs->atodgain[3], ""))
		return (status);
	if (PutKeyFlt (x.globalhdr, "READNSEA", acs->readnoise[0], ""))
		return (status);
	if (PutKeyFlt (x.globalhdr, "READNSEB", acs->readnoise[1], ""))
		return (status);
	if (PutKeyFlt (x.globalhdr, "READNSEC", acs->readnoise[2], ""))
		return (status);
	if (PutKeyFlt (x.globalhdr, "READNSED", acs->readnoise[3], ""))
		return (status);
	trlmessage ("\n");

	PrRefInfo ("ccdtab", acs->ccdpar.name, acs->ccdpar.pedigree,
	acs->ccdpar.descrip, acs->ccdpar.descrip2);

	if (extver == 1) {
		if (CCDHistory (acs, x.globalhdr))
			return (status);	
	}

	/* Get overscan region information from OSCNTAB */
	if( FindOverscan (acs, x.sci.data.nx, x.sci.data.ny, &overscan))
        return (status);

	/* Data quality initialization and (for the CCDs) check saturation. */
	dqiMsg (acs, extver);
	if (acs->dqicorr == PERFORM || acs->dqicorr == DUMMY) {
	    if (doDQI (acs, &x))
		return (status);
	    PrSwitch ("dqicorr", COMPLETE);
	    if (acs->printtime)
		TimeStamp ("DQICORR complete", acs->rootname);
	}
	if (extver == 1 && !OmitStep (acs->dqicorr))
	    if (dqiHistory (acs, x.globalhdr))
		return (status);

	/* Analog to digital correction. */
	AtoDMsg (acs, extver);
	if (acs->atodcorr == PERFORM) {
	    if (doAtoD (acs, &x))
		return (status);
	    PrSwitch ("atodcorr", COMPLETE);
	    if (acs->printtime)
		TimeStamp ("ATODCORR complete", acs->rootname);
	}
	if (extver == 1 && !OmitStep (acs->atodcorr))
	    if (atodHistory (acs, x.globalhdr))
			return (status);

	/* Subtract bias from overscan. */
	BlevMsg (acs, extver);
    blevcorr = PERFORM;
	if (acs->blevcorr == PERFORM) {
	
		/* This is set based on results of FindOver */
		done = overscan;
	    
	    if (doBlev (acs, &x, acs->chip, &meanblev, &done, &driftcorr))
		return (status);
        
	    if (done) {
            
			trlmessage ("Bias level from overscan has been subtracted;");

	    	if (PutKeyFlt (&x.sci.hdr, "MEANBLEV", meanblev,
				"mean of bias levels subtracted")) {
				return (status);
			}

			sprintf (MsgText, "     mean of bias levels subtracted was %.6g.",
				meanblev);
			trlmessage (MsgText);
	    
        } else {
			trlmessage ("Default bias level from CCDTAB was subtracted.");
	    }

        /* Provide immediate feedback to the user on the values computed for each AMP,
            but only for those amps used for the chip being processed. 
        */
        for (i = 0; i < 4; i++){
            if (acs->blev[i] != 0.) {
		        sprintf (MsgText, "     bias level of %.6g was subtracted for AMP %c.", acs->blev[i], AMPSORDER[i]);
		        trlmessage (MsgText); 
            }
        }
        
	    PrSwitch ("blevcorr", COMPLETE);

		/* Set this to complete so overscan-trimmed image will be
			written out.
		*/      
		blevcorr = COMPLETE;
	    
		if (acs->printtime)
		TimeStamp ("BLEVCORR complete", acs->rootname);
	}
	if (extver == 1 && !OmitStep (acs->blevcorr))
    if (blevHistory (acs, x.globalhdr, done, driftcorr))
      return (status);

  /* Fill in the error array, if it initially contains all zeros. */
	if (acs->noisecorr == PERFORM) {
    if (doNoise (acs, &x, &done))
      return (status);
    if (done) {
      if (extver == 1) {
		    if (noiseHistory (x.globalhdr))
          return (status);
      }
      trlmessage ("    Uncertainty array initialized,");
      buff[0] = '\0';
      
      sprintf(MsgText, "    readnoise =");
      for (i=0; i < NAMPS-1; i++) {
        if (acs->readnoise[i] > 0) {
          sprintf (buff, "%.5g,",acs->readnoise[i]);
          strcat (MsgText, buff);
        }
      }
      if (acs->readnoise[NAMPS-1] > 0) {
        sprintf(buff, "%.5g",acs->readnoise[NAMPS-1]);
        strcat (MsgText, buff);
      }
      trlmessage (MsgText);
      
      sprintf(MsgText, "    gain =");
      for (i=0; i < NAMPS-1; i++) {
        if (acs->atodgain[i] > 0) { 
          sprintf(buff, "%.5g,",acs->atodgain[i]);
          strcat (MsgText, buff);
        }
      }
      if (acs->atodgain[NAMPS-1] > 0) { 
        sprintf(buff, "%.5g",acs->atodgain[NAMPS-1]);
        strcat (MsgText, buff);
      }
      trlmessage (MsgText);
      sprintf(MsgText, "   default bias levels =");
      for (i=0; i < NAMPS-1; i++) {
        if (acs->ccdbias[i] > 0) { 
          sprintf(buff, "%.5g,",acs->ccdbias[i]);
          strcat (MsgText, buff);
        }
      }
      if (acs->ccdbias[NAMPS-1] > 0) { 
        sprintf(buff, "%.5g",acs->ccdbias[NAMPS-1]);
        strcat (MsgText, buff);
      }
      trlmessage (MsgText);
      
      if (acs->printtime)
		    TimeStamp ("Uncertainty array initialized", acs->rootname);
		}
	}
  
	/* Subtract bias image. */
	BiasMsg (acs, extver);
	if (acs->biascorr == PERFORM) {
	    if (doBias (acs, &x))
		return (status);
	    PrSwitch ("biascorr", COMPLETE);
	    if (acs->printtime)
		TimeStamp ("BIASCORR complete", acs->rootname);
	}
	if (extver == 1 && !OmitStep (acs->biascorr))
	    if (biasHistory (acs, x.globalhdr))
		return (status);

	/* Subtract post-flash image. */
	FlashMsg (acs, extver);
	if (acs->flashcorr == PERFORM) {
		/* Initialize this to a set value... */
		meanflash = 0.0;
		
	    if (doFlash (acs, &x, &meanflash))
		return (status);
	
		/* Report mean of post-flash image subtracted from science image,
			if it was performed...*/
		if (meanflash > 0.){
			sprintf(MsgText,"Mean of post-flash image (MEANFLSH) = %g",meanflash);
			trlmessage(MsgText);

			/* If they want to add this keyword, we can uncomment this code. */
			if (PutKeyFlt (&x.sci.hdr, "MEANFLSH", meanflash,
				"mean of post-flash values subtracted"))
				return (status);	    

		    PrSwitch ("flshcorr", COMPLETE);
		} else {
		    PrSwitch ("flshcorr", SKIPPED);
	 	}
				
	    if (acs->printtime) TimeStamp ("FLSHCORR complete", acs->rootname);
	}
	if (extver == 1 && !OmitStep (acs->flashcorr)) {
	    if (flashHistory (acs, x.globalhdr))
			return (status);
	}
	
	/* Write this imset to the output file.  If extver is one, the
	   CAL_VER and FILENAME keywords will be updated, and the primary
	   header will be written.
	*/
	if (extver == 1) {
	    UCalVer (x.globalhdr);
	    UFilename (acs->output, x.globalhdr);
	}
	
	/* If BLEVCORR was performed, then output trimmed image,
		otherwise output full image... 
	*/
	if (blevcorr == COMPLETE) {
		/* BLEVCORR was completed, so overscan regions can be trimmed... */
		if (acs->verbose) {
			sprintf(MsgText,"Writing out image with trimx = %d,%d, trimy = %d,%d",acs->trimx[0],acs->trimx[1], acs->trimy[0],acs->trimy[1]);
			trlmessage(MsgText);

			sprintf(MsgText,"Outputting chip %d",acs->chip);
			trlmessage(MsgText);
		}

		/* Output overscan-trimmed, bias-subtracted image 
		**  using the routine 'putSingleGroupSect' from HSTIO V2.0 (or later) 
		**  to write it directly from the full image in memory.
		*/
		sizex = x.sci.data.nx - (acs->trimx[0] + acs->trimx[1]);
		sizey = x.sci.data.ny - (acs->trimy[0] + acs->trimy[1]);

        /* 
        Update SIZAXIS keywords to reflect the new trimmed image size.
        */
		if (PutKeyInt (&x.sci.hdr, "SIZAXIS1", sizex,""))
			return (status);	    
		if (PutKeyInt (&x.sci.hdr, "SIZAXIS2", sizey,""))
			return (status);	    

        /* Now, write out updated trimmed data to disk... */
		putSingleGroupSect (acs->output, extver, &x, acs->trimx[0],
		acs->trimy[0], sizex, sizey, option);
		
	} else {
		/* BLEVCORR was not completed, so keep overscan regions... */
		if (acs->verbose) {
			sprintf(MsgText,"Writing out FULL image with overscan regions");
			trlmessage(MsgText);

			sprintf(MsgText,"Outputting chip %d",acs->chip);
			trlmessage(MsgText);
		}
		
		putSingleGroup (acs->output, extver, &x, option);		
	
	}
		
	if (hstio_err()) {
	    sprintf (MsgText, "Couldn't write imset %d.", extver);
		trlerror (MsgText);
	    return (status = WRITE_FAILED);
	}
	if (acs->printtime)
	    TimeStamp ("Output written to disk", acs->rootname);

	/* Free x, which still has memory allocated. */
	freeSingleGroup (&x);

	return (status);
}

static void AtoDMsg (ACSInfo *acs, int extver) {

	int OmitStep (int);
	void PrSwitch (char *, int);
	void PrRefInfo (char *, char *, char *, char *, char *);

	trlmessage ("\n");
	PrSwitch ("atodcorr", acs->atodcorr);

	if (extver == 1 && !OmitStep (acs->atodcorr)) {

	    PrRefInfo ("atodtab", acs->atod.name, acs->atod.pedigree,
			acs->atod.descrip, acs->atod.descrip2);
	}
}

static void BiasMsg (ACSInfo *acs, int extver) {

	int OmitStep (int);
	void PrSwitch (char *, int);
	void PrRefInfo (char *, char *, char *, char *, char *);

	trlmessage ("\n");
	PrSwitch ("biascorr", acs->biascorr);

	if (extver == 1 && !OmitStep (acs->biascorr)) {

	    PrRefInfo ("biasfile", acs->bias.name, acs->bias.pedigree,
			acs->bias.descrip, "");
	}
}

static void FlashMsg (ACSInfo *acs, int extver) {

	int OmitStep (int);
	void PrSwitch (char *, int);
	void PrRefInfo (char *, char *, char *, char *, char *);

	trlmessage ("\n");
	PrSwitch ("flshcorr", acs->flashcorr);

	if (extver == 1 && !OmitStep (acs->flashcorr)) {

	    PrRefInfo ("flshfile", acs->flash.name, acs->flash.pedigree,
			acs->flash.descrip, "");
	}
}

static void BlevMsg (ACSInfo *acs, int extver) {

	void PrSwitch (char *, int);
	void PrRefInfo (char *, char *, char *, char *, char *);
	int OmitStep (int);

	trlmessage ("\n");
	PrSwitch ("blevcorr", acs->blevcorr);

	if (extver == 1 && !OmitStep (acs->blevcorr)) {

	    PrRefInfo ("oscntab", acs->oscn.name, acs->oscn.pedigree,
			acs->oscn.descrip, acs->oscn.descrip2);
	}

}

static void dqiMsg (ACSInfo *acs, int extver) {

	int OmitStep (int);
	void PrSwitch (char *, int);
	void PrRefInfo (char *, char *, char *, char *, char *);

	trlmessage ("\n");
	PrSwitch ("dqicorr", acs->dqicorr);

	if (extver == 1 && !OmitStep (acs->dqicorr)) {

	    PrRefInfo ("dqitab", acs->bpix.name, acs->bpix.pedigree,
			acs->bpix.descrip, acs->bpix.descrip2);
	}
}
