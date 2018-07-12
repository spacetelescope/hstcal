# include <stdio.h>
# include <string.h>

#include "hstcal.h"
# include "hstio.h"
# include "wf3.h"
# include "wf3info.h"
# include "hstcalerr.h"

/* This routine scales CHIP2 in the UVIS data so that the flux correction over both chips
   is uniform. It uses the ratio of PHTFLAM2 / PHTFLAM1 and saves that value for reference
   in the header as PHTRATIO. This will be performed by default in the pipeline and the
   previous PHOTFLAM keyword will be valid for both chips after the correction.

   If users don't wish to perform the correction then they set FLUXCORR OMIT and can use
   the values of PHTFLAM1 and PHTFLAM2 to correct the flux in the respective chips

    MLS Dec 6, 2013

    Update: Nov 2014
    This step takes place after all the WF32d steps and after the FLT file has been initially
    written to disk. This is necessary because of how the pipeline is setup to process the chips,
    which chip is getting corrected here, and where the data is - we need information from the header
    of both imsets to process chip2 correctly.

    MLS June 24, 2015:

    Added more code to deal with correct scaling of subarray images

    M. Sosey: 27 March 2017
    Update keyword descriptions   
*/

/*prototypes*/

void FluxMsg(WF3Info *);
int doFlux(WF3Info *);
int doStat(SingleGroup *, short);
int TimeStamp(char *, char *);
int PutKeyDbl(Hdr *, char *, double, char *);
int UpdateSwitch(char *, int, Hdr*, int*);
int GetKeyDbl (Hdr *, char *, int, double, double *);

void FluxMsg (WF3Info *wf32d) {

	int OmitStep (int);
	void PrSwitch (char *, int);

	trlmessage ("");

	if (!OmitStep (wf32d->fluxcorr)) {
		PrSwitch ("fluxcorr", wf32d->fluxcorr);
	}
}


int doFlux (WF3Info *wf32d){

	/* arguments:
	   WF3Info *wf3     i: calibration switches, etc
	 */

	extern int status;
	SingleGroup chip2;
	SingleGroup chip1;
	double ratio;
	double phtflam1;
	double phtflam2;
	int logit;
	int multk1d (SingleGroupLine *a, float k);
	int i,j;


	if (wf32d->fluxcorr != PERFORM)
		return (status=CAL_STEP_NOT_DONE);

	if (wf32d->verbose){
		sprintf (MsgText, "Starting FLUXCORR");
		trlmessage(MsgText);
	}

	initSingleGroup (&chip2);
    if (wf32d->subarray == NO){
    	initSingleGroup (&chip1);
    }

	/* Open the input image. I'm reading in the entire thing here because
	   I can't get  the overwrite mode to work in hstio, and it might be a bug

	   so in order to check the output is correct, I'm writing a completely new image
	   with the updated group information

	 */

    /*If a subarray in Chip 2 is passed it needs to be scaled*/
    if (wf32d->subarray == YES){
        getSingleGroup (wf32d->output, 1, &chip2); /*chip2 is in sci,1*/
    	if (hstio_err())
		    return (status = OPEN_FAILED);

        memcpy(&phtflam1,&wf32d->chip1_flam,sizeof(double));
        sprintf(MsgText,"Copied wf332d->chip1_flam: %g to phtflam1: %g",wf32d->chip1_flam,phtflam1);
        trlmessage(MsgText);

    } else {
    	getSingleGroup (wf32d->output, 1, &chip2); /*chip2 is in sci,1*/
	    getSingleGroup (wf32d->output, 2, &chip1); /*chip1 is in sci,2*/
    	if (hstio_err())
		    return (status = OPEN_FAILED);
    	if (GetKeyDbl (&chip1.sci.hdr, "PHTFLAM1", USE_DEFAULT, 1., &phtflam1 ))
		    return (status);

    }

	if (GetKeyDbl (&chip2.sci.hdr, "PHTFLAM2", USE_DEFAULT, 1., &phtflam2 ))
		return (status);

	ratio = phtflam2/phtflam1;

	if (wf32d->verbose){
		sprintf (MsgText, "flam1 %g, flam2 %g, ratio %g", phtflam1, phtflam2, ratio);
		trlmessage(MsgText);
    	sprintf(MsgText,"Using PHTRATIO: %g ",ratio);
	    trlmessage(MsgText);
    }

	for (i=0; i < chip2.sci.data.nx ; i++) {
		for (j=0; j < chip2.sci.data.ny ; j++) {
			Pix(chip2.sci.data,i,j) *= ratio;
			Pix(chip2.err.data,i,j) *= ratio;
		}
	}


	/* Recompute min, max, mean, etc. of good science data. */
	if (doStat (&chip2, wf32d->sdqflags))
		return (status);
	TimeStamp ("Image statistics recomputed", wf32d->rootname);

    /*doing it the long way, calling all the same functions that
    putSingleGroup uses to get around the overwrite issue with hstio */

    IODescPtr out;
    Hdr phdr;
    initHdr (&phdr);

    out = openUpdateImage(wf32d->output,"",0,&phdr);
	if (PutKeyDbl (&phdr, "PHTRATIO", ratio,
				"PHTFLAM2/PHTFLAM1 ratio"))
		return (status);
	if (PutKeyDbl (&phdr, "PHTFLAM2", phtflam2,
				"Chip2 Inv Sens, use when FLUXCORR=OMIT"))
		return (status);
	if (PutKeyDbl (&phdr, "PHTFLAM1", phtflam1,
				"Chip1 Inv Sens, same as PHOTFLAM"))
		return (status);
    putHeader(out);
    closeImage(out);
    freeHdr(&phdr);


    chip2.sci.iodesc = openUpdateImage(wf32d->output, "SCI", chip2.group_num, &chip2.sci.hdr);
    putFloatData(chip2.sci.iodesc,&chip2.sci.data);
	if (PutKeyDbl (&chip2.sci.hdr, "PHTRATIO", ratio,
				"PHTFLAM2/PHTFLAM1 ratio"))
		return (status);

	if (PutKeyDbl (&chip2.sci.hdr, "PHTFLAM1", phtflam1,
				"Chip1 Inv Sens, same as PHOTFLAM"))
		return (status);
    closeImage(chip2.sci.iodesc);

    chip2.err.iodesc = openUpdateImage(wf32d->output, "ERR", chip2.group_num, &chip2.err.hdr);
    putFloatData(chip2.err.iodesc,&chip2.err.data);
    closeImage(chip2.err.iodesc);



	if (hstio_err()) {
		sprintf (MsgText, "Couldn't write imset %d.", 1);
		trlerror (MsgText);
		return (status = 1001);
	}

	if (wf32d->printtime)
		TimeStamp ("Output written to disk", wf32d->rootname);

	logit=0; /*means dont log reference file names*/
	if (UpdateSwitch ("FLUXCORR", wf32d->fluxcorr, chip2.globalhdr, &logit))
		return (status);

    /*clear out memory*/
    freeSingleGroup(&chip2);
    if (wf32d->subarray == NO)
        freeSingleGroup(&chip1);

	return (status);
}
