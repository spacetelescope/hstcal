# include <stdio.h>
# include <string.h>

# include "hstio.h"
# include "wf3.h"
# include "wf3info.h"
# include "wf3err.h"

/* This routine scales CHIP2 in the UVIS data so that the flux correction over both chips
   is uniform. It uses the ratio of PHTFLAM1 / PHTFLAM2 and saves that value for reference
   in the header as PHTRATIO. This will be performed by default in the pipeline and the
   previous PHOTFLAM keyword will be valid for both chips after the correction.
   
   If users don't wish to perform the correction then they set FLUXCORR OMIT and can use
   the values of PHTFLAM1 and PHTFLAM2 to correct the flux in the respective chips
   
   MLS Dec 6, 2013
   
 */

int doFlux (WF3Info *wf32d, SingleGroup *x) {

    /* arguments:
       WF3Info *wf3     i: calibration switches, etc
       SingleGroup *x	io: image to be calibrated; written to in-place
     */

    extern int status;

    float ratio;
    int multk1d (SingleGroupLine *a, float k);
	int GetKeyFlt (Hdr *, char *, int, float, float *);
    int i,j;
	SingleGroupLine y, z;	/* y and z are scratch space */
    
    if (wf32d->fluxcorr != PERFORM)
        return (status);

   /* Get the values of PHTFLAM* from header which was already updated in dophot*/
	if (GetKeyFlt (&x->sci.hdr, "PHTRATIO",  USE_DEFAULT, 1., &ratio))
	    return (status); 
        
    sprintf(MsgText,"Using PHTRATIO from header: %f",ratio);
    trlmessage(MsgText);

    
    for (i=0; i < x->sci.data.nx ; i++) {   
         for (j=0; j < x->sci.data.ny ; j++) { 
            Pix(x->sci.data,i,j) *= ratio;
            Pix(x->err.data,i,j) *= ratio;
        }
    }

    return (status);
}

