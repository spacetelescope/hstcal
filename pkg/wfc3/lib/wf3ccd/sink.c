/* WFC3 -- Detect and mark SINK pixels in the DQ array

    This is done  for the CTE  AND  non-CTE corrected data.
    See WFC3-2014-22 ISR for more information.
    Ignore the example code in that  paper, it doesn't match what the text
    says to do and it's wrong.
    
    Go through the reference image pixel by pixel, if the reference image has
    a pixel value > 999 then that indicates that this is a sink pixel with 
    its value representing the turn on date. If the turn on date of the pixel is AFTER 
    the exposure date, then we ignore the pixel in this science exposure and move on to the next
    sink pixel.

    if the sink pixels turn on date is BEFORE the exposure date, then the pixel was compromised 
    at the time of the observation. 

    The charge-trap flag 1024 will be used in the DQ file, and is designated with TRAP.

    If the pixel beneath the sink pixel in the reference file has a value of -1, then
    we activate the charge trap flag for that pixel as well. We then proceed vertically up
    from the sink pixel and compare each pixel in the reference file to the value of the sink pixel
    in the science expsoure at hand, not the pixel which corresponds to the trail pixel.  If the 
    value of the sink pixel in the exposure is below the value of the upstead pixel in the reference 
    image we flag that pixel wth 1024 in the DQ file of the exposure.

    We continue to flag until there are no more flagable pixels for this sink pixel (until the 
    value of the pixel in the reference image is zero), or until the value of the sink pixel in the exposure
    is greater than the value of the upstream pixel in the reference image

    The pixel mask is saved to the group DQ image in the raz which is passed

    As long as the un-raz'ed DQ information makes it to the BLC_TMP file that wf3cte saves
    we should be good since the rest of the code uses a logical OR for all the DQ 
    math.

    I'm running this code after wf3ccd has been executed for both the non-cte and cte data.
    It will take the input science filename and work on BOTH chips at the same time, upating
    the DQ arrays separately.
    
    MLS, May 2015

    MLS, June 1, 2015. oops. The reference file the team is giving me is in the raw pre-overscan
    trimmed format, so this has to be performed AFTER doDQI, before BLEVCORR and after BIASCORR.

*/

# include <time.h>
# include <string.h>
#include "hstcal.h"
# include "wf3dq.h"
# include "hstio.h"
# include "wf3.h"
# include "wf3info.h"
# include "wf3corr.h"		/* calibration switch names */

int makedqRAZ(SingleGroup *, SingleGroup *);
int makeSciSingleRAZ(SingleGroup *, SingleGroup *);
int undodqRAZ(SingleGroup *, SingleGroup *);
int makeFloatRaz(FloatTwoDArray *, FloatTwoDArray  *, int);
int getFloatHD(char *, char *, int , FloatHdrData *);

int SinkDetect(WF3Info *wf3, SingleGroup *x){

    extern int status;
    int i,j, jj;
    short dqval=0;
    float scipix; /*to save the value of the science pixel*/
    float refdate=50000.;
    int keep_going=1;
    
    sprintf(MsgText,"\nPerforming SINK pixel detection for imset %i",x->group_num);
    trlmessage(MsgText);
    

    /*THE SCIENCE IMAGE*/
    SingleGroup raz; /*quad rotated image to work with*/    

    /* INIT THE SCIENCE INPUT  */
    initSingleGroup (&raz);
    allocSingleGroup (&raz,RAZ_COLS/2, RAZ_ROWS, True);
        
    /*CONVERT DQ DATA TO RAZ FORMAT FOR SCIENCE FILE*/
    makedqRAZ(x, &raz);
    makeSciSingleRAZ(x, &raz);

	/* GET THE SINK FILE REFERENCE IMAGE FROM SINKFILE AND INITIALIZE */
    FloatHdrData sinkref;
    initFloatHdrData(&sinkref);
    getFloatHD(wf3->sink.name,"SCI",x->group_num,&sinkref);
                 
    /*NOW TURN THE SINK REFERENCE IMAGES INTO RAZ FORMAT*/
    FloatTwoDArray sinkraz;    
    initFloatData(&sinkraz); /*float 2d arrays*/
    allocFloatData(&sinkraz,RAZ_COLS/2, RAZ_ROWS, True);

    makeFloatRaz(&sinkref.data,&sinkraz,x->group_num);

       
    /*THE MJD OF THE SCIENCE EXPOSURE IS THE COMPARISON DATE
     THE FOLLOWING TRANSLATION TAKEN FROM ISR WFC3-2014-22.PDF */
    
    scipix=0.;
    for (i=0;i<(RAZ_COLS/2);i++){
        for (j=0; j<RAZ_ROWS; j++){
        
            if (  (PPix(&sinkraz,i,j) > refdate)  &&  ( wf3->expstart > PPix(&sinkraz,i,j))  ){
                keep_going=1;
                
                /*FLAG THE PRIMARY SINK PIXEL*/
                dqval = TRAP | DQPix (raz.dq.data, i, j);
                DQSetPix (raz.dq.data, i, j, dqval);
                scipix = Pix(raz.sci.data,i,j);
                
                /*FLAG THE DOWNSTREAM PIXEL*/
                if (j>0 && PPix(&sinkraz,i,j-1) < 0 ){
                    dqval = TRAP | DQPix (raz.dq.data, i, j-1);
                    DQSetPix (raz.dq.data, i, j-1, dqval);
                }

                /*FLAG THE UPSTREAM PIXELS*/
                for (jj=j+1; jj<RAZ_ROWS; jj++){
                    if ((int) PPix(&sinkraz,i,jj) == 0)
                        keep_going=0;
                    if ( PPix(&sinkraz,i,jj) > refdate)
                        keep_going=0;
                    if ( 0. < PPix(&sinkraz,i,jj) &&  PPix(&sinkraz,i,jj) < 1000. && keep_going){
                        if (scipix <= PPix(&sinkraz,i,jj) ){
                           dqval = TRAP | DQPix (raz.dq.data, i, jj);
                           DQSetPix (raz.dq.data, i, jj, dqval);
                        }                
                    } else {
                        keep_going=0;
                    }
                }                
            } /*end if*/ 
        } /*end j*/
    }/*end i*/   

    /*format the dq data back to expected orientation*/
    undodqRAZ(x,&raz);

    freeSingleGroup(&raz);
    freeFloatData(&sinkraz);
    freeFloatHdrData(&sinkref);
    trlmessage("Sink pixel flagging complete");
    return(status);
}



