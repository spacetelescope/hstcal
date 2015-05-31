/* WFC3 -- Detect and mark SINK pixels in the DQ array

    This is done  for the CTE  AND  non-CTE corrected data.

    go through the reference image pixel by pixel, if the reference image has
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
    in the science expsoure at hand. If the value of the sink pixel in the exposure is below the value
    of the upstead pixel in the reference image we flag that pixel wth 1024 in the DQ file of the exposure.

    We continue to flag until there are no more flagable pixels for this sink pixel (until the 
    value of the pixel in the reference image is zero), or until the value of the sink pixel in the exposure
    is greater than the value of the upstream pixel in the reference image

    The pixel mask is saved to the group DQ image in the raz which is passed

    As long as the un-raz'ed DQ information makes it to the BLC_TMP file that wf3cte saves
    we should be good since the rest of the code uses a logical or for all the DQ 
    math.

    I'm running this code after wf3ccd has been executed for both the non-cte and cte data.
    It will take the input science filename and work on BOTH chips at the same time, upating
    the DQ arrays separately.
    
    MLS, May 2015

*/

# include <time.h>
# include <string.h>
# include "wf3dq.h"
# include "hstio.h"
# include "wf3.h"
# include "wf3info.h"
# include "wf3corr.h"		/* calibration switch names */

int makedqRAZ(SingleGroup *, SingleGroup *, SingleGroup *);
int undodqRAZ(SingleGroup *, SingleGroup *, SingleGroup *);
int makeFloatRaz(FloatTwoDArray *, FloatTwoDArray *, FloatTwoDArray  *);
int getFloatHD(char *, char *, int , FloatHdrData *);

int SinkDetect(WF3Info *wf3){

    extern int status;
    int i,j, jj;
    short dqval;
    float refdate;

    trlmessage("\nPerforming SINK pixel detection and flagging");
    
    refdate=51544.; /*Year 2000*/
    if (wf3->verbose){
        sprintf(MsgText,"Reference date for Sink Pixels: %f",refdate);
        trlmessage(MsgText);
    }

    /*THE SCIENCE IMAGE*/
    SingleGroup raz; /*quad rotated image to work with*/
    SingleGroup cd; /* input image */
    SingleGroup ab; /*input image */
    

    /* INIT THE SCIENCE INPUT  */
    initSingleGroup (&cd);
    initSingleGroup (&ab);
    initSingleGroup (&raz);
    allocSingleGroup (&raz,RAZ_COLS, RAZ_ROWS);
    
    trlmessage("Finished allocating sink pixel image memory");
    
    getSingleGroup (wf3->input, 1, &cd);
    if (hstio_err())
        return (status = OPEN_FAILED);

    getSingleGroup (wf3->input, 2, &ab);
    if (hstio_err())
        return (status = OPEN_FAILED);
        
    /*CONVERT DQ DATA TO RAZ FORMAT FOR SCIENCE FILE*/
    makedqRAZ(&cd, &ab, &raz);
    trlmessage("Converted science image to RAZ format");


	/* GET THE SINK FILE REFERENCE IMAGE FROM SINKFILE AND INITIALIZE */
    
    FloatHdrData sinkrefab;
    FloatHdrData sinkrefcd;
    initFloatHdrData(&sinkrefab);
    initFloatHdrData(&sinkrefcd);

    getFloatHD(wf3->sink.name,"SCI",1,&sinkrefcd);
    getFloatHD(wf3->sink.name,"SCI",2,&sinkrefab);    
    
    sprintf(MsgText,"Finished reading Sinkpixel file: %s",wf3->sink.name);
    trlmessage(MsgText);
        
     
    /*NOW TURN THE SINK REFERENCE IMAGES INTO RAZ FORMAT*/
    FloatTwoDArray sinkraz;
    initFloatData(&sinkraz); /*float 2d arrays*/
    allocFloatData(&sinkraz,RAZ_COLS, RAZ_ROWS);     
    makeFloatRaz(&sinkrefcd.data,&sinkrefab.data,&sinkraz);

    trlmessage("Turned sink pixel mask into RAZ format");
    
    /*THE MJD OF THE SCIENCE EXPOSURE IS THE COMPARISON DATE
     THE FOLLOWING TRANSLATION TAKEN FROM ISR WFC3-2014-22.PDF */
    
    
    for (i=1;i<RAZ_COLS-1;i++){
        for (j=1; j<RAZ_ROWS-1; j++){
            if ( PPix(&sinkraz,i,j) >  refdate && wf3->expstart > PPix(&sinkraz,i,j) ){
                /*FLAG THE PRIMARY SINK PIXEL*/
                dqval = TRAP | DQPix (raz.dq.data, i, j);
                DQSetPix (raz.dq.data, i, j, dqval);
            
                if (PPix(&sinkraz,i,j-1) != 0){
                    /*FLAG THE DOWNSTREAM PIXEL*/
                    dqval = TRAP | DQPix (raz.dq.data, i, j);
                    DQSetPix (raz.dq.data, i, j, dqval);
                }

                for (jj=j+1; jj<RAZ_ROWS; jj++){
                    if (PPix(&sinkraz,i,j) != 0){
                        if (Pix(raz.dq.data,i,j) <= PPix(&sinkraz,i,jj) ){
                        /*FLAG THIS UPSTREAM PIXEL*/
                        dqval = TRAP | DQPix (raz.dq.data, i, jj);
                        DQSetPix (raz.dq.data, i, jj, dqval);
                        }
                
                    }
                }
                
            } /*end if*/ 
        } /*end j*/
    }/*end i*/   
    
    trlmessage("Undoing RAZ for science image");
    undodqRAZ(&cd,&ab,&raz);

    freeSingleGroup(&ab);
    freeSingleGroup(&cd);
    freeSingleGroup(&raz);
    freeFloatData(&sinkraz);
    freeFloatHdrData(&sinkrefab);
    freeFloatHdrData(&sinkrefcd);
}



