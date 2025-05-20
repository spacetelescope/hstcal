# include <string.h>
#include "hstcal.h"
# include "hstio.h"
# include "wf3.h"
# include "wf3info.h"
# include "hstcalerr.h"

/* These are subroutines which can be used to go between subarray
   and full-frame chip images

   Megan Sosey, July 2016
*/

int PutKeyInt(Hdr *, char *, int, char *);
int PutKeyDbl(Hdr *, char *, double , char *);
int PutKeyStr(Hdr *, char *, char *, char *);
int PutKeyBool(Hdr *, char *, int, char *);
int GetCorner (Hdr *, int, int *, int *);


int CreateEmptyChip(WF3Info *wf3, SingleGroup *full){
  /* Create a full size, but empty group, which contains necessary meta data
     alloc and init the full data with the group and sizes you want
  */

  int row, col;

  if (full->group_num == 1){
      if (PutKeyDbl(&full->sci.hdr, "LTV2", 0.0, "offset in X to light start")) {
        trlmessage("Error putting LTV1 keyword in header");
        return (status=HEADER_PROBLEM);
      }
 } else {
     if (PutKeyDbl(&full->sci.hdr, "LTV2", 19.0, "offset in X to light start")) {
       trlmessage("Error putting LTV1 keyword in header");
       return (status=HEADER_PROBLEM);
     }
 }
  if (PutKeyDbl(&full->sci.hdr, "LTV1", 25.0, "offset in X to light start")) {
    trlmessage("Error putting LTV1 keyword in header");
    return (status=HEADER_PROBLEM);
  }
  if (PutKeyDbl(&full->sci.hdr, "LTM1_1", 1.0, "reciprocal of sampling rate in X")) {
    trlmessage("Error putting LTV2 keyword in header");
    return (status=HEADER_PROBLEM);
  }
  if (PutKeyDbl(&full->sci.hdr, "LTM2_2", 1.0, "reciprocal of sampling rate in Y")) {
    trlmessage("Error putting LTV2 keyword in header");
    return (status=HEADER_PROBLEM);
  }
  if (PutKeyStr(&full->sci.hdr, "CCDAMP", wf3->ccdamp, "CCD amplifier")){
    trlmessage("Error updating CCDAMP keyword in image header");
    return (status=HEADER_PROBLEM);
  }

  /*Zero out the large arrays*/
  for(row=0; row < full->sci.data.nx; row++){
    for(col=0; col < full->sci.data.ny; col++){
      PPix(&full->sci.data,row,col) = 0.;
      PPix(&full->dq.data,row,col) = 0;
    }
  }

return(status);
}


int Sub2Full(WF3Info *wf3, SingleGroup *x, SingleGroup *full, int real_dq, int flag, int virtual) {
  /* Place a subarray into a full frame image, this only  works on one chip

  flag is used to populate the DQ array when real_dq is false
  the chip is decided by group num
  */

  int sci_bin[2];			/* bin size of science image */
  int sci_corner[2];		/* science image corner location */
  int ref_bin[2];			/* bin size of full image */
  int ref_corner[2];		/* full image corner location */
  int rsize = 1;       /* reference pixel size */
  int col = 0;
  int row = 0;
  int scix, sciy;


  /*zero dq boolean says whether to use the actual dq values or another values
  flag can be used to initialize the dq array to something other than 0
  */

  if (!wf3->subarray){
    sprintf(MsgText,"Image is not a subarray, check image...");
    trlmessage(MsgText);
    return(status = HEADER_PROBLEM);
  }

  if (wf3->detector != CCD_DETECTOR){
    sprintf(MsgText,"Sub2Full only valid for UVIS subarrays");
    trlmessage(MsgText);
  }

  /* Find the subarray location*/
  if (GetCorner (&x->sci.hdr, rsize, sci_bin, sci_corner))
    return (status);

  if (GetCorner (&full->sci.hdr, rsize, ref_bin, ref_corner))
      return (status);

  scix=sci_corner[0] - ref_corner[0] ;
  sciy=sci_corner[1] - ref_corner[1] ;

  sprintf(MsgText,"(subtools) Sci corner at x0,y0 = %i, %i", scix,sciy);
  trlmessage(MsgText);

  /*Zero out the large arrays*/
  for(row=0; row < full->sci.data.nx; row++){
    for(col=0; col < full->sci.data.ny; col++){
      PPix(&full->sci.data,row,col) = 0.;
      PPix(&full->dq.data,row,col) = 0;
    }
  }

  if(virtual){
      if (scix >= 2072){ /*image starts in B or D regions and we can just shift the starting pixel*/
          sprintf(MsgText,"Subarray starts in B or D region, moved from (%d,%d) to ",scix,sciy);
          trlmessage(MsgText);
          scix += 60;
          sprintf(MsgText,"(%d,%d) to avoid virtual overscan",scix,sciy);
          trlmessage(MsgText);
      }
  }

  /* Copy the data from the sub-array to the full-array */
  for(row=0; row < x->sci.data.nx; row++){
    for(col=0; col < x->sci.data.ny; col++){
      PPix(&full->sci.data,row+scix,col+sciy) = PPix(&x->sci.data,row,col);
      if (real_dq){
        /* if real_dq is true(one), then the actual dq values are used*/
        PPix(&full->dq.data,row+scix, col+sciy) = PPix(&x->dq.data,row,col);
      } else {
        PPix(&full->dq.data,row+scix, col+sciy) = flag;
      }
    }
  }

  /*Save scix, sciy and the extent into the full frame header?*/
  return (status);
}

int Full2Sub(WF3Info *wf3, SingleGroup *x, SingleGroup *full, int dq, int sci, int virtual) {
  /* Take a full frame image and cut out the original subarray, this only
  works on one chip
  */
  int sci_bin[2];			/* bin size of science image */
  int sci_corner[2];		/* science image corner location */
  int ref_bin[2];			/* bin size of full image */
  int ref_corner[2];		/* full image corner location */
  int rsize = 1;       /* reference pixel size */
  int col, row;
  int scix=0;
  int sciy=0;


  /*dq boolean says to  copy back the DQ values to the full array
    sci boolean says to copy the science pixels*/

  if (!wf3->subarray){
    sprintf(MsgText,"Original image is not a subarray, check image...");
    trlmessage(MsgText);
    return(status = HEADER_PROBLEM);
  }

  if (wf3->detector != CCD_DETECTOR){
    sprintf(MsgText,"Full2Sub only valid for UVIS subarrays");
    trlmessage(MsgText);
  }

  /* Find the subarray location, visible corner*/
  if (GetCorner (&x->sci.hdr, rsize, sci_bin, sci_corner))
    return (status);

  if (GetCorner (&full->sci.hdr, rsize, ref_bin, ref_corner))
      return (status);

  scix=sci_corner[0] - ref_corner[0];
  sciy=sci_corner[1] - ref_corner[1];

  sprintf(MsgText,"(subtools) Sci corner at x0,y0 = %i, %i", scix,sciy);
  trlmessage(MsgText);

  if (virtual){
      if (scix >= 2072){ /*image starts in B or D regions and we can just shift the starting pixel*/
          sprintf(MsgText,"Subarray starts in B or D region, moved from (%d,%d) to ",scix,sciy);
          trlmessage(MsgText);
          scix += 60;
          sprintf(MsgText,"(%d,%d) to avoid virtual overscan",scix,sciy);
          trlmessage(MsgText);
      }
  }
  /* Copy the data from the sub-array to the full-array */
  for(row=0; row< x->sci.data.tot_nx; row++){
    for(col=0; col< x->sci.data.tot_ny; col++){
      if (dq) {
        PPix(&x->dq.data,row,col) = PPix(&full->dq.data,row+scix,col+sciy);
      }
      if (sci) {
        PPix(&x->sci.data,row,col) = PPix(&full->sci.data,row+scix,col+sciy);
      }
    }
  }
  return (status); /*returns the reference, full is deallocated outside of function*/
}
