/* This file contains:
 doDark
 AvgSciVal
 */


# include <stdio.h>
# include <stdlib.h>		/* calloc */
# include <math.h>		/* fabs */

#include "hstcal_memory.h"
#include "hstcal.h"
# include "hstio.h"
# include "acs.h"
# include "acsinfo.h"
# include "hstcalerr.h"

/* This routine subtracts the dark image from x (in-place).
 For CCD data, the dark image is multiplied by the dark time
 before subtracting.  The dark time is the exposure time + FLASHDUR;
 for post-SM4 non-BIAS WFC, it also INCLUDES the idle time since the last
 flushing of the chip or the readout time.
 
 For MAMA data, the dark image is just multiplied by the exposure time
 before subtracting.
 
 Reference image should have been selected to have
 the same binning factor as the science image, so
 assume ratio of bin factors to be 1.
 
 The value of MEANDARK is calculated based on the weighted average
 of the image where the weighting is based on the percent of 
 good pixels in the image.  Only pixels not flagged BAD (in some way)
 will contribute to the average.
 
 Warren Hack, 1998 June 11:
 Initial ACS Version.
 Warren Hack, 1999 July 27:
 Added weighting to calculation of MEANDARK.  Use same code for all bin
 cases.
 Warren Hack, 2000 April 14:
 Fixed the code which applies the gain to each line so that if
 there is only 1 AMP used per line, it doesn't try to apply a 
 second amp with gain = 0.
 Warren Hack, 2001 April 16:
 Revised calling sequence for 'multgn1d' to be more general
 and accept arbitrary scaling factor instead of being hardwired
 to 'exptime'.
 Warren Hack, 2001 December 4:
 Compute darktime from exposure time keywords EXPSTART,EXPEND
 and use it instead of EXPTIME for scaling dark image.
 Warren Hack, 2002 Apr 24:
 Reverted back to using EXPTIME for scaling dark image.
 Pey Lian Lim, 2012 Dec 11:
 Now use DARKTIME for scaling dark image.
 M.D. De La Pena, 2018 June 05:
 Dark correction processing is now done on the entire image instead of 
 line by line. Also, performed some clean up and new commentary.
 M.D. De La Pena, 2019 November 26:
 Removed hard-coded darktime scaling value and read new 
 post-flashed and unflashed columns from updated CCDTAB reference 
 file to use for the offset to the DARKTIME FITS keyword value 
 under appropropriate date and supported subarray criteria.  The 
 DARKTIME keyword value is now the default scaling factor for the 
 superdarks, with the offset being an additive correction to DARKTIME 
 under appropriate circumstances.
 
 */
static const char *subApertures[] = {"WFC1A-2K", "WFC1B-2K", "WFC2C-2K", "WFC2D-2K",
                                     "WFC1A-1K", "WFC1B-1K", "WFC2C-1K", "WFC2D-1K",
                                     "WFC1A-512", "WFC1B-512", "WFC2C-512", "WFC2D-512",
                                     "WFC1-IRAMPQ", "WFC1-MRAMPQ", "WFC2-ORAMPQ",
                                     "WFC1-POL0UV", "WFC1-POL0V", 
                                     "WFC1-SMFL"};
static const int numSupportedSubApertures = sizeof(subApertures) / sizeof(subApertures[0]);

int doDark (ACSInfo *acs2d, SingleGroup *x, float *meandark) {
  
  /* arguments:
   ACSInfo *acs     i: calibration switches, etc
   SingleGroup *x  io: image to be calibrated; written to in-place
   float *meandark  o: mean of dark image values subtracted
  */
  
  extern int status;

  int extver = 1;	/* get this imset from dark image */

  /* Assumption is the science and dark images are the same size and bin factor */
  int rx = 1, ry = 1;  /* for binning dark image down to size of x */
  int x0 = 0, y0 = 0;  /* offsets of science image */
  int same_size = 1;   /* true if no binning of dark image required */

  int scicols;      /* number of columns in science image */
  int scirows;      /* number of rows in science image */
  double mean;      /* mean value for the image */
  double weight;    /* weight value for the image (No. good pixels/all pixels) */
  int update;
  float darktimeFromHeader;  /* base darktime value based upon FITS keyword DARKTIME from image header */
  float darktime;        /* final darktime value after applicable offset, if any */
  float darktimeOffset;  /* overhead offset time (s) based upon full-frame/subarray and post-flashed/unflashed */
  
  int DetCCDChip (char *, int, int, int *);

  /* 2D routines based on 1D counterparts - process in 2D instead of line by line */
  int trim2d (SingleGroup *, int, int, int, int, int, SingleGroup *);
  int sub2d (SingleGroup *, SingleGroup *);
  int multk2d (SingleGroup *, float);
  void AvgSciVal (SingleGroup *, short, double *, double *);
  int FindBin (SingleGroup *, SingleGroup *, int *, int *, int *, int *, int *);
  
  /* Get the dimensions of the science data */
  scicols = x->sci.data.nx;
  scirows = x->sci.data.ny; 

  /* Get the DARKTIME FITS keyword value stored in the ACSInfo structure */
  darktimeFromHeader = (float)acs2d->darktime;

  /*
     The overhead offset time is a function of full-frame vs subarray and post-flash vs unflashed 
     and has been determined empirically for CCD data.  The appropriate overhead for full-frame or
     subarray was extracted from the calibration file during the table read.
  */

  /* Unflashed observation */
  darktimeOffset = acs2d->overhead_unflashed;
 
  /* Post-flashed observation */
  if ((acs2d->flashdur > 0.0) && (strcmp(acs2d->flashstatus, "SUCCESSFUL") == 0)) {
        darktimeOffset = acs2d->overhead_postflashed;
  }

  sprintf(MsgText, "Darktime from header: %f Exp: %f FlashDur: %f darktimeOffset: %f\n", darktimeFromHeader, acs2d->exptime,
      acs2d->flashdur, darktimeOffset);
  trlmessage(MsgText);

  /* 
     Compute the final darktime based upon the date of full-frame or subarray data.
     The full-frame overhead offset is applicable to all data post-SM4.  The subarray 
     overhead offset is applicable to all data post-CYCLE24 and only for supported
     subarrays. 
  */
  darktime = darktimeFromHeader;  /* Default */
  if (acs2d->detector != MAMA_DETECTOR) {

     /* Full-frame data */
     if (acs2d->subarray == NO) {
         if (acs2d->expstart >= SM4MJD)
              darktime = darktimeFromHeader + darktimeOffset;

         sprintf(MsgText, "Full Frame darktime: %f\n", darktime);
         trlmessage(MsgText);

     /* Subarray data */
     } else {
         if (acs2d->expstart >= CYCLE24) {
             for (unsigned int i = 0; i < numSupportedSubApertures; i++) {
                 if (strcmp(acs2d->aperture, subApertures[i]) == 0) {
                     darktime = darktimeFromHeader + darktimeOffset;
                     sprintf(MsgText, "Supported Subarray darktime: %f aperture: %s\n", darktime, subApertures[i]);
                     trlmessage(MsgText);
                     break;
                 }
                 sprintf(MsgText, "Supported Subarray aperture: %s\n", subApertures[i]);
                 trlmessage(MsgText);
             }
         }
     }
  }
        
  sprintf(MsgText, "Final darktime: %f", darktime);
  trlmessage(MsgText);

  /* Compute correct extension version number to extract from
     reference image to correspond to CHIP in science data.
  */
  if (acs2d->pctecorr == PERFORM) {
     if (DetCCDChip (acs2d->darkcte.name, acs2d->chip, acs2d->nimsets, &extver) )
        return (status);
  } else {
     if (DetCCDChip (acs2d->dark.name, acs2d->chip, acs2d->nimsets, &extver) )
        return (status);
  }
	
  if (acs2d->verbose) {
     sprintf(MsgText,"Performing dark subtraction on chip %d in imset %d",acs2d->chip, extver);
     trlmessage(MsgText);
  }
  
  /* Create a pointer register for bookkeeping purposes and to ease cleanup */
  PtrRegister ptrReg;
  initPtrRegister(&ptrReg);

  /* Initialize and get the dark image data */
  SingleGroup y; /* scratch space */
  initSingleGroup (&y);

  addPtr(&ptrReg, &y, &freeSingleGroup);

  if (acs2d->pctecorr == PERFORM) {
      getSingleGroup (acs2d->darkcte.name, extver, &y);
  } else {
      getSingleGroup (acs2d->dark.name, extver, &y);
  }

  if (hstio_err()) {
     freeOnExit(&ptrReg);
     return (status = OPEN_FAILED);
  }

  /* Compare binning of science image and reference image; get same_size 
     and high_res flags, and get info about binning and offset for 
     use by bin2d.
  */
  if (FindBin (x, &y, &same_size, &rx, &ry, &x0, &y0)) {
     freeOnExit(&ptrReg);
     return (status);
  }
  
  if (rx != 1 || ry != 1) {
     sprintf(MsgText,"Reference image and input are not binned to the same pixel size!");
     trlmessage(MsgText);
  }
  if (acs2d->verbose){
     sprintf(MsgText,"Image has an offset of %d,%d",x0,y0);
     trlmessage(MsgText);
  }
  
  /* Trim the dark image (y->z) down to the actual size of the science image (x). */
  SingleGroup z; /* scratch space */
  initSingleGroup (&z);
  allocSingleGroup (&z, scicols, scirows, True);

  addPtr(&ptrReg, &z, &freeSingleGroup);

  if (hstio_err()) {
     freeOnExit(&ptrReg);
     return (status = ALLOCATION_PROBLEM);
  }
    
  /* We are working with a sub-array and need to apply the
     proper section from the reference image to the science image.
     update = NO = do not update the header of the scratch image
  */
  update = NO;
  if (trim2d (&y, x0, y0, rx, ry, update, &z)) {
     trlerror ("(darkcorr) size mismatch.");
     freeOnExit(&ptrReg);
     return (status);
  }
    
  /* Multipy the dark image (SCI and ERR extensions) by the darktime constant. */
  if (multk2d(&z, darktime)) {
     freeOnExit(&ptrReg);
     return (status);
  }
    
  /* Compute the average of all the good pixels in the science image, as well as 
     a weight = (number of good pixels / total number of pixels).
  */
  mean   = 0.0;
  weight = 0.0;
  AvgSciVal (&z, acs2d->sdqflags, &mean, &weight);

  /* Subtract the dark data from the science data */
  if (sub2d (x, &z)) {
     freeOnExit(&ptrReg);
     return (status);
  }

  /* Compute the weighted mean for the image */	
  /* *** MDD FIX - Not sure this is really needed anymore */
  /*
  *meandark = 0.0;
  if ( (weight > 0.0) && (scirows > 0) )
     *meandark = (float) (mean / weight); 
  */
  /* This is to force a compatibility match to the previous version of the
     code.  
  */
  *meandark = mean;
	
  /* Free up the scratch data */
  freeOnExit(&ptrReg);

  return (status);
}
