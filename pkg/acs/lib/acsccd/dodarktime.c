#include <string.h>
#include "acs.h"
#include "acsinfo.h"


/* This routine updates the DARKTIME value such that every BLV_TMP file will
   have the correct DARKTIME keyword.

   DARKTIME has to be calculated in ACSCDD step and not where DARKCORR happens in ACS2D.
   See 24-May-2021 MDD note in doccd.c file.
*/
void computeDarktime(ACSInfo *acs, float *darktime) {
    const char *darkSubApertures[] = {
        "WFC1A-2K", "WFC1B-2K", "WFC2C-2K", "WFC2D-2K",
        "WFC1A-1K", "WFC1B-1K", "WFC2C-1K", "WFC2D-1K",
        "WFC1A-512", "WFC1B-512", "WFC2C-512", "WFC2D-512",
        "WFC1-IRAMPQ", "WFC1-MRAMPQ", "WFC2-ORAMPQ",
        "WFC1-POL0UV", "WFC1-POL0V", "WFC1-POL60UV", "WFC1-POL60V",
        "WFC1-POL120UV", "WFC1-POL120V", "WFC1-SMFL"};
    const size_t numSupportedDarkSubApertures = sizeof(darkSubApertures) / sizeof(darkSubApertures[0]);

    float darktimeFromHeader;  /* base darktime value from FITS keyword DARKTIME in image header */
    float darktimeOffset;      /* overhead offset time (s) based upon full-frame/subarray and post-flashed/unflashed */

    /* Get the DARKTIME FITS keyword value stored in the ACSInfo structure */
    darktimeFromHeader = (float)acs->darktime;

    /*
       The overhead offset time is a function of full-frame vs subarray and post-flash vs unflashed
       and has been determined empirically for CCD data.  Both the unflashed and post-flash
       overhead values for full-frame or subarray were extracted from the calibration file
       during the table read.  Now it is just necessary to determine which offset actually
       applies.
    */

    /* Unflashed observation */
    darktimeOffset = acs->overhead_unflashed;

    /* Post-flashed observation */
    if ((acs->flashdur > 0.0) && (strcmp(acs->flashstatus, "SUCCESSFUL") == 0)) {
          darktimeOffset = acs->overhead_postflashed;
    }

    /*
       Compute the final darktime based upon the date of full-frame or subarray data.
       The full-frame overhead offset is applicable to all data post-SM4.  The subarray
       overhead offset is applicable to all data post-CYCLE24 and ONLY for supported
       subarrays.

       Effectively the additive factor only applies to ACS/WFC as the HRC was no longer
       operational by SM4MJD or CYCLE24, and SBC is a MAMA detector.
    */
    *darktime = darktimeFromHeader;  /* Default */

    /* Full-frame data */
    if (acs->subarray == NO) {
      if (acs->expstart >= SM4MJD) {
        *darktime = darktimeFromHeader + darktimeOffset;
      }
      trlmessage("Full Frame adjusted Darktime: %f", *darktime);

    /* Subarray data */
    } else {
      if (acs->expstart >= CYCLE24) {
        for (unsigned int i = 0; i < numSupportedDarkSubApertures; i++) {
            if (strcmp(acs->aperture, darkSubApertures[i]) == 0) {
                *darktime = darktimeFromHeader + darktimeOffset;
                trlmessage("Supported Subarray adjusted Darktime: %f for aperture: %s", *darktime, darkSubApertures[i]);
                break;
            }
        }
      }
    }
    trlmessage("DARKTIME from SCI header: %f  Offset from CCDTAB: %f  Final DARKTIME: %f", darktimeFromHeader, darktimeOffset, *darktime);
}