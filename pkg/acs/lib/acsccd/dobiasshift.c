#include <string.h>
#include "acs.h"
#include "acsinfo.h"


/* If this is subarray data, determine if it is supported for the bias shift
   correction.  Compare the aperture name against the array of supported
   apertures.  Use the aperture name and not the date as the discriminant
   as the date is fuzzy - some supported subarrays were acquired before the
   official name change, and some old style subarrays were acquired after the
   change.  This scheme works fine for non-polarization and non-ramp data as the
   "*-2K" apertures have both physical and virtual overscan.  These are needed
   to compute the "fat zero" (aka magic square) value needed for the bias
   shift correction.  The "*-1K and *-512" subarrays do not have virtual
   overscan which is needed at this time for fat zero.

   For polarization and ramp subarrays, the aperture names did not change.  For
   this data, the dimensions of the image also have to be checked to determine
   if there are both physical and virtual overscan sections so the bias shift
   correction can be applied.
*/
int isValidBiasShiftSubArrWithVirtOscn(int is_subarray, char *aperture, int virt_oscn_value) {
    if (is_subarray != YES)
        return NO;

    /* This variable contains the subarray aperture names for the officially supported
       subarrays since ~Cycle 24.  The polarizer and ramp subarray apertures did NOT
       change aperture names, even though the dimensions of the data changed to accommodate
       overscan.  In addition to checking the aperture name, the dimensions will also need
       to be checked for pol/ramp data to determine if the bias shift correction
       can be applied.

       The non-pol and non-ramp data adopted new aperture names for easy identification.
       These apertures contain physical overscan.  The new "-2K" subarrays also have virtual
       overscan so the bias shift correction can be applied to this data. At this time,
       bias shift correction can only be applied to data with both physical and virtual
       overscan so the "fat zero" value can be computed.
    */
    const char *subApertures[] = {
        "WFC1A-2K", "WFC1B-2K", "WFC2C-2K", "WFC2D-2K",
        "WFC1-IRAMPQ", "WFC1-MRAMPQ", "WFC2-MRAMPQ", "WFC2-ORAMPQ",
        "WFC1-POL0UV", "WFC1-POL0V", "WFC2-POL0UV", "WFC2-POL0V",
        "WFC1-POL60UV", "WFC1-POL60V", "WFC2-POL60UV", "WFC2-POL60V",
        "WFC1-POL120UV", "WFC1-POL120V", "WFC2-POL120UV", "WFC2-POL120V"
    };
    const size_t numSupportedSubApertures = sizeof(subApertures) / sizeof(subApertures[0]);
    Bool isSupportedSubAperture = NO;

    for (size_t i = 0; i < numSupportedSubApertures; i++) {
        /* If this is polarization or ramp data, then the size of the science
           array must also be checked to make sure there is physical and
           virtual overscan data.  This does not need to be checked explicitly
           for non-polarization or non-ramp data, but the comparison works for
           this data too.  Only need to check one virtOverscan value in the array.

           NOTE: Just because the aperture is supported does not mean the bias
           shift correction will be applied.  The DS_int and gain criteria also
           have to be satisfied.
        */
        if ((strcmp(aperture, subApertures[i]) == 0) && virt_oscn_value) {
            isSupportedSubAperture = YES;
            trlmessage("This subarray data is supported for the bias shift correction.");
            break;
        }
    }

    return isSupportedSubAperture;
}
