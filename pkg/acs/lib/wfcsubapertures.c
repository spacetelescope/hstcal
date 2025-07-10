#include <stdio.h>
#include "wfcsubapertures.h"

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

const char *darkSubApertures[] = {
    "WFC1A-2K", "WFC1B-2K", "WFC2C-2K", "WFC2D-2K",
    "WFC1A-1K", "WFC1B-1K", "WFC2C-1K", "WFC2D-1K",
    "WFC1A-512", "WFC1B-512", "WFC2C-512", "WFC2D-512",
    "WFC1-IRAMPQ", "WFC1-MRAMPQ", "WFC2-ORAMPQ",
    "WFC1-POL0UV", "WFC1-POL0V", "WFC1-POL60UV", "WFC1-POL60V",
    "WFC1-POL120UV", "WFC1-POL120V", "WFC1-SMFL"};
const size_t numSupportedDarkSubApertures = sizeof(darkSubApertures) / sizeof(darkSubApertures[0]);
