# include <stdio.h>
# include "xtables.h"


# include "stis.h"
# include "calstis6.h"


/* 
   This routine adds together the offsets:  MAMA "dither", wavecal,
   and aperture. The aperture offset will be converted from arcseconds
   to pixels. Both the total and aperture offsets will be saved in sts.



   Revision history:
   ----------------
   20 Feb 97  -  Adapted from a similar routine in calstis7 (I.Busko)
   24 Feb 97  -  Renamed routine to avoid conflict with cs7 (IB)
   18 Mar 97  -  Deleted all offset computations for the dispersion
                 axis ([0]), they are never used (IB)
*/

void AddOffsets6 (StisInfo6 *sts, ApInfo *slit) {

/* arguments:
StisInfo6 *sts      i: calibration switches and info
ApInfo *slit        i: description of slit (aperture offset)
*/
	sts->ap_offset[1] = slit->ap_offset[1] / (sts->cd[1] * 3600.);
	sts->total_offset[1] = sts->dither_offset[1] +
                               sts->msm_offset[1] + sts->ap_offset[1];
}
