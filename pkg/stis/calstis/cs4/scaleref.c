# include <stdio.h>
# include "stis.h"
# include "calstis4.h"

/* This routine multiplies the shift for each axis by the scale factor
   to convert the shift from image pixels to reference pixels.

   If the shift in either axis is undefined, it will not be modified.
   The dispersion axis must be either one or two; if not, the shifts
   will not be changed.
*/

void ScaleRef (StisInfo4 *sts, double *w_shift, double *s_shift) {

/* arguments:
StisInfo4 *sts    i: calibration switches and info
double *w_shift   io: the shift in the dispersion direction
double *s_shift   io: the shift in the spatial direction
*/

	if (sts->dispaxis == 1) {

	    if (*w_shift != UNDEFINED_SHIFT)
		*w_shift *= sts->scale[0];

	    if (*s_shift != UNDEFINED_SHIFT)
		*s_shift *= sts->scale[1];

	} else if (sts->dispaxis == 2) {

	    if (*w_shift != UNDEFINED_SHIFT)
		*w_shift *= sts->scale[1];

	    if (*s_shift != UNDEFINED_SHIFT)
		*s_shift *= sts->scale[0];
	}
}
