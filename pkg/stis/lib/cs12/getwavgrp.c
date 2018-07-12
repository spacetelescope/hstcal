# include <stdio.h>
# include "c_iraf.h"
# include "hstio.h"

# include "stis.h"
# include "calstis12.h"
# include "stisdef.h"

/* This routine gets the time of the middle of the exposure and the
   shift for each axis from an extension header.
*/

int GetWavGrp (StisInfo12 *wav, Hdr *hdr, int extver) {

/* arguments:
StisInfo12 *wav   io: info for wavecal
Hdr *hdr          i: header of current imset
int extver        i: used to find index for storing info in wav
*/

	int status;

	double expstart, expend;
	int index;			/* index for arrays, extver - 1 */
	int use_def = 0;		/* missing keyword is fatal error */

	index = extver - 1;

	if ((status = Get_KeyD (hdr, "EXPSTART", use_def, 0., &expstart)))
	    return (status);
	if ((status = Get_KeyD (hdr, "EXPEND", use_def, 0., &expend)))
	    return (status);

	wav->midpt[index] = (expstart + expend) / 2.;

	if ((status = Get_KeyD (hdr, "SHIFTA1",
                                use_def, 0., &wav->shift1[index])))
	    return (status);
	if ((status = Get_KeyD (hdr, "SHIFTA2",
                                use_def, 0., &wav->shift2[index])))
	    return (status);

	return (0);
}
