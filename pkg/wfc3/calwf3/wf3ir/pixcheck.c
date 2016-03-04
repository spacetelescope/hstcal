# include "hstio.h"	/* defines HST I/O functions */
# include "wf3.h"

/* PIXOK: Function to check whether or not an image pixel is OK, based
** on its DQ flag settings. Returns boolean False if pixel is bad; True
** if pixel is OK.
**
** Revision history:
** H.Bushouse	Oct. 2000	Initial CALNICA to CALWF3 port.
** H.Bushouse	30-Sep-2010	Updated to WFC3 IR DQ assignments.
**				(PR #66080, Trac #607)
*/

/* The following bitmask value corresponds to these DQ settings:
**  Reed-Solomon error  =     1
**  Replaced by fill    =     2
**  Bad/dead pixel      =     4
**  Deviant zero read   =     8
**  Hot pixel           =    16 
**  Unstable pixel      =    32
**  Warm pixel          =    64
**  Bad reference pixel =   128
**  Saturated           =   256
**  Bad in flat (blobs) =   512
**  unassigned          =  1024 
**  Zero-read signal    =  2048 
**  MultiDrizzle CR     =  4096
**  calwf3 CR           =  8192
**  Crosstalk           = 16384
**                        -----
*/
# define BITMASK 260

Bool pixOK (SingleNicmosGroup *im, int i, int j) {

/* Arguments:
** 	im	i: input image
**	i	i: pixel coord
**	j	i: pixel coord
*/


	/* Test the pixel DQ value against the BitMask */
	if (DQPix(im->dq.data,i,j) & BITMASK)
	    return (False);

	/* No bad bits set; return True */
	return (True);
}

