# include "hstio.h"     /* defines HST I/O functions */
# include "wf3.h"
# include "wf3dq.h"

/* SATCHECK: Flag pixels as saturated in a MultiAccum group if they're
** flagged as such in the preceding group.
**
** Revision history:
** H.Bushouse	Oct. 2000	Initial CALNICA to CALWF3 port.
** H.Bushouse	08-May-2002	Added use of wf3dq.h and changed DQ flag macro
**				"SATURATED" to "SATPIXEL".
*/

void satcheck (SingleNicmosGroup *group1, SingleNicmosGroup *group2) {

/* Arguments:
**	group1	 i: first image group
**	group2	io: second image group
*/

	/* Local variables */
	int i, j;		/* loop indexes */

	/* Loop through the DQ image of group 1 */
	for (j=0; j<group1->dq.data.ny; j++) {
	     for (i=0; i<group1->dq.data.nx; i++) {

		  /* If a pixel has a saturation flag, make sure the
		  ** flag is also set in the next group */

		  if (DQPix(group1->dq.data,i,j) & SATPIXEL)
		      DQSetPix(group2->dq.data,i,j,
			 DQPix(group2->dq.data,i,j) | SATPIXEL);
	     }
	}

}
