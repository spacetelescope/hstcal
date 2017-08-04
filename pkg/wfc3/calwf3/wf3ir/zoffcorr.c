# include "hstio.h"    /* defines HST I/O functions */
# include "wf3.h"
# include "wf3info.h"
# include "trlbuf.h"

extern int status;

static int zoffcorr (WF3Info *, SingleNicmosGroup *, SingleNicmosGroup *);

/* DOZOFF: Call the ZOFFCORR step for all MULTIACCUM groups.
**
** Revision history:
** H.Bushouse	Oct. 2000	Initial CALNICA to CALWF3 port.
** H.Bushouse	08-May-2002	Modified to use trlkwerr.
*/

int doZoffIR (WF3Info *wf3, MultiNicmosGroup *input, SingleNicmosGroup *zoff) {

/* Arguments:
**	wf3	 i: WFC3 info structure
**	input	io: MULTIACCUM science image
**	zoff	 i: zeroth read image
*/

	/* Local variables */

	/* Function definitions */
	void PrSwitch (char *, int);

	/* Do the MultiAccum zero-read subtraction for each group */
	if (wf3->zoffcorr == PERFORM) {

	    for (wf3->group=wf3->ngroups; wf3->group >= 1; wf3->group--) {
		 if (zoffcorr (wf3, &(input->group[wf3->group-1]), zoff))
		     return (status);
	    }

	    PrSwitch ("zoffcorr", COMPLETE);
	}

	/* Successful return */
	return (status = 0);
}

/* ZOFFCORR: Perform zero-read subtraction for MULTIACCUM data
** groups. The science images are subtracted and the DQ arrays 
** are combined. The ERR, SAMP, arrays are unchanged.
** The exposure time for the group being corrected is reduced
** by an amount equal to the exposure time of the zero-read.
**
** Revision history:
** H.Bushouse	Oct. 2000	Initial CALNICA to CALWF3 port.
*/

static int zoffcorr (WF3Info *wf3, SingleNicmosGroup *input,
		     SingleNicmosGroup *zoff) {

/* Arguments:
**	wf3	 i: WFC3 info structure
**	input	io: image to be zero-subtracted
**	zoff	 i: zero-read image
*/

	/* Local variables */
	int i, j;		/* loop indexes */
	double ztime;		/* zero-read exposure time */

	/* Function definitions */
	void aor (SingleNicmosGroup *, SingleNicmosGroup *);

	/* Subtract the science arrays from one another */
	for (j=0; j < input->sci.data.ny; j++) {
	     for (i=0; i < input->sci.data.nx; i++) {
		  Pix(input->sci.data,i,j) -= Pix(zoff->sci.data,i,j);
	     }
	}

	/* Combine (i.e. logical "OR") the data quality arrays */
	aor (input, zoff);

	/* Subtract the exposure time of the zero-read image
	** from the exposure time of the science image being processed */
	ztime = 0;
	if (getKeyD (&(zoff->sci.hdr), "SAMPTIME", &ztime)) {
	    trlkwerr ("SAMPTIME", wf3->zoff.name);
	    return (status = 1);
	}
	wf3->exptime[wf3->group-1] -= ztime;

	/* Subtract the time arrays from one another */
	for (j=0; j < input->intg.data.ny; j++) {
	     for (i=0; i < input->intg.data.nx; i++) {
		  Pix(input->intg.data,i,j) -= Pix(zoff->intg.data,i,j);
	     }
	}

	/* Successful return */
	return (status = 0);
}

