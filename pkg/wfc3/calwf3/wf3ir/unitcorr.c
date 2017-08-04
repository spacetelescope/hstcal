# include <stdio.h>

#include "hstcal.h"
# include "hstio.h"	/* defines HST I/O functions */
# include "wf3.h"
# include "wf3info.h"
# include "trlbuf.h"

extern int status;

/* DOUNIT: Call UNITCORR for all readouts of a MultiAccum.
**
** Revision history:
** H.Bushouse	Oct. 2000	Initial CALNICA to CALWF3 port.
*/

int doUnitIR (WF3Info *wf3, MultiNicmosGroup *input) {

/* Arguments:
**	wf3	 i: WFC3 info structure
**	input	io: input image
*/

	/* Local variables */

	/* Function definitions */
	int unitcorr (WF3Info *, SingleNicmosGroup *);
	void PrSwitch (char *, int);

	/* Do the units correction for each group */
	if (wf3->unitcorr == PERFORM) {

	    for (wf3->group=wf3->ngroups; wf3->group >= 1; wf3->group--) {
		 if (unitcorr (wf3, &(input->group[wf3->group-1])))
		     return (status);
	    }

	    PrSwitch ("unitcorr", COMPLETE);
	}

	/* Successful return */
	return (status = 0);
}

/* UNITCORR: Convert data from units of counts to countrates.
** The input SCI and ERR arrays are divided by the exposure time.
** The DQ, SAMP, and TIME arrays are unchanged.
** The BUNIT keyword in the SCI and ERR image headers are updated.
**
** Revision history:
** H.Bushouse	Oct. 2000	Initial CALNICA to CALWF3 port.
** H.Bushouse	10-Apr-2002	Modified to skip IR reference pixels by using
**				amulk_noref routine and basing loop limits on
**				trim values from OSCNTAB.
** H.Bushouse	14-Jan-2010	No longer need to check ZSIGCORR before using
**				sampzero, because it's computed in doIR now.
**				(calwf3 v2.0)
** H.Bushouse	21-Oct-2010	Upgraded to check flatcorr status to decide
**				proper units for BUNIT keyword update, to
**				support re-entrant processing. (PR 66081)
*/

int unitcorr (WF3Info *wf3, SingleNicmosGroup *input) {

/* Arguments:
**	wf3	 i: WFC3 info structure
**	input	io: input image
*/

	/* Local variables */
	int i, j;		/* pixel indexes */
	int ibeg, iend;		/* loop limits */
	int jbeg, jend;		/* loop limits */
	float time;		/* exposure time */

	/* Function definitions */
	void amulk_noref (WF3Info *, SingleNicmosGroup *, float);
	int  PutKeyStr (Hdr *, char *, char *, char *);

	/* Skip conversion if units are already countrate */
	if (wf3->bunit[wf3->group-1] == COUNTRATE) {
	    sprintf (MsgText,
	       "Data already in units of countrates; UNITCORR will be skipped");
	    trlwarn (MsgText);
	    wf3->unitcorr = SKIP;
	    return (status = 0);
	}

	/* If we're processing a MultiAccum zeroth read, use the value of
	** wf3->sampzero for the exposure time */
	if (wf3->group == wf3->ngroups) {
	    time = 1.0 / wf3->sampzero;
	    amulk_noref (wf3, input, time);

	/* Otherwise, divide the input SCI and ERR arrays by the TIME array */
	} else {

	    ibeg = wf3->trimx[0]; iend = input->sci.data.nx - wf3->trimx[1];
	    jbeg = wf3->trimy[0]; jend = input->sci.data.ny - wf3->trimy[1];

	    for (j = jbeg; j < jend; j++) {
		 for (i = ibeg; i < iend; i++) {
		      time = Pix(input->intg.data,i,j);
		      if (time != 0) {
			  Pix(input->sci.data,i,j) /= time;
			  Pix(input->err.data,i,j) /= time;
		      } else {
			  Pix(input->sci.data,i,j) = 0.0;
			  Pix(input->err.data,i,j) = 0.0;
		      }
		 }
	    }
	}

	/* Update the units keyword in the SCI and ERR headers */
	if (wf3->flatcorr == COMPLETE) {
	    if (PutKeyStr (&input->sci.hdr, "BUNIT", "ELECTRONS/S", ""))
		return (status);
	    if (PutKeyStr (&input->err.hdr, "BUNIT", "ELECTRONS/S", ""))
		return (status);
	} else {
	    if (PutKeyStr (&input->sci.hdr, "BUNIT", "COUNTS/S", ""))
		return (status);
	    if (PutKeyStr (&input->err.hdr, "BUNIT", "COUNTS/S", ""))
		return (status);
	}
	wf3->bunit[wf3->group-1] = COUNTRATE;

	/* Successful return */
	return (status = 0);
}

