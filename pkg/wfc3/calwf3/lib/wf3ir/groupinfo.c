# include <stdio.h>
# include <string.h>

#include "hstcal.h"
# include "hstio.h"/* defines HST I/O functions */
# include "wf3.h"
# include "wf3info.h"
# include "trlbuf.h"

extern int status;

/* GETGROUPINFO: Get group-specific information, such as exposure time,
** data units, and TDF status.
**
** Revision history:
** H.Bushouse	Oct. 2000	Initial CALNICA to CALWF3 port.
** H.Bushouse	08-May-2002	Modified to use trlkwerr.
** H.Bushouse	21-Oct-2010	Upgraded getDataUnits to recognize BUNIT
**				values of ELECTRONS, to support re-entrant
**				processing. (PR 66081)
*/

int getGroupInfo (WF3Info *wf3, SingleNicmosGroup *in) {

/* Arguments:
**	wf3	io: WFC3 info structure
**	in	 i: input image data
*/

	/* Local variables */

	/* Function definitions */
	int getExpTime   (WF3Info *, Hdr *);
	int getDataUnits (WF3Info *, Hdr *);
	/*int getTDFTrans  (WF3Info *, Hdr *);*/

	/* Get the exposure time */
	if (getExpTime (wf3, &(in->sci.hdr)))
	    return (status);

	/* Get the data units */
	if (getDataUnits (wf3, &(in->sci.hdr)))
	    return (status);

	/* Get the TDF status --see #980
	if (getTDFTrans (wf3, &(in->sci.hdr)))
	    return (status);
    */
    
	/* Successful return */
	return (status = 0);
}

/* GETEXPTIME: Read the exposure time keyword from group header */

int getExpTime (WF3Info *wf3, Hdr *header) {

/* Arguments:
**	wf3	io: WFC3 info structure
**	header	 i: image header
*/

	/* Initialize the exposure time */
	wf3->exptime[wf3->group-1] = 0;

	/* Read the exposure time keyword from the header */
	if (getKeyD (header, "SAMPTIME", &(wf3->exptime[wf3->group-1]))) {
	    trlkwerr ("SAMPTIME", wf3->input);
	    return (status = 1);
	}

        /* Check the exposure time value for validity */
	if (wf3->exptime[wf3->group-1] < 0) {
	    sprintf (MsgText, "SAMPTIME keyword value \"%g\" not valid in %s",
			      wf3->exptime[wf3->group-1], wf3->input);
	    trlerror (MsgText);
	    return (status = 1);
	} else if (wf3->exptime[wf3->group-1] == 0) {
	    if (wf3->group != wf3->ngroups) {
		sprintf (MsgText, "SAMPTIME equal to zero in %s",
			 wf3->input);
		trlwarn (MsgText);
	    }
	}

	/* Successful return */
	return (status = 0);
}

/* GETDATUNITS: Read the BUNIT keyword from group header */

int getDataUnits (WF3Info *wf3, Hdr *header) {

/* Arguments:
**	wf3	io: WFC3 info structure
**	header	 i: image header
*/

	/* Local variables */
	char units[12];				/* BUNIT keyword value */

	/* Read the BUNIT keyword */
	units[0] = '\0';
	if (getKeyS (header, "BUNIT", units)) {
	    trlkwerr ("BUNIT", wf3->input);
	    return (status = 1);
	}

	/* Check the BUNIT keyword value for validity */
	if (strncmp (units, "COUNTS/S", 8) == 0 ||
	    strncmp (units, "ELECTRONS/S", 11) == 0)
	    wf3->bunit[wf3->group-1] = COUNTRATE;
	else if (strncmp (units, "COUNTS", 6) == 0 ||
		 strncmp (units, "ELECTRONS", 9) == 0)
	    wf3->bunit[wf3->group-1] = COUNTS;
	else {
	    sprintf (MsgText, "Unrecognized BUNIT value \"%s\" in %s\n", units,
		     wf3->input);
	    trlerror (MsgText);
	    return (status = 1);
	}

	/* Successful return */
	return (status = 0);
}

/* GETTDFTRANS: Read the TDFTRANS keyword from group header */

int getTDFTrans (WF3Info *wf3, Hdr *header) {

/* Arguments:
**	wf3	io: WFC3 info structure
**	header	 i: image header
*/

	/* Read the TDFTRANS keyword */
	wf3->tdftrans[wf3->group-1] = -1;
	if (getKeyI (header, "TDFTRANS", &(wf3->tdftrans[wf3->group-1]))) {
	    trlkwerr ("TDFTRANS", wf3->input);
	    return (status = 1);
	}

	/* Check the TDFTRANS keyword value for validity */
	if (wf3->tdftrans[wf3->group-1] < 0) {
	    sprintf (MsgText, "TDFTRANS value \"%d\" in %s is illegal\n",
		     wf3->tdftrans[wf3->group-1], wf3->input);
	    trlerror (MsgText);
	    return (status = 1);
	}

	/* Successful return */
	return (status = 0);
}

