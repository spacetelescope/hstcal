# include <stdlib.h>
# include <string.h>
# include <ctype.h>

# include "hstio.h"	/* defines HST I/O functions */
# include "wf3.h"
# include "wf3info.h"
# include "hstcalerr.h"


/* PHOTCALC: Store the photometry parameter values in the global header
** keywords PHOTMODE, PHOTFLAM, PHOTFNU, PHOTZPT, PHOTPLAM, and PHOTBW.
** 
** The input data are NOT modified (only the global header is modified).
**
** Revision history:
** H.Bushouse	Oct. 2000	Initial CALNICA to CALWF3 port.
** H.Bushouse	19-Nov-2001	Modified to track CALACS changes: Call
**				Synphot library functions to compute
**				photometry parameters, instead of reading
**				from the PHOTTAB ref table. Copied from
**				CALWF3 UVIS routine.
** H.Bushouse	 7-Sep-2011	Removed calls to synphot routines and replaced
**				with interface to new imphttab ref table.
*/

static void Phot2Obs (char *, char *);

int photcalc (WF3Info *wf3, MultiNicmosGroup *input) {
    
/* Arguments:
**	wf3	 i: WFC3 info structure
**	input	io: input image
*/

	/* Local variables */
	PhotPar obs;
	float photfnu;
	char  photmode[SZ_LINE+1], obsmode[SZ_LINE+1];
    int status;

	/* Function definitions */
	int GetKeyStr (Hdr *, char *, int, char *, char *, int);
	int PutKeyFlt (Hdr *, char *, float, char *);
	void PrSwitch (char *, int);

	if (wf3->photcorr == PERFORM) {

	/* Extract photmode from sci extension header */
	if ( (status=GetKeyStr (input->group[0].globalhdr, "PHOTMODE", USE_DEFAULT, "",
		       photmode, SZ_LINE)))
	    return (status);

	/* Convert PHOTMODE string into synphot OBSMODE syntax */
	Phot2Obs (photmode, obsmode);
	if (wf3->verbose) {
	    sprintf (MsgText, "Created obsmode of: %s", obsmode);
	    trlmessage (MsgText);
	}

	/* Initialize PhotPar structure */
	InitPhotPar (&obs, wf3->phot.name, wf3->phot.pedigree);

	/* Get the photometry values from the IMPHTTAB table */
	if (GetPhotTab (&obs, obsmode)) {
	    trlerror ("Error return GetPhotTab.");
	}

	if (wf3->verbose) {
	    sprintf (MsgText, "Computed PHOTFLAM value of %g", obs.photflam);
	}

	/* Update the photometry keyword values in the header */
	if (PutKeyFlt (input->group[0].globalhdr, "PHOTFLAM", obs.photflam, ""))
	    return (status);

	if (PutKeyFlt (input->group[0].globalhdr, "PHOTZPT", obs.photzpt, ""))
	    return (status);

	if (PutKeyFlt (input->group[0].globalhdr, "PHOTPLAM", obs.photplam, ""))
	    return (status);

	if (PutKeyFlt (input->group[0].globalhdr, "PHOTBW", obs.photbw, ""))
	    return (status);

	photfnu = 3.33564e+4 * obs.photflam * obs.photplam*obs.photplam;

	if (PutKeyFlt (input->group[0].globalhdr, "PHOTFNU", photfnu, ""))
	    return (status);

	FreePhotPar (&obs);

	PrSwitch ("photcorr", COMPLETE);
	}

	/* Successful return */
	return (status = 0);
}

/* This function converts the PHOTMODE string into an OBSMODE string
   suitable for use with synphot functions.
        PHOTMODE - all upper case component names separated by blanks
        OBSMODE - all lower case names separated by commas
   NOTE: The operation occurs **in-place** on the photmode string.
*/
static void Phot2Obs (char *photmode, char *obsmode) {

        char blank = 32, comma = 44;
        int i, len, c1;

        len = strlen (photmode);

        for (i=0; i < len; i++) {
             c1 = photmode[i];
             if (c1 == blank) {
                 obsmode[i] = comma;
             } else
                 obsmode[i] = tolower(c1);
        }
        obsmode[len] = '\0';
}

