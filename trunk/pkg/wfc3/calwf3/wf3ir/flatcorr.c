# include <string.h>	/* for strncmp */

# include "hstio.h"	/* defines HST I/O functions */
# include "wf3.h"
# include "wf3info.h"

extern int status;

static int flatcorr (WF3Info *, SingleNicmosGroup *, SingleNicmosGroup *);

/* DOFLAT: Call FLATCORR for each readout of a MultiAccum.
**
** Revision history:
** H.Bushouse	Oct. 2000	Initial CALNICA to CALWF3 port.
** H.Bushouse	28-Mar-2002	Upgraded to handle 3 types of flats and
**				subarrays - added use of getFlats and
**				RebinRef routines.
*/

int doFlatIR (WF3Info *wf3, MultiNicmosGroup *input,
	      SingleNicmosGroup *crimage) {

/* Arguments:
**	wf3	 i: WFC3 info structure
**	input	io: image to be flat fielded
**	crimage io: cr-rej image to be flat fielded
*/

	/* Local variables */
	int avg = 1;
	SingleNicmosGroup flat;

	/* Function definitions */
	int getFlats (WF3Info *, SingleNicmosGroup *, SingleNicmosGroup *);
	int RebinRef (SingleNicmosGroup *, SingleNicmosGroup *, int);
	void PrSwitch (char *, int);

	/* Apply the flatfield correction */
	if (wf3->flatcorr == PERFORM) {

	    /* Load the flat field ref image(s) */
	    initSingleNicmosGroup (&flat);
	    if (getFlats (wf3, &(input->group[0]), &flat))
		return (status);

	    /* Rebin or extract subarray from flat, if necessary */
	    if (RebinRef (&(input->group[0]), &flat, avg))
		return (status);

	    /* Apply the flat to each MultiAccum group */
	    for (wf3->group=wf3->ngroups; wf3->group >= 1; wf3->group--) {
		 if (flatcorr (wf3, &(input->group[wf3->group-1]), &flat))
		     return (status);
	    }

	    /* Also apply the flat to the cr-combined image */
	    if (wf3->crcorr == PERFORM) {
		if (flatcorr (wf3, crimage, &flat))
		    return (status);
	    }

	    PrSwitch ("flatcorr", COMPLETE);

	    freeSingleNicmosGroup (&flat);
	}

	/* Successful return */
	return (status = 0);
}

/* FLATCORR: Flat field an image. The science image
** is divided in-place by the flat field image.
** Errors and DQ flags from the flat field are combined with
** the science data errors and flags. The input SAMP and TIME
** arrays are unchanged.
**
** Revision history:
** H.Bushouse	Oct. 2000	Initial CALNICA to CALWF3 port.
** H.Bushouse	10-Apr-2002	Upgraded to skip reference pixels by using
**				adiv_noref routine.
** H.Bushouse	 9-Oct-2008	Upgraded to convert data to units of electrons
**				by multiplying by the gain after flat field
**				has been applied. Uses new function "mult_gain".
** H.Bushouse   10-Dec-2008	Fixed mult_gain to apply offsets to ampx/ampy
**				values for proper operation with subarrays.
** H.Bushouse	13-Aug-2009	Set BUNIT to electrons after gain is applied.
** H.Bushouse	21-Oct-2009	Updated mult_gain to apply the mean gain
**				from all amps to the entire image, except in
**				the case of grism images.
** H.Bushouse	24-Nov-2009	Updated mult_gain to use mean gain for all
**				images, including grism. (calwf3 v2.0)
*/

int flatcorr (WF3Info *wf3, SingleNicmosGroup *input,
	      SingleNicmosGroup *flat) {

/* Arguments:
**	wf3	 i: WFC3 info structure
**	input	io: image to be flat fielded
**	flat	 i: flat field image
*/

	/* Function definitions */
	void adiv_noref (WF3Info *, SingleNicmosGroup *, SingleNicmosGroup *);
	void mult_gain  (WF3Info *, SingleNicmosGroup *);
	int  PutKeyStr (Hdr *, char *, char *, char *);

	/* Do the flat fielding in-place in input image data */
	adiv_noref (wf3, input, flat);

	/* Now convert the data from DN to electrons by multiplying by gain */
	mult_gain (wf3, input);

	/* Update the BUNIT keyword to reflect gain correction */
	if (wf3->unitcorr == PERFORM || wf3->unitcorr == COMPLETE) {
	    if (PutKeyStr (&input->sci.hdr, "BUNIT", "ELECTRONS/S", ""))
		return (status);
	    if (PutKeyStr (&input->err.hdr, "BUNIT", "ELECTRONS/S", ""))
		return (status);
	} else {
	    if (PutKeyStr (&input->sci.hdr, "BUNIT", "ELECTRONS", ""))
		return (status);
	    if (PutKeyStr (&input->err.hdr, "BUNIT", "ELECTRONS", ""))
		return (status);
	}

	/* Successful return */
	return (status = 0);
}

void mult_gain (WF3Info *wf3, SingleNicmosGroup *in) {

/* Multiply a SingleNicmosGroup by the gain value.

   (*in) *= gain

   The science and error arrays are multiplied by the gain. The data
   quality array is unchanged.
*/

	/* Local variables */
        int i, j;       /* array indexes */

	/* Apply the mean gain value to the entire image */
	for (j=0; j<in->sci.data.ny; j++) {
	     for (i=0; i<in->sci.data.nx; i++) {
		  Pix(in->sci.data,i,j) *= wf3->mean_gain;
		  Pix(in->err.data,i,j) *= wf3->mean_gain;
	     }
	}

}

