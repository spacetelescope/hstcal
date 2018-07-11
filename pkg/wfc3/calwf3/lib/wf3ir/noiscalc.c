# include <math.h>
# include <stdio.h>
# include <string.h>

# include "hstio.h"	/* defines HST I/O functions */
# include "wf3.h"
# include "wf3info.h"
# include "trlbuf.h"

extern int status;

/* DoNoisIR: Call NOISCALC for each readout of a MultiAccum.
**
** Revision history:
** H.Bushouse	Oct. 2000	Initial CALNICA to CALWF3 port.
** H.Bushouse	Mar. 2002	Modified to eliminated use of NICMOS-style
**				mask reference image.
** H.Bushouse	21-Oct-2010	Modified doNoisIR to print trailer message
**				and noiscorr switch values and also indicate
**				if noiscalc is skipped, to support re-entrant
**				processing. (PR 66081)
*/

int doNoisIR (WF3Info *wf3, MultiNicmosGroup *input) {

/* Arguments:
**	wf3	 i: WFC3 info structure
**	input	io: input image
*/

	/* Function definitions */
	int noiscalc (WF3Info *, SingleNicmosGroup *);
	int OmitStep (int);
	void PrSwitch (char *, int);

	/* Do the noise calculation for each group */
	if (wf3->noiscorr == PERFORM) {

	    for (wf3->group=wf3->ngroups; wf3->group >= 1; wf3->group--) {
		 if (noiscalc (wf3, &(input->group[wf3->group-1]))) {
		     PrSwitch ("noiscorr", SKIPPED);
		     return (status=0);
		 }
	    }
	}

	/* Print status to trailer */
	if (!OmitStep (wf3->noiscorr))
	    trlmessage ("Uncertainty array initialized.");
	PrSwitch ("noiscorr", COMPLETE);

	/* Successful return */
	return (status = 0);
}

/* NOISCALC: Perform noise calculation (populate ERR image).
** This is done by combining the detector readnoise and the Poisson noise
** in the science image.
**
** The input ERR image is modified in-place.
** The input SCI, DQ, SAMP, and TIME arrays are unchanged.
**
** NOTE: Routine assumes that the readnoise values are in units of
** ELECTRONS (not DNs). The noise calculation is performed in units of
** electrons and then converted back to DNs.
**
** Revision history:
** H.Bushouse	10-Apr-2002	Modified to skip WFC3 reference pixels by
**				setting	loop limits based on OSCNTAB trim
**				values.	This automatically also handles
**				subarray images because the loop limits are
**				based on the size of the input science image.
**
** H.Bushouse	27-May-2005	Fixed bug in noise computation by adding
**				variable "noise2" to hold readnoise value.
**
** H.Bushouse	27-Aug-2008	Modified to apply individual gain and readnoise
**				values to each amp region of the images.
**
** H.Bushouse	20-Oct-2010	Updated to check for uninitialized ERR array
**				before computing new values, to support
**				re-entrant processing. (PR 66081)
*/

int noiscalc (WF3Info *wf3, SingleNicmosGroup *input) {

/* Arguments:
**	wf3	 i: WFC3 info structure
**	input	io: input image
*/

	/* Local variables */
	int i, j;		/* loop indexes */
	int ibeg, iend;		/* loop limits */
	int jbeg, jend;		/* loop limits */
	int ampx, ampy;		/* AMP readout boundaries */
	float gain;		/* gain value */
	float rn2;		/* read noise squared */
	float noise;		/* noise value */
	float signal;		/* science image value */

	/* Loop over all pixels in science image, computing err
	** estimate pixel-by-pixel. Note that loop limits exclude
	** the IR detector reference pixels. */

	ibeg = wf3->trimx[0]; iend = input->err.data.nx - wf3->trimx[1];
	jbeg = wf3->trimy[0]; jend = input->err.data.ny - wf3->trimy[1];

	/* First check to see if the ERR array has been populated before.
	** If it has, then just return without doing anything. */
	for (j = jbeg; j < jend; j++) {
	     for (i = ibeg; i < iend; i++) {
		  if (Pix (input->err.data,i,j) != 0.) {
		      wf3->noiscorr = SKIPPED;
		      return (status=1);  /* not an empty array */
		  }
	     }
	}

	/* Correct AMP readout boundaries for subarray offsets */
	ampx = wf3->ampx + wf3->offsetx - wf3->trimx[0];
	ampy = wf3->ampy + wf3->offsety - wf3->trimy[0];

	/* Quad 2 = Amp B */
	gain = wf3->atodgain[1];
	rn2  = wf3->readnoise[1] * wf3->readnoise[1];
	for (j = jbeg; j < ampy; j++) {
	     for (i = ibeg; i < ampx; i++) {

		  /* Combine (in quadrature) the detector readnoise and the
		  ** photon noise for each pixel. */
		  signal = Pix(input->sci.data,i,j);  /* photon noise (in DN) */
		  noise = sqrt (rn2 + fabs(signal*gain));
		  Pix(input->err.data,i,j) = noise / gain;
	     }
	}

	/* Quad 3 = Amp C */
	gain = wf3->atodgain[2];
	rn2  = wf3->readnoise[2] * wf3->readnoise[2];
	for (j = jbeg; j < ampy; j++) {
	     for (i = ampx; i < iend; i++) {

		  /* Combine (in quadrature) the detector readnoise and the
		  ** photon noise for each pixel. */
		  signal = Pix(input->sci.data,i,j);  /* photon noise (in DN) */
		  noise = sqrt (rn2 + fabs(signal*gain));
		  Pix(input->err.data,i,j) = noise / gain;
	     }
	}

	/* Quad 1 = Amp A */
	gain = wf3->atodgain[0];
	rn2  = wf3->readnoise[0] * wf3->readnoise[0];
	for (j = ampy; j < jend; j++) {
	     for (i = ibeg; i < ampx; i++) {

		  /* Combine (in quadrature) the detector readnoise and the
		  ** photon noise for each pixel. */
		  signal = Pix(input->sci.data,i,j);  /* photon noise (in DN) */
		  noise = sqrt (rn2 + fabs(signal*gain));
		  Pix(input->err.data,i,j) = noise / gain;
	     }
	}

	/* Quad 4 = Amp D */
	gain = wf3->atodgain[3];
	rn2  = wf3->readnoise[3] * wf3->readnoise[3];
	for (j = ampy; j < jend; j++) {
	     for (i = ampx; i < iend; i++) {

		  /* Combine (in quadrature) the detector readnoise and the
		  ** photon noise for each pixel. */
		  signal = Pix(input->sci.data,i,j);  /* photon noise (in DN) */
		  noise = sqrt (rn2 + fabs(signal*gain));
		  Pix(input->err.data,i,j) = noise / gain;
	     }
	}

	/* Successful return */
	return (status = 0);
}

