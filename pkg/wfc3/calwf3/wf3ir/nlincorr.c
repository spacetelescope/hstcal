# include <stdio.h>
# include <math.h>

#include "hstcal.h"
# include "hstio.h"	/* defines HST I/O functions */
# include "wf3.h"
# include "wf3info.h"
# include "wf3dq.h"
# include "trlbuf.h"

extern int status;

static int nlincorr (WF3Info *, SingleNicmosGroup *, NlinData *,
		     SingleNicmosGroup *);

/* DONLIN: Call NLINCORR for all readouts of a MultiAccum.
**
**	   After each MultiAccum group is corrected the routine also
**	   sets saturation flags in the next group for those pixels
**	   that are flagged as saturated in the current group. This
**	   is necessary because the SCI image value of a saturated
**	   pixel will sometimes start to go back down in subsequent
**	   reads after saturation occurs, which means it may fall
**	   below the nominal saturation threshold value and therefore
**	   not get flagged by the normal threshold checking technique.
**
** Revision history:
** H.Bushouse	Oct. 2000	Initial CALNICA to CALWF3 port.
*/

int doNlinIR (WF3Info *wf3, MultiNicmosGroup *input,
	      SingleNicmosGroup *zsig) {

/* Arguments:
**	wf3	 i: WFC3 info structure
**	input	io: input image to be corrected
**	zsig	 i: MULTIACCUM zero-read signal image
*/

	/* Local variables */
	NlinData nlin;

	/* Function definitions */
	int getNlinData (WF3Info *, NlinData *);
	void freeNlinData (NlinData *);
	void satcheck (SingleNicmosGroup *, SingleNicmosGroup *);
	void PrSwitch (char *, int);

	/* Do the non-linearity correction for each group */
	if (wf3->nlincorr == PERFORM) {

	    /* Load the nlin reference file data */
	    if (getNlinData (wf3, &nlin))
		return (status);

	    /* Loop over MultiAccum groups */
	    for (wf3->group=wf3->ngroups; wf3->group >= 1; wf3->group--) {

		 if (nlincorr (wf3, &(input->group[wf3->group-1]), &nlin,
			       zsig))
		     return (status);

		 /* Flag pixels in the next group as saturated if they're
		 ** flagged as saturated in the current group */
		 if (wf3->group-1 > 0)
		     satcheck (&(input->group[wf3->group-1]),
			       &(input->group[wf3->group-2]));
	    }

	    freeNlinData (&nlin);

	    PrSwitch ("nlincorr", COMPLETE);
	}

	/* Successful return */
	return (status = 0);
}

/* NLINCORR: Correct science image values for detector nonlinearity.
** This routine assumes that the correction is represented by a nth-order
** polynomial for each pixel in the detector. N+1 arrays of coefficients
** and (co)variances define the polynomial function.
** It is further assumed that there is one "node" or transition value array
** in the reference data that defines the signal level (in units of DNs)
** above which a given pixel saturates. Science image values below this node
** receive the non-linearity correction; values above the node are deemed 
** saturated, are flagged as such, and receive no correction.
** There is also a DQ flag array in the reference data, which is
** propagated into the science data.
**
** A correction is made for estimated signal in the zeroth read of a
** MultiAccum observation. Without the correction, the zero-subtracted pixel
** values can fall into the wrong correction regime. This can lead to
** pixels that are truly saturated going unflagged because their values
** are below the saturation limits in the NLINFILE, or, for those that fall
** within the linear correction (middle) regime it can lead to the wrong
** amount of correction. The correction for signal in the zeroth read is
** accomplished by adding the estimated zero-read signal to the incoming
** SCI image values before applying the linearity correction, and then
** removing it again afterwards.
**
** The SCI and ERR arrays are updated, and the DQ values are propagated.
** The SAMP and TIME arrays are not modified.
**
** Revision history:
** H.Bushouse	Oct. 2000	Initial CALNICA to CALWF3 port.
** H.Bushouse	21-Mar-2002	Modified to support WFC3 IR subarrays.
** H.Bushouse	10-Apr-2002	Modified to skip WFC3 IR reference pixels
**				by setting loop limits based on trim values
**				from OSCNTAB.
** H.Bushouse	08-May-2002	Added use of wf3dq.h and changed DQ flag
**				macro "SATURATED" to "SATPIXEL".
** H.Bushouse	14 Aug 2003	Modified to use only 1 node array from ref
**				table; assumption is that we will have only
**				a saturation node for WFC3, with no lower
**				limit for linearity correction.
** H.Bushouse	26-Jul-2005	Fixed bug in calculation of nlin ref image
**				pixel indexes.
** H.Bushouse	29-Aug-2008	Modified to use WFC3/IR third-order correction
**				and use ncoeff and nerr from NlinData struct.
**				No update is done to the input ERR values at
**				this time. This will be added in the future.
*/

static int nlincorr (WF3Info *wf3, SingleNicmosGroup *input, NlinData *nlin,
		     SingleNicmosGroup *zsig) {

/* Arguments:
**	wf3	 i: WFC3 info structure
**	input	io: input image to be corrected
**	nlin	 i: nonlinearity reference data
**	zsig	 i: MULTIACCUM zero-read signal image
*/

	/* Local variables */
	int i, j, li, lj, k;	/* pixel indexes */
	int ibeg, iend;		/* loop limits */
	int jbeg, jend;		/* loop limits */
	int li_beg, lj_beg;	/* loop limits */
	int rsize = 1;		/* for use by GetCorner */
	int sci_bin[2];		/* bin size of science image */
	int sci_corner[2];	/* science image corner location */
	int ref_bin[2];		/* bin size of reference image */
	int ref_corner[2];	/* ref image corner location */
	int nsatpix;		/* number of saturated pixels */
	float sval, eval;	/* science and err image values */
	float corr;		/* correction value */
	float n1;		/* node value */

	/* Function definitions */
	int GetCorner (Hdr *, int, int *, int*);

	/* Compute subarray offsets, if any, between ref data
	** science data. */
	if ( (status = GetCorner(&input->sci.hdr, rsize, sci_bin, sci_corner)))
	    return (status);
	if ( (status = GetCorner(&nlin->coeff[0].hdr, rsize, ref_bin, ref_corner)))
	    return (status);

	/* Initialize saturated pixel counter */
	nsatpix = 0;

	/* Loop through science image */
	ibeg = wf3->trimx[0]; iend = input->sci.data.nx - wf3->trimx[1];
	jbeg = wf3->trimy[0]; jend = input->sci.data.ny - wf3->trimy[1];
	li_beg = (sci_corner[0] - ref_corner[0]) + ibeg;
	lj_beg = (sci_corner[1] - ref_corner[1]) + jbeg;

	for (j = jbeg, lj = lj_beg; j < jend; j++, lj++) {
	for (i = ibeg, li = li_beg; i < iend; i++, li++) {

	     /* Get the science and error image values */
	     sval = Pix(input->sci.data,i,j);
	     eval = Pix(input->err.data,i,j);

	     /* Temporarily add the MULTIACCUM zero-read signal back into the
	     ** the pixel value, but only if ZSIG step is turned on and only 
	     ** for groups other than the zeroth-read itself */
	     if (wf3->zsigcorr == PERFORM && wf3->group != wf3->ngroups) {
		 sval += Pix(zsig->sci.data,i,j);
		 if (DQPix(zsig->dq.data,i,j) & ZEROSIG) {
		     DQSetPix(input->dq.data,i,j,
			DQPix(input->dq.data,i,j) | ZEROSIG);
		 }
	     }

	     /* Get the node values for this pixel */
	     n1 = Pix(nlin->nodes[0].data,li,lj);

	     /* Propagate the DQ value from the NLIN ref data */
	     DQSetPix(input->dq.data,i,j,
		DQPix(input->dq.data,i,j) | DQPix(nlin->dqual[0].data,li,lj));

	     /* If it's already flagged as saturated,
	     ** skip the correction */
	     if (DQPix(input->dq.data,i,j) & SATPIXEL) {
		 nsatpix++;

	     /* Apply the correction for the non-linear region */
	     /*} else if (sval >= n1 && sval <= n2) {*/
	     } else if (sval <= n1) {

	       /* Compute the new science image pixel value */
	       corr = 1.0;
	       for (k=0; k < nlin->ncoeff; k++)
		    corr += Pix(nlin->coeff[k].data,li,lj) * (pow(sval,k));
	       Pix(input->sci.data,i,j) = sval * corr;

	       /* Remove the MULTIACCUM zero-read signal that was added in
	       ** above, but only if ZSIG step is turned on and only for 
	       ** groups other than the zeroth-read itself */
	       if (wf3->zsigcorr == PERFORM && wf3->group != wf3->ngroups)
		   Pix(input->sci.data,i,j) -= Pix(zsig->sci.data,i,j);

	     /* Above the saturation node, just mark the pixel as saturated */
	     } else if (sval > n1) {
	       nsatpix++;
	       DQSetPix(input->dq.data,i,j,
		  DQPix(input->dq.data,i,j) | SATPIXEL);
	     }
	}}

	/* Report the number of saturated pixels */
	sprintf (MsgText, "NLINCORR detected %d saturated pixels in imset %d",
		 nsatpix, wf3->group);
	trlmessage (MsgText);

	/* Successful return */
	return (status = 0);
}

