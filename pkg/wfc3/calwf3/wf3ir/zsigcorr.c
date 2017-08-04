# include <math.h>
# include <stdio.h>

#include "hstcal.h"
# include "hstio.h"    /* defines HST I/O functions */
# include "wf3.h"
# include "wf3info.h"
# include "wf3dq.h"
# include "trlbuf.h"

extern int status;

static int zsigcorr (WF3Info *, MultiNicmosGroup *, NlinData *, 
		     SingleNicmosGroup *);


/* DOZSIG: Calls ZSIGCORR calibration step for a MultiAccum science image.
**
** Revision history:
** H.Bushouse	Oct. 2000	Initial CALNICA to CALWF3 port.
** H.Bushouse	24-May-2001	Eliminated use of mask image for WFC3.
** H.Bushouse	14-Aug-2003	Use the first (and only) nlin node array as
				the saturation value at each pixel.
*/

int doZsigIR (WF3Info *wf3, MultiNicmosGroup *input, SingleNicmosGroup *zsig) {

/* Arguments:
**	wf3	i: WFC3 info structure
**	input	i: input MultiAccum image
**	zsig	o: zero-read signal image
*/

	/* Local variables */
	NlinData nlin;

	/* Function definitions */
	int getNlinData  (WF3Info *, NlinData *);
	void freeNlinData (NlinData *);
	void PrSwitch (char *, int);
	
	if (wf3->zsigcorr == PERFORM) {

	    /* Get nonlinearity reference data;
	    ** we need this for the saturation values */
	    if (getNlinData (wf3, &nlin))
		return (status);

	    /* Compute zero-read signal correction image */
            if (zsigcorr (wf3, input, &nlin, zsig))
        	return (status);
	
	    /* Free the reference data */
	    freeNlinData (&nlin);

	    PrSwitch ("zsigcorr", COMPLETE);
	}
 
	/* Successful return */
	return (status = 0);

}

# define ZTHRESH 4.0	/* default signal clipping threshold (sigma) */

/* ZSIGCORR: Calculate the amount of signal from sources in a MultiAccum
** zeroth-read image. This is done by computing the amount of signal that
** arrived between the zeroth and first reads, and then scaling that to the
** exposure time of the zeroth read in order to estimate how much signal
** probably came in between the reset and the zeroth read. Note that this
** does not work well if the source flux is so high that it's beginning to
** go non-linear already in either the zeroth or first read. The raw data
** values in the zeroth and first readouts are also differenced with the
** super zero read reference image and checked for saturation in those
** readouts.
**
** Revision History:
** H.Bushouse	Oct. 2000	Initial CALNICA to CALWF3 port.
** H.Bushouse	24-May-2001	Eliminated use of mask image for WFC3.
** H.Bushouse	20-Mar-2002	Upgraded to handle subarrays.
** H.Bushouse	11-Apr-2002	Upgraded to skip reference pixels in all
**				computations.
** H.Bushouse	08-May-2002	Added use of wf3dq.h and changed DQ flag macro
**				"SATURATED" to "SATPIXEL".
** H.Bushouse	26-Jul-2005	Fixed bug in calculation of nlin ref image
**				pixel indexes.
** H.Bushouse	01-Dec-2008	Fixed the calculation of the zeroth read
**				exposure time to handle WFC3/IR subarray
**				readouts.
** H.Bushouse	13-Jan-2010	Modified to: set ZEROSIG DQ values in ZSIG
** image along with SATPIXEL values; set and count pixels as saturated in the
** first read if they're saturated in the zeroth read; only check for saturation
** in first read if not	already saturated in zeroth read; for pixels saturated
** in the zeroth or first read, recompute zsig from difference of zeroth read
** and super-zero zsci. Moved computation of sampzero into doIR. (calwf3 v2.0)
**
** H.Bushouse	30-Sep-2010	Modified to no longer call pixOK function
** before operating on a pixel. Always do the calculation regardless of DQ
** values. (PR #66080, Trac #607)
**
** H.Bushouse	11-Mar-2011	Modified whole algorithm to form the zsig image
** by always just taking the difference of the science zeroth-read with the
** super-zero read in the nlinfile, instead of doing the backwards
** extrapolation of the first-zero difference. The new method does not use the
** effective exposure time of the zeroth-read (sampzero), which avoids the
** problems created by the variable zeroth-read exposure time in IR subarrays.
** Also no longer need the dark ref file. (PR #67728, Trac #681)
*/

static int zsigcorr (WF3Info *wf3, MultiNicmosGroup *input, NlinData *nlin,
		     SingleNicmosGroup *zsig) {

/* Arguments:
**	wf3	i: WFC3 info structure
**	input	i: input images
**	nlin	i: non-linearity correction data
**	zsig	o: zero-read signal image
*/

	/* Local variables */
	int i, j, li, lj;		/* pixel indexes */
	int ibeg, iend, jbeg, jend;	/* loop limits */
	int li_beg, lj_beg;		/* loop limits */
	int rsize = 1;			/* for use by GetCorner */
	int sci_bin[2];			/* bin size of science image */
	int sci_corner[2];		/* science image corner location */
	int ref_bin[2];			/* bin size of reference image */
	int ref_corner[2];		/* ref image corner location */
	int nsat0, nsat1;		/* number of saturated pixels */
	float noise;			/* total noise in zsig+zerr */

	/* Function definitions */
	int  copyGroup (SingleNicmosGroup *, SingleNicmosGroup *);
	void asub (SingleNicmosGroup *, SingleNicmosGroup *);
	void asub_noref (WF3Info *, SingleNicmosGroup *, SingleNicmosGroup *);
	int  noiscalc (WF3Info *, SingleNicmosGroup *);
	int  GetCorner (Hdr *, int, int *, int *);

	/* Initialize counters */
	nsat0 = 0; nsat1 = 0;

	/* Set proper threshold value; use default if no user input */
	if (wf3->zsthresh == 0)
	    wf3->zsthresh = ZTHRESH;

	/* Copy zeroth read to zsig image */
	if (copyGroup (zsig, &(input->group[wf3->ngroups-1])))
	    return (status);

	/* Compute subarray offsets, if any, between ref data
	** and science data. */
	if ((status = GetCorner(&zsig->sci.hdr, rsize, sci_bin, sci_corner)))
	    return (status);
	if ( (status = GetCorner(&nlin->coeff[0].hdr, rsize, ref_bin, ref_corner)))
	    return (status);

	/* Set loop limits; this automatically handles subarrays because
	** it is based on the size of the input science image. */
	ibeg = wf3->trimx[0]; iend = zsig->sci.data.nx - wf3->trimx[1];
	jbeg = wf3->trimy[0]; jend = zsig->sci.data.ny - wf3->trimy[1];
	li_beg = (sci_corner[0] - ref_corner[0]) + ibeg;
	lj_beg = (sci_corner[1] - ref_corner[1]) + jbeg;

	/* Subtract super-zero read ref image from science zero read */
	for (j = jbeg, lj = lj_beg; j < jend; j++, lj++) {
	for (i = ibeg, li = li_beg; i < iend; i++, li++) {

	     Pix(zsig->sci.data,i,j) -= Pix(nlin->zsci[0].data,li,lj);
	}}

	/* Compute noise in the zsig image */
	if (noiscalc (wf3, zsig))
	    return (status);

	/* Loop over the zsig image, skipping reference pixels */
	for (j = jbeg, lj = lj_beg; j < jend; j++, lj++) {
	for (i = ibeg, li = li_beg; i < iend; i++, li++) {

	     /* Compute total noise in pixel, which includes noise in the
	     ** science zero read and noise in the super-zero read. */
	     noise = sqrt (Pix(zsig->err.data,i,j)*Pix(zsig->err.data,i,j) +
		   Pix(nlin->zerr[0].data,li,lj)*Pix(nlin->zerr[0].data,li,lj));

	     /* Flag pixels with greater than zsthresh*err signal as
	     ** having signal in the zeroth read */
	     if (Pix(zsig->sci.data,i,j) >= wf3->zsthresh*noise) {

		 DQSetPix(zsig->dq.data,i,j, DQPix(zsig->dq.data,i,j) |
			  ZEROSIG);
		 DQSetPix(input->group[wf3->ngroups-1].dq.data,i,j,
			  DQPix(input->group[wf3->ngroups-1].dq.data,i,j) |
			  ZEROSIG);

	     /* Mask out low signal pixels by setting them to zero */
	     } else
		 Pix(zsig->sci.data,i,j) = 0.0;

	     /* Flag pixels that are saturated in the zeroth read */    
	     if (Pix(zsig->sci.data,i,j) > Pix(nlin->nodes[0].data,li,lj)) {

		 DQSetPix(zsig->dq.data,i,j, DQPix(zsig->dq.data,i,j) |
			  (SATPIXEL+ZEROSIG));
		 DQSetPix(input->group[wf3->ngroups-1].dq.data,i,j,
			  DQPix(input->group[wf3->ngroups-1].dq.data,i,j) |
			  (SATPIXEL+ZEROSIG));
		 DQSetPix(input->group[wf3->ngroups-2].dq.data,i,j,
			  DQPix(input->group[wf3->ngroups-2].dq.data,i,j) |
			  SATPIXEL);
		 nsat0++;
		 nsat1++;

	     /* Flag pixels that are saturated in the first read */
	     } else if (Pix(input->group[wf3->ngroups-2].sci.data,i,j) -
			Pix(nlin->zsci[0].data,li,lj) >
			Pix(nlin->nodes[0].data,li,lj)) {

		 DQSetPix(zsig->dq.data,i,j, DQPix(zsig->dq.data,i,j) |
			  (SATPIXEL+ZEROSIG));
		 DQSetPix(input->group[wf3->ngroups-2].dq.data,i,j,
			  DQPix(input->group[wf3->ngroups-2].dq.data,i,j) |
			  SATPIXEL);
		 nsat1++;
	     }
	}}

	/* Report the number of saturated pixels detected */
	sprintf (MsgText,
		 "ZSIGCORR detected %d saturated pixels in 0th read", nsat0);
	trlmessage (MsgText);
	sprintf (MsgText,
		 "ZSIGCORR detected %d saturated pixels in 1st read", nsat1);
	trlmessage (MsgText);

	/* Add zsig image values to zeroth-read image of input data */
	for (j = jbeg; j < jend; j++) {
	     for (i = ibeg; i < iend; i++) {
		  Pix(input->group[wf3->ngroups-1].sci.data,i,j) +=
			Pix(zsig->sci.data,i,j);
	     }
	}
	
	/* Successful return */
	return (status = 0);

}

