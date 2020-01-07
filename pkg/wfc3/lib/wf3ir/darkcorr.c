# include <stdio.h>
# include <float.h>

#include "hstcal.h"
# include "hstio.h"	/* defines HST I/O functions */
# include "wf3.h"
# include "wf3info.h"
# include "trlbuf.h"

extern int status;

static int darkcorr (WF3Info *, SingleNicmosGroup *, SingleNicmosGroup *);

/* DoDarkIR: Call DARKCORR routine for each readout of a MultiAccum.
**	     The appropriate dark image is loaded for each readout
**	     from the DARKFILE reference file.
**
** Revision history:
** H.Bushouse	Oct. 2000	Initial CALNICA to CALWF3 port.
*/

int doDarkIR (WF3Info *wf3, MultiNicmosGroup *input) {

/* Arguments:
**	wf3	 i: WFC3 info structure
**	input	io: image to be dark subtracted
*/

	/* Local variables */
	SingleNicmosGroup dark;

	/* Function definitions */
	int getDarkImage (WF3Info *, SingleNicmosGroup *, int);
	void PrSwitch (char *, int);

	/* Do the dark current subtraction for each group */
	if (wf3->darkcorr == PERFORM) {

	    for (wf3->group=wf3->ngroups; wf3->group >= 1; wf3->group--) {

		 /* Load the appropriate dark image for this group */
	    	 initSingleNicmosGroup (&dark);
		 if (getDarkImage (wf3, &dark, wf3->group))
		     return (status);

		 /* Do dark subtraction for this group */
		 if (darkcorr (wf3, &(input->group[wf3->group-1]), &dark))
		     return (status);

		 /* Free the dark image */
	    	 freeSingleNicmosGroup (&dark);
	    }

	    /* Print status to trailer */
	    PrSwitch ("darkcorr", COMPLETE);
	}

	/* Successful return */
	return (status = 0);
}

/* DARKCORR: Subtract dark current from a single science image.
** The dark current reference image is subtracted in-place from the
** science image. The dark image errors and DQ flags are combined
** with the science data errors and DQ flags. The input SAMP and TIME
** arrays are unchanged.
**
** Revision history:
** H.Bushouse	Oct. 2000	Initial CALNICA to CALWF3 port.
** H.Bushouse	20-Mar-2002	Added use of RebinRef to extract subarray
**				from dark ref image, if necessary.
** H.Bushouse	10-Apr-2002	Upgraded to eliminate reference pixels from
**				the calculation by using asub_noref.
** H.Bushouse	14-Feb-2007	Eliminated use of RebinRef, because we don't
**				want to extract subarrays from a full-frame
**				dark ref image, we want to instead use a
**				matching subarray dark ref image.
** H.Bushouse	14-May-2010	Added computation of MEANDARK.
*/

static int darkcorr (WF3Info *wf3, SingleNicmosGroup *input,
		     SingleNicmosGroup *dark) {

/* Arguments:
**	wf3	 i: WFC3 info structure
**	dark	 i: dark current image
**	input	io: image to be dark subtracted
*/

        /* Local variables */
        int i1, i2, j1, j2;             /* stats pixel limits */
        short dqmask;                   /* mask for stat rejs */
        float low, high;                /* stat rej limits */
        float mean, median, mode;       /* stats */
        float stdv, min, max;           /* more stats */

	/* Function definitions */
	void asub_noref (WF3Info *, SingleNicmosGroup *, SingleNicmosGroup *);
        int stats (SingleNicmosGroup *, int, int, int, int, float, float, short,
                   float *, float *, float *, float *, float *, float *);
	int PutKeyFlt (Hdr *, char *, float, char *);

	/* Report which ref file frames are being used */
	if (wf3->DarkType == MATCH) {
	    sprintf (MsgText,
	   "DARKCORR using dark imset %2d for imset %2d with exptime=%8.6g",
		     wf3->darkframe1, wf3->group, wf3->exptime[wf3->group-1]);
	    trlmessage (MsgText);

	} else if (wf3->DarkType == INTERP) {
	    sprintf (MsgText,
		     "DARKCORR using dark imsets %d and %d for imset %d",
		     wf3->darkframe1, wf3->darkframe2, wf3->group);
	    trlwarn (MsgText);
	    sprintf (MsgText,
		     "         interpolated to exptime=%g",
		     wf3->exptime[wf3->group-1]);
	    trlwarn (MsgText);

	} else if (wf3->DarkType == EXTRAP && wf3->darkframe1 != 0) {
	    sprintf (MsgText, "DARKCORR using dark imset %d for imset %d",
		     wf3->darkframe1, wf3->group);
	    trlwarn (MsgText);
	    sprintf (MsgText, "         extrapolated to exptime=%g",
		     wf3->exptime[wf3->group-1]);
	    trlwarn (MsgText);

	} else if (wf3->DarkType == EXTRAP && wf3->darkframe2 != 0) {
	    sprintf (MsgText, "DARKCORR using dark imset %d for imset %d",
		     wf3->darkframe2, wf3->group);
	    trlwarn (MsgText);
	    sprintf (MsgText, "         extrapolated to exptime=%g",
		     wf3->exptime[wf3->group-1]);
	    trlwarn (MsgText);
	}

	/* Do the dark subtraction in-place in input;
	** this subtraction does NOT include the reference pixels. */
	asub_noref (wf3, input, dark);

	/* Compute the mean of the dark image and populate MEANDARK */
	dqmask = 4+8+16+32+128+256+512;
	low = -FLT_MAX; high = FLT_MAX;
	i1 = wf3->trimx[0]; i2 = input->sci.data.nx - wf3->trimx[1] - 1;
	j1 = wf3->trimy[0]; j2 = input->sci.data.ny - wf3->trimy[1] - 1;
	if (stats (dark,i1,i2,j1,j2,low,high,dqmask,&mean,&median,&mode,
		   &stdv,&min,&max))
	    return (status = -1);

	if (PutKeyFlt (&input->sci.hdr, "MEANDARK", mean, ""))
	    return (status);
	
	/* Successful return */
	return (status = 0);
}

