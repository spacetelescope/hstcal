# include <float.h> 
# include <math.h>
# include <stdio.h>

#include "hstcal.h"
# include "hstio.h"	/* defines HST I/O functions */
# include "wf3.h"
# include "wf3info.h"
# include "trlbuf.h"

extern int status;



/* DOBLEV: Apply blevcorr to each readout of a MultiAccum.
**
** Revision history:
** H.Bushouse	Oct. 2000	Initial CALNICA to CALWF3 port.
** H.Bushouse	13-Jan-2010	Added zoff image as I/O argument and call to
**				process it with blevcorr. This is necessary
**				due to reordering of zsig and blev steps.
**				(calwf3 v2.0)
*/

int doBlevIR (WF3Info *wf3, MultiNicmosGroup *input, SingleNicmosGroup *zoff) {

/* Arguments:
**	wf3	 i: WFC3 info structure
**	input	io: image to be bias corrected
*/

	/* Local variables */

	/* Function definitions */
	int blevcorr (WF3Info *, SingleNicmosGroup *);
	void PrSwitch (char *, int);
    
	/* Do the bias correction for each group */
        if (wf3->blevcorr == PERFORM) {

	    for (wf3->group=wf3->ngroups; wf3->group >= 1; wf3->group--) {
		 if (blevcorr (wf3, &(input->group[wf3->group-1])))
		     return (status);
            }
	    if (blevcorr (wf3, zoff))
		return (status);

	    PrSwitch ("blevcorr", COMPLETE);
	}

	/* Successful return */
	return (status = 0);
}

/* BLEVCORR: Bias level correction. Subtract the mean reference pixel value
** from each pixel in the image. The science image is modified
** in-place. No other image extensions are modified.
**
** Revision history:
** H.Bushouse	Oct. 2000	Initial CALNICA to CALWF3 port. Just a dummy
**				routine for now.
** H.Bushouse	02-May-2002	Initial attempt at a real algorithm.
** H.Bushouse	14-Feb-2005	Set reference pixel region boundaries from
**				biassect instead of trimx values.
** H.Bushouse	21-Apr-2005	Fixed bug in calculation of j2 values for
**				quads 3 and 4.
** H.Bushouse	01-Aug-2006	Enhancement to handle transposition of raw
**				images in OPUS. Swap the ordering of the
**				quads based on OPUS processing date.
** H.Bushouse	20-Aug-2007	Remove enhancements implemented in previous
**				version, which handled different image
**				rotations, because all old data have been
**				reprocessed to new orientation, and reset all
**				quad boundaries for new orientation.
** M. Sosey	August 2008 	The reference pixel subtraction has been 
**				changed so that	the mean of all the reference 
**				pixels in the image is subtracted from the
**				the whole image, rather than computing and
**				subtracting a separate mean for each quad.
**				A method similar to the resistant_mean routine
**				in IDL was implemented in a new statistics 
**				routine that just accepts an array of values,
**				rather than a SingleNicmosGroup pointer. The 
**				new statistics routine is "resistmean".
*/

int blevcorr (WF3Info *wf3, SingleNicmosGroup *input) {

/* Arguments:
**	wf3	 i: WFC3 info structure
**	input	io: image to be bias corrected
*/

	/* Local variables */
	int i, j, q;		/* loop limits */
	int i1, i2, j1, j2;	/* pixel indexes */
	float mean;		/* ref pixel stats */
	float stdv, min, max;	/* more stats */
	float *refpix;		/* array to hold all the reference pixels */
	int arrsize;		/* size of ref pixel array */        
	float sigrej;		/* rejection limit for statistics */
	int pixcount;		/* counter for reference pixels */
           
	/* Function definitions */
	int resistmean(float *, int, float, float *,float *,float *,float *);

	/* Set defaults */
	sigrej = 3.0; /*sigma rejection limit for mean calculation*/
    
	/* Allocate memory for the temporary reference pixel array:
	   there are 5 pixels on each end of each row, but 1 is ignored
	   on each side for a total of 8 being used per row */
	arrsize= (input->sci.data.ny) * 8;
	refpix = (float *) calloc(arrsize, sizeof(float));
	if (refpix == NULL) {
	    sprintf (MsgText, "Memory allocation failure in blevcorr");
	    trlerror (MsgText);
	    return (status = 1);
	}
    
	/* zero out the memory here just to be sure */
	for (i=0; i<arrsize; i++)
	     refpix[i]=0.0;
        
	/* Loop over the 4 quads of the image to gather the reference 
	   pixels into a new array that can be passed to the resistant 
	   mean stats routine */
	pixcount = 0;
	for (q = 1; q <= 4; q++) {      
        
	     /* Set the bounds of the ref pixels in each quadrant;
	     ** note that these are zero-indexed. */
	     switch (q) {
		case 1: i1 = wf3->biassecta[0];
			i2 = wf3->biassecta[1];
			j1 = (input->sci.data.ny / 2);
			j2 = input->sci.data.ny - wf3->trimy[1] - 1;
			for (j=j1; j<=j2; j++) {
			     for (i=i1; i<=i2; i++) {
				  refpix[pixcount]=Pix(input->sci.data,i,j);
				  pixcount++;
			     }
			}
			break;

		case 2: i1 = wf3->biassecta[0];
			i2 = wf3->biassecta[1];
			j1 = wf3->trimy[0];
			j2 = (input->sci.data.ny / 2) - 1;
			for (j=j1; j<=j2; j++) {
			     for (i=i1; i<=i2; i++) {
				  refpix[pixcount]=Pix(input->sci.data,i,j);
				  pixcount++;
			     }
			}
			break;

		case 3: i1 = wf3->biassectb[0];
			i2 = wf3->biassectb[1];
			j1 = wf3->trimy[0];
			j2 = (input->sci.data.ny / 2) - 1;
			for (j=j1; j<=j2; j++) {
			     for (i=i1; i<=i2; i++) {
				  refpix[pixcount]=Pix(input->sci.data,i,j);
				  pixcount++;
			     }
			}
			break;

		case 4: i1 = wf3->biassectb[0];
			i2 = wf3->biassectb[1];
			j1 = (input->sci.data.ny / 2);
			j2 = input->sci.data.ny - wf3->trimy[1] - 1;
			for (j=j1; j<=j2; j++) {
			     for (i=i1; i<=i2; i++) {
				  refpix[pixcount]=Pix(input->sci.data,i,j);
				  pixcount++;
			     }
			}
	     }
	                  
	}
     
	/* Compute stats of the ref pixels */
	if (resistmean(refpix, pixcount, sigrej, &mean, &stdv, &min, &max))
	    return (status = 1);
         
	/* Record the computed mean in the header for the science extension*/
	if (putKeyF (&input->sci.hdr, "MEANBLEV", mean, ""))
	    return (status = 1);
     
	/* Subtract the mean value from entire image */
	for (j=0; j < input->sci.data.ny; j++) {
	     for (i=0; i < input->sci.data.nx; i++) {
		  Pix(input->sci.data,i,j) -= mean;
	     }
	}

	/* Free local memory */
	free (refpix);

	/* Successful return */
	return (status = 0);
}

