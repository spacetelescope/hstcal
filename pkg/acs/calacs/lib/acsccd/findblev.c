# include <stdio.h>
# include <stdlib.h>		/* for calloc and qsort */

#include "hstcal.h"
# include "hstio.h"
# include "hstcalerr.h"	/* for NO_GOOD_DATA */
# include "acsdq.h"		/* for GOODPIXEL */
# include "acs.h"        /* for MsgText */

/* This routine determines the bias level for one line of an image
 by taking the median of the values in the overscan region.
 ** NO MAJOR MODIFICATIONS for ACS...
 */

int FindBlev (SingleGroup *x, int j, int *biassect, short sdqflags,
              double *biaslevel, int *npix) {
  
  /* arguments:
   SingleGroup *x      i: needed for science data and data quality
   int j               i: line of x data to use to get overscan
   int biassect[4]     i: beginning and end of region(s) to use for overscan
   double *biaslevel   o: median bias level for current (j) line
   int *npix           o: number of pixels used to compute bias level
   */
  
	extern int status;
  
	double *over;	/* values extracted from overscan region */
	int nvals;	/* number of good pixels extracted from overscan */
  int nbias;  /* number of pixels from both overscans to be used*/
	int i;
	int inplace = 1;	/* sort the array in-place */
  
	double MedianDouble (double *, int, int);
  
	/* Allocate space for the overscan, and copy out good data. */
  /* If no trailing region is used, biassect[2] and [3] will be set
   to zero, so it will have no effect on nbias. */
  nbias = (biassect[1]-biassect[0]) + (biassect[3]-biassect[2]) + 2;
	over = calloc (nbias, sizeof (double));
  
	nvals = 0;
  
  /* check to see if first overscan region is used */
  if (biassect[1] != 0) {
    for (i = biassect[0];  i <= biassect[1];  i++) {
      if (DQPix (x->dq.data, i, j) == GOODPIXEL || !(DQPix (x->dq.data, i, j) & sdqflags)) {
        over[nvals] = Pix (x->sci.data, i, j);
        nvals++;
      }
    }
  }
  
  /* Check to see if second overscan region is used...*/
  if (biassect[2] != 0){
    /* Include the pixels from 2nd overscan region in fit here... */
    for (i = biassect[2];  i <= biassect[3];  i++) {
      if (DQPix (x->dq.data, i, j) == GOODPIXEL || !(DQPix (x->dq.data, i, j) & sdqflags)) {
		    over[nvals] = Pix (x->sci.data, i, j);
		    nvals++;
      }
    }        
  } /* End of 2nd overscan pixels */
  
	*npix = nvals;
  
	if (nvals < 1) {
    free (over);
    return (status = NO_GOOD_DATA);
	}
  
	/* Find the median. */
  *biaslevel = MedianDouble (over, nvals, inplace);
  
	free (over);
  
	return (status);
}
