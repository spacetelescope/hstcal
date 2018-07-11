# include   <stdio.h>
# include   <string.h>
# include   <stdlib.h>
# include   <limits.h>
# include   <math.h>

#include "hstcal.h"
# include "hstio.h"
# include   "wf3.h"    /* for message output */
# include   "hstcalerr.h"
# include   "wf3info.h"
# include   "rej.h"

# define    MINVAL      -15000
# define    BIN_WIDTH   1
# define    MIN_BINS    1000
# define    MAX_BINS	10000
# define    SIGREJ	4.0  /* resistmean sigma rejection threshold */

/* rej_sky -- Calculate the sky for an image. */

void rej_sky (char *sky, IODescPtr ipsci[], IODescPtr ipdq[], int nimgs,
	      short badinpdq, float efac[], DataUnits bunit[], float skyval[]) {

/* Revision history:
**
** H. Bushouse	18-Jun-2002	Made bin width checking more robust. Also
**				allocate and free histogram for each image
**				(following CALACS changes).
** H. Bushouse	06-Dec-2007	Added calls to getHeader before each call to
**				getShortLine to prevent getShortLine from
**				crashing on null DQ arrays.
** H. Bushouse	22-May-2008	Avoid arithmetic overflow in binning
**				calculations.
** H. Bushouse	08-Oct-2008	Added capabilities for "mean" sky calculation
**				mode, using resistmean function.
** H. Bushouse	14-Dec-2011	Upgraded to rescale input data that are in
**				units of count rates. (PR 69969; Trac #814)
*/

    extern int status;

    int         *histgrm;   /* pointer to the histogram */
    int         nbins;      /* number of bins */
    int         min, max;
    float       data_min;
    float       hwidth;     /* histogram resolution */
    float       hmin;       /* minimum histogram value */
    int         i, k, h;
    float       *a;
    short       *b;
    int         line, npt;
    int         dimx, dimy;
    float       sum, mean;
    Hdr		dqhdr;

    Bool   mode, rmean;	    /* sky calculation mode flags */
    float *skyarr;	    /* pointer to sky values array */
    float  ssig, smin, smax; /* sky resistantmean values */

    float cr_mode (int *, int, float, float);
    int	  resistmean (float *, int, float, float *, float *, float *, float *);

    /* -------------------------------- begin ------------------------------- */

    nbins=0;
    min=0;
    max=0;
    hwidth=0.0f;
    hmin=0.0f;
    npt=0;
    histgrm=NULL;
    skyarr=NULL;

    /* decide what to do according to sky */
    if (strcmp (sky, "none") == 0) {
        for (k = 0; k < nimgs; ++k) {
            skyval[k] = 0.;
        }
        return;
    } else if (strcmp (sky, "mode") == 0) {
	mode = True;
	rmean = False;
    } else if (strcmp (sky, "mean") == 0) {
	rmean = True;
	mode = False;
    } else {
        trlerror ("illegal sky value");
        status = INVALID_VALUE;
        return;
    }

    dimx = getNaxis1(ipsci[0]);
    dimy = getNaxis2(ipsci[0]);

    a = (float *) calloc (dimx, sizeof(float));
    b = (short *) calloc (dimx, sizeof(short));
    if (a == NULL || b == NULL) {
        trlerror ("Couldn't allocate memory for sky arrays");
        status = OUT_OF_MEMORY;
        return; 
    }

    if (mode) {

	/* compute MIN and MAX values for data */
	/* use the minimum and twice of the mean to determine the data range */
	data_min = INT_MAX;
	sum = 0.;
	npt = 0;

	initHdr (&dqhdr);
	getHeader(ipdq[0], &dqhdr);

	for (line = 0; line < dimy; line++) {

	     /* read the data in */
	     getFloatLine (ipsci[0], line, a);
	     getShortLine (ipdq[0], line, b);

	     /* Rescale data to counts, if necessary */
	     if (bunit[0] == COUNTRATE) {
		 for (i = 0; i < dimx; ++i) {
		      a[i] *= efac[0];
		 }
	     }

	     for (i = 0; i < dimx; ++i) {
		  if ( (b[i] & badinpdq) == WF3_OK) {
			data_min = (a[i] < data_min) ? a[i] : data_min;
			sum += a[i];
			npt++;
		  }
	     }
	} /* End of loop over lines */

	freeHdr(&dqhdr);

	/* Compute min and max for histogram.
	MIN is min of good data or MINVAL, which ever is greater
	DELTA is difference between mean of data and MIN
	MAX is mean plus DELTA
	This insures that the mean falls in the center of the range
	between min and max. */
	if (npt == 0) npt = 1;
	min = (data_min < MINVAL) ? MINVAL : data_min;
	mean = (sum > 0.) ? (int) ( (sum / (float)npt) + 1) : 1;
	max = 2 * mean - min;
    
	/* use the mode as sky value, use the bin size of 1 (DN) */
	nbins = max - min + 1; 

	/* Insure that there are at least MIN_BINS in the histogram
	and reset the width accordingly. */ 
	if (nbins < MIN_BINS) {
	    nbins = MIN_BINS;
	    hwidth = (float) nbins / (float)MIN_BINS;
	} else if (nbins > MAX_BINS) {
	    hwidth = (float) nbins / (float)MAX_BINS;
	    nbins = MAX_BINS;
	} else {
	    hwidth = 1.;
	}
	hmin = (float) min;

	/*
	sprintf (MsgText, 
	"sky computed using min %d, max %d, mean %g, bins %d, and width %g",
	min, max, mean, nbins, hwidth);
	trlmessage (MsgText);
	*/
    }

    /* Now loop over the input images, computing the sky value for each
    ** image, using either the mode or the resistant mean */
    for (k = 0; k < nimgs; ++k) {

	if (mode) {
	    /*  set up a new histogram array for each image */
	    histgrm = (int *) calloc (nbins, sizeof(int));
	    if (histgrm == NULL){
		    trlerror ("Couldn't allocate memory for sky histogram array");
		    status = OUT_OF_MEMORY;
		return;   
	    }
	} else if (rmean) {
	    skyarr = (float *) calloc (dimx*dimy, sizeof(float));
	    if (skyarr == NULL){
		trlerror ("Couldn't allocate memory for sky array");
		status = OUT_OF_MEMORY;
		return;   
	    }
	    npt = 0;
	}

        initHdr (&dqhdr);
        getHeader (ipdq[k], &dqhdr);
        for (line = 0; line < dimy; line++) {

            /* read the data in */
            getFloatLine (ipsci[k], line, a);
            getShortLine (ipdq[k], line, b);

	    /* Rescale data to counts, if necessary */
	    if (bunit[k] == COUNTRATE) {
		for (i = 0; i < dimx; ++i) {
		     a[i] *= efac[k];
		}
	    }

            /* construct the histogram for the mode calculation */
	    if (mode) {
		for (i = 0; i < dimx; ++i) {
		     if (b[i] == WF3_OK) {

			 /* Adjust the bin position by the width of each bin */
			 if (fabs((a[i]-min)/hwidth) < INT_MAX) {
			     h = (int)((a[i] - min) / hwidth);
			     if (h >= 0 && h < nbins)
				 histgrm[h]++;
			 }
		     }
		}

	    /* load the sky array for calculating the mean */
	    } else if (rmean) {
		for (i = 0; i < dimx; ++i) {
		     if (b[i] == WF3_OK) {
			 skyarr[npt] = a[i];
			 npt++;
		     }
		}
	    }
        } /* End of loop over lines */
        freeHdr(&dqhdr);

        /* calculate the mode from the histogram */
	if (mode) {
            skyval[k] = cr_mode (histgrm, nbins, hwidth, hmin);
            free (histgrm);

	/* calculate the resistant mean */
	} else if (rmean) {
	    resistmean (skyarr, npt, SIGREJ, &skyval[k], &ssig, &smin, &smax);
	    free (skyarr);
	}
    }

    free(a);
    free(b);
}
