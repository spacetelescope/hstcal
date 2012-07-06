# include	<stdio.h>
# include	<string.h>
# include	<stdlib.h>
# include	"hstio.h"

# define	MINVAL	-1000

/* crrej_sky -- Calculate the sky for an image.

  Description:
  ------------
  Computes the sky level which will be subtracted during cosmic ray
  rejection and added back at the end of processing.
  
  Date		Author		Description
  ----		------		-----------
  	       J.-C. Hsu	Original
  10-Feb-2000  Phil Hodge	Return int instead of void,
				and replace exit with return;
				check for out of memory.
*/

int crrej_sky (char *sky, IODescPtr ipsci[], IODescPtr ipdq[], 
		int nimgs,    float skyval[])
{
	int	*histgrm;	/* pointer to the histogram */
	int	nbins;		/* number of bins */
	int	min, max, npt;
	float	amin, amax;
	float	hwidth;		/* histogram resolution */
	float	hmin;		/* minimum histogram value */
	float	sum;
	int	i, j, k, h;
	short	OK = 0;
	FloatTwoDArray	a;
	ShortTwoDArray	b;

	float 	cr_mode (int *, int, float, float);

/* -------------------------------- begin ---------------------------------- */

	/* decide what to do according to sky */
	if (strcmp (sky, "none") == 0) {
	    for (k = 0; k < nimgs; ++k) {
	        skyval[k] = 0.;
	    }
	    return (0);
	}
	if (strcmp (sky, "mode") != 0) {
	    printf ("illegal sky value\n");
	    return (2);
	}

	/* loop all images */
	initFloatData (&a);
	initShortData (&b);

	for (k = 0; k < nimgs; ++k) {

	    /* read the data in */
	    getFloatData (ipsci[k], &a);
	    if (hstio_err()) {
		printf ("ERROR    %s\n", hstio_errmsg());
		return (2);
	    }
	    getShortData (ipdq[k], &b);
	    if (hstio_err()) {
		printf ("ERROR    %s\n", hstio_errmsg());
		return (2);
	    }

	    amin = 0.;
	    npt = 0;
	    sum = 0;
	    for (j = 0; j < a.ny; ++j) {
	        for (i = 0; i < a.nx; ++i) {
		    if (DQPix(b,i,j) == OK) {
		    	sum += Pix(a,i,j);
			npt++;
			if (Pix(a,i,j) < amin) amin = Pix(a,i,j);
		    }
		}
	    }
	    
	    /* use the minimum and twice of the mean to determine the data 
		range */
		
	    if (amin < MINVAL) 
		min = MINVAL;
	    else
  	        min = (int) amin - 1;
	    amax = sum/(float)npt;
	    if (amax <= 0.) 
		max = 1;
	    else
		max = (int) (amax + 1.) * 2;
	    nbins = max - min + 1;

	    /* use the mode as sky value, use the bin size of 1 (DN) */
	    hmin = (float) min;
	    hwidth = 1.;

	    /*  set up the histogram array */
	    histgrm = calloc (nbins, sizeof(int));
	    if (histgrm == NULL) {
		printf ("ERROR    out of memory in crrej_sky\n");
		return (2);
	    }

	    /* construct the histogram */
	    for (j = 0; j < a.ny; ++j) {
	        for (i = 0; i < a.nx; ++i) {
		    if (DQPix(b,i,j) == OK) {
		    	h = (int) Pix(a,i,j) - min;
			if (h >= 0 && h < nbins)
			    histgrm[h]++;
		    }
		}
	    }
	
	    /* calculate the mode from the histogram */
	    skyval[k] = cr_mode (histgrm, nbins, hwidth, hmin);

	    freeFloatData(&a);
	    freeShortData(&b);
	    free (histgrm);
	}

	return (0);
}
