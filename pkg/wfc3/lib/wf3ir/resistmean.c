# include <math.h>
# include <float.h>
# include <stdio.h>
# include <stdlib.h>
#include "hstcal.h"
# include "msg.h"
# include "trlbuf.h"
 
/* RESISTMEAN: Compute the clipped mean and stddev of  
   the values that are passed in the array. This was modeled
   on the resistant_mean procedure in IDL, and initially used
   to calculate the statistics for the reference pixels in 
   BLEVCORR for WFC3 IR images.

   Since this was written to deal with the bias pixels in the
   wfc3 arrays, I'm going to assume they wont exceed FLOAT
   type numbers....the science array has the possibility to
   go to float I believe, but the bias pixels should stay
   lower than that. It complicates the sorting and arrays
   to move stuff to float, but it can be done if necessary.

   Megan Sosey  August 2008

Processing outline:
  -get the median of the input array
  -subtract the median from every element in the array 
  -get the absolute values of the resulting differences
  -get the median of the absolute deviation array and divide by a constant 
   with some logic attached
  -compute the cutoff value in terms of the mean absolute deviation

  -make a new array that contains just the values from the original
   array where the absolute deviation is less than the cutoff value
     
  -get the average of the new array of "good" values
  -compute the standard deviation of the new array
  -correct the stddev to take into account the fact you left out values 
   of the original sample
  -recompute the cutoff using the corrected stddev
  
  -get a new array of good values using the new cutoff value
  -recompute the mean of the remaining good values
  -recompute the standard deviation of the remaining good values
  -correct the stddev to take into account the fact you left out values
   of the original sample
  
Revision history:

 M. Sosey	Aug. 2008	Created initial version
 H. Bushouse	10-Apr-2009	Modified findMean and findSigma routines to
				use double precision in sums. Also some general
				code cleanup throughout.
*/

#define PFLAG -9999. /*for rejecting pixels in the array*/

int resistmean (float *in, int npix, float sigrej, float *mean, 
		float *sigmean, float *min, float *max) {

/* Arguments:
**	in	i: pointer to array of input values 
**	npix	i: size of the input array
**	sigrej	i: the sigma level for pixel rejection
**	mean	o: mean value of unrejected pixels
**	sigmean	o: standard deviation of unrejected pixels
**		   (NOT the uncertainty or sigma of the mean)
**      min	o: minimum of input values
**      max	o: maximum of input values
*/

	/* Local variables */
	int i, j;	  /* loop indexes */
	float median;     /* median value of input array */
	float *absdev; 	  /* absolute deviation array */
	float medabsdev;  /* median of the absdev array */
	float *tempdata;  /* array of good pixel values */
	float *tempdata2; /* array of good pixel values */
	double sum;       /* sum of pixel values */
	float cutoff;	  /* value limitation */
	int   npix1, npix2;  /* number of pixels in temp arrays */
    
	/* Function definitions */
	float findRMedian (float *, int );        /* compute median of array */
   	float findMean(float *, int );              /* compute mean of array */
	float findSigma(float *, int , float *);  /* compute stddev of array */
	void arrayDiff(float *, float *, float , int );   /* subtract median */
	void arrayAbs(float *,int);       /* compute absolute value of array */
	void wheregood(float *, float , int * );   /* flag outliers in array */
    
	/* Initialize the counters and results */
	npix1 = npix;
	*mean = 0.;
	*sigmean = 0.;
	*min = 0;
	*max = 0;
	median = 0.;
	medabsdev = 0.;
	cutoff = 0.;

	/* allocate temp arrays to store and manipulate the input array */
	if (npix != 0) {
        
	    /* an array to store a copy of the data*/
	    tempdata = (float *) calloc(npix, sizeof(float));
	    if (tempdata == NULL) {
		sprintf (MsgText, "Memory allocation failure in resistmean");
        	trlmessage (MsgText);
		return (1);
	    }

	    /* and an array to store absolute deviation */
	    absdev = (float *) calloc(npix, sizeof(float));
	    if (absdev == NULL) {
		sprintf (MsgText, "Memory allocation failure in resistmean");
        	trlmessage (MsgText);
		return (1);
	    }

	} else {
	    sprintf (MsgText, "Zero size array passed to resistmean");
	    trlerror (MsgText);
	    return (1);
	}
              
	/* Copy the input array, computing the initial sum, min, max */
	sum = 0.;
	*min = tempdata[0];
	*max = tempdata[0];
	for (i=0; i<npix1; i++) {
	     tempdata[i] = in[i];
	     sum  += tempdata[i];
	     if (tempdata[i] < *min) *min = tempdata[i];
	     if (tempdata[i] > *max) *max = tempdata[i];
	}

	/* Compute the mean and median of the unrejected values */
	*mean = (float) (sum/(double)npix1);
	median = findRMedian (tempdata, npix1);

	/* Subtract the median from every element in the array */
	arrayDiff (tempdata, absdev, median, npix1);

	/* Get the absolute value for each element in the resulting array */
	arrayAbs (absdev, npix1);

	/* Get the median of the absolute deviation array and divide by a 
	   constant with some logic attached */
	medabsdev = findRMedian (absdev, npix1) / 0.6745; 
	if (medabsdev < 1.0E-24)
	    medabsdev = findMean(absdev,npix1) / 0.8;

	/* Compute the cutoff value in terms of the median absolute deviation */
	cutoff = (sigrej) * medabsdev;

	/* Flag values in the array that are beyond the cutoff */
	wheregood (absdev, cutoff, &npix1);

	/* Start a new temp array to populate with the good values;
	   npix1 should be the count of unrejected pixels here */
	tempdata2 = (float *) calloc(npix1, sizeof(float));
	if (tempdata2 == NULL) {
	    sprintf (MsgText, "Memory allocation failure in resistmean");
	    trlmessage (MsgText);
	    return (1);
	}

	/* Copy only the good values into the new temp array.
	   npix is still the number of original input values, while npix1
	   is now the number of remaining good values. */
	j=0;
	for (i=0; i<npix; i++) {
	     if (absdev[i] != PFLAG) {	
		 tempdata2[j]=tempdata[i];
		 j++;
	     }
	}

	/************ ROUND 2 ***********/

	/* Compute the mean and stddev of the values in the new array */
	*mean = findMean(tempdata2, npix1);    
	*sigmean = findSigma(tempdata2, npix1, mean); 

	/* Compensate sigma for truncation and compute new cutoff */
 	if (sigrej <= 4.5) {
	    *sigmean = *sigmean / (-0.15405+0.90723*sigrej - 0.23584*sigrej*sigrej+0.020142*sigrej*sigrej*sigrej);
	}
 	cutoff = sigrej * (*sigmean); 

	/* Reinitialize the absdev array to hold new set of values */
	free(absdev); /* clear the old array */
	absdev = (float *) calloc(npix1, sizeof(float));
	if (absdev == NULL) {
	    sprintf (MsgText, "Memory allocation failure in resistmean");
	    return (1);
	}
    
	/* Find the median of the good values. */
	median = findRMedian(tempdata2, npix1);

	/* Subtract the median from every element in the array */
	arrayDiff (tempdata2, absdev, median, npix1);

	/* Get the absolute value for each element in the resulting array */
	arrayAbs (absdev, npix1);

	/* Get the median of the absolute deviation array */
	medabsdev = findRMedian(absdev, npix1) / 0.6745; 
	if (medabsdev < 1.0E-24)
	    medabsdev = findMean(absdev, npix1) / 0.8;

	/* Flag values in the array that are beyond the new cutoff */
	npix2 = npix1;
	wheregood(absdev, cutoff, &npix2);

	/* Start a new temp array to populate with just the final good values */
	free (tempdata); /* free the previous copy */
	tempdata = (float *) calloc(npix2, sizeof(float));
	if (tempdata == NULL) {
	    sprintf (MsgText, "Memory allocation failure in resistmean");
	    return (1);
	}
    
	/* Copy only the good values to the new array */
	j=0;
	for (i=0; i<npix1; i++) {
	     if (absdev[i] != PFLAG) {	
		 tempdata[j] = tempdata2[i];
		 j++;
	     }
	}

	/* Compute the mean and stddev of the latest array of good values */
	*mean    = findMean (tempdata, npix2);
	*sigmean = findSigma(tempdata, npix2, mean); 

	/* Compensate sigma for truncation :*/
 	if (sigrej <= 4.5) {
	    *sigmean = *sigmean / (-0.15405+0.90723*sigrej-0.23584*sigrej*sigrej+0.020142*sigrej*sigrej*sigrej);
	}

	/* Free all local arrays */
	free(absdev);
	free(tempdata);
	free(tempdata2);

	return (0);    	/* Successful return */
}


/* Find the median vaue in the array. Handles odd and even number of inputs. */
float findRMedian (float *arr, int npts) {

	int i;
	float *tarr;
	float median;
    
	int sort (float *, int);  /* sort the given array*/

	median = 0.;

	/* Check for trivial cases */
	if (npts == 0) {
	    return(0.0);

	} else if (npts == 1) {
	    return(arr[0]);

	} else {
        
	    /*create a temporary array to sort the values*/
	    tarr = (float *) calloc(npts, sizeof(float));
	    if (tarr == NULL) {
		sprintf (MsgText, "Memory allocation failure in resistmean");
		trlmessage (MsgText);
		return (0.0);
	    }
    	    for (i=0; i<npts; i++)
        	tarr[i]=arr[i];
        
	    /* Sort the array of values */
	    if (sort(tarr-1, npts))
		return (0.0);

	    /* Find the median */
	    if ((npts % 2) == 0)
		median = 0.5 * (tarr[npts/2-1] + tarr[npts/2]);
	    else
		median = tarr[npts/2];

	    free(tarr);
	}

	return (median);
}

/*find the standard deviation of the given array */
float findSigma (float *arr, int npts, float *mean) {

	int i;
	double stdv, dmean, sum2;

	dmean = (double)(*mean);

	if (npts <= 1) {
	    stdv = 0.0;
	} else {
	    sum2 = 0.0; 
	    for (i = 0; i < npts; i++)
		 sum2 += arr[i]*arr[i];
        
	    stdv = (double)(npts/(npts-1.))*((sum2/(double)npts) - dmean*dmean);
	    if (stdv >= 0.)
		stdv = sqrt ((float)stdv);
	    else
		stdv = 0.0;
	}

	return ((float)stdv);
}

/* Take the difference of an array and a scalar */
void arrayDiff (float *in, float *out, float scalar, int npts) {

	int i;
    
	for (i=0; i<npts; i++)
	     out[i]=in[i]-scalar;   

	return;
}

/* Find the mean of the input array */
float findMean (float *in, int npts) {

	int i;
	double sum;

	sum = 0.;
	for (i=0; i<npts; i++)
	     sum+= in[i];

	return ((float)(sum/(double)npts));
}   

/* Take the absolute value of each item in the array */
void arrayAbs (float *arr, int npts) {

	int i;
    
	for (i=0; i<npts; i++)
	     arr[i] = fabsf(arr[i]);
    
	return;
}


/* This function will set the values in the input array to FLAG value
 which dont satisfly the criteria. You can then count them up inside
 a loop in the main program and make a new array. */

void wheregood (float *arr, float limit, int *npix) {
	int i;
	int nrej;
    
	nrej = 0;
	for (i=0; i< *npix; i++) {
	     if (arr[i] > limit) {
		 arr[i]=PFLAG;
		 nrej++;
	     }
	}
	*npix -= nrej;  /*now npix is the number of good pixels*/ 
        
	return;
}

#undef PFLAG
