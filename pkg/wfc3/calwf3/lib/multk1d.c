# include "hstio.h"
# include "wf3.h"
# include "wf3info.h"
# include "hstcalerr.h"
# include <string.h>

/* The science data and the error array values are multiplied by k;
   the data quality array is not modified.

   Howard Bushouse, 2001 Nov 16:
	Revised to track CALACS changes - Modified calling sequence and
	internal operation of multgn1d.
    
   Megan Sosey, 2013 Sep 6
   Added new routine to deal with subarrays which span amp regions
*/

int multk1d (SingleGroupLine *a, float k) {

/* arguments:
SingleGroupLine *a	io: input data; output product
float k			 i: multiply a by this constant
*/

	int i;
	int dimx;

	if (k == 1.)
	    return (0);
       
	/* science data */
	/* error array */
	dimx = a->sci.tot_nx;

	for (i = 0;  i < dimx;  i++) {
	     a->sci.line[i] = k * a->sci.line[i];
	     a->err.line[i] = k * a->err.line[i];
	}
    
	return (0);
}

int multkline (float *a, float k, int dimx) {

/* arguments:
float   *a              io: input data; output product
float   k            	 i: multiply a by this constant
int     dimx             i: size of array
*/
	int i;

	if (k == 1.)
	    return (0);

	for (i = 0;  i < dimx;  i++)
	     a[i] = k * a[i];

	return (0);
}

/* This routine computes the average of the science data values
** from a single line of the image
** that are not flagged as bad by the data quality array.
*/

void AvgSciValLine (SingleGroupLine *y, short sdqflags, float *mean,
		    float *weight) {

	double sum;
	int numgood;			/* number of good pixels */
	int i;
	short flagval;			/* data quality flag value */
	int dimx;

	numgood = 0;
	sum = 0.;
	dimx = y->sci.tot_nx;
    
	for (i = 0;  i < dimx;  i++) {
	     flagval = y->dq.line[i];
	     if ( ! (sdqflags & flagval) ) {
		 /* no serious flag bit set */
		 sum += y->sci.line[i];
		 numgood++;
	     }
	}

	if (numgood > 0) {
	    *mean = sum / (double) numgood;
	    *weight = (float) numgood / (float)y->sci.tot_nx;
	} else { 
	    *mean = 0.;
	    *weight = 0.;
	}
}

void multgn1dsub(SingleGroupLine *a, int line, float *gain, float k0, char *ccdamp) {

/* arguments:
SingleGroupLine *a	io: input data; output product
int line		i: line number of input line
int ampx, ampy		i: amp parameters
float *gain		i: atogain parameter
float k0		i: value for scaling image data
			(exptime or flashdur or ...)
            
This uses logic which is good for subarrays    
            
*/
	int i;
	float k;
	int dimx;


	/* Determine k */
	k = 0;
	dimx = a->sci.tot_nx;
            
    /* Since both the science and error arrays operate on the 
        same line of data, we will combine them into one loop over
        X values.
    */
    
    if (strcmp (ccdamp, "A")==0) {   
            k = k0 / gain[AMP_A];
    } else if (strcmp (ccdamp, "B")==0) {   
             k = k0 / gain[AMP_B];
    } else if  (strcmp (ccdamp, "C")==0) {   
             k = k0 / gain[AMP_C];
    } else if  (strcmp (ccdamp, "D")==0) {   
             k = k0 / gain[AMP_D];
    } else {
        trlerror("Bad AMP assignment in multgn1dsub");
    }
             
    /* Apply scale factors  */
    for (i = 0;  i < dimx;  i++) {
	     a->sci.line[i] = k * a->sci.line[i];
	     a->err.line[i] = k * a->err.line[i];        
    }
 
}


void multgn1d (SingleGroupLine *a, int line, int ampx, int ampy, float *gain,
	       float k0) {

/* arguments:
SingleGroupLine *a	io: input data; output product
int line		i: line number of input line
int ampx, ampy		i: amp parameters
float *gain		i: atogain parameter
float k0		i: value for scaling image data
			(exptime or flashdur or ...)
            
None of this logic is taking into account the subarray cases where only 1 amp
is used to read out the entire array, no mater where it is taken from            
            
*/
	int i;
	float k;
	int dimx;


	/* Determine k */
	k = 0;
	dimx = a->sci.tot_nx;
    
    /* Since both the science and error arrays operate on the 
        same line of data, we will combine them into one loop over
        X values.
    */
	/* science data and error array */
	if (line < ampy) {
        if (gain[AMP_C] > 0.) {
            k = k0 / gain[AMP_C];
        } else {
            k = 0.;
        } 
    
        /* This line has 2-AMP readout */
        /* Apply gain for first amp to first part of line */
        for (i = 0;  i < ampx;  i++) {
	     a->sci.line[i] = k * a->sci.line[i];
	     a->err.line[i] = k * a->err.line[i];        
        }
        /* Apply gain for second amp over remainder of line...
            only if there is a second amp which needs applying.
            WJH 14 Apr 2000
        */
        if (ampx < dimx) {
            k = k0 / gain[AMP_D];
            for (i = ampx;  i < dimx;  i++) {
		        a->sci.line[i] = k * a->sci.line[i];
		        a->err.line[i] = k * a->err.line[i];        
            }
        }

    } else {
        if (gain[AMP_A] > 0.) {
            k = k0 / gain[AMP_A];
        } else {
            k = 0.;
        }
        /* This line has 2-AMP readout */
        /* Apply gain for first amp to first part of line */
        for (i = 0;  i < ampx;  i++) {
	     a->sci.line[i] = k * a->sci.line[i];
	     a->err.line[i] = k * a->err.line[i];        
        }
        /* Apply gain for second amp over remainder of line... 
            only if there is a second amp which needs applying.
            WJH 14 Apr 2000
        */
        if (ampx < dimx) {
            k = k0 / gain[AMP_B];
            for (i = ampx;  i < dimx;  i++) {
		 a->sci.line[i] = k * a->sci.line[i];
		 a->err.line[i] = k * a->err.line[i];        
            }
        }
    }

}

