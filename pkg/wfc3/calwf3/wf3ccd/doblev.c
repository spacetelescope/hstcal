/* This file contains:
	doBlev
	FitToOverscan
*/

# include <stdio.h>
# include <math.h>		/* sqrt */
# include <string.h>

# include "hstio.h"

# include "wf3.h"
# include "wf3info.h"
# include "hstcalerr.h"
# include "wf3dq.h"		/* for BADFLAT */

static void FitToOverscan (SingleGroup *, int, int, int *, float, short, float);

/* This routine subtracts the bias determined from the overscan region
   on each line.  Only a portion of the trailing overscan region will be
   used to find the bias level.

   If an output text file for bias levels (third command-line argument)
   was specified, this routine also creates this file, writes the bias
   level for each image line, and then closes the file.

   The CRPIX and LTV keyword values in the output extension headers will be
   updated.

   NOTE:  This task now performs the subtractions in-place on x, resulting
   in output which still has the overscan regions.  

   Warren Hack, 1998 June 2:
   	Revised for ACS data. 
   Howard Bushouse, 2000 Aug 29:
	Revised for WFC3 data.
   H.Bushouse, 2001 May 7:
	Added check to prevent computed dimensions from extending past edge
	of chip (same as calacs change).
   H.Bushouse, 2001 Nov 16:
	Changes to track CALACS updates - Revised to use CCDOFST[A,B,C,D] to
	select default bias levels from CCDBIAS[A,B,C,D] (which are new to
	CCDTAB) in this version. Updated to use VY values of <=0 to detect
	lack of virtual overscan use and NOT apply BlevDrift().
   H.Bushouse, 2002 March 1:
	Upgraded to handle WFC3 CCD serial virtual overscan region: changed
	the computations of begx,endx to take extra trimx3,trimx4 regions
	into account. 
	Upgraded to handle splitting the parallel virtual overscan region, as
	defined by VX,VY vectors, into two separate regions; one region for
	each amp. Changed the call to BlevDrift to pass the appropriate set of
	VX,VY values for each amp.
	Fixed a bug in the loop that subtracts the bias and drift levels,
	which was causing the "di" column index to be incremented twice.
	Modified the FitToOverscan routine so that the fit is computed as a
	function of input image line number, not output (trimmed) line number.
   H.Bushouse, 2003 Oct. 23:
	Modified for WFC3 to use the serial physical overscan regions to
	measure bias if neither of the serial virtual overscan regions are
	specified (e.g. for subarray images that don't have virtual overscan
	but do contain some physical overscan). Also added cleanBiasFit
	routine and associated changes for removing outliers from array of
	bias values (CALACS changes).
   H.Bushouse, 2004 Feb 19:
	Modified CleanBiasFit to use readnoise in more robust estimate of
	outliers in bias region. Modified calls to all higher-level routines
	to pass readnoise value into CleanBiasFit (following CALACS changes).
   H.Bushouse, 2006 Jul 26:
	Modified doBlev to send the readnoise value for the correct amp
	to the FitToOverscan routine and also convert the readnoise value
	to units of DN's so that it's the same as the science data.
   H.Bushouse, 2008 Sep 17:
	Pass readnoise to BlevDrift routine for use in cleanDriftFit. 
	Modified cleanBiasFit to use different clip values on each pass.
   H.Bushouse, 2009 Jan 16:
	Upgraded the methods used in cleanBiasFit to compute the mean and
	standard deviation of unrejected values so that they use only the
	good values returned by FindBlev. Also added checks for potential
	divide-by-zero conditions.
   H.Bushouse, 2009 Apr 24:
	Added print statements for verbose mode to give the overscan column
	limits being used.
*/

int doBlev (WF3Info *wf3, SingleGroup *x, int chip, float *meanblev, 
	    int *overscan, int *driftcorr) {

/* arguments:
WF3Info *wf3     i: calibration switches, etc
SingleGroup *x	io: image to be calibrated
int chip         i: chip number
float meanblev   o: mean value of bias levels that were subtracted
int *overscan   io: true if bias determined from an overscan region in x
int *driftcorr   o: true means correction for drift along lines was applied
*/

	extern int status;

	double biaslevel;	 /* value to subtract from a line of image */
	double averagedrift;	 /* drift averaged over all (output) columns */
	double sumbias, sumblev; /* sum of bias levels (for getting getting */
	int binx, biny;		 /* bin factors */
	int trimx1, trimx2;	/* width to trim off ends of each line */
	int trimx3, trimx4;	/* width to trim off middle of each line */
	int trimy1, trimy2;	/* amount to trim off ends of each column */
	int biassect[4];	/* section(s) to use for finding bias level */
	int vx[2], vy[2];	/* section(s) to use for finding bias drift */
	int i, j;
	double di, dj;
	int endx, endy;		/* boundaries of actual image */
	int begx, begy;		/* boundaries of actual image */
	int sizex, sizey;
	int biasnum;		/* array position for single amp used */
	int amp;		/* Counter for amps used */
	int numamps;		/* Number of AMPS used to readout chip */
	char *ccdamp;		/* Amps which are used for this chip */
	int dodrift = YES;	/* use virtual overscan to get blevdrift? */
	char *amploc; 
	int bias_loc, amp_indx;
	int bias_ampx, bias_ampy;
	int bias_orderx[4] = {0,1,0,1};
	int bias_ordery[4] = {0,0,1,1};
	double deval;
	double ampslope, ampintercept; /* value from the bias level fit */
	float ccdbias;		/* default bias level from CCDTAB */
	float rn;		/* readnoise value */

	/* Function definitions */
	int BlevDrift (SingleGroup *, int *, int *, int, int *, int *, short,
		       float);
	double DriftEval (double);
	double DriftMean (double);
	double BlevEval (double);
	void BlevResults (double *, double *);
	void parseWFCamps (char *, int, char *);
	int selectBias (char *);

	/* Allocate space for ccdamp... */
	ccdamp = (char *) calloc (NAMPS+1, sizeof(char));

	/* initial values */
	*driftcorr = 0;
	*meanblev = 0.;
	
	binx = wf3->bin[0];
	biny = wf3->bin[1];

	if ((binx != 1 && binx != 2 && binx != 3) ||
	    (biny != 1 && biny != 2 && biny != 3)) {
	    trlerror ("(doBlev) bin size must be 1, 2, or 3.");
	    free (ccdamp);
	    return (status = 1001);
	}

	/* Get biassect, the location of the region to use for determining
	   the bias level, and vx & vy, the location in the virtual
	   overscan to use for determining the drift with column number.
	   Also get the widths of the trim regions, the amount to remove
	   from each axis of the input when copying to output.
	*/
	/* ASSUMPTION: Only have no overscan with SINGLE AMP readouts
	   Only use FIRST AMP listed in CCDAMP to determine which
	   CCDOFST[A,B,C,D] value to use.
	*/

	biasnum = selectBias (wf3->ccdamp);
	if (strlen(wf3->ccdamp) > 1) {
	    if (wf3->detector == CCD_DETECTOR && wf3->chip == 2)
		biasnum += 2;
	}
	ccdbias = wf3->ccdbias[biasnum];

	if (*overscan != YES) {

	    /* If no overscan region, subtract BIAS level obtained from CCDTAB*/
	    trlwarn ("Overscan region is too small to do BLEVCORR; ");
	    sprintf (MsgText,"     bias from CCDTAB of %g will be subtracted.",
		     ccdbias);
	    trlmessage (MsgText);

	    for (j = 0;  j < x->sci.data.ny;  j++) {
		 for (i = 0;  i < x->sci.data.nx;  i++)
		      Pix (x->sci.data,i,j) -= ccdbias;
	    }
	    *meanblev = ccdbias;
	    wf3->blev[biasnum] = ccdbias;
	    free (ccdamp);
	    return (status);
	}

	/* Determine bias levels from BIASSECT columns specified in OSCNTAB */
	
	/* Copy out overscan info for ease of reference in this function*/
	trimx1 = wf3->trimx[0];
	trimx2 = wf3->trimx[1];
	trimx3 = wf3->trimx[2];
	trimx4 = wf3->trimx[3];
	trimy1 = wf3->trimy[0];
	trimy2 = wf3->trimy[1];

	/* Establish which amps from ccdamp string are appropriate 
	** for this chip */
	ccdamp[0] = '\0';
	
	/* Set up the 2 AMPS used per line */ 
	if (wf3->detector == CCD_DETECTOR) {
	    parseWFCamps (wf3->ccdamp, chip, ccdamp);
	}

	/* How many amps are used for this chip */	
	numamps = strlen (ccdamp);
	
	/* Are we going to calculate drift in bias from virtual overscan? */
	/* If the end points of vx and vy are zero, no section was specified */
	if (wf3->vy[0] <= 0 && wf3->vy[1] <= 0) {
	    dodrift = NO;	
	    trlmessage("(blevcorr) No virtual overscan region specified.");
	    trlmessage("(blevcorr) Bias drift correction will not be applied.");
	}

	/* For each amp used, determine bias level and subtract from 
	** appropriate region of chip */
	sumblev = 0.;
	for (amp = 0; amp < numamps; amp++) {
	     bias_loc = 0;

	     /* determine which section goes with the amp */
	     /* bias_amp = 0 for BIASSECTA, = 1 for BIASSECTB */
	     amploc = strchr (AMPSORDER, ccdamp[amp]);
	     bias_loc = *amploc - *ccdamp;
	     amp_indx = *amploc - *AMPSORDER;
	     bias_ampx = bias_orderx[bias_loc];
	     bias_ampy = bias_ordery[bias_loc];

	     /* Requirement: at least 1 section should be specified!	 */
	     /* If both bias sections are specified,... */
	     if (wf3->biassectc[1] > 0 && wf3->biassectd[1] > 0) {

		 /* select section nearest the amp based on bias_amp */
		 biassect[0] = (bias_ampx == 0) ? wf3->biassectc[0] :
						  wf3->biassectd[0];
		 biassect[1] = (bias_ampx == 0) ? wf3->biassectc[1] :
						  wf3->biassectd[1];
		 /* If only 1 AMP is used for each line, use both sections */
		 if (numamps == 1) {
		     biassect[2] = (bias_ampx == 0) ? wf3->biassectd[0] :
						      wf3->biassectc[0];
		     biassect[3] = (bias_ampx == 0) ? wf3->biassectd[1] :
						      wf3->biassectc[1];
		 } else {
		     biassect[2] = 0;
		     biassect[3] = 0;
		 }

	     /* if neither of the virtual overscan regions are specified
	     ** use the physical overscan instead */
	     } else if (wf3->biassectc[1] <= 0 && wf3->biassectd[1] <= 0) {

		 biassect[0] = (wf3->biassecta[1] == 0) ? wf3->biassectb[0] :
							  wf3->biassecta[0];
		 biassect[1] = (wf3->biassecta[1] == 0) ? wf3->biassectb[1] :
							  wf3->biassecta[1];
		 /* then set trailing section to zero indicating it's not used*/
		 biassect[2] = 0;
		 biassect[3] = 0;
		 
	     } else {

		 /* otherwise, select the non-zero bias section */
		 biassect[0] = (wf3->biassectc[1] == 0) ? wf3->biassectd[0] :
							  wf3->biassectc[0];
		 biassect[1] = (wf3->biassectc[1] == 0) ? wf3->biassectd[1] :
							  wf3->biassectc[1];
		 /* then set trailing section to zero indicating it's not used*/
		 biassect[2] = 0;
		 biassect[3] = 0;
	     }
	     if (wf3->verbose) {
		 sprintf (MsgText, "Using overscan columns %d to %d",
			  biassect[0]+1, biassect[1]+1);
		 trlmessage (MsgText);
		 if (biassect[2] != 0) {
		 sprintf (MsgText, "           and columns %d to %d",
		 	  biassect[2]+1, biassect[3]+1);
		 trlmessage (MsgText); }
	     }

	     /* Compute range of pixels affected by each amp */	
	     begx = trimx1 + (wf3->ampx + trimx3 + trimx4) * bias_ampx;
	     endx = (bias_ampx == 0 && wf3->ampx != 0) ? wf3->ampx + trimx1 :
							x->sci.data.nx - trimx2;
	     begy = trimy1 + wf3->ampy* bias_ampy;
	     endy = (bias_ampy == 0 && wf3->ampy != 0) ? wf3->ampy +trimy1 :
							x->sci.data.ny - trimy2;

	     /* Make sure that endx and endy do not extend beyond the bounds of
	     ** the image... HAB 7 May 2001 (same as calacs change). */
	     if (endx > x->sci.data.nx) endx = x->sci.data.nx;
	     if (endy > x->sci.data.ny) endy = x->sci.data.ny;
	     sizex = endx - begx;
	     sizey = endy - begy;		

	     /* Compute the actual readnoise for this amp, in units
	     ** of DN's */
	     if (wf3->atodgain[amp_indx] != 0)
	         rn = wf3->readnoise[amp_indx] / wf3->atodgain[amp_indx];
	     else
		 rn = wf3->readnoise[amp_indx];

	     /* At this point, we decide how we are to determine the bias level.
	     ** Average two sections across the line for 1 amp, or use only
	     ** one section. */

	     if (dodrift == YES) {

		 /* HAB 2002 Mar 1: This section added to accomodate splitting
		 ** the parallel overscan region into two; one for each amp. */
		 if (bias_ampx == 0) {
		     vx[0] = wf3->vx[0];
		     vx[1] = wf3->vx[1];
		     vy[0] = wf3->vy[0];
		     vy[1] = wf3->vy[1];
		 } else {
		     vx[0] = wf3->vx[2];
		     vx[1] = wf3->vx[3];
		     vy[0] = wf3->vy[2];
		     vy[1] = wf3->vy[3];
		 }

		 /* Fit a line to the virtual overscan region as a function of
		 ** column number.  */
		 if (BlevDrift (x, vx, vy, trimx1, biassect,
				driftcorr, wf3->sdqflags, rn)) {
        	     free (ccdamp);
		     return (status);
        	 }

		 /* Evaluate the fit for each line and subtract from the data.*/
		 averagedrift = DriftMean ((double)sizex);

	     } else {

		 averagedrift = 0.;
	     }
	     dj = (double)begy;
	     sumbias = 0.;

	     /* For each image line, determine the bias level from the
	     ** overscan in x, and fit to these values as a function of
	     ** line number in the input image.  */
	     FitToOverscan (x, sizey, trimy1, biassect, ccdbias, wf3->sdqflags,
			    rn);

	     /* Report the slope and intercept from the BLEV fit here */
	     BlevResults (&ampslope, &ampintercept);

	     /* Optional results for possible future use ...
	     bslope[amp_indx] = ampslope;
	     bzero[amp_indx]  = ampintercept;
	     */

	     for (j = begy;  j < endy;  j++) {
		  biaslevel = BlevEval (dj);
		  dj = dj + 1.0;

		  /* bias for current line plus average drift (constant) */
		  sumbias += (biaslevel + averagedrift);

		  di = (double)begx;
		  for (i = begx;  i < endx;  i++) {
		       if (dodrift == YES) {
			   deval = DriftEval(di);
		       } else {
			   deval = 0;
		       }

		       di = di + 1.0;
		       Pix (x->sci.data,i,j) = Pix (x->sci.data,i,j) -
					       biaslevel - deval;
		  }
	     }

	     wf3->blev[amp_indx] = sumbias / (double)sizey;
	     sumblev += (sumbias/(double)sizey);

	} /* End loop over AMPs */

	*overscan = YES;
	
	/* This is the mean value of all the bias levels subtracted. */
	*meanblev = sumblev / numamps;
	
	/* free ccdamp space allocated here... */
	free (ccdamp);

	return (status);
}

/* This function determines the bias level from the overscan in the input
   image x, and it fits a straight line to these values as a function of
   line number in the input image.
*/

static void FitToOverscan (SingleGroup *x, int ny, int trimy1, int *biassect,
			   float ccdbias, short sdqflags, float rn) {

/* arguments:
SingleGroup *x    i: input image
int ny            i: size of second axis after overscan has been removed
int trimy1        i: offset between input and output line numbers
int biassect[4]   i: section to use for bias level determination
float ccdbias     i: bias level to subtract if we can't get it from overscan
short sdqflags    i: dq flag value representing a bad pixel
floar rn          i: calibrated readnoise level for this amp
*/

	extern int status;

	double biaslevel;	/* bias level in one line of input image */
	int too_few = 0;	/* number of lines with too few good pixels */
	int npix;		/* number of pixels used to compute bias */
	double *biasvals;	/* intermediate array for biaslevel values */
	int *biasmask;		/* mask array for biasvals: 0 means don't use */
	int j;
	void BlevInit (int);
	void BlevAccum (int, double);
	int BlevFit (void);
	void BlevSet (double);
	int FindBlev (SingleGroup *, int, int *, short, double *, int *);
	void cleanBiasFit (double *, int *, int, float);

	BlevInit (ny / 2);			/* initialize for fitting */

	/* Allocate space for biaslevel array */
	biasvals = (double *) calloc (ny+1, sizeof(double));
	biasmask = (int *) calloc (ny+1, sizeof(int));

	/* For each line, determine the bias level from the overscan in x.
	   Note that we loop over the number of pixels in the output image.
	   The argument j+trimy1 to FindBlev is the line number in the
	   input image x, with trimy1 the offset to the illuminated region.
	*/
	for (j = 0;  j < ny;  j++) {
	    if (FindBlev (x, j+trimy1, biassect, sdqflags, &biaslevel, &npix)) {
		too_few++;
		status = 0;			/* not fatal */
	    } else {
		/* add biaslevel value to intermediate array;
		** this array will be analyzed later for outliers. */
		biasvals[j] = biaslevel;
		biasmask[j] = 1;
	    }
	}

	/* Analyze biasvals for outliers and mask them by setting the
	** corresponding value in the biasmask array to zero. */
	cleanBiasFit (biasvals, biasmask, ny, rn);

	/* Now that we have read in all bias level values from the overscan
	** regions and thrown out any outliers (cosmic-ray hits), we can
	** accumulate the fitting statistics now. */
	for (j=0; j < ny; j++) {
	     if (biasmask[j] == 1) {
		 /* Note: j is line number in output image. */
		 BlevAccum (j, biasvals[j]);   /* increment sums */
	     }
	}

	if (too_few > 0) {
	    sprintf (MsgText, "(blevcorr) %d image line", too_few);
	    if (too_few == 1)
		strcat (MsgText, " has");
	    else
		strcat (MsgText, "s have");
	    strcat (MsgText, " too few usable overscan pixels.");
	    trlwarn (MsgText);
	}

	/* Fit a curve to the bias levels found. */
	if (BlevFit()) {
	    trlwarn ("No bias level data, or singular fit; ");
	    trlmessage ("            bias from CCDTAB will be subtracted.");
	    BlevSet (ccdbias);		/* assign the default value */
	}

	/* Free up local memory */
	free (biasvals);
	free (biasmask);

}

/* This routine uses iterative sigma clipping to reject outliers from
   the array of bias values.
*/

void cleanBiasFit (double *barray, int *bmask, int ny, float rn) {

	int j;
	double bsum, bmean, sdev, svar;
	double s;
	int nsum;
	float clip;
	int nrej=0;

	bsum = bmean = sdev = svar = 0.0;
	nsum = 0;

	for (j=0; j < ny; j++) {
	     if (bmask[j] == 1) {
		 bsum += barray[j];
		 nsum++;
	     }
	}
	if (nsum == 0) return;

	bmean = bsum / nsum;

	for (j=0; j < ny; j++) {
	     if (bmask[j] == 1) {
		 s = barray[j] - bmean;
		 svar += (s*s);
	     }
	}
	sdev = sqrt (svar/nsum);

	/* Reset stddev from mean to stddev of poisson dist centered
	** on mean. This will keep cosmic-ray hits or bleeding from
	** bright sources into the overscan columns from biasing
	** the fit. HAB 19-Feb-2004
	*/
	if (sdev > sqrt(bmean)) sdev = sqrt(bmean);
	clip = 3.5;

	/* With statistics in hand, ID and flag outliers */
	for (j=0; j < ny; j++) {
	     if (barray[j] > clip*sdev+bmean) {
		 bmask[j] = 0;
		 nrej++;
	     }
	}

	/* Recompute the mean based on clipped values. */
	bsum = bmean = 0.0;
	nsum = 0;
	for (j=0; j < ny; j++) {
	     /* if value has not already been thrown out, use it
	     ** to compute new mean ... */
	     if (bmask[j] != 0) {
		 bsum += barray[j];
		 nsum++;
	     }
	}
	if (nsum == 0) return;
	bmean = bsum / nsum;
	clip = 2.0;

	/* With statistics in hand, ID and flag outliers based on
	** readnoise as sigma to further refine the value ... */
	for (j=0; j < ny; j++) {
	     if (barray[j] > clip*rn+bmean) {
		 if (bmask[j] != 0) {
		     bmask[j] = 0;
		     nrej++;
		 }
	     }
	}

	sprintf (MsgText,
		"(blevcorr) Rejected %d bias values from serial fit.",nrej);
	trlmessage (MsgText);
}

int selectBias (char *ccdamp) {

	int i;
	char *ampstr = "ABCD";
	char *loc;

	/* Search for CCDAMP */
	loc = strchr (ampstr, ccdamp[0]);

	/* Find out what array position corresponds to that amp */
	i = loc - ampstr;

	return i;
}

