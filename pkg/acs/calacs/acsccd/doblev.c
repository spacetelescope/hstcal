/* This file contains:
 doBlev
 FitToOverscan
 */

# include <stdio.h>
# include <math.h>		/* sqrt */
# include <string.h>

# include "hstio.h"

# include "acs.h"
# include "acsinfo.h"
# include "hstcalerr.h"
# include "acsdq.h"		/* for CALIBDEFECT */

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

 2-June-1998 WJH: Initial Version. Revised for ACS data...

 29-Oct-2001 WJH: Revised to use CCDOFST[A,B,C,D] to select default
 bias levels from CCDBIAS[A,B,C,D] (which are new to CCDTAB) in
 this version. Updated to use VY values of <=0 to detect lack of
 virtual overscan use and NOT apply BlevDrift().
 27-Feb-2002 WJH/RB: Fixed a bug in evaluating virtual overscan drift
 where the column id was being incremented twice.
 */

int doBlev (ACSInfo *acs, SingleGroup *x, int chip, float *meanblev,
            int *overscan, int *driftcorr) {

  /* arguments:
   ACSInfo *acs    i: calibration switches, etc
   SingleGroup *x    io: image to be calibrated
   int chip        i: chip number
   float meanblev   o: mean value of bias levels that were subtracted
   int *overscan     io: true if bias determined from an overscan region in x
   int *driftcorr    o: true means correction for drift along lines was applied
   */

  extern int status;

  double biaslevel;           /* value to subtract from a line of image */
  double averagedrift;        /* drift averaged over all (output) columns */
  double sumbias, sumblev;    /* sum of bias levels (for getting mean) */
  int binx, biny;             /* bin factors */
  int trimx1, trimx2;	/* width to trim off ends of each line */
  int trimy1, trimy2;	/* amount to trim off ends of each column */
  int biassect[4];	/* section(s) to use for finding bias level */
  int i, j;
  double di, dj;
  int endx, endy;		/* boundaries of actual image */
  int begx, begy;		/* boundaries of actual image */
  int sizex, sizey;
  int biasnum;        /* array position for single amp used */
  int amp;			/* Counter for amps used */
  int numamps;		/* Number of AMPS used to readout chip */
  char *ccdamp; /* Amps which are used for this chip */
  int dodrift = YES;		/* use virtual overscan to get blevdrift? */
  char *amploc;
  int bias_loc, amp_indx;
  int bias_ampx, bias_ampy;
  int bias_orderx[4] = {0,1,0,1};
  int bias_ordery[4] = {0,0,1,1};
  double deval;
  double rn;
  double ampslope, ampintercept;  /* values from the bias level fit*/
  float ccdbias;      /* default bias level from CCDTAB */

  /* Function definitions */
  int BlevDrift (SingleGroup *, int *, int *, int, int *, int *, short);
  double DriftEval (double);
  double DriftMean (double);
  double BlevEval (double);
  void BlevResults(double *, double *);
  int FindBlev (SingleGroup *, int, int *, short, double *, int *);
  void parseWFCamps (char *, int, char *);
  int selectBias(char *);

  /* Allocate space for ccdamp... */
  ccdamp = (char *) calloc (NAMPS+1, sizeof(char));

  /* initial values */
  *driftcorr = 0;
  *meanblev = 0.;

  binx = acs->bin[0];
  biny = acs->bin[1];

  if ((binx != 1 && binx != 2 && binx != 4) ||
	    (biny != 1 && biny != 2 && biny != 4)) {
    trlerror ("(doBlev) bin size must be 1, 2, or 4.");
    free (ccdamp);
    return (status = 1001);
  }

  /* Get biassect, the location of the region to use for determining
   the bias level, and vx & vy, the location in the virtual
   overscan to use for determining the drift with column number.
   Also get the widths of the trim regions, the amount to remove
   from each axis of the input when copying to output.
   ** This routine is new to ACS...  It reads in overscan and trim
   **	region specifications from an OSCNTAB, rather than having
   **	it hardwired into the code as it was in CALSTIS.
   */
  /*
   ASSUMPTION: Only have no overscan with SINGLE AMP readouts
   Only use FIRST AMP listed in CCDAMP to determine which
   CCDOFST[A,B,C,D] value to use.
   */
  biasnum = selectBias(acs->ccdamp);
  if (strlen(acs->ccdamp) > 1){
    if (acs->detector == WFC_CCD_DETECTOR && acs->chip == 2) biasnum += 2;
  }
  ccdbias = acs->ccdbias[biasnum] * acs->atodgain[biasnum];

  if (*overscan != YES) {
    /* If no overscan region, subtract default BIAS level
     obtained from CCDTAB. */
    trlwarn ("Overscan region is too small to do BLEVCORR; ");
    sprintf (MsgText,"     bias from CCDTAB of %g will be subtracted.",ccdbias);
    trlmessage(MsgText);

    for (j = 0;  j < x->sci.data.ny;  j++) {
      for (i = 0;  i < x->sci.data.nx;  i++)
		    Pix (x->sci.data,i,j) = Pix (x->sci.data,i,j) - ccdbias;
    }

    *meanblev = ccdbias;
    acs->blev[biasnum] = ccdbias;
    free (ccdamp);
    return (status);
  }

  /* Determine bias levels from BIASSECT columns specified in OSCNTAB */

  /* Copy out overscan info for ease of reference in this function*/
  trimx1 = acs->trimx[0];
  trimx2 = acs->trimx[1];
  trimy1 = acs->trimy[0];
  trimy2 = acs->trimy[1];

  /* Establish which amps from ccdamp string are appropriate
   for this chip */
  ccdamp[0] = '\0';

  if (acs->detector == WFC_CCD_DETECTOR) {
    /* Set up the 2 AMPS used per line */
    parseWFCamps (acs->ccdamp, chip, ccdamp);

  } else {
    /* HRC observation; use full CCDAMP string */
    strcpy (ccdamp, acs->ccdamp);
  }

  /* How many amps are used for this chip */
  numamps = strlen (ccdamp);


  /* Are we going to calculate drift in bias from virtual overscan? */
  /* If the end points of vx and vy are zero, no section was specified */
  if (acs->vy[0] <= 0 && acs->vy[1] <= 0 ){
    dodrift = NO;
    trlmessage ("(blevcorr) No virtual overscan region specified.");
    trlmessage ("(blevcorr) Bias drift correction will not be applied.");
  }

  sumblev = 0.;

  /* For each amp used, determine bias level and subtract from
   appropriate region of chip */
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
    if (acs->biassecta[1] > 0 && acs->biassectb[1] > 0) {
      /* select section nearest the amp based on bias_amp */
      biassect[0] = (bias_ampx == 0) ? acs->biassecta[0] : acs->biassectb[0];
      biassect[1] = (bias_ampx == 0) ? acs->biassecta[1] : acs->biassectb[1];
      /* If only 1 AMP is used for each line, then use both sections*/
      if (numamps == 1) {
        biassect[2] = (bias_ampx == 0) ? acs->biassectb[0] : acs->biassecta[0];
        biassect[3] = (bias_ampx == 0) ? acs->biassectb[1] : acs->biassecta[1];
      } else {
        biassect[2] = 0;
        biassect[3] = 0;
      }
    } else {
      /* otherwise, select the non-zero bias section */
      biassect[0] = (acs->biassecta[1] == 0) ? acs->biassectb[0] : acs->biassecta[0];
      biassect[1] = (acs->biassecta[1] == 0) ? acs->biassectb[1] : acs->biassecta[1];
      /* then set trailing section to zero, indicating it is not used.*/
      biassect[2] = 0;
      biassect[3] = 0;
    }

    /* make sure that the biassect doesn't go beyond the edge of the chip */
    if (biassect[0] < 0)
      biassect[0] = 0;
    if (biassect[2] < 0)
      biassect[2] = 0;
    if (biassect[1] >= x->sci.data.nx)
      biassect[1] = x->sci.data.nx - 1;
    if (biassect[3] >= x->sci.data.nx)
      biassect[3] = x->sci.data.nx - 1;

    /* Compute range of pixels affected by each amp */
    begx = trimx1 + acs->ampx * bias_ampx;
    endx = (bias_ampx == 0 && acs->ampx != 0) ? acs->ampx + trimx1 : x->sci.data.nx - trimx2;
    begy = trimy1 + acs->ampy* bias_ampy;
    endy = (bias_ampy == 0 && acs->ampy != 0) ? acs->ampy +trimy1 : x->sci.data.ny - trimy2;

    /* Make sure that endx and endy do not extend beyond the bounds of the
     image... WJH 8 Sept 2000
     */
    if (endx > x->sci.data.nx) endx = x->sci.data.nx;
    if (endy > x->sci.data.ny) endy = x->sci.data.ny;
    sizex = endx - begx;
    sizey = endy - begy;

    /* At this point, we decide how we are to determine the bias level...
     average two sections across the line for 1 amp, or use only
     one section...
     */

    if (dodrift == YES) {
      /* Fit a line to the virtual overscan region as a function of
       column number.
       */
      if (BlevDrift (x, acs->vx, acs->vy, trimx1, biassect, driftcorr,
                     acs->sdqflags)) {
        free (ccdamp);
        return (status);
      }

      /* Evaluate the fit for each line, and subtract from the data. */
      averagedrift = DriftMean ((double)sizex);

    } else {
      averagedrift = 0.;
    }
    dj = (double)begy;
    sumbias = 0.;

    /* For each image line, determine the bias level from the overscan
     in x, and fit to these values as a function of line number in
     the output image.
     The readnoise which gets used here needs to be specific to the amp
     used for the readout region. (H. Bushouse, 28 July 2006).
     */
    rn = acs->readnoise[amp_indx];
    FitToOverscan (x, sizey, trimy1, biassect, ccdbias, acs->sdqflags, rn);

    /* Report the slope and intercept from the BLEV fit here. */
    BlevResults (&ampslope, &ampintercept);

    /* Optional results for later use...
     bslope[amp_indx] = ampslope;
     bzero[amp_indx] = ampintercept;
     */
    for (j = begy;  j < endy;  j++) {
      biaslevel = BlevEval (dj);

      dj = dj + 1.0;
      /* bias for current line plus average drift (constant) */
      sumbias += (biaslevel + averagedrift);

      di = (double)begx;
      for (i = begx;  i < endx;  i++){
        if (dodrift == YES) {
          deval = DriftEval(di);
        } else {
          deval = 0;
        }
        di = di + 1.0;
        Pix (x->sci.data,i,j) = Pix (x->sci.data,i,j) - biaslevel - deval;
      }
    }
    acs->blev[amp_indx] = sumbias / (double)sizey;
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
 line number in the output image.
 */

static void FitToOverscan (SingleGroup *x, int ny, int trimy1,
                           int *biassect, float ccdbias, short sdqflags, float rn) {

  /* arguments:
   SingleGroup *x    i: input image
   int ny            i: size of second axis after overscan has been removed
   int trimy1        i: offset between input and output line numbers
   int biassect[2]   i: section to use for bias level determination
   float ccdbias  i: bias level to subtract if we can't get it from overscan
   float rn         i: calibrated readnoise level for this amp
   */

	extern int status;

	double biaslevel;	/* bias level in one line of input image */
	int too_few = 0;	/* number of lines with too few good pixels */
	int npix;		/* number of pixels used to compute bias */
  double *biasvals;  /* intermediate array for biaslevel values */
  int *biasmask;     /* mask array for biasvals- 0 means do not use */
	int j;
	void BlevInit (int);
	void BlevAccum (int, double);
	int BlevFit (void);
	void BlevSet (double);
	int FindBlev (SingleGroup *, int, int *, short, double *, int *);
  void cleanBiasFit(double *, int *, int, float);

	BlevInit (ny / 2);			/* initialize for fitting */

  /* Allocate space for biaslevel arrays... */
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
      /* add biaslevel value to intermediate array
       This array will then be analyzed for outliers and
       with those outliers being removed.
       */
      biasvals[j] = biaslevel;
      biasmask[j] = 1;
      if (j == 0) {
          sprintf (MsgText,
              "(FitToOverscan) biassecta=(%d,%d) biassectb=(%d,%d) npix=%d",
              biassect[0], biassect[1], biassect[2], biassect[3], npix);
          trlmessage(MsgText);
      }
    }
	}
  /* Analyze biasvals for outliers and mask them
   by setting the corresponding value in the biasweight
   array to 0.
   */
  cleanBiasFit(biasvals, biasmask, ny, rn);

  /*
   Now that we have read in all the bias level values from the overscan
   regions and thrown out any outliers (cosmic-ray hits), we can accumulate
   the fitting statistics now.
   */
	for (j = 0;  j < ny;  j++) {
    if (biasmask[j] == 1){
      /* Note:  j is line number in output image. */
      BlevAccum (j, biasvals[j]);	/* increment sums */
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

  /* Clean up memory */
  free(biasvals);
  free(biasmask);
}

void cleanBiasFit(double *barray, int *bmask, int ny, float rn){
  int j;
  double bsum, bmean, sdev, svar;
  double s;
  int nsum;
  float clip = 3.0;
  int nrej=0;

  bsum = bmean = sdev = svar = 0.0;
  nsum = 0;
	for (j = 0;  j < ny;  j++) {
    if (bmask[j] == 1){
      bsum += barray[j];
      nsum += 1;
    }
  }

  if (nsum == 0)
    return;

  bmean = bsum / nsum;
	for (j = 0;  j < ny;  j++) {
    if (bmask[j] == 1) {
      s = barray[j] - bmean;
      svar += (s*s);
    }
  }
  sdev = sqrt(svar/(nsum-1));
  /*
   Reset stddev from mean to stddev of poisson dist
   centered on mean.  This will keep cosmic-ray hits or
   bleeding from bright sources into the overscan columns
   from biassing the fit.  WJH  5-Nov-2003
   */
  if (sdev > sqrt(bmean)) sdev = sqrt(bmean);

  /* With statistics in hand, ID and flag outliers*/
	for (j = 0;  j < ny;  j++) {
    if (barray[j] > abs((clip*sdev)+bmean)) {
      bmask[j] = 0;
      nrej++;
    }
  }

  /* Recompute the mean based on clipped values. */
  bsum = bmean = 0.0;
  nsum = 0;
	for (j = 0;  j < ny;  j++) {
    /* If value has not already been thrown out,
     use it to compute new mean... */
    if ( bmask[j] != 0 ){
      bsum += barray[j];
      nsum ++;
    }
  }
  bmean = bsum / nsum;

  /* With statistics in hand, ID and flag outliers based on
   readnoise as sigma to further refine the value... */
	for (j = 0;  j < ny;  j++) {
    if (barray[j] > abs((clip*rn)+bmean)) {
      bmask[j] = 0;
      nrej++;
    }
  }

  sprintf(MsgText,"(blevcorr) Rejected %d bias values from fit.",nrej);
  trlmessage(MsgText);
}


int selectBias (char *ccdamp) {
  int i;
  char *ampstr = "ABCD";
  char *loc;

  /* Search for CCDAMP */
  loc = strchr(ampstr, ccdamp[0]);
  /* Find out what array position corresponds to that amp */
  i = loc - ampstr;

  return i;
}
