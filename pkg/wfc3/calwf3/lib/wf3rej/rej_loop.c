# include   <stdio.h>
# include   <string.h>
# include   <stdlib.h>
# include   <math.h>

#include "hstcal.h"
# include "hstio.h"

# include   "wf3.h"
# include   "wf3rej.h"
# include   "rej.h"
# include   "hstcalerr.h"
# include   "wf3info.h"

/* local mask values */
# define    OK          (short)0
# define    HIT         (short)8
# define    SPILL       (short)2
# define    EXCLUDE     (short)4
# define    BAD         (short)~OK
# define    ABS(x)      ((x>0.) ? x : -x)

/* Local routines for dealing with arrays of pointers to
**  arrays.  */

static float **allocFloatBuff (int, int);
static short **allocShortBuff (int, int);
static void freeFloatBuff (float **, int);
static void freeShortBuff (short **, int);
static void scrollFloatBuff (float *, int, int, int, int, float **, float *);
static void scrollShortBuff (short *, int, int, int, int, short **, short *);
static void InitShortSect (short **, short *, IODescPtr *, int, int, int);
static void InitFloatSect (float **, float *, IODescPtr *, int, int, int);

static Byte ***allocBitBuff (int, int, int);
static void freeBitBuff (Byte ***, int, int);
static void readBitLine (Byte ***, int, int, int, short, short, short *);
static void writeBitLine (short *, int, int, int, short, Byte ***); 

/* For debugging purposes:
static void printBitLine (Byte ***, int, int, int);
*/
static int  initShad (Hdr *, int, char *, int, int *, int *, int *, int *,
		      int *);
static void getShadLine (float **, int, int, int, float *);
static void getShadBuff (IODescPtr *, int, int, int, int, int , int , int ,
			 int , float **);
static void getShadcorr (float **, int, int, int, float, int, float *);

/*  crrej_loop -- Perform cosmic ray rejection

Description:
------------
This routine performs many operations involving multiple input images
and scrolling buffers.  Cosmic rays are detected through iteration using
different sigma sizes on each iteration with neighboring pixels being
rejected as well.  

Code Outline:
-------------
    - Initialize all data structures, 
         malloc'ing space for scrolling buffers.
         initialize shadcorr buffers 
         determine expansion factors, if applying a shadfile.
         setup variables based on input values
    - Determine gain and amp values used in image
    - Start the rejection iteration
        reset single line buffers to all zero
        reset detection mask
        
        - Loop over each line in the image
            Setup gain and noise values for this line
            Get shading correction for this line from the reference image
            Initialize or scroll sliding buffers
                - make sure all data in pic buffer is normalized the same

            - Loop over each input image to be combined:
                Calculate the threshold for each pixel above which a cosmic
		ray would stand out
                    - Use either 'minimum' or 'median' method here
                Mark SPILL pixels as EXCLUDE
                For each pixel in the line, compare pixel value with thresshold
                    - for each detected pixel, mark as HIT
                    - mark surrounding pixels out to RADIUS as SPILL if above a
                    more stringent threshhold value
                Accumulate the counts for each good pixel in this line from
                    this image into the output image buffer.
                On the last iteration:
                    - calculate the variance(error) for each pixel
                    - output the buffer for the DQ array complete with
                        which pixels were masked as HITs from this image.
            - For this combined line, 
                perform shading correction
                calculate the new average value for each pixel based on 
                    how many images contributed to the output value.
     
    - Write out CR-hit information to all the input images if par->mask was set
    
    - Free memory used by all the buffers                    
            
                           
  
  Date          Author          Description
  ----          ------          -----------
  26-Apr-1996   J.-C. Hsu       adapt from the SPP code crrej_loop.x
  28-Aug-1998   W.J. Hack       Modified for section-by-section detection
  20-Nov-1998   W.J. Hack       revised to handle trailer file comments
                                and exit more gracefully (see goto statements)
  13-Sep-1999   W.J. Hack       Heavily revised to more directly incorporate 
                                shadcorr calculation (rather than using if()),
                                and only use a single set of noise and gain
				values.
                                Other revisions were also made for speed.
  18-Oct-1999   W.J. Hack       Corrected buffer management, SPILL-pixel 
                                radius computation, dqpat usage.  Buffer now
                                contains all normalized values.
   8-Dec-1999   W.J. Hack       Added explanatory comments and added buffer
                                for shading correction in addition to
                                shading reference data from shadfile.
                                These corrections were based on re-review.
  14-Apr-2000   W.J. Hack       Revised to work with sub-arrays and 1-AMP
				readouts
  29-Aug-2000	H.A. Bushouse	Revised for WFC3 use.
  16-Nov-2001   H.A. Bushouse	Made allocation of zl in getShadBuff more robust
  18-Jun-2002   H. Bushouse	Added a check to not include bad pixels in
				combined error image (following CALACS changes).
  20-Feb-2004   H. Bushouse	Removed use of SCALENSE in computing errors.
				Fixed memory allocation bug in AllocBitBuff to
				handle arbitrary buffer lengths, not just those
				divisible by 8 (following CALACS changes).
  06-Dec-2007   H. Bushouse     Added calls to getHeader before each call to
                                getShortLine to prevent getShortLine from
                                crashing on null DQ arrays.
  28-Jan-2009   H. Bushouse	Added use of dq2 internal array to hold DQ
				values intended to be written back into blv_tmp
				files so that they're separate from the values
				in the dq array, which get written to the
				output crj file. This allows the crj dq value
				to not contain CR or SPILL values if at least
				1 good pixel exists. Also fixed the error
				inherited from acsrej that was propagating CR
				DQ values from 1 blv_tmp file to all following
				blv_tmp files.
  18-Feb-2009	H. Bushouse	Fixed bug in test to exclude flagged pixels from
				being tested for CR hits so that pixels marked
				previously as SPILL still get tested for CR.
  14-Dec-2011	H. Bushouse	Upgraded to rescale input data that are in
				units of count rates. (PR 69969; Trac #814)
*/

int rej_loop (IODescPtr ipsci[], IODescPtr ipdq[], char imgname[][CHAR_FNAME_LENGTH+1],
	      int grp [], int nimgs, clpar *par, int niter, int dim_x,
	      int dim_y, float sigma[], multiamp noise, multiamp gain, 
	      float efac[], float skyval[], DataUnits bunit[],
	      FloatTwoDArray *ave, FloatTwoDArray *avevar, float *efacsum, 
              ShortTwoDArray *dq, int *nrej, char *shadfile)
{
    extern int status;

    Hdr     dqhdr;              /* data quality header structure */
    int     width;
    int     i, j, k, n, jndx;   /* loop indices */
    int     iter;
    int     ii, jj, j2;
    float   sig2, psig2, rej2;  /* square of something */
    float   *exp2;              /* square of exptime per image*/
    float   scale, val, dum, pixsky;
    short   sval, crflag, nocr, nospill, dqpat;
    short   maskdq;

    float   efacn, skyvaln, exp2n;

    ShortTwoDArray dq2;		/* local array for storing output blv DQ vals */

    /* local data arrays */
    float   ***pic;
    float   ***thresh, ***spthresh;
    short   ***mask;

    float   *sum;
    float   *sumvar;
    float   *buf;
    short   *bufdq;
    
    Byte    ***crmask;          /* Compressed CR HIT mask for all images */

    /* local variables for sections */
    int     buffheight, line, bufftop;
    IODescPtr	ipdqn;

    /* Noise and Gain value for a pixel */
    float       rog2[NAMPS];
    float       noise2[NAMPS], gain2[NAMPS];
    float       gn[2], nse[2];
    int         detector, chip, ampy, ampx;
    int         numpix;

    /* Parameters for applying SHADCORR */
    IODescPtr   ipshad;
    float       **shadbuff;
    float       *shadline, *shadcorr;
    int         shad_dimy;
    Hdr         scihdr;
    int         rx, ry, x0, y0;
    float       pixexp;
    int         shadf_x;        /* Number of x pixels in SHADFILE */
    float       *zerofbuf;       /* line buffer of all FLOAT zeroes  */
    short       *zerosbuf;      /* line buffer of all SHORT zeroes */

    /* Functions for dealing with MULTIAMP values of gain and noise */
/*  int       LoadHdr (char *, Hdr *);*/
    void      WhichError (int);
    void      get_nsegn (int, int, int, int, float *, float*, float *, float *);
    void      TimeStamp (char *, char *);

    /********************************** Begin Code ****************************/
    /* Initialization */
    ipshad=NULL;
    shadf_x=0;
    crflag = par->crval;
    dqpat = par->badinpdq;
    scale = par->scalense/100.;
    nocr = ~crflag;
    nospill = ~SPILL;
    numpix = dim_x * dim_y;
    
    /* Set up mask for detecting CR-affected pixels */
    maskdq = OK | EXCLUDE;
    maskdq = maskdq | HIT;
    maskdq = maskdq | SPILL;
    
    /* Define the buffer size for scrolling  
        up the image.
    */
    width = (int) ceil(par->radius);
    buffheight = 1+2*width;
    
    /* allocate data arrays */
    pic = (float ***) calloc(nimgs, sizeof(float **));
    if (pic == NULL) {
        trlerror ("Couldn't allocate memory for scratch array in REJ_LOOP.");
        return (status = OUT_OF_MEMORY);
    }
    thresh = (float ***) calloc (nimgs, sizeof(float **));
    if (thresh == NULL) {
        trlerror ("Couldn't allocate memory for scratch array in REJ_LOOP.");
        return (status = OUT_OF_MEMORY);
    }
    spthresh = (float ***) calloc (nimgs, sizeof(float **));
    if (spthresh == NULL) {
        trlerror ("Couldn't allocate memory for scratch array in REJ_LOOP.");
        return (status = OUT_OF_MEMORY);
    }
    mask = (short ***) calloc (nimgs, sizeof(short **));
    if (mask == NULL) {
        trlerror ("Couldn't allocate memory for scratch array in REJ_LOOP.");
        return (status = OUT_OF_MEMORY);
    }
    exp2 = (float *) calloc (nimgs, sizeof(float));
    if (exp2 == NULL) {
        trlerror ("Couldn't allocate memory for scratch array in REJ_LOOP.");
        return (status = OUT_OF_MEMORY);
    }

    for (k = 0; k < nimgs; k++) {
        mask[k] = allocShortBuff (buffheight, dim_x);
        pic[k] = allocFloatBuff(buffheight, dim_x);
        thresh[k] = allocFloatBuff (buffheight, dim_x);
        spthresh[k] = allocFloatBuff (buffheight, dim_x);
    }

    zerofbuf = calloc (dim_x, sizeof(float));
    zerosbuf = calloc (dim_x, sizeof(short));

    sum = calloc(dim_x, sizeof(float));
    sumvar = calloc (dim_x, sizeof(float));

    buf = calloc (dim_x, sizeof(float));
    bufdq = calloc (dim_x, sizeof(short));

    /* 
        If we want to perform SHADCORR here, ...
    */
    if (par->shadcorr == PERFORM) {


        /* Read in the SCI header here... */
        initHdr (&scihdr);	
	    if (getHeader (ipsci[0], &scihdr) )
		    status = HEADER_PROBLEM;	
	    if (hstio_err() || status) {
		    freeHdr (&scihdr);
	        return (status = OPEN_FAILED);
	    }

        /* Read in SHADFILE header
        if (LoadHdr (shadfile, &scihdr))
            return (status = HEADER_PROBLEM);
        */
        
        /* Determine expansion factors dimensions, and offsets for SHADFILE */
        if (initShad (&scihdr, dim_x, shadfile, grp[0], &shadf_x, &rx, &ry,
                      &x0, &y0) )
            WhichError (status);
        ipshad = openInputImage (shadfile, "SCI", grp[0]);

    } else {

        rx = 1;
        ry = 1;
        x0 = 0;
        y0 = 0;    
        /* Use this value to flag that no shading correction 
            will be applied, since this refers to the number of 
            pixels in each line of the shading correction file...
        */
        shadf_x = 0;
    }

    /* Buffer for use with SHADCORR */
    shad_dimy = SECTLINES * ry;
    shadbuff = allocFloatBuff (shad_dimy, dim_x);
    
    /* We always want this to be defined, 
        with it all ZERO when NOT performing SHADCORR.
    */
    shadline = calloc (dim_x, sizeof(float));
    
    /* This buffer is used for the scaled shadcorr data.
        It defaults to ONE.    
    */
    shadcorr = calloc (dim_x, sizeof(float));

    /* Allocate space for the CR-hit mask */	
    crmask = allocBitBuff (nimgs, dim_y, dim_x);

    /* readout is in DN */
    rej2 = SQ(par->radius);
    psig2 = SQ(par->thresh);
    if (par->thresh <= 0.)
        psig2 = -1.;

    /* initialize the output DQ arrays: dq2 will hold values to be written
       back into the input DQ arrays at the end of processing, while dq
       will hold values to write to the output CRJ DQ extension */
    initShortData (&dq2);
    allocShortData (&dq2, dim_x, dim_y, True);
    for (j = 0; j < dim_y; j++) {
        for (i = 0; i < dim_x; i++) {
            PDQSetPix(dq,i,j,crflag);
            PDQSetPix(&dq2,i,j,OK);
	}
    } 

    /* All Observations will have the same CCDAMP */
    ampx = gain.colx;
    ampy = gain.coly;
    detector = gain.detector;
    chip = gain.chip;

    /* Assumption: ALL images have the same noise/gain values */
    for (k = 0; k < NAMPS; k++){
        gain2[k] = 0.;
        noise2[k] = 0.;
        rog2[k] = SQ(noise.val[k]);
    }

    for (n = 0; n < nimgs; n++) {
        exp2[n] = SQ(efac[n]);
    }

    /* Set up gain and values used for each image */
    get_nsegn (detector, chip, ampx, ampy, gain.val, rog2, gain2, noise2); 

    /* start the rejection iteration */
    for (iter = 0; iter < niter; iter++) {
        if (par->verbose) { 
            sprintf (MsgText, "iteration %d", iter+1);
            trlmessage (MsgText);
        }
        
        sig2 = SQ(sigma[iter]);

	if (iter > 0) {
            /* Re-initialize the arrays... */
            for (j = 0; j < numpix; j++) {
                *(efacsum+j) = 0.;
            }
             
            memcpy (sum, zerofbuf, dim_x * sizeof(float));
            memcpy (sumvar, zerofbuf, dim_x * sizeof(float));
            memcpy (buf, zerofbuf, dim_x * sizeof(float));
            memcpy (bufdq, zerosbuf, dim_x * sizeof(short));

            *nrej = 0;

            /* ...and reset the mask  */
            for (k = 0; k < nimgs; k++) {
                for (j = 0; j < buffheight; j++) {
                    memcpy (mask[k][j], zerosbuf, dim_x * sizeof(short));
                }
            }
        } /* End initialization section */
        
        /* Start loop over lines in image */
        for (line =0; line < dim_y; line++) { 
            
            /* Zero out this buffer for the next line */
            memcpy (sum, zerofbuf, dim_x * sizeof(float));
            memcpy (sumvar, zerofbuf, dim_x * sizeof(float));

            /* Set up the gain and noise values used for this line in
	    ** ALL images */
            if (line < ampy ) {
                gn[0] = gain2[AMP_C];
                gn[1] = gain2[AMP_D];
                nse[0] = noise2[AMP_C];
                nse[1] = noise2[AMP_D];
            } else {
                gn[0] = gain2[AMP_A];
                gn[1] = gain2[AMP_B];
                nse[0] = noise2[AMP_A];
                nse[1] = noise2[AMP_B];            
            }

            /* If we are doing SHADCORR, then fill buffer and get
                a single SHADLINE to be applied. 
            
             Get (binned?) line from shadfile reference image
                and put into 'shadbuff'.  However, we only need
                to get a line every time we reach the end of
                the buffer. 
                If no shading correction is being performed, it
                will return a buffer of all ZEROES.
            */
            if ((line % shad_dimy) == 0 || line == 0) {
                getShadBuff (ipshad, line, shad_dimy, dim_x, shadf_x, rx,
                ry, x0, y0, shadbuff);
            }
           
           /* Manage scrolling buffers for each line here... */
            if (line > 0) {								
                /* Scroll buffers so new line can be inserted into middle row
		 of buffer add BLANK line to bottom of buffer 
                */
                bufftop = line+width;
                for (k = 0; k < nimgs; k++) {
                    
                    if (bufftop < dim_y) {
                        initHdr(&dqhdr);
                        getHeader(ipdq[k],&dqhdr);
                        getFloatLine (ipsci[k], bufftop, buf);
                        getShortLine (ipdq[k], bufftop, bufdq);
                        freeHdr(&dqhdr);

			/* Rescale the inputs by exposure time, if needed */
			if (bunit[k] == COUNTRATE) {
			    for (i = 0; i < dim_x; i++) {
				 buf[i] *= efac[k];
			    }
			}

                         /* Scale the input values by the sky and exposure time
                            for comparison to the detection threshhold.
                        */ 
                        for (i = 0; i < dim_x; i++){
			    if (efac[k] > 0.) {
                                buf[i] = (buf[i] - skyval[k]) / efac[k];
			    } else {
				buf[i] = 0.;
			    }
                        }   
                    } else {                       
                        memcpy (buf, zerofbuf, dim_x * sizeof(float));
                        memcpy (bufdq, zerosbuf, dim_x * sizeof(short));
                    }   
                    
                    scrollShortBuff (bufdq, line, dim_y, buffheight, dim_x,
				     mask[k], zerosbuf);
                    scrollFloatBuff (buf, line, dim_y, buffheight, dim_x,
				     pic[k], zerofbuf);                   
                    scrollFloatBuff (buf, line, dim_y, buffheight, dim_x,
				     thresh[k], zerofbuf);
                    scrollFloatBuff (buf, line, dim_y, buffheight, dim_x,
				     spthresh[k], zerofbuf);
                }
            } else {
                /* Working with first line, so we need to initialize the
                    scrolling buffers properly...
                */
                for (k = 0; k < nimgs; k++) {
                    /* Put initial lines of data into scrolling buffers here. */
                    InitFloatSect (pic[k], buf, ipsci[k], line, width, dim_x);
                    InitFloatSect (thresh[k], buf, ipsci[k], line, width,dim_x);
                    InitFloatSect (spthresh[k], buf, ipsci[k],line,width,dim_x);
                    InitShortSect (mask[k], bufdq, ipdq[k], line, width, dim_x);

		    /* Rescale input data by exposure time, if necessary */
		    if (bunit[k] == COUNTRATE) {
			for (ii = 0; ii < buffheight; ii++) {
			     for (i=0; i < dim_x; i++) {
				  pic[k][ii][i] *= efac[k];
			     }
			}
		    }

                    /* Scale the pic value by the sky value and exposure time
                        for comparison to the detection threshhold.
                    */
                    
                    for (ii = 0; ii < buffheight; ii++){
                        for (i=0; i < dim_x; i++){
			    if (efac[k] > 0.) {
                                pic[k][ii][i] = (pic[k][ii][i] - skyval[k]) /
					         efac[k];
			    } else {
				pic[k][ii][i] = 0.;
			    }
                        }
                    }
                    
                } /* End loop over images */
            } /* End if...else line > 0 */

            for (n = 0; n < nimgs; n++) {
                            
                /* 
                    Select line of data to process 
                    This data has already been read into the buffer.
                    However, this data gets modified 
                */
                efacn = efac[n];
                skyvaln = skyval[n];
                exp2n = exp2[n];
                
		/* Only process an image if it has an exposure time > 0. */
		if (efacn > 0.0) {
                    memcpy (bufdq, mask[n][width], dim_x * sizeof(short));
                
                    /* If no shading correction, this will return all ONEs for
                        shadcorr, to avoid divide by ZERO errors.
                    */
                    getShadcorr (shadbuff, line, shad_dimy, dim_x, efacn,
				 shadf_x, shadcorr);

                    /* calculate the threshold for each pixel
                    If initgues is set to minimum, calculate threshhold based
                        on the sigma for this iteration times the variance
                    Otherwise, compute the threshhold based on the pixel
                        values directly corrected by amp gain/noise and 
                        shading correction.
                    */
                    if (strncmp(par->initgues,"minimum",3) == 0 && iter == 0) {
                        for (i = 0; i < dim_x; i++) {
	                    thresh[n][width][i] = sig2 * PPix(avevar,i,line);
	                    spthresh[n][width][i] = sig2 * PPix(avevar,i,line);
                        }
                    } else {
                        for (i = 0; i < ampx; i++) {
                            /* APPLY SHADCORR correction here, as necessary 
                                SHADCORR buffer defaults to ONE if SHADCORR is 
                                not performed.
                            */
                            dum = PPix(ave,i,line)*efacn/shadcorr[i] + skyvaln;

                            /* clip the data at zero */
                            val = (dum > 0.) ? dum : 0.;
			    /* compute sky subtracted pixel value for use with 
			    ** SCALENSE */
			    pixsky = (dum-skyvaln > 0.) ? dum-skyvaln : 0.;

                            /* Apply noise and gain appropriate for AMP used for
			    ** this pixel */
                            thresh[n][width][i] = sig2 *
			    ((nse[0] + val/gn[0] + SQ(scale * pixsky))) / exp2n;
			    /* Compute threshold without SCALENSE for use with
			    ** SPILL pixels */
			    spthresh[n][width][i] = sig2 *
			        ((nse[0] + val/gn[0])) / exp2n;
                        
                        } /* End of loop over first amp used for line */                    
                        for (i = ampx; i < dim_x; i++) {

                            dum = PPix(ave,i,line)*efacn/shadcorr[i] + skyvaln;

                            /* clip the data at zero */
                            val = (dum > 0.) ? dum : 0.;
			    /* compute sky subtracted pixel value for use with
			    * SCALENSE */
			    pixsky = (dum-skyvaln > 0.) ? dum-skyvaln : 0.;

                            /* Apply noise and gain appropriate for AMP used
			    ** for this pixel */
                            thresh[n][width][i] = sig2 *
			    ((nse[1] + val/gn[1] + SQ(scale * pixsky))) / exp2n;
			    /* Compute threshold without SCALENSE for use with
			    ** SPILL pixels */
			    spthresh[n][width][i] = sig2 *
			        ((nse[1] + val/gn[1])) / exp2n;

                        } /* End of loop over second amp used for this line */                    
                    }
 
                    /* exclude points: pixels marked with SPILL will 
                        not propagate the flagging to its neighbors.
                        Also, if pixels are marked with a SERIOUS DQ flag
                        in the input, reject it as well.  */
                    for (i = 0; i < dim_x; i++) {
                        if (((bufdq[i] & dqpat) != OK) && (bufdq[i]!=SPILL)) {
                            mask[n][width][i] = EXCLUDE;
                        }
                    }

                    for (i = 0; i < dim_x; i++){
            
                        /* find the CR by using statistical rejection */
                        if (SQ(pic[n][width][i]-PPix(ave,i,line)) > 
                          thresh[n][width][i] && mask[n][width][i] != EXCLUDE) {
                            mask[n][width][i] = HIT;
                        
                            if (width == 0) continue;
                            /* mark the surrounding pixels also as CR */
                            for (jj = 0; jj < buffheight; jj++) {
                                jndx = line - width + jj;
                                if (jndx < 0 || jndx >= dim_y) continue;
                            
                                /* Distance from buffer center */
                                j2 = SQ(width - jj);

                                for (ii = i-width; ii <= i+width; ii++) {
                                    if ((float)(SQ(ii-i)+j2) > rej2) continue;
                                    if (ii >= dim_x || ii < 0) continue;

                                   if (SQ(pic[n][jj][ii]-PPix(ave,ii,jndx)) <=
				        psig2*spthresh[n][jj][ii]) continue;
                                    if (mask[n][jj][ii] != HIT) {
				        mask[n][jj][ii] = SPILL;
				    }
                                }
                            }
                        }
                    } /* End of loop over i */

                    /* accumulate the total counts in each good pixel */
                    for (i = 0; i < dim_x; i++)  {                                       
                        if ( (mask[n][width][i] & maskdq) == OK ) {
                            /* add the sky-subtracted but UN-scaled counts */
                            sum[i] += pic[n][width][i] * efacn;
                            PIX(efacsum,i,line,dim_x) += efacn;
                        } 
                    } /* End loop over i (X) */

                    /* On the last iteration accumulate variance and DQ vals */
                    if (iter == (niter-1)) {
                
                        /* accumulate the variance only during the last
                           iteration and ONLY for non-HIT or non-SPILL pixels.
                        */
                        for (i = 0; i < ampx; i++){

			    if ((mask[n][width][i] & maskdq) == OK) {
                                 dum = pic[n][width][i]*efacn + skyvaln;
                                 /* clip the data at zero */
                                 val = (dum > 0.) ? dum : 0.;
				 /* Removed SCALENSE from error. 20-Feb-2004 HAB
                                 sumvar[i]+= nse[0] + val/gn[0] + SQ(scale*val);
				 */
                                 sumvar[i]+= nse[0] + val/gn[0];
			    }
                        }

                        for (i = ampx; i < dim_x; i++){

			    if ((mask[n][width][i] & maskdq) == OK) {
                                 dum = pic[n][width][i]*efacn + skyvaln;
                                 /* clip the data at zero */
                                 val = (dum > 0.) ? dum : 0.;
				 /*
                                 sumvar[i]+= nse[1] + val/gn[1] + SQ(scale*val);
				 */
                                 sumvar[i]+= nse[1] + val/gn[1];
			    }
                        }
                    
                        for (i = 0; i < dim_x; i++) {

                            /* output DQF is just the logical OR of all 
                            (original) input DQF */
                            bufdq[i] = bufdq[i] | PDQPix(&dq2,i,line);

			    /* remove any CR or SPILL flags that were in the
			    input files from a previous run of wf3rej or were
			    set from a previous image */
                            bufdq[i] = bufdq[i] & nocr;
                            bufdq[i] = bufdq[i] & nospill;

			    /* sval - which is used to set the DQ value in the
			    output crj file - will initially pick up a value of
			    CR from the initialization of the dq array, but as
			    soon as 1 good input is encountered the CR value
			    will be removed */
                            sval = bufdq[i] | PDQPix(dq,i,line);
                        
                            /* A pixel marked as CR or SPILL will have a CR
			       value in the output images */
                            if (mask[n][width][i] == HIT ||
			        mask[n][width][i] == SPILL) {
                                bufdq[i] = bufdq[i] | crflag;
                                (*nrej)++;

			    /* if the input is good, remove the CR value from
			    the output image array */
			    } else
				sval = sval & nocr;
                        
			    /* Store the values arrived at so far in the
			    output arrays */
                            PDQSetPix(&dq2,i,line,bufdq[i]);
                            PDQSetPix(dq,i,line,sval);

                        } /* End loop over x position */
                    
                        /* compress bufdq into byte mask, line-by-line */
			/* this will be uncompressed later to be written back
			   into the DQ arrays of the input files */
                        writeBitLine (bufdq, n, line, dim_x, crflag, crmask);

                    } /* End of last iteration block */
		} /* End if(efacn > 0.) block */
            } /* End loop over images */

            /* If no shading correction is done, this will return all zeroes.
                getShadCorr can't be called here because each pixel now has
                a different exposure time, so the raw reference data will be
                applied directly.            
             */
            getShadLine (shadbuff, line, shad_dimy, dim_x/rx, shadline);

            pixexp = 0.;
            /* calculate the new average after the rejection */
            /* Then, on last iteration, calculate final ERROR arrays */
            for (i = 0; i < dim_x; i++) {
                pixexp = PIX(efacsum,i,line,dim_x);
                if (pixexp > 0.) {
                    PPix(ave,i,line) = (sum[i] / pixexp)/(1+shadline[i]/pixexp);
                    if (iter == (niter-1)) {
                        PPix(avevar,i,line) = sqrt(sumvar[i])/pixexp;
		    }
                } else {
                    if (iter == (niter-1)){
                        PPix(ave,i,line) = par->fillval;
                        PPix(avevar,i,line) = par->fillval;
                    }
                }
            } 
        } /* End of loop over lines */	 
    } /* End loop for each iteration */
    
    if (par->verbose) {
        trlmessage("Finished all iterations, now writing out results...");
        TimeStamp("Finished all iterations...","");
    }

    /* Write out CR hit information to input images,
       if par->mask was set...
    */
    if (par->mask) {
        /* Close all references to the images so we can open the
           data quality ones as read/write. */
        for (n=0; n<nimgs; n++) {
            closeImage (ipsci[n]);
            closeImage (ipdq[n]);
        }

        for (n=0; n<nimgs; n++) {
            /*  reopen DQ as read/write*/
            initHdr (&dqhdr);

            ipdqn = openUpdateImage (imgname[n], "dq", grp[n], &dqhdr);

            for (line = 0; line < dim_y; line++) {
                getShortLine (ipdqn, line, bufdq);
                readBitLine (crmask, n, line, dim_x, crflag, nocr, bufdq);
                putShortLine (ipdqn, line, bufdq);
            } /* End loop over lines in each image */

            /* close images... This should be done by the calling routine! */
            closeImage (ipdqn);
            freeHdr (&dqhdr);
        } /* End loop over images */

        /* Reopen all images in readonly mode. */
        for (n=0; n<nimgs; n++) {
            ipsci[n] = openInputImage (imgname[n], "sci", grp[n]);
            ipdq[n] = openInputImage (imgname[n], "dq", grp[n]);
        }
    } /* End if */

    /* Use this marker to allow easier clean-up after an error condition 
    An error condition will have already set status to something else
    that should be passed on...
    
    cleanup: ;*/


    /* free memory */
    free (sum);
    free (sumvar);
    free (zerofbuf);
    free (zerosbuf);
    
    for (k=0;k<nimgs;k++) {
        freeFloatBuff (pic[k], buffheight);
        freeFloatBuff (thresh[k], buffheight);
        freeFloatBuff (spthresh[k], buffheight);
        freeShortBuff (mask[k], buffheight);
    }
    free (mask);
    free (pic);	
    free (thresh);
    free (spthresh);
    free (buf);
    free (bufdq);
    free (exp2);
    freeBitBuff (crmask, nimgs, dim_y);	
    free (shadline);
    free (shadcorr);
    freeFloatBuff(shadbuff, shad_dimy);
    freeShortData (&dq2);

    if (par->shadcorr == PERFORM) {
        freeHdr (&scihdr);
        closeImage (ipshad);
    }

    return (status);
}

/* ------------------------------------------------------------------*/
/*                          allocFloatBuff                           */
/* ------------------------------------------------------------------*/


static float **allocFloatBuff (int lines, int numpix) {
    float **sect;
    int i;

    sect = (float **) calloc (lines, sizeof(float *));

    for (i=0; i<lines; i++){
        sect[i] = (float *) calloc (numpix, sizeof(float));
    }

    return (sect);
}

/* ------------------------------------------------------------------*/
/*                          allocShortBuff                           */
/* ------------------------------------------------------------------*/


static short **allocShortBuff (int lines, int numpix){
    short **sect;
    int i;

    sect = (short **) calloc (lines, sizeof(short *));

    for (i=0; i<lines; i++){
        sect[i] = (short *) calloc (numpix, sizeof(short));
    }

    return (sect);
}


/* ------------------------------------------------------------------*/
/*                          freeFloatBuff                            */
/* ------------------------------------------------------------------*/


static void freeFloatBuff (float **sect, int lines){
    int i;

    for (i=0; i<lines; i++) free(sect[i]) ;
    free (sect);
}

/* ------------------------------------------------------------------*/
/*                          freeShortBuff                            */
/* ------------------------------------------------------------------*/


static void freeShortBuff (short **sect, int lines){
    int i;

    for (i=0; i<lines; i++) free(sect[i]) ;
    free (sect);
}


/* ------------------------------------------------------------------*/
/*                          allocBitBuff                             */
/* ------------------------------------------------------------------*/


static Byte ***allocBitBuff (int nimgs, int lines, int numbits){
    Byte ***mask;
    int y,n;
    float remainder;
    int nsize, nrem;

    mask = (Byte ***) calloc (nimgs, sizeof(Byte **));
    nsize = numbits/SIZE_BYTE;
    /* Compute remainder, if any, for nsize. */
    remainder = ((float) numbits/SIZE_BYTE) - nsize;
    /* Convert remainder to either 1 or 0 extra bytes that needs to
    ** be added to buffer. */
    nrem = (remainder > 0.) ? 1 : 0;

    for (n=0; n<nimgs; n++){
        /* Scale the buffer to 1 BIT per pix */
        mask[n] = (Byte **) calloc ( lines, sizeof(Byte *));
        for (y=0; y<lines; y++){
            /* Scale the buffer to 1 BIT per pix
                'numbits/SIZE_BYTE' should always be even, such that
                there is never any remainder...
		However, when it is not even, we need to add 1 for the
		remainder.
            */
            mask[n][y] = (Byte *) calloc ( (nsize+1*nrem), sizeof(Byte));
        }
    }

    return (mask);
}

/* ------------------------------------------------------------------*/
/*                          freeBitBuff                              */
/* ------------------------------------------------------------------*/


static void freeBitBuff (Byte ***crmask, int nimgs, int lines){
    int y,n;

    for (n=0; n<nimgs; n++) {
        for (y=0; y<lines; y++) {
            free(crmask[n][y]) ;
        }
        free(crmask[n]) ;
    }
    free (crmask);
}

/* ------------------------------------------------------------------*/
/*                          printBitLine                             */
/* ------------------------------------------------------------------*/


static void printBitLine (Byte ***crmask, int img, int line, int nx){

    int     x,i;
    Byte    bit;
    int     pix;

    /* Print out values in crbuff to STDOUT using 'X' and '.' */

    for (x = 0; x < (nx/SIZE_BYTE); x++) {
        /* Set each bit in compressed buffer */
        for (bit = 0x80,i=0; bit > 0; bit=(bit>>1),i++) {
            if ( (crmask[img][line][x] & bit) > 0) { 
                pix = x * SIZE_BYTE + i;
                sprintf(MsgText, "Compressed hit at %d,%d", pix, line);
                trlmessage (MsgText); 
            } /*	printf("X");
            else
                printf("."); 
            */
        }		
    }
    /*printf("\n"); */

}


/* ------------------------------------------------------------------*/
/*                          readBitLine                              */
/* ------------------------------------------------------------------*/


static void readBitLine (Byte ***crmask, int img, int line, int numbits, short
crflag, short nocr, short *bufdq) {

    int     x;
    Byte    bit;
    /* int pix, i; */
    int     buffx=0;

    /* Print out values in mask to line-by-line short int buffer */
    for (x = 0; x < (numbits/SIZE_BYTE); x++) {
        /* Decompress buffer here */
        for (bit = 0x80; bit > 0; bit=(bit>>1)) {
            /* pix = x * SIZE_BYTE + i; */

            if ( (crmask[img][line][x] & bit) != 0 ) {
                /*
                ** printf("Compressed hit at %d,%d  with bit = %x ",pix,line,bit);
                ** printf(" and mask = %x\n",crmask[img][line][x]);
                */
                bufdq[buffx] = bufdq[buffx] | crflag;					
            } else
                bufdq[buffx] = bufdq[buffx] & nocr;

            buffx++;
        }		
    }

}

/* ------------------------------------------------------------------*/
/*                          writeBitLine                             */
/* ------------------------------------------------------------------*/


static void writeBitLine (short *bufdq, int img, int line, int numbits, short
crflag, Byte ***crmask) 
{
    int x;					/* Loop counters */
    int buffx;

    /* Set values in mask based on bufdq  */
    /* Each CRHIT was set using : bufdq[i] | crflag */
    for (x = 0; x < numbits; x++) {
        /* Set each bit in compressed buffer */
        if ( ( (bufdq[x]) & crflag) == crflag)  {
            buffx=x/SIZE_BYTE;				
            crmask[img][line][buffx] |=  (0x80 >> ((x)%(SIZE_BYTE)));
        } 
    }
}

/* ------------------------------------------------------------------*/
/*                          scrollFloatBuff                          */
/* ------------------------------------------------------------------*/


/* This function will scroll a subsection up a 2-d array */
static void scrollFloatBuff (float *sect, int line, int nlines, int bufflines, int numpix, float **subsect, float *zero) {
/* Parameters:
**	float *sect         i: line from input image (2-d array)
**	int line            i: line number from image to add to sub-section
**	int nlines          i: number of lines in original image
**	int bufflines       i: number of lines in subsection
**	int numpix          i: number of pixels per line
**	float **subsect     o: scrolled subsection
*/

    int i;	
    float *begptr; /* Use for first line in buffer */

    /* If there is only one line in the buffer, simply copy it out */
    if (bufflines == 1) {
        memcpy (subsect[0], sect, numpix * sizeof(float));
        return;
    }
    
    begptr = *subsect; /* Save first line in buffer */
        
    /* Shift lines in the buffer up, 
    **	moving subsect[1] into subssect[0], and so on..
    **
    ** NOTE: This could probably be optimized to use pointers
    **		more directly to get faster code.
    */ 
    for (i=0; i < bufflines-1; i++){
        *(subsect+i) = *(subsect+i+1); 
    }
    
    /* Now, put pointer from first line back as last line
    **  This recycles the pointer; the data will be overwritten.
    */
    *(subsect + bufflines -1) = begptr;

    /* Finally, copy in new data into last line of buffer.
    **	If the line we want to write to buffer is valid,...
    */
    if (line < nlines && bufflines > 1) { 
        /* Copy new line into last line of buffer. */
        memcpy (subsect[bufflines-1], sect, numpix * sizeof(float));	 
    } else {
        /* Otherwise, set buffer values to default values */
        memcpy (subsect[bufflines-1], zero, numpix * sizeof(float));
    }

}

/* ------------------------------------------------------------------*/
/*                          scrollShortBuff                          */
/* ------------------------------------------------------------------*/


/* This function will scroll a subsection up a 2-d array */
static void scrollShortBuff (short *ssect, int line, int nlines, int bufflines, int numpix, short **subssect, short *szero) {
/* Parameters:
**	short *ssect         i: line from input DQ image (2-d short array)
**	int line            i: line number from sect to add to sub-section
**	int nlines          i: number of lines in original section
**	int bufflines       i: number of lines in subsection
**	int numpix          i: number of pixels per line
**	short **ssubsect     o: scrolled subsection
*/

    int i;	
    short *begsptr; /* Use for first line in buffer */

    /* If there is only one line in the buffer, simply copy it out */
    if (bufflines == 1) {
        memcpy (subssect[0], ssect, numpix * sizeof(short));
        return;
    }
    
    begsptr = *subssect; /* Save first line in buffer */
    
    /* Shift lines in the buffer up, 
    **	moving subsect[1] into subssect[0], and so on..
    **
    ** NOTE: This could probably be optimized to use pointers
    **		more directly to get faster code.
    */ 
    for (i=0; i < bufflines-1; i++){
        subssect[i] = subssect[i+1]; 
    }

    /* Now, put pointer from first line back as last line */
    *(subssect + bufflines - 1) = begsptr;

    /* Finally, copy in new data into last line of buffer.
    **	If the line we want to read in is valid,...
    */
    if (line < nlines && bufflines > 1) {
        /* Copy new line into last line of buffer. */
        memcpy (subssect[bufflines-1], ssect, numpix * sizeof(short));
    } else {
        /* Otherwise, set buffer values to default values */
        memcpy (subssect[bufflines-1], szero, numpix * sizeof(short));
    }

}

/* ------------------------------------------------------------------*/
/*                          initShad                                 */
/* ------------------------------------------------------------------*/


static int initShad (Hdr *scihdr, int dimx, char *shadname, int chipext,
		     int *refx, int *rx, int *ry, int *x0, int *y0) {

    SingleGroupLine y;
    extern int status;
    int same_size = 0;

    int FindLineHdr (Hdr *, Hdr *, int, int, int *, int *, int *, int *, int *);

    initSingleGroupLine (&y);

    /* Get the shutter shading image data. */
    openSingleGroupLine (shadname, chipext, &y);
    if (hstio_err())
        return (status = OPEN_FAILED);

    *refx = y.sci.tot_nx;
    /* Compare binning of science image and reference image;
        get the same_size flag, and get info about binning and offset
        for use by bin2d.
    */
    if(FindLineHdr(scihdr, &y.sci.hdr, dimx, *refx, &same_size, rx, ry, x0, y0))
       return (status);

    closeSingleGroupLine (&y);
    freeSingleGroupLine (&y);

    return (status);
}

/* ------------------------------------------------------------------*/
/*                          getShadLine                              */
/* ------------------------------------------------------------------*/


static void getShadLine (float **shad, int line, int nlines, int dimx,
			 float *shadline) {

/* Parameters:
	shad                i: shadfile buffer
	line                i: line number from input image we are working with
	nlines              i: number of lines in shadfile buffer
	dimx                i: size of shad buffer
	shadline            o: line of shadfile data to be applied to image
*/	
    int offset;

    offset = line % nlines;

    /* Copy appropriate line from buffer into shadline */
    memcpy (shadline, shad[offset], dimx * sizeof(float));
	
}

/* ------------------------------------------------------------------*/
/*                          getShadcorr                              */
/* ------------------------------------------------------------------*/


static void getShadcorr (float **shad, int line, int nlines, int dimx,
			 float efacn, int shadf_x, float *shadline) {

/* Parameters:
    shad                i: shadfile buffer
    line                i: line number from input image we are working with
    nlines              i: number of lines in shadfile buffer
    efacn               i: fraction of exposure time for image
    dimx                i: size of shad buffer
    shadf_x             i: number of pixels in input shadfile line
                         [This will be ZERO if no shading correction is applied]
    shadline            o: line of shadfile data to be applied to image
*/	
    int offset;

    int         multkline (float *, float, int);
    int         addkline (float *, float, int);

    offset = line % nlines;

    /* Copy appropriate line from buffer into shadline */
    memcpy (shadline, shad[offset], dimx * sizeof(float));
    
    /* Only when doing shading correction, perform this expensive
        math operation. 
    */
    if (shadf_x > 0)  multkline (shadline, 1./efacn, dimx);
    
    /* Always make sure that shadline is either 1 or 1+1/efacn,
        to avoid divide by zero errors.
    */
    addkline (shadline, 1., dimx);
	
}

/* ------------------------------------------------------------------*/
/*                          getShadBuff                              */
/* ------------------------------------------------------------------*/


static void getShadBuff (IODescPtr *ipshad, int line, int shad_dimy, int dim_x,
		 int shad_x, int rx, int ry, int x0, int y0, float **shad) {

    int i,j,k,m,zline;
    float **ysect, **zsect;
    float *zl;
    int shad_dimx;

    void unbinArray (float **, int, int, int, int, float **);

    /* Compute the final expanded size in X for the shadfile */
    shad_dimx = shad_x * rx;
    if (shad_x > 0)
	zl = (float *) calloc (shad_dimx, sizeof(float));
    else
	zl = (float *) calloc (dim_x, sizeof(float));
    
    /* If we are NOT performing SHADCORR, then populate the buffer 'shad'
        with ZEROES.
    */
    if (shad_x == 0) {
        for (j = 0; j < shad_dimy; j++) {
            memcpy (shad[j], zl, dim_x * sizeof(float));     
        }
        /* We don't need to do anything else once this is set... */
        free (zl);
        return;   
    }

    /* Apply the shutter shading image. */

    /* 
    Shading correction image may be binned down more than image, and
	    might need to be expanded to match the image size.
        If rx and ry are 1, then it simply uses the shadfile as is.
    */
    /* Set up scratch spaces to expand shadfile		*/
    /* Input file section */
    ysect = allocFloatBuff (SECTLINES, shad_x);
    /* Output (fully expanded) file section */
    zsect = allocFloatBuff (shad_dimy, shad_dimx);

    /* Initialize row counter for shadfile, based on offset
        calculated by FindLine.					
        X and Y offsets are input in binned coordinates,
        so we need to convert X to unbinned coordinates
        for use below. 
    */
    j = y0 + (int)(line/ry);
    k = x0 * rx;
    
    /* Read in shadfile here */
    for (i=0; i < SECTLINES; i++)
        getFloatLine (ipshad, j+i, ysect[i]);

    /* Expand binned reference data, if necessary.	*/
    unbinArray (ysect, shad_x, SECTLINES, rx, ry, zsect);

    /* For each line in expanded section, copy out a 
        single line and trim it to the science image size */
    for (zline=0; zline < shad_dimy; zline++) {

        /* Copy out individual expanded lines from reference data */ 
        memcpy (zl, zsect[zline], shad_dimx * sizeof(float));

        /* Trim each line down here, as necessary (k could be 0)... */
        for (m = 0, i = k;  m < dim_x;  m++, i++)
	        shad[zline][m] = zl[i];
    }	

    /* Clean up scratch space... */
    freeFloatBuff (zsect,shad_dimy);
    freeFloatBuff (ysect, SECTLINES);

    free(zl);           /* done with z */
}

/* ------------------------------------------------------------------*/
/*                          unbinArray                               */
/* ------------------------------------------------------------------*/


/* This routine will expand a simple 2-D float array by binx and biny, 
	without offsets. Furthermore, this routine does NOT update header 
	information.
*/
void unbinArray (float **a, int inx, int iny, int binx, int biny, float **b){

    int     m, n;               /* pixel index in output array */
    int     i, j;               /* pixel index in input array */
    float   p, q, r, s;         /* for interpolating */
    float   xoffset, yoffset;   /* for getting location of output in input */
    float   ai, aj;             /* location in input corresponding to m,n */
    float   value;              /* interpolated value */
    int     onx, ony;           /* size of output array */

    void InterpInfo (float, int, int *, float *, float *);

    onx = inx * binx;
    ony = iny * biny;
    xoffset = 0.0;
    yoffset = 0.0;

    if (binx == 1 && biny == 1) {

        /* Same size, so just copy. */

        /* Copy the science data. */
        for (n = 0;  n < ony;  n++)
            memcpy (b[n], a[n], onx * sizeof(float));

    } else {

        /* Science data array. */
        for (n = 0;  n < ony;  n++) {
            aj = ((float)n - yoffset) / (float)biny;
            InterpInfo (aj, iny, &j, &r, &s);
            for (m = 0;  m < onx;  m++) {
                ai = ((float)m - xoffset) / (float)binx;
                InterpInfo (ai, inx, &i, &p, &q);
                value = p * r * a[j][i] +
                    q * r * a[j][i+1] +
                    p * s * a[j+1][i] +
                    q * s * a[j+1][i+1];
                b[n][m] = value;
            }
        }
    }	

}


/* ------------------------------------------------------------------*/
/*                          InitFloatSect                            */
/* ------------------------------------------------------------------*/

/* 
    This routine performs all the initial bookkeeping for the scrolling
    data buffers. 
        - Initializes the scrolling buffer by populating it with the first
            lines of data from the image.
*/
static void InitFloatSect (float **sect, float *buf, IODescPtr *ipdat, int line, int width, int dimx){
/* Parameters:
float     **sect  i/o: scrolling buffer
float     *buf    i: single line of scratch space
IODescPtr *ipdat  i: file handle for working image (image being processed)
int       line    i: number of current line from working image
int       width   i: number of lines in buffer on either side of current line
int       dimx    i: number of pixels in each line
*/	
    int     l;
        
    /* Fill initial buffer with first lines of image */
    for (l = 0; l <= width; l++){
        getFloatLine (ipdat, line+l, buf);
        /* Copy new line into last line of buffer. */
        memcpy (sect[width+l], buf, dimx * sizeof(float));	 
    }
}


/* ------------------------------------------------------------------*/
/*                          InitShortSect                            */
/* ------------------------------------------------------------------*/

/* 
    This routine performs all the initial bookkeeping for the scrolling
    data buffers. 
        - Initializes the scrolling buffer by populating it with the first
            lines of data from the image.
*/
static void InitShortSect (short **sect, short *sbuf, IODescPtr *ipdat,
			   int line, int width, int dimx) {
/* Parameters:
short       **sect   i/o: scrolling buffer
short       *sbuf    i: single line of scratch space
IODescPtr   *ipdat   i: file handle for working image (image being processed)
int         line     i: number of current line from working image
int         width    i: number of lines in buffer on either side of current line
int         dimx     i: number of pixels in each line
*/	
    int     l;
    Hdr     dqhdr;

    /* Fill initial buffer with first lines of image */
    initHdr(&dqhdr);
    getHeader(ipdat,&dqhdr);
    for (l = 0; l <= width; l++){
        getShortLine (ipdat, line+l, sbuf);
        /* Copy new line into last line of buffer. */
        memcpy (sect[width+l], sbuf, dimx * sizeof(short));	 
    }
    freeHdr(&dqhdr);

}
