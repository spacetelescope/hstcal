# include	<stdio.h>
# include	<string.h>
# include	<stdlib.h>
# include	<math.h>
# include	"hstio.h"

# include	"stis.h"
# include	"cs2.h"
# include	"calstis2.h"

/*  crrej_loop -- Perform cosmic ray rejection

Description:
------------
Iterate the rejection using different sigma sizes.  Also reject pixels 
next to the cosmic ray hit.
  
  Date		Author		Description
  ----		------		-----------
  26-Apr-1996  J.-C. Hsu	adapt from the SPP code crrej_loop.x
  17-Jun-1999  J.-C. Hsu	unset the CR bit from the input DQ's
  10-Feb-2000  Phil Hodge	Check for out of memory; return int.
   8-May-2000  Phil Hodge	Line number in get/putTypeLine is zero-indexed.
   7-Oct-2002  Phil Hodge	Copy Warren Hack's changes to acsrej_loop.c
				to this file, to handle scalense correctly;
				i.e. subtract the sky before applying scalense.
  22-May-2012  Phil Hodge	Change the declaration of imgname.
*/

int crrej_loop (IODescPtr ipsci[], IODescPtr ipdq[], 
			char *imgname[], int grp [], int nimgs, 
			clpar *par, int niter, int dim_x, int dim_y, 
			float sigma[], float noise[], float gain[], 
			float efac[], float skyval[], FloatTwoDArray *ave, 
			FloatTwoDArray *avevar, float *efacsum, 
			ShortTwoDArray *dq, int *nrej)
{
	Hdr	dqhdr;			/* data quality header structure */
	int	width;
	int	npts;
	int	i, j, k, n, indx, jndx;		/* loop indices */
	int	iter;
	int	ii, jj, j2;
	float	rog2, sig2, psig2, rej2, exp2;	/* square of something */
	float	scale, val, dum, pixsky;
	short	sval, crflag, nocr, dqpat;
	
	/* local mask values */
	short	OK = 0;
	short	HIT = 8;
	short	SPILL = 2;
	short	EXCLUDE = 4;

	/* local data arrays */
	float	*pic;
	float	*sum;
	float	*sumvar;
	float	*thresh, *spthresh;
	short	*mask;
	float	*buf;
	short	*bufdq;

/* -------------------------------- begin ---------------------------------- */

	initHdr (&dqhdr);

	crflag = par->crval;
	dqpat = par->badbits;
	scale = par->scalenoise/100.;
	nocr = ~crflag;

	npts = dim_x * dim_y;

	/* allocate data arrays */
	pic = calloc (npts, sizeof(float));
	thresh = calloc (npts, sizeof(float));
	spthresh = calloc (npts, sizeof(float));
	sum = calloc (npts, sizeof(float));
	sumvar = calloc (npts, sizeof(float));
	mask = calloc (npts, sizeof(short));
	buf = calloc (dim_x, sizeof(float));
	bufdq = calloc (dim_x, sizeof(short));
	if (pic == NULL || thresh == NULL || spthresh == NULL ||
	    sum == NULL || sumvar == NULL ||
	    mask == NULL || buf == NULL || bufdq == NULL) {
	    printf ("ERROR    out of memory in crrej_loop\n");
	    return (2);
	}
	
	/* readout is in DN */
	rej2 = SQ(par->rej);
	psig2 = SQ(par->psigma);
	if (par->psigma <= 0.)
	    psig2 = -1.;

	width = (int) (par->rej+1.e-5);

	/* reset the (ouput) DQF */
	/* set to crflag so it is easier to do the logical AND and OR later */
	for (j = 0; j < dim_y; ++j) {
	    for (i = 0; i < dim_x; ++i)
	        PDQSetPix(dq,i,j,crflag);
	} 

	/* start the rejection iteration */
	for (iter = 0; iter < niter; ++iter) {
	    if (par->verbose) 
		printf ("iteration %d\n", iter+1);

	    sig2 = SQ(sigma[iter]);

	    /* reset the sum/file count arrays */
	    for (indx = 0; indx < dim_x*dim_y; ++indx) 
		sum[indx] = 0.;
	    for (j = 0; j < dim_y; ++j) {
		for (i = 0; i < dim_x; ++i)
		    PIX(efacsum,i,j,dim_x) = 0.;
	    }

	    /* at the last iteration, write out individual masks */
	    if (iter == (niter-1)) {

		/* must close all images first and then reopen only one group*/
		for (k = 0; k < nimgs; ++k) {
		    closeImage (ipsci[k]);
		    closeImage (ipdq[k]);
		}
	    }

	    *nrej = 0;

	    for (n = 0; n < nimgs; ++n) {
		rog2 = SQ(noise[n]);
		exp2 = SQ(efac[n]);

 		if (iter == (niter-1)) {
		    ipsci[n] = openInputImage (imgname[n], "SCI", grp[n]);
		    ipdq[n] = openInputImage (imgname[n], "DQ", grp[n]);
            	    getHeader (ipdq[n], &dqhdr);
		}

	        /* reset the mask */
		for (indx = 0; indx < npts; ++indx)
		    mask[indx] = OK;

                /* calculate the threshold for each pixel */
	        if (strncmp(par->initial, "minimum", 3) == 0 && iter == 0) {
		    for (j = 0; j < dim_y; ++j) {
                        indx = j*dim_x;
			for (i = 0; i < dim_x; ++i) {
		            thresh[indx+i] = sig2 * PPix(avevar,i,j);
		            spthresh[indx+i] = sig2 * PPix(avevar,i,j);
			}
                    }
	        } else {
		    for (j = 0; j < dim_y; ++j) {
                        indx = j*dim_x;
			for (i = 0; i < dim_x; ++i) {
                            dum = PPix(ave,i,j)*efac[n]+skyval[n];

			    /* clip the data at zero */
                            val = (dum > 0.) ? dum : 0.;
			    /* Compute sky subracted pixel value for use
				with SCALENSE. */
			    pixsky = (dum-skyval[n] > 0.) ? dum-skyval[n] : 0.;
		            thresh[indx+i] = sig2 * ((rog2 + val/gain[n] + 
					     SQ(scale * pixsky))) / exp2;
			    /* Compute threshold without SCALENSE for use
				with SPILL pixels. */
		            spthresh[indx+i] = sig2 * (rog2 + val/gain[n]) /
					     exp2;
			}
                    }
                }

		/* read data, subtract sky and scale by the exposure factor */
		getHeader (ipdq[n], &dqhdr);

		for (j = 0; j < dim_y; ++j) {
		    getFloatLine (ipsci[n], j, buf);
		    getShortLine (ipdq[n], j, bufdq);

		    indx = j*dim_x;
		    for (i = 0; i < dim_x; ++i) {

			/* exclude points: pixels marked with SPILL will 
			   not propagate the flagging to its neighbors */
			if ((bufdq[i] & dqpat) != OK)
			    mask[indx+i] = EXCLUDE;
			pic[indx+i] = (buf[i] - skyval[n]) / efac[n];
		    }
		}

		freeHdr (&dqhdr);

		for (j = 0; j < dim_y; ++j) {
		    indx = j*dim_x;
		    for (i = 0; i < dim_x; ++i) {

			/* find the CR by using statistical rejection */
			if (SQ(pic[indx+i]-PPix(ave,i,j)) > 
				thresh[indx+i] && mask[indx+i] != EXCLUDE) {
			    mask[indx+i] = HIT;

			    /* mark the surrounding pixels also as CR */
			    for (jj = j-width; jj <= j+width; ++jj) {
				if (jj >= dim_y || jj < 0) continue;
				j2 = (jj-j) * (jj-j);
				jndx = jj*dim_x;
				for (ii = i-width; ii <= i+width; ++ii) {
				    if ((float)((ii-i)*(ii-i)+j2) > rej2) 
					continue;
				    if (ii >= dim_x || ii < 0) continue;

				    if (SQ(pic[jndx+ii]-PPix(ave,ii,jj)) <=
					    psig2*spthresh[jndx+ii]) continue;
				    if (mask[jndx+ii] != HIT)
					mask[jndx+ii] = SPILL;
				}
			    }
			}

/* for debugging */
/*
if (i==676 && j==642) {
if (n==0) printf("At pixel (%d,%d) ave=%0.3g sigma=%0.3g\n", i+1, j+1, PPix(ave,i,j), sqrt(thresh[indx+i]/sig2));
printf("file=%d scaled DN/sec=%0.3g mask = %d\n", n, pic[indx+i], mask[indx+i]);
}
*/
		    }
		}

		/* accumulate the total counts in each pixel */
		for (j = 0; j < dim_y; ++j) {
		    indx = j*dim_x;
		    for (i = 0; i < dim_x; ++i) {
			if (mask[indx+i] == OK) {

			    /* add the sky-subtracted but UN-scaled counts */
			    sum[indx+i] += pic[indx+i] * efac[n];
			    PIX(efacsum,i,j,dim_x) += efac[n];
			
			    /* accumulate the variance only during the
				last iteration */
			    if (iter == (niter-1)) {
                                dum = pic[indx+i]*efac[n]+skyval[n];

			        /* clip the data at zero */
                                val = (dum > 0.) ? dum : 0.;
				sumvar[indx+i] += rog2 + val/gain[n] + 
							SQ(scale * val);
			    }
			}
		    }
		}

		/* at the last iteration, write out individual masks */
		if (iter == (niter-1)) {

		    /* close images then reopen only DQ as read/write*/
		    closeImage (ipsci[n]);

		    if (par->mask) {
		        closeImage (ipdq[n]);
		        ipdq[n] = openUpdateImage (imgname[n], "dq", 
							grp[n], &dqhdr);
		    }

		    for (j = 0; j < dim_y; ++j) {
		        getShortLine (ipdq[n], j, bufdq);
		        indx = j*dim_x;
			for (i = 0; i < dim_x; ++i) {

			    /* need to unset the CR bit first, in case the 
				task is rerun, otherwise the output DQ will
				always be the same as the last image's DQ */
			    bufdq[i] = bufdq[i] & nocr;

			    /* output DQF is just the logical OR of all 
				(original) input DQF */
			    sval = PDQPix(dq,i,j) | bufdq[i];
			    if (mask[indx+i]==HIT || mask[indx+i] == SPILL) {
				(*nrej)++;
			    } else
				sval = sval & nocr;
			    PDQSetPix(dq,i,j,sval);

			    /* All masked pixels have the same flag value and
			       is combined with what was in the mask, if any. */
		            if (par->mask) {
			        if (mask[indx+i]==HIT || mask[indx+i] == SPILL)
				    bufdq[i] = bufdq[i] | crflag;
				else
				    bufdq[i] = bufdq[i] & nocr;
			    }
		        }
		        if (par->mask)
		            putShortLine (ipdq[n], j, bufdq);
		    }
		    closeImage (ipdq[n]);
	        }
	    }

	    /* calculate the new average after the rejection */
	    for (j = 0; j < dim_y; ++j) {
		indx = j*dim_x;
		for (i = 0; i < dim_x; ++i) {
		    if (PIX(efacsum,i,j,dim_x) > 0.) {
		        PPix(ave,i,j) = sum[indx+i] / PIX(efacsum,i,j,dim_x);
		    } else {
			if (iter == (niter-1)) {
			    PPix(ave,i,j) = par->fillval;
			}
		    }
		}
	    } 
	}

	/* calculate the (unnormalized) error array */
	for (j = 0; j < dim_y; ++j) {
	    indx = j * dim_x;
	    for (i = 0; i < dim_x; ++i) {
		if (PIX(efacsum,i,j,dim_x) > 0.) {
	            PPix(avevar,i,j) = sqrt(sumvar[indx+i]) / 
					PIX(efacsum,i,j,dim_x);
		} else {
		    PPix(avevar,i,j) = par->fillval;
		}
	    }
	} 

	/* free memories */
	free (pic);
	free (thresh);
	free (spthresh);
	free (sum);
	free (sumvar);
	free (mask);
	free (buf);
	free (bufdq);

	return (0);
}
