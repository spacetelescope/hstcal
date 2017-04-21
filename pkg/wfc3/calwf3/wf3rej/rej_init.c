# include   <stdio.h>
# include   <string.h>
# include   <stdlib.h>
# include "hstio.h"

# include   "wf3.h"
# include   "rej.h"
# include   "wf3rej.h"
# include   "hstcalerr.h"
# include   "wf3info.h"

# define    OK          (short)0

/*  rej_init -- get the initial average pixel values

  Description:
  ------------
  Get the initial average according to the specified scheme
  This function 
  
  Date          Author          Description
  ----          ------          -----------
  22-Sep-1998   W.J. Hack       initial version, uses multiamp noise,gain
  20-Mar-2000   W.J. Hack       corrected problem with rog2, needed to be array
  29-Aug-2000   H.A. Bushouse	revised for WFC3 use
  18-Jun-2002   H. Bushouse	correctly zero'd out npts array, added DQ
				checking for BOTH median and minimum images,
				also added ERR array for median image (following
				CALACS changes).
  26-Nov-2002	H. Bushouse	Modified threshold for minimum to use SCALENSE
				only with sky-subtracted pixels (following
				CALACS changes).
  27-Oct-2003	H. Bushouse	Various upgrades to handle 1 or more inputs
				with EXPTIME = 0 (CALACS changes).
  06-Dec-2007	H. Bushouse	Added calls to getHeader before each call to
				getShortLine to prevent getShortLine from
				crashing on null DQ arrays.
  16-Jun-2011	H. Bushouse	Added missing call to free(ipts) at end.
  14-Dec-2011   H. Bushouse	Upgraded to rescale input data that are in
				units of count rates. (PR 69969; Trac #814)
*/

int rej_init (IODescPtr ipsci[], IODescPtr ipdq[], clpar *par, int nimgs,
	      int dim_x, int dim_y, multiamp noise, multiamp gain, float efac[],
	      float skyval[], DataUnits bunit[], SingleGroup *sg, float *work) {

    extern int status;

    float  scale, val, raw, raw0, signal0;
    float  *buf;
    short  *bufdq;
    float  *exp2;
    int    i, j, n;
    int    dum;
    int    *npts, *ipts;
    float  noise2[NAMPS], rog2[NAMPS];
    float  gain2[NAMPS];
    float  nse[2], gn[2];
    int    ampx, ampy, detector, chip;
    int    k, p;
    short  dqpat;
    float  exp2n, expn;
    int    non_zero;
    Hdr    dqhdr;

    void ipiksrt (float [], int, int[]);
    void get_nsegn (int, int, int, int, float *, float*, float *, float *);

    /* -------------------------------- begin ------------------------------ */
    expn=0.0f;
    scale = par->scalense / 100.;
    ampx = gain.colx;
    ampy = gain.coly;
    detector = gain.detector;
    chip = gain.chip;
    dqpat = par->badinpdq;

    ipts = calloc (nimgs, sizeof(int));
    npts = calloc (dim_x, sizeof(int));
    buf = calloc (dim_x, sizeof(float));
    exp2 = (float *) calloc (nimgs, sizeof(float));
    bufdq = calloc (dim_x, sizeof(short));

    for (k = 0; k < NAMPS; k++) {
        gain2[k] = 0.;
        noise2[k] = 0.;
        /* Assumption: ALL images have the same noise/gain values */
        rog2[k] = SQ(noise.val[k]);
    }

    non_zero = 0;
    for (n = 0; n < nimgs; n++) {
        exp2[n] = SQ(efac[n]);
	if (efac[n] > 0.) non_zero++;
    }
    get_nsegn (detector, chip, ampx, ampy, gain.val, rog2, gain2, noise2); 

    /* Use the stack median to construct the initial average */
    if (strncmp(par->initgues, "median", 3) == 0) {
        for (j = 0; j < dim_y; j++) {
            memset (npts, 0, dim_x*sizeof(int));

            /* Set up the gain and noise values used for this line
	    ** in ALL images */
            if (j < ampy ) {
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

            for (n = 0; n < nimgs; n++) {
                initHdr(&dqhdr);
                getHeader(ipdq[n],&dqhdr);
                getFloatLine (ipsci[n], j, buf);
                getShortLine (ipdq[n], j, bufdq);
                freeHdr(&dqhdr);

		/* Rescale SCI data, if necessary */
		if (bunit[n] == COUNTRATE) {
		    for (i = 0; i < dim_x; i++) {
			 buf[i] *= efac[n];
		    }
		}

                for (i = 0; i < dim_x; i++) {
		     if (efac[n] > 0.) {
                         /* Only use GOOD pixels to build initial image */
                         if ((bufdq[i] & dqpat) == OK) {
                             PIX(work,npts[i],i,nimgs) =
						(buf[i] - skyval[n]) / efac[n];
                             npts[i] += 1;
                         }
		     } else {
			 PIX(work,npts[i],i,nimgs) = 0.;
		     }
                }
            }
 
            for (i = 0; i < ampx; i++) {
                dum = npts[i];

                if (dum == 0)
                    Pix(sg->sci.data,i,j) = 0.0F;

                else {

		    /* Setup index array for sorting... */
		    for (p=0; p < nimgs; p++) ipts[p] = p;
		    /* Sort pixel stack and corresponding index array. */
                    ipiksrt (&PIX(work,0,i,nimgs), dum, ipts);

		    /* Use sorted index array to match proper exptimes to
		       selected pixels for use in ERR array calculation. */
                    if ((dum/2)*2 == dum) {
                        /* Even number of input images for this pixel */
                        Pix(sg->sci.data,i,j) = (PIX(work,dum/2-1,i,nimgs) +
						 PIX(work,dum/2,i,nimgs)) / 2.;
			expn = (exp2[ipts[dum/2-1]] + exp2[ipts[dum/2]]) / 2.;
                    } else {
			/* Odd number of input images for this pixel */
                        Pix(sg->sci.data,i,j) = PIX(work,dum/2,i,nimgs);
			expn = exp2[ipts[dum/2]];
		    }
                }
                
                raw0 = Pix(sg->sci.data,i,j);
		exp2n = (expn > 0.) ? expn : 1.;
                Pix(sg->err.data,i,j) = (nse[0]+ raw0/gn[0] + SQ(scale*raw0)) /
					exp2n;
            } /* End loop over FIRST AMP used on pixels in the line */
             
            for (i = ampx; i < dim_x; i++) {
                dum = npts[i];
                if (dum == 0)
                    Pix(sg->sci.data,i,j) = 0.0F;
                else {
		    for (p=0; p < nimgs; p++) ipts[p] = p;
                    ipiksrt (&PIX(work,0,i,nimgs), dum, ipts);
                    if ((dum/2)*2 == dum) {
                        /* Even number of input images for this pixel */
                        Pix(sg->sci.data,i,j) = (PIX(work,dum/2-1,i,nimgs) +
						 PIX(work,dum/2,i,nimgs)) / 2.;
			expn = (exp2[ipts[dum/2-1]] + exp2[ipts[dum/2]]) / 2.;
                    } else {
                        Pix(sg->sci.data,i,j) = PIX(work,dum/2,i,nimgs);
			expn = exp2[ipts[dum/2]];
		    }
                }
                
                raw0 = Pix(sg->sci.data,i,j);
		exp2n = (expn > 0.) ? expn : 1.;
                Pix(sg->err.data,i,j) = (nse[1]+ raw0/gn[1] + SQ(scale*raw0)) /
					exp2n;
            } /* End loop over SECOND AMP used on pixels in the line */
        } /* End loop over lines */

    } else {

        /* use the minimum to construct the initial average */
        if (strncmp(par->initgues, "minimum", 3) != 0) {
            sprintf (MsgText,"Invalid INITGUES value %s, reset it to 'minimum'",
		     par->initgues);
            trlwarn (MsgText);
            strcpy (par->initgues, "minimum");
        }

        for (n = 0; n < nimgs; n++) {
            initHdr(&dqhdr);
            getHeader(ipdq[n],&dqhdr);
            for (j = 0; j < dim_y; j++) { 
                /* Set up the gain and noise values used for this line
		** in ALL images */
                if (j < ampy ) {
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

                getFloatLine (ipsci[n], j, buf);
		getShortLine (ipdq[n],  j, bufdq);

		/* Rescale SCI data, if necessary */
		if (bunit[n] == COUNTRATE) {
		    for (i = 0; i < dim_x; i++) {
			 buf[i] *= efac[n];
		    }
		}

                /* AMPS C and D */
                for (i = 0; i < ampx; i++) {
                    raw = buf[i];
                    raw0 = (raw > 0.)? raw : 0.;
		    signal0 = ((raw - skyval[n]) > 0.) ? (raw - skyval[n]) : 0.;

		    if (efac[n] > 0.) {
                        val = (raw - skyval[n]) / efac[n];
		    } else {
			val = 0.;
		    }
                    if ((n == 0) || (val < Pix(sg->sci.data,i,j)) ) {
			if ((bufdq[i] & dqpat) == OK && (efac[n] > 0.)) {
                             Pix(sg->sci.data,i,j) = val;
                             /*Pix(sg->err.data,i,j) =
			      (nse[0]+ raw0/gn[0] + SQ(scale*raw0)) / exp2[n];*/
			     Pix(sg->err.data,i,j) =
			    (nse[0] + raw0/gn[0] + SQ(scale*signal0)) / exp2[n];
			} else {
			     Pix(sg->sci.data,i,j) = 0.;
			     Pix(sg->err.data,i,j) = 0.;
			}
                    } 
                } /* End of loop over FIRST AMP for this line in each image */

                for (i = ampx; i < dim_x; i++) {
                    raw = buf[i];
                    raw0 = (raw > 0.)? raw : 0.;
		    signal0 = ((raw - skyval[n]) > 0.) ? (raw - skyval[n]) : 0.;

		    if (efac[n] > 0.) {
                        val = (raw - skyval[n]) / efac[n];
		    } else {
			val = 0.;
		    }
                    if ((n == 0) ||
		       (val<Pix(sg->sci.data,i,j) && ((bufdq[i]&dqpat)==OK))) {
                        Pix(sg->sci.data,i,j) = val;
			if (efac[n] > 0.) {
                            Pix(sg->err.data,i,j) = 
			     (nse[1]+ raw0/gn[1] + SQ(scale*signal0)) / exp2[n];
			} else {
			    Pix(sg->err.data,i,j) = 0.;
			}
                    } 
                } /* End of loop over SECOND AMP for this line in each image */

            } /* End of loop over lines in image (y) */
            freeHdr(&dqhdr);
        } /* End of loop over images in set */
    }

    /* free the memory */
    free (ipts);
    free (npts);
    free (buf);
    free (exp2);
    free (bufdq);

    return (status);
}

