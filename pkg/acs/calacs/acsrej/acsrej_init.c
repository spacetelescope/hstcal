# include   <stdio.h>
# include   <string.h>
# include   <stdlib.h>
#include "hstcal.h"
# include   "hstio.h"

# include   "acs.h"
# include   "rej.h"
# include   "acsrej.h"
# include   "hstcalerr.h"

# define    OK          (short)0

/*  acsrej_init -- get the initial average pixel values

  Description:
  ------------
  Get the initial average according to the specified scheme

  Date          Author          Description
  ----          ------          -----------
  22-Sep-1998   W.J. Hack       initial version, uses multiamp noise,gain
  20-Mar-2000   W.J. Hack       corrected problem with rog2, needed to be array
  15-Apr-2002   W.J. Hack       correctly zero'd out npts array, added
                                DQ checking for BOTH median and minimum images,
                                also added ERR array for median image.
  26-Aug-2002   J. Blakeslee    Modified threshhold for minimum to use SCALENSE
                                only with sky-subtracted pixels.
*/

int acsrej_init (IODescPtr ipsci[], IODescPtr ipdq[], clpar *par, int nimgs,
                 int dim_x, int dim_y,
                 multiamp noise, multiamp gain, float efac[], float skyval[],
                 SingleGroup *sg, float *work) {
    /*
      Parameters:

      ipsci   i: Array of pointers to SCI extension of the given EXTVER,
                 each pointer is an input image. Unit now in electrons.
      ipdq    i: Array of pointers to DQ extension of the given EXTVER,
                 each pointer is an input image. Unit now in electrons.
      par     i: User specified parameters.
      nimgs   i: Number of input images.
      dim_x, dim_y  i: Image dimension taken from the first input image.
                       All images must have the same dimension.
      noise   i: Calibrated readnoise converted to DN in acsrej_check.c.
                 It is squared in this function as rog2 and then copied to
                 noise2 (using get_nsegn) and subsequently nse.
      gain    i: Calibrated gain in e/DN from primary header keyword
                 ATODGN[AMP], where [AMP] can be A, B, C, or D.
      efac    i: Array of EXPTIME for each image. If all inputs have
                 EXPTIME=0 (all biases), then the values are all set to 1.
      skyval  i: Array of sky values for each input image.
                 Unit now in electrons.
      sg      o: Average image for comparison during CR rejection.
      work    o: Intermediate result used to calculate sg but not used
                 outside this function.
    */

    extern int status;

    float     scale, val, raw, raw0, signal0;
    float     *buf;
    short     *bufdq;
    float     *exp2;
    int       i, j, n;
    int       dum;
    int       *npts, *ipts;
    float     noise2[NAMPS], rog2[NAMPS];
    float     gain2[NAMPS];
    float     nse[2], gn[2];
    int       ampx, ampy, detector, chip;
    int       k, p;
    short     dqpat;
    float     exp2n, expn;
    int       non_zero;
    Hdr       dqhdr;
    int       newbias;

    void      ipiksrt (float [], int, int[]);
    void      get_nsegn (int, int, int, int, float *, float*, float *, float *);

    /* -------------------------------- begin ------------------------------- */

    scale = par->scalense / 100.;
    ampx = gain.colx;
    ampy = gain.coly;
    detector = gain.detector;
    chip = gain.chip;
    dqpat = par->badinpdq;
    newbias = par->newbias;

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


    /* use the stack median to construct the initial average */
    if (strncmp(par->initgues, "median", 3) == 0) {
        for (j = 0; j < dim_y; j++) {
            memset (npts, 0, dim_x*sizeof(int));
            /* Set up the gain and noise values used for this line in
               ALL images */
            if (j < ampy ) {
                gn[0] = gain2[AMP_C];
                gn[1] = gain2[AMP_D];
                nse[0] = noise2[AMP_C];  /* DN^2 */
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
                getFloatLine (ipsci[n], j, buf);  /* electrons */
                getShortLine (ipdq[n], j, bufdq);
                freeHdr(&dqhdr);

                for (i = 0; i < ampx; i++) {
                    if (efac[n] > 0.) {
                        /* Only use GOOD pixels to build initial image.
                           work = DN/s */
                        if ((bufdq[i] & dqpat) == OK) {
                            PIX(work, npts[i], i, nimgs) =
                                (buf[i] - skyval[n]) / efac[n] / gn[0];
                            npts[i] += 1;
                        }
                    } else {
                        PIX(work, npts[i], i, nimgs) = 0.;
                    }
                }
                for (i = ampx; i < dim_x; i++) {
                    if (efac[n] > 0.){
                        /* Only use GOOD pixels to build initial image.
                           work = DN/s */
                        if ((bufdq[i] & dqpat) == OK) {
                            PIX(work, npts[i], i, nimgs) =
                                (buf[i] - skyval[n]) / efac[n] / gn[1];
                            npts[i] += 1;
                        }
                    } else {
                        PIX(work, npts[i], i, nimgs) = 0.;
                    }
                }
            }  /* End of nimgs loop */

            /* AMPS C and D */
            for (i = 0; i < ampx; i++) {
                dum = npts[i];  /* Number of good data points */
                if (dum == 0)
                    Pix(sg->sci.data, i, j) = 0.0F;
                else {
                    /* Setup index array for sorting... */
                    for (p=0; p < nimgs; p++) ipts[p] = p;
                    /* Sort pixel stack and corresponding index array. */
                    ipiksrt (&PIX(work, 0, i, nimgs), dum, ipts);
                    /*
                        Use sorted index array to match proper exptimes to
                        selected pixels for use in ERR array calculation.
                    */
                    /* Even number of input images for this pixel */
                    if ((dum / 2) * 2 == dum) {
                        Pix(sg->sci.data, i, j) =
                            (PIX(work, dum / 2 - 1, i, nimgs) +
                             PIX(work, dum / 2, i, nimgs)) / 2.;
                        expn = (exp2[ipts[dum / 2 - 1]] +
                                exp2[ipts[dum / 2]]) / 2.;
                    } else {
                        Pix(sg->sci.data, i, j) = PIX(work, dum / 2, i, nimgs);
                        expn = exp2[ipts[dum / 2]];
                    }
                }
                raw0 = Pix(sg->sci.data, i, j);   /* DN/s */
                exp2n = (expn > 0.) ? expn : 1.;  /* s^2 */
                if (newbias == 0) {
                    /* This is the 2nd term in DHB formula.
                       ERR = (READNSE^2 + VAL/GN? + (SCALE * VAL)^2) / EXPTIME^2
                       nse is DN^2
                       raw0 is DN/s - why divided by gn?
                                      is /s supposed to be there?
                       exp2n is s^2
                    */
                    Pix(sg->err.data, i, j) =
                        (nse[0] + raw0 / gn[0] + SQ(scale * raw0)) / exp2n;
                } else {
                    /* Only consider readnoise for biases.
                       ERR = READNSE^2 / EXPTIME^2 */
                    Pix(sg->err.data,i,j) = (nse[0]) / exp2n;
                }
            } /* End loop over FIRST AMP used on pixels in the line */

            for (i = ampx; i < dim_x; i++) {
                dum = npts[i];  /* Number of good data points */
                if (dum == 0)
                    Pix(sg->sci.data, i, j) = 0.0F;
                else {
                    /* Setup index array for sorting... */
                    for (p=0; p < nimgs; p++) ipts[p] = p;
                    /* Sort pixel stack and corresponding index array. */
                    ipiksrt (&PIX(work, 0, i, nimgs), dum, ipts);
                    /*
                        Use sorted index array to match proper exptimes to
                        selected pixels for use in ERR array calculation.
                    */
                    /* Even number of input images for this pixel */
                    if ((dum / 2) * 2 == dum) {
                        Pix(sg->sci.data,i,j) =
                            (PIX(work, dum / 2 - 1, i, nimgs) +
                             PIX(work, dum / 2, i, nimgs)) / 2.;
                        expn = (exp2[ipts[dum / 2 - 1]] +
                                exp2[ipts[dum / 2]]) / 2.;
                    } else {
                        Pix(sg->sci.data, i, j) = PIX(work, dum / 2, i, nimgs);
                        expn = exp2[ipts[dum / 2]];
                    }
                }
                raw0 = Pix(sg->sci.data, i, j);  /* DN/s */
                exp2n = (expn > 0.) ? expn : 1.;  /* s^2 */
                if (newbias == 0){
                    /* This is the 2nd term in DHB formula.
                       ERR = (READNSE^2 + VAL/GN? + (SCALE * VAL)^2) / EXPTIME^2
                       nse is DN^2
                       raw0 is DN/s - why divided by gn?
                                      is /s supposed to be there?
                       exp2n is s^2
                    */
                    Pix(sg->err.data, i, j) =
                        (nse[1] + raw0 / gn[1] + SQ(scale * raw0)) / exp2n;
                } else {
                    /* Only consider readnoise for biases.
                       ERR = READNSE^2 / EXPTIME^2 */
                    Pix(sg->err.data, i, j) = (nse[1]) / exp2n;
                }
            } /* End loop over SECOND AMP used on pixels in the line */
        } /* End loop over lines */

    /* use the minimum to construct the initial average */
    } else {
        if (strncmp(par->initgues, "minimum", 3) != 0) {
            sprintf (MsgText,
                     "Invalid INITGUES value %s, reset it to 'minimum'",
                     par->initgues);
            trlwarn (MsgText);
            strcpy (par->initgues, "minimum");
        }

        for (n = 0; n < nimgs; n++) {
            initHdr(&dqhdr);
            getHeader(ipdq[n],&dqhdr);
            for (j = 0; j < dim_y; j++) {
                /* Set up the gain and noise values used for this line in
                   ALL images */
                if (j < ampy ) {
                    gn[0] = gain2[AMP_C];
                    gn[1] = gain2[AMP_D];
                    nse[0] = noise2[AMP_C];  /* DN^2 */
                    nse[1] = noise2[AMP_D];
                } else {
                    gn[0] = gain2[AMP_A];
                    gn[1] = gain2[AMP_B];
                    nse[0] = noise2[AMP_A];
                    nse[1] = noise2[AMP_B];
                }

                getFloatLine (ipsci[n], j, buf);  /* electrons */
                getShortLine (ipdq[n], j, bufdq);

                /* AMPS C and D */
                for (i = 0; i < ampx; i++) {
                    raw = buf[i] / gn[0];  /* DN */
                    raw0 = (raw > 0.) ? raw : 0.;  /* DN */
                    signal0 = ((raw - (skyval[n] / gn[0])) > 0.) ?
                        (raw - (skyval[n] / gn[0])) : 0.;  /* DN */

                    if (efac[n] > 0.) {
                        /* Why not just use signal0? */
                        val = (raw - (skyval[n] / gn[0])) / efac[n];  /* DN/s */
                    } else {
                        val = 0.;
                    }

                    if ( (n == 0) || (val < Pix(sg->sci.data, i, j)) ) {
                        if ((bufdq[i] & dqpat) == OK && (efac[n] > 0.)) {
                            Pix(sg->sci.data, i, j) = val;  /* DN/s */
                            /*Pix(sg->err.data,i,j) = (nse[0]+ raw0/gn[0] +
                              SQ(scale*raw0)) / exp2[n];*/  /* NOT USED */
                            if (newbias == 0){
                                /* This is the 2nd term in DHB formula.
                       ERR = (READNSE^2 + VAL/GN? + (SCALE * VAL)^2) / EXPTIME^2
                       nse is DN^2
                       raw0 is DN - with sky; why divided by gn?
                       signal0 is DN - sky subtracted
                       exp2n is s^2
                                */
                                Pix(sg->err.data, i, j) =
                                    (nse[0] + (raw0 / gn[0]) +
                                     SQ(scale * signal0)) / exp2[n];
                            } else {
                                /* Only consider readnoise for biases.
                                   ERR = READNSE^2 / EXPTIME^2 */
                                Pix(sg->err.data, i, j) = (nse[0]) / exp2[n];
                            }
                        } else {
                            Pix(sg->sci.data, i, j) = 0.;
                            Pix(sg->err.data, i, j) = 0.;
                        }
                    }
                } /* End of loop over FIRST AMP for this line in each image */

                for (i = ampx; i < dim_x; i++) {
                    raw = buf[i] / gn[1];  /* DN */
                    raw0 = (raw > 0.) ? raw : 0.;  /* DN */
                    signal0 = ((raw - (skyval[n] / gn[1])) > 0.) ?
                        (raw - (skyval[n] / gn[1])) : 0.;  /* DN */

                    if (efac[n] > 0.) {
                        /* Why not just use signal0? */
                        val = (raw - (skyval[n] / gn[1])) / efac[n];  /* DN/s */
                    } else {
                        val = 0.;
                    }

                    if ( (n == 0) || (val < Pix(sg->sci.data, i, j) &&
                                      ((bufdq[i] & dqpat) == OK)) ) {
                        Pix(sg->sci.data, i, j) = val;  /* DN/s */
                        if (efac[n] > 0.) {
                            if (newbias == 0){
                                /* This is the 2nd term in DHB formula.
                       ERR = (READNSE^2 + VAL/GN? + (SCALE * VAL)^2) / EXPTIME^2
                       nse is DN^2
                       raw0 is DN - with sky; why divided by gn?
                       signal0 is DN - sky subtracted
                       exp2n is s^2
                                */
                                Pix(sg->err.data, i, j) =
                                    (nse[1] + (raw0 / gn[1]) +
                                     SQ(scale * signal0)) / exp2[n];
                            } else {
                                /* Only consider readnoise for biases.
                                   ERR = READNSE^2 / EXPTIME^2 */
                                Pix(sg->err.data, i, j) = (nse[1]) / exp2[n];
                            }
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
    free (exp2);
    free (ipts);
    free (npts);
    free (buf);
    free (bufdq);
    return (status);
}
