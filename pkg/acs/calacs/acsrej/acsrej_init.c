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


/* acsrej_init -- get the initial average pixel values

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
  01-Dec-2015   P.L. Lim        Calculations now entirely in electrons.
  13-Jan-2016   P.L. Lim        Removed variance init and cleaned up function.
*/
int acsrej_init (IODescPtr ipsci[], IODescPtr ipdq[], clpar *par, int nimgs,
                 int dim_x, int dim_y, float efac[], float skyval[],
                 SingleGroup *sg, float *work) {
    /*
      Parameters:

      ipsci   i: Array of pointers to SCI extension of the given EXTVER,
                 each pointer is an input image. Unit now in electrons.
      ipdq    i: Array of pointers to DQ extension of the given EXTVER,
                 each pointer is an input image.
      par     i: User specified parameters.
      nimgs   i: Number of input images.
      dim_x, dim_y  i: Image dimension taken from the first input image.
                       All images must have the same dimension.
      efac    i: Array of EXPTIME for each image. If all inputs have
                 EXPTIME=0 (all biases), then the values are all set to 1.
      skyval  i: Array of sky values for each input image.
                 Unit now in electrons.
      sg      o: Its "sci" component is the average image used for
                 comparison during CR rejection. Unit is e/s.
                 This is really the median or minimum depending on
                 "initgues" provided by the user.
      work    o: Intermediate result used to calculate sg but not used
                 outside this function.
    */
    extern int status;

    float     val, raw, dumf;
    int       i, j, n, p, dum;
    float     *buf;
    short     *bufdq;
    int       *npts, *ipts;
    short     dqpat;
    Hdr       dqhdr;

    void      ipiksrt (float [], int, int[]);

    /* -------------------------------- begin ------------------------------- */

    dqpat = par->badinpdq;

    ipts = calloc (nimgs, sizeof(int));
    npts = calloc (dim_x, sizeof(int));
    buf = calloc (dim_x, sizeof(float));
    bufdq = calloc (dim_x, sizeof(short));

    /* Use the stack median to construct the initial average. */
    if (strncmp(par->initgues, "median", 3) == 0) {
        for (j = 0; j < dim_y; j++) {
            memset (npts, 0, dim_x*sizeof(int));

            for (n = 0; n < nimgs; n++) {
                initHdr(&dqhdr);
                getHeader(ipdq[n], &dqhdr);
                getFloatLine (ipsci[n], j, buf);  /* electrons */
                getShortLine (ipdq[n], j, bufdq);
                freeHdr(&dqhdr);

                /* Only use GOOD pixels to build initial image.
                   work array is already initialized to zeroes in acsrej_do.c */
                if (efac[n] > 0.) {
                    for (i = 0; i < dim_x; i++) {
                        if ((bufdq[i] & dqpat) == OK) {
                            PIX(work, npts[i], i, nimgs) =
                                (buf[i] - skyval[n]) / efac[n];  /* e/s */
                            npts[i] += 1;
                        }
                    }
                }
            }  /* End of nimgs loop */

            /* ALL AMPS */
            for (i = 0; i < dim_x; i++) {
                dum = npts[i];  /* Number of good data points */
                if (dum == 0)
                    Pix(sg->sci.data, i, j) = 0.0F;
                else {
                    /* Setup index array for sorting... */
                    for (p=0; p < nimgs; p++) ipts[p] = p;
                    /* Sort pixel stack and corresponding index array. */
                    ipiksrt (&PIX(work, 0, i, nimgs), dum, ipts);
                     /* Even number of input images for this pixel */
                    if ((dum / 2) * 2 == dum) {
                        Pix(sg->sci.data, i, j) =
                            (PIX(work, dum / 2 - 1, i, nimgs) +
                             PIX(work, dum / 2, i, nimgs)) / 2.;
                    } else {
                        Pix(sg->sci.data, i, j) = PIX(work, dum / 2, i, nimgs);
                    }
                }
            } /* End loop over ALL AMPS used on pixels in the line */
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
            getHeader(ipdq[n], &dqhdr);

            for (j = 0; j < dim_y; j++) {
                getFloatLine (ipsci[n], j, buf);  /* electrons */
                getShortLine (ipdq[n], j, bufdq);

                /* ALL AMPS */
                for (i = 0; i < dim_x; i++) {
                    raw = buf[i];  /* e */
                    dumf = raw - skyval[n];  /* e */

                    if (efac[n] > 0.) {
                        /* Can be negative */
                        val = dumf / efac[n];  /* e/s */
                    } else {
                        val = 0.;
                    }

                    if ( (n == 0) || (val < Pix(sg->sci.data, i, j)) ) {
                        /* If this pixel is bad in the first image,
                           then the min is automatically set to 0. As a
                           result, only negative val is going to be stored and
                           valid positive min is ignored.
                           SLIGHTLY BUGGY HERE??? */
                        if ((bufdq[i] & dqpat) == OK && (efac[n] > 0.)) {
                            Pix(sg->sci.data, i, j) = val;  /* e/s */
                        } else {
                            Pix(sg->sci.data, i, j) = 0.;
                        }
                    }
                } /* End of loop over ALL AMPS for this line in each image */
            } /* End of loop over lines in image (y) */

            freeHdr(&dqhdr);
        } /* End of loop over images in set */
    }

    /* free the memory */
    free (ipts);
    free (npts);
    free (buf);
    free (bufdq);

    return (status);
}
