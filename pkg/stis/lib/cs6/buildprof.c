# include <stdio.h>
# include <string.h>
# include <math.h>
# include <float.h>

# include "xtables.h"
# include "hstio.h"

# include "stis.h"
# include "calstis6.h"
# include "hstcalerr.h"

# define DEFAULT_PSTEP  256  /* default pixel range */
# define MAX_ITER       20   /* maximum iterations in sigma-clip */

static void Drizzle (StisInfo6 *, int, int, int);

/*
   Builds compressed profile array from the raw array created
   by the spectrum extraction routine.

   This code works based on pixel ranges only. At initialization,
   it verifies if there is a wavelength range in place, and if
   there is one, it converts it into a pixel range using a linear
   approximation to the dispersion relation. This should be
   accurate enough for the purposes of profile building, since no
   abrupt changes in the profile dependency with wavelength are
   expected.

   On 16Jun00 a major revision took place. Instead of compressing the
   profile data along the dispersion direction within each pixel range,
   the contents of the 2-D profile array are NOT compressed anymore.
   Instead, the information stored in the range index array is later
   used by the output routine to output the entire 2-D contents of each
   pixel range into a single row of the output table. This enables the
   calstis6 extraction code to access the entire raw profile data for
   each range, complete with individual Y positions for each extracted pixel.
   This in turn enables the extraction code to rebuild local profiles with
   better resolution due to the intrinsic dithering that exists in the
   profile data. This dithering comes from the fact the spectrum traces
   are not parallel to the image lines.

   Profiles are still sigma-cleaned and normalized.

   Another major revision took place on 01Dec00. The raw profile data is
   again beign compressed into bins, but now instead of just coadding
   along the dispersion direction, we "drizzle" (shift-and-add with
   subsampling) into a finer sampled array.

   Profiles are *not* sigma-cleaned for now.




   Revision history:
   ----------------
   25 Jun 98  -  Implemented (I.Busko)
   29 Jun 98  -  Rejection ranges (IB)
   04 Sep 98  -  Normalize by data's total flux, set bad regions to -1 (IB)
   19 Oct 98  -  Fix bug in profile array initialization (IB)
   27 Oct 98  -  Warn if any profile value < 0 or extremes > 0.1 flux (IB)
   01 Dec 99  -  Clean cosmic rays (IB)
   07 Dec 99  -  Rejection flag in profile builder (IB)
   11 Apr 00  -  Interpolated profile is one pixel smaller (IB)
   16 Jun 00  -  Keep entire 2-D profile array (IB)
   01 Dec 00  -  Shift-and-add (drizzle) in profile building (IB)
*/

void BuildOutputProfile (StisInfo6 *sts, RowContents *row) {

        int i, j, k, i1, i2, iclip;
        int npoints;
        double wstep, pix;
        double sum, sumsq, lupper, llower, avg, stddev;

        /* If no ranges were specified, adopt some sensible default. */
        if (sts->profile_wstep == 0.0 && sts->profile_pstep == 0)
            sts->profile_pstep = DEFAULT_PSTEP;
        /* Translate wavelength range into pixel range. */
        if (sts->profile_wstep != 0.0) {
            sts->profile_pstep = (int)(sts->profile_x /
                                 ((fabs (row->wave[(sts->profile_x)-1] -
                                 row->wave[0]) / sts->profile_wstep)));
        }

        /* If pixel range results in nonsense, adopt default. */
        if (sts->profile_pstep < 1 || sts->profile_pstep > sts->profile_x)
            sts->profile_pstep = DEFAULT_PSTEP;

        /* Compute how many pixel ranges are needed to encompass the
           current spectral order. The last range might be shorter
           that the others.
        */
        sts->profile_msize = sts->profile_x / sts->profile_pstep +
            ((sts->profile_x % sts->profile_pstep != 0) ? 1 : 0);

        /* Compute linear dispersion coefficient. */
        wstep = (row->wave[(sts->profile_x)-1] - row->wave[0]) /
                 sts->profile_x;

        /* Alloc memory for output wavelength and pixel ranges. */
        sts->profile_sn   = (double *) calloc (sts->profile_msize,
                                               sizeof (double));
        sts->profile_minw = (double *) calloc (sts->profile_msize,
                                               sizeof (double));
        sts->profile_maxw = (double *) calloc (sts->profile_msize,
                                               sizeof (double));
        sts->profile_minp = (short *) calloc (sts->profile_msize,
                                              sizeof (short));
        sts->profile_maxp = (short *) calloc (sts->profile_msize,
                                              sizeof (short));

        /* Initialize range end point indexes to point to first range. */
        i1 = 0;
        i2 = sts->profile_pstep - 1;

        /* Loop over ranges. */
        for (k = 0; k < sts->profile_msize; k++) {

            /* Clean cosmic rays. Bin pixels in the A1 direction,
               determining the mean and stddev for each bin. Flag
               deviant pixels and iterate MAX_ITER times. The dq
               array comes initialized with values set to NO_GOOD_DATA
               where pixels were excluded by the extraction routine in
               X1DSpec based on DQ flags.
            */

            for (j = 0; j < sts->profile_y; j++) {
                lupper =  DBL_MAX;
                llower = -DBL_MAX;
                for (iclip = 0; iclip <= MAX_ITER; iclip++) {
                    sum   = 0.0;
                    sumsq = 0.0;
                    npoints = 0;
                    for (i = i1; i <= i2; i++) {
                        if (sts->profile_dq[i][j] != NO_GOOD_DATA &&
                            sts->profile[i][j]    <  lupper       &&
                            sts->profile[i][j]    >  llower) {
                            pix = sts->profile[i][j];
                            sum   += pix;
                            sumsq += pix * pix;
                            npoints++;
                        } else
                            sts->profile_dq[i][j] = NO_GOOD_DATA;
                    }
                    if (npoints > 2) {
                        avg = sum / (double)npoints;
                        stddev = (double)sqrt(((sumsq / (double)npoints) -
                                 (avg * avg)));
                        lupper = avg + stddev * sts->psclip;
                        llower = avg - stddev * sts->psclip;
                    }
                }
            }

            /* Compute total counts from good pixels in current range. */
            /*
            sts->profile_sn[k] = 0.0;
            empty = 1;
            for (j = 0; j < sts->profile_y; j++) {
                for (i = i1; i <= i2; i++) {
                    if (!sts->profile_rej[i] &&
                        sts->profile_dq[i][j] != NO_GOOD_DATA) {
                        sts->profile_sn[k] += sts->profile[i][j];
                        empty = 0;
                    }
                }
            }
            */
            /* S/N. If lower than threshold, flag profile. */
            /*
            if (sts->profile_sn[k] > 0.0)
                sts->profile_sn[k] = sqrt (sts->profile_sn[k]);
            else
                sts->profile_sn[k] = 0.0;
            if (sts->profile_sn[k] < sts->profile_minsn) {
                for (j = 0; j < (sts->profile_y)-1; sts->profile[k][j++] = -1.);
            }
            */

            /* Normalize. */

            for (i = i1; i <= i2; i++) {
                sum = 0.0;
                for (j = 0; j < sts->profile_y; j++) {
                    if (!sts->profile_rej[i] &&
                        sts->profile_dq[i][j] != NO_GOOD_DATA) {
                        sum += sts->profile[i][j];
                    } else
                        sts->profile_rej[i] = 1;
                }
                if (sum != 0.0)
                    for (j = 0; j < sts->profile_y;
                         sts->profile[i][j++] /= sum);
            }

            /* Build subsampled profile */
            Drizzle (sts, i1, i2, k);

            /* Update min/max wavelength/pixel for this range. */
            sts->profile_minw[k] = i1 * wstep +row->wave[0];
            sts->profile_maxw[k] = i2 * wstep +row->wave[0];
            sts->profile_minp[k] = i1 + 1;
            sts->profile_maxp[k] = i2 + 1;

            /* Update range end point indexes. */
            i1 += sts->profile_pstep;
            i2 += sts->profile_pstep;
            if (i2 >= row->npts)
                i2 = row->npts -1;
        }
}



static void Drizzle (StisInfo6 *sts, int i1, int i2, int k) {

        int i;                  /* scan columns inside bin                */
        int j;                  /* scan pixels in column                  */
        int ix;                 /* scan subpixels                         */
        double x1, x2;          /* edges of large pixel in subpixel array */
        int ix1, ix2;           /* and their integer counterparts         */
        double sum, npt;        /* normalization */

        /* Loop over columns in bin. */

        for (i = i1; i < i2; i++) {
            if (!sts->profile_rej[i]) {

                /* Loop over input pixels in column. */

                for (j = 0; j < sts->profile_y; j++) {
                    if (sts->profile_dq[i][j] != NO_GOOD_DATA) {

                        /* Compute where in subpixel array the edges of
                           the current pixel are.
                        */
                        x1 = (j - (int)(sts->profile_y / 2) - 0.5 +
                              sts->profile_offset[i]) * sts->subscale +
                              (int)(sts->subprof_size / 2);

                        x2 = (j - (int)(sts->profile_y / 2) + 0.5 +
                              sts->profile_offset[i]) * sts->subscale +
                              (int)(sts->subprof_size / 2);

                        ix1 = (int)x1;
                        ix2 = (int)x2;

                        if (ix1 >= 0 && ix2 < sts->subprof_size) {

                            /* Add contributions from edge areas. */

                            sts->subprof[k][ix1] += sts->profile[i][j] *
                                                (1.0 - x1 + ix1);
                            sts->subprof[k][ix2] += sts->profile[i][j] *
                                                (x2 - ix2);

                            /* Scan subpixels in between. */

                            if ((ix2 - ix1) > 1) {
                                for (ix = ix1 + 1; ix < ix2; ix++) {
                                    sts->subprof[k][ix] += sts->profile[i][j];
                                }
                            }
                        }
                    }
                }
            }
        }

        /* Normalize subpixelized profile. */

        sum = 0.0;
        npt = 0.0;
        for (i = 0; i < sts->subprof_size; i++) {
            sum += sts->subprof[k][i];
            npt += 1.0;
        }
        for (i = 0; i < sts->subprof_size; i++) {
            sts->subprof[k][i] /= sum;
        }
}



void FreeOutputProfile (StisInfo6 *sts) {

        free (sts->profile_minw);
        free (sts->profile_maxw);
        free (sts->profile_minp);
        free (sts->profile_maxp);
        free (sts->profile_sn);
        sts->profile_msize = 0;
}



