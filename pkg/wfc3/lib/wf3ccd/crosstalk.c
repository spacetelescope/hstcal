# include <math.h>		/* for sqrt */

#include "hstio.h"
#include "wf3.h"
#include "wf3info.h"


/* remove amplifier cross-talk */
int cross_talk_corr(WF3Info *wf3, SingleGroup *x) {
    extern int status;
    const int arr_rows = x->sci.data.ny;
    const int arr_cols = x->sci.data.nx;
    const int half_nx = x->sci.data.nx / 2;
    const int tot_pix = arr_rows * arr_cols;
    int i, j, cur_amp;  /* iteration variables */
    int n_skipped[NAMPS] = {0, 0, 0, 0};
    float *y;
    double corr_fac, cur_err;

    /* Correction coefficients (ABCD) from WFC3 ISR 2012-02 */
    double intercept[NAMPS] = {0.0180206, 0.15501201, -0.038376406, 0.19124641};
    const double slope[NAMPS] = {-6.0494304e-5, -2.0746221e-4, -7.9701178e-5, -2.3177171e-4};

    if ((y = calloc(arr_cols, sizeof(float))) == NULL) {
        return (status = OUT_OF_MEMORY);
    }

    /* Factor in gain into intercept to avoid data unit conversion.
       Crosstalk correction should have been in electrons but data in DN.
    */
    for (i = 0; i < NAMPS; i++) {
        intercept[i] /= wf3->atodgain[i];
    }

    /* Crosstalk correction. */
    for (i = 0; i < arr_rows; i++) {
        /* Copy out original row for corr_fac calculation. */
        for (j = 0; j < arr_cols; j++) {
            y[j] = Pix(x->sci.data, j, i);
        }

        for (j = 0; j < arr_cols; j++) {
            if (x->group_num == 1) {
                if (j < half_nx) {
                    cur_amp = AMP_C;
                } else {
                    cur_amp = AMP_D;
                }
            } else {
                if (j < half_nx) {
                    cur_amp = AMP_A;
                } else {
                    cur_amp = AMP_B;
                }
            }
            corr_fac = -1.0 * (intercept[cur_amp] + y[arr_cols - j - 1] * slope[cur_amp]);

            /* Only fix when we can recover the signal, not removing more signal */
            if (corr_fac > 0) {
                Pix(x->sci.data, j, i) = y[j] + corr_fac;

                /* Propagate error; assume ERR of correction is sqrt(corr_fac) */
                cur_err = Pix(x->err.data, j, i);
                Pix(x->err.data, j, i) = sqrt(cur_err * cur_err + corr_fac);
            } else {
                n_skipped[cur_amp]++;
            }
        }
    }

    if (x->group_num == 1) {
        for (i=AMP_C; i<NAMPS; i++) {
            trlmessage("    amp=%d slope=%lf intercept=%lf n_skipped=%d/%d", i, slope[i], intercept[i], n_skipped[i], tot_pix);
        }
    } else {
        for (i=0; i<AMP_C; i++) {
            trlmessage("    amp=%d slope=%lf intercept=%lf n_skipped=%d/%d", i, slope[i], intercept[i], n_skipped[i], tot_pix);
        }
    }

    free(y);
    return status;
}
