# include <math.h>		/* for sqrt */

#include "hstio.h"
#include "wf3.h"
#include "wf3info.h"


/* remove amplifier cross-talk */
int cross_talk_corr(WF3Info *wf3, SingleGroup *x) {
    extern int status;
    const int arr_rows = x->sci.data.ny;
    const int arr_cols = x->sci.data.nx;
    int i, j, cur_amp;  /* iteration variables */
    float *y;
    double corr_fac, cur_err;

    /* Correction coefficients (ABCD) from WFC3 ISR 2012-02 */
    const double intercept[NAMPS] = {0, 0, 0, 0};
    const double slope[NAMPS] = {-6.0494304e-5, -2.0746221e-4, -7.9701178e-5, -2.3177171e-4};

    /* Crosstalk correction in electrons but then converted back to DN */
    for (i = 0; i < arr_rows; i++) {
        /* Copy out original row for corr_fac calculation. */
        if ((y = calloc(arr_cols, sizeof(float))) == NULL) {
            return (status = OUT_OF_MEMORY);
        }
        for (j = 0; j < arr_cols; j++) {
            y[j] = Pix(x->sci.data, j, i);
        }

        for (j = 0; j < arr_cols; j++) {
            if (x->group_num == 1) {
                if (j < wf3->ampx) {
                    cur_amp = AMP_C;
                } else {
                    cur_amp = AMP_D;
                }
            } else {
                if (j < wf3->ampx) {
                    cur_amp = AMP_A;
                } else {
                    cur_amp = AMP_B;
                }
            }
            corr_fac = intercept[cur_amp] + y[arr_cols - j - 1] * wf3->atodgain[cur_amp] * slope[cur_amp];
            Pix(x->sci.data, j, i) = (y[j] * wf3->atodgain[cur_amp] - corr_fac) / wf3->atodgain[cur_amp];

            /* Propagate error; assume ERR of correction is sqrt(corr_fac) */
            cur_err = Pix(x->err.data, j, i) * wf3->atodgain[cur_amp];
            Pix(x->err.data, j, i) = sqrt(cur_err * cur_err + fabs(corr_fac)) / wf3->atodgain[cur_amp];
        }
    }

    return status;
}
