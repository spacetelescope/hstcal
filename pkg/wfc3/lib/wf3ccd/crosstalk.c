# include <math.h>		/* for sqrt */

#include "hstio.h"
#include "wf3.h"
#include "wf3info.h"

static void dn_to_e (SingleGroup *, int, int, int, float *);
static void e_to_dn (SingleGroup *, int, int, int, float *);


/* remove amplifier cross-talk */
int cross_talk_corr(WF3Info *wf3, SingleGroup *x) {
    extern int status;
    const int arr_rows = x->sci.data.ny;
    const int arr_cols = x->sci.data.nx;
    int i, j, cur_amp;  /* iteration variables */
    float corr_fac, cur_err, *y;

    /* Correction coefficients (ABCD) from WFC3 ISR 2012-02 */
    const float intercept[NAMPS] = {0.0180206, 0.15501201,-0.038376406, 0.19124641};
    const float slope[NAMPS] = {-0.060494304, -0.20746221, -0.079701178, -0.23177171};

    /* Convert to electrons */
    for (i = 0; i < arr_rows; i++) {
        dn_to_e(x, i, wf3->ampx, wf3->ampy, wf3->atodgain);
    }

    /* Crosstalk correction */
    for (i = 0; i < arr_rows; i++) {
        /* Copy out original row for corr_fac calculation. */
        if ((y = calloc(arr_cols, sizeof(float))) == NULL) {
            return (status = OUT_OF_MEMORY);
        }
        for (j = 0; j < arr_cols; j++) {
            y[j] = Pix(x->sci.data, j, i);
        }

        for (j = 0; j < arr_cols; j++) {
            if (i < wf3->ampy) {
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
            corr_fac = intercept[cur_amp] + y[arr_cols - j - 1] * slope[cur_amp];
            Pix(x->sci.data, j, i) = y[j] - corr_fac;

            /* Propagate error; assume ERR of correction is sqrt(corr_fac) */
            cur_err = Pix(x->err.data, j, i);
            Pix(x->err.data, j, i) = (float) sqrt(cur_err * cur_err + fabs(corr_fac));
        }
    }

    /* Convert back to DN */
    for (i = 0; i < arr_rows; i++) {
        e_to_dn(x, i, wf3->ampx, wf3->ampy, wf3->atodgain);
    }

    return status;
}


static void dn_to_e (SingleGroup *a, int line, int ampx, int ampy, float *gain) {

/* arguments:
SingleGroupLine *a	io: input data; output product
int line		i: line number of input line
int ampx, ampy		i: amp parameters
float *gain		i: atogain parameter

Adapted from multgn1d for fullframe only.
*/
	int i;
	int dimx=a->sci.data.nx;

    /* Since both the science and error arrays operate on the
       same line of data, we will combine them into one loop over
       X values.
    */
	if (line < ampy) {
        /* This line has 2-AMP readout */
        /* Apply gain for first amp to first part of line */
        for (i = 0;  i < ampx;  i++) {
            Pix(a->sci.data, i, line) = gain[AMP_C] * Pix(a->sci.data, i, line);
	        Pix(a->err.data, i, line) = gain[AMP_C] * Pix(a->err.data, i, line);
        }
        /* Apply gain for second amp over remainder of line */
        for (i = ampx;  i < dimx;  i++) {
            Pix(a->sci.data, i, line) = gain[AMP_D] * Pix(a->sci.data, i, line);
	        Pix(a->err.data, i, line) = gain[AMP_D] * Pix(a->err.data, i, line);
        }

    } else {
        /* This line has 2-AMP readout */
        /* Apply gain for first amp to first part of line */
        for (i = 0;  i < ampx;  i++) {
            Pix(a->sci.data, i, line) = gain[AMP_A] * Pix(a->sci.data, i, line);
	        Pix(a->err.data, i, line) = gain[AMP_A] * Pix(a->err.data, i, line);
        }
        /* Apply gain for second amp over remainder of line */
        for (i = ampx;  i < dimx;  i++) {
            Pix(a->sci.data, i, line) = gain[AMP_B] * Pix(a->sci.data, i, line);
            Pix(a->err.data, i, line) = gain[AMP_B] * Pix(a->err.data, i, line);
        }
    }
}


static void e_to_dn (SingleGroup *a, int line, int ampx, int ampy, float *gain) {

    /* arguments:
    SingleGroupLine *a	io: input data; output product
    int line		i: line number of input line
    int ampx, ampy		i: amp parameters
    float *gain		i: atogain parameter

    Adapted from multgn1d for fullframe only.
    */
    int i;
    int dimx=a->sci.data.nx;

    /* Since both the science and error arrays operate on the
       same line of data, we will combine them into one loop over
       X values.
    */
    if (line < ampy) {
        /* This line has 2-AMP readout */
        /* Apply gain for first amp to first part of line */
        for (i = 0;  i < ampx;  i++) {
            Pix(a->sci.data, i, line) = Pix(a->sci.data, i, line) / gain[AMP_C];
            Pix(a->err.data, i, line) = Pix(a->err.data, i, line) / gain[AMP_C];
        }
        /* Apply gain for second amp over remainder of line */
        for (i = ampx;  i < dimx;  i++) {
            Pix(a->sci.data, i, line) = Pix(a->sci.data, i, line) / gain[AMP_D];
            Pix(a->err.data, i, line) = Pix(a->err.data, i, line) / gain[AMP_D];
        }

    } else {
        /* This line has 2-AMP readout */
        /* Apply gain for first amp to first part of line */
        for (i = 0;  i < ampx;  i++) {
            Pix(a->sci.data, i, line) = Pix(a->sci.data, i, line) / gain[AMP_A];
            Pix(a->err.data, i, line) = Pix(a->err.data, i, line) / gain[AMP_A];
        }
        /* Apply gain for second amp over remainder of line */
        for (i = ampx;  i < dimx;  i++) {
            Pix(a->sci.data, i, line) = Pix(a->sci.data, i, line) / gain[AMP_B];
            Pix(a->err.data, i, line) = Pix(a->err.data, i, line) / gain[AMP_B];
        }
    }
}
