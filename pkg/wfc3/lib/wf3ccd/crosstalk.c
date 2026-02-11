# include <math.h>		/* for sqrt */

#include "hstio.h"
#include "wf3.h"
#include "wf3info.h"

static void dn_to_e (SingleGroup *, int, int, int, float *);
static void e_to_dn (SingleGroup *, int, int, int, float *);


/* remove amplifier cross-talk */
void cross_talk_corr(WF3Info *wf3, SingleGroup *x) {
    const int arr_rows = x->sci.data.ny;
    const int arr_cols = x->sci.data.nx;
    int i, j;  /* iteration variables */

    /* Convert to electrons */
    for (i = 0; i < arr_rows; i++) {
        dn_to_e(x, i, wf3->ampx, wf3->ampy, wf3->atodgain);
    }

    /* WFC3 ISR 2012-02 */
    for (i = 0; i < arr_rows; i++) {
        if (i < wf3->ampy) {
            for (j = 0; j < wf3->ampx; j++) {
                // c_corr = c - (-0.038376406 + np.flip(d, 1) * -0.079701178)
                Pix(x->sci.data, j, i) -= -0.038376406 + Pix(x->sci.data, arr_cols - j - 1, i) * -0.079701178;
            }
            for (j = wf3->ampx; j < arr_cols; j++) {
                // d_corr = d - (0.19124641 + np.flip(c, 1) * -0.23177171)
                Pix(x->sci.data, j, i) -= 0.19124641 + Pix(x->sci.data, arr_cols - j - 1, i) * -0.23177171;
            }
        } else {
            for (j = 0; j < wf3->ampx; j++) {
                // a_corr = a - (0.0180206 + np.flip(b, 1) * -0.060494304)
                Pix(x->sci.data, j, i) -= 0.0180206 + Pix(x->sci.data, arr_cols - j - 1, i) * -0.060494304;
            }
            for (j = wf3->ampx; j < arr_cols; j++) {
                // b_corr = b - (0.15501201 + np.flip(a, 1) * -0.20746221)
                Pix(x->sci.data, j, i) -= 0.15501201 + Pix(x->sci.data, arr_cols - j - 1, i) * -0.20746221;
            }
        }
    }

    /* Propagate error */
    for (i = 0; i < arr_rows; i++) {
        for (j = 0; j < arr_cols; j++) {
            Pix(x->err.data, j, i) = (float) sqrt(fabs(Pix(x->sci.data, j, i)));
        }
    }

    /* Convert back to DN */
    for (i = 0; i < arr_rows; i++) {
        e_to_dn(x, i, wf3->ampx, wf3->ampy, wf3->atodgain);
    }
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
