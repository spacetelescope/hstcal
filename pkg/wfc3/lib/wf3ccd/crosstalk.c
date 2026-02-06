#include "hstio.h"
#include "wf3info.h"

/* remove amplifier cross-talk */
int cross_talk_corr(WF3Info *wf3, SingleGroup *x) {
    extern int status;

    /* iteration variables */
    int i, j;

    /* cross talk scaling constant */
    double cross_scale = 9.1e-5;

    double temp;

    const int arr_rows = x->sci.data.ny;
    const int arr_cols = x->sci.data.nx;

    /* TODO: add Norman requirements, GH 701; add change log etc */

    for (i = 0; i < arr_rows; i++) {
        for (j = 0; j < arr_cols; j++) {
            temp = Pix(x->sci.data, arr_cols-j-1, i) * cross_scale;

            Pix(x->sci.data, j, i) += (float) temp;
        }
    }

    return (status);
}
