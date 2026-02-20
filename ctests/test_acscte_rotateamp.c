#include "hstio.h"
#include "acs.h"
#include "pcte_gen3_funcs.h"

static int setup_input_array(FloatTwoDArray *, int, int);
static int compare_arrays(FloatTwoDArray *, float *, int, int);
static int worker();

static int setup_input_array(FloatTwoDArray *da, int nx, int ny) {
    int iy, ix;

    initFloatData(da);
    if (allocFloatData(da, nx, ny, False)) {
        return OUT_OF_MEMORY;
    }
    for (iy=0; iy<ny; iy++) {
        for (ix=0; ix<nx; ix++) {
            Pix(*da, ix, iy) = iy * nx + ix;
        }
    }

    return ACS_OK;
}

static int compare_arrays(FloatTwoDArray *da, float *truth, int truth_nx, int truth_ny) {
    int iy, ix, ii, test_status=ACS_OK;

    if (da->nx != truth_nx || da->ny != truth_ny) {
        printf("ERROR: Dimension mismatch!\n\n  Expected: nx=%d ny=%d\n       Got: nx=%d ny=%d\n\n",
               truth_nx, truth_ny, da->nx, da->ny);
        test_status = SIZE_MISMATCH;
    }

    for (iy=0; iy<truth_ny; iy++) {
        for (ix=0; ix<truth_nx; ix++) {
            ii = iy * truth_nx + ix;
            if (Pix(*da, ix, iy) != truth[ii]) {
                test_status = ERROR_RETURN;
                break;
            }
        }
        if (test_status) {
            break;
        }
    }

    if (test_status) {
        printf("ERROR: Arrays are different!\n\n  Expected:\n    ");
        for (iy=0; iy<truth_ny; iy++) {
            for (ix=0; ix<truth_nx; ix++) {
                ii = iy * truth_nx + ix;
                printf("%.1f ", truth[ii]);
            }
            printf("\n    ");
        }
        printf("\n  Got:\n    ");
        for (iy=0; iy<da->ny; iy++) {
            for (ix=0; ix<da->nx; ix++) {
                printf("%.1f ", Pix(*da, ix, iy));
            }
            printf("\n    ");
        }
        printf("\n");
    }

    return test_status;
}

static int worker() {
    const char amps[5] = "ABCD\0";
    int amp_id=AMP_C, test_status;
    FloatTwoDArray da;
    int nx=3, ny=4;

    /* Input array before rotation:

       0 1 2
       3 4 5
       6 7 8
       9 10 11
     */
    if ((test_status = setup_input_array(&da, nx, ny))) {
        freeFloatData(&da);
        return test_status;
    }

    // TODO: Need to test all amps
    printf("==== rotateAmpData_acscte AMP %c ====\n", amps[amp_id]);
    if ((test_status = rotateAmpData_acscte(&da, amp_id))) {
        freeFloatData(&da);
        return test_status;
    }

    /* AMP_B or AMP_C
        9 6 3 0
       10 7 4 1
       11 8 5 2
    */
    int truth_nx=ny;
    int truth_ny=nx;
    float *truth = malloc(sizeof(float) * truth_nx * truth_ny);
    truth[0] = 9; truth[1] = 6; truth[2] = 3; truth[3] = 0;
    truth[4] = 10; truth[5] = 7; truth[6] = 4; truth[7] = 1;
    truth[8] = 11; truth[9] = 8; truth[10] = 5; truth[11] = 2;

    test_status = compare_arrays(&da, truth, truth_nx, truth_ny);
    free(truth);
    freeFloatData(&da);
    return test_status;
}

int main(int argc, char **argv) {
    int test_status;

    test_status = worker();

    return test_status;
}
