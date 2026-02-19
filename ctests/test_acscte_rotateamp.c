#include "hstio.h"
#include "acs.h"
#include "pcte_gen3_funcs.h"

int main(int argc, char **argv) {
    extern int status;
    int iy, ix, amp_id, test_status=0;
    FloatTwoDArray da;

    /* AMP_B or AMP_C, [ny][nx] */
    float expected[3][4] = {{9, 6, 3, 0},
                            {10, 7, 4, 1},
                            {11, 8, 5, 2}};

    /* Input array before rotation (nx=3, ny=4):

       0 1 2
       3 4 5
       6 7 8
       9 10 11
     */
    initFloatData(&da);
    if (allocFloatData(&da, 3, 4, False)) {
        return -1;
    }
    for (iy=0; iy<4; iy++) {
        for (ix=0; ix<3; ix++) {
            Pix(da, ix, iy) = iy * 3 + ix;
        }
    }

    // TODO: Need to test all amps
    amp_id = AMP_C;
    if ((status = rotateAmpData_acscte(&da, amp_id))) {
        return status;
    }

    // TODO: check against truth
    for (iy=0; iy<3; iy++) {
        for (ix=0; ix<4; ix++) {
            if (Pix(da, ix, iy) != expected[iy][ix]) {
                if (!test_status) {
                    printf("test_acscte_rotateamp failed for amp=%d\n", amp_id);
                    test_status = 1;
                }
                printf("ix=%d, iy=%d\n  expected=%.1f\n  got=%.1f\n", ix, iy, expected[iy][ix], Pix(da, ix, iy));
            }
        }
    }

    freeFloatData(&da);
    return test_status;
}
