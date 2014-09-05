#include <assert.h>
#include <stdio.h>
#include "hstio.h"

int
main(int argc, char **argv) {
    ShortTwoDArray da;
    IODescPtr iodesc;
    Hdr hdr;
    int x, y;

    initShortData(&da);
    initHdr(&hdr);

    iodesc = openUpdateImage(argv[argc-1], "TST", 1, &hdr);
    getShortData(iodesc, &da);

    for (y = 0; y < da.ny; ++y) {
        for (x = 0; x < da.nx; ++x) {
            assert(PPix(&da, x, y) == 42);
            if (PPix(&da, x, y) != 42) {
                return 1;
            }
        }
    }

    /* Now make it a non-constant array */
    PPix(&da, 12, 34) = 43;

    if (putShortData(iodesc, &da)) {
        return 1;
    }

    closeImage(iodesc);
    freeShortData(&da);

    return 0;
}
