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

    iodesc = openUpdateImage(argv[argc-1], "TST", 1, &hdr);
    getShortData(iodesc, &da);

    for (y = 0; y < da.ny; ++y) {
        for (x = 0; x < da.nx; ++x) {
            assert(PPix(&da, x, y) == 42);
        }
    }

    /* Now make it a non-constant array */
    PPix(&da, 12, 34) = 43;

    putShortData(iodesc, &da);

    closeImage(iodesc);
    freeShortData(&da);

    return 0;
}
