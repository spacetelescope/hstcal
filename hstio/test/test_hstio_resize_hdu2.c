#include <assert.h>
#include <stdio.h>
#include "hstio.h"

int
main(int argc, char **argv) {
    ShortTwoDArray da;
    IODescPtr iodesc;
    short buffer[2048];
    Hdr hdr;
    int x, y;

    initShortData(&da);
    initHdr(&hdr);

    iodesc = openUpdateImage(argv[argc-1], "TST", 1, &hdr);
    getShortData(iodesc, &da);

    for (y = 0; y < da.ny; ++y) {
        for (x = 0; x < da.nx; ++x) {
            assert(PPix(&da, x, y) == 42);
        }
    }

    /* Now make it a non-constant array */
    for (x = 0; x < 2048; ++x) {
        buffer[x] = x;
    }

    if (putShortLine(iodesc, 1000, buffer)) {
        return 1;
    }

    closeImage(iodesc);
    freeShortData(&da);

    return 0;
}
