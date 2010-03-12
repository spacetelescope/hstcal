# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include "ctables.h"

/* This version is just for testing the file name template functions. */

int main (int argc, char **argv) {

        IRAFPointer imt;
        int done, lenfilename;
        char comment[81];

        if (argc != 2) {
            printf ("syntax:  a.out filename\n");
            exit (2);
        }

        printf ("file name template functions:\n");
        imt = c_imtopen (argv[1]);
        printf ("  number of names in list = %d\n", c_imtlen (imt));
        done = 0;
        while (!done) {
            lenfilename = c_imtgetim (imt, comment, 80);
            done = (lenfilename == 0);
            if (!done)
                printf ("  file name = %s, len = %d\n", comment, lenfilename);
        }
        c_imtclose (imt);
        printf ("file name template closed\n");

        exit (0);
}
