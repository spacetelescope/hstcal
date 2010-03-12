# include <stdio.h>
# include <stdlib.h>
# include "ctables.h"

int main (int argc, char **argv) {

	IRAFPointer tp;

	if (argc != 2) {
	    printf ("syntax:  a.out filename\n");
	    exit (2);
	}
	tp = c_tbtopn (argv[1], IRAF_READ_ONLY, 0);
	printf ("debug:  file opened\n");
	c_tbtclo (tp);
	printf ("debug:  file closed\n");

	exit (0);
}
