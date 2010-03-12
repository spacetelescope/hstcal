/* cs2 -- Reject cosmic rays */

# include <stdio.h>
# include <stdlib.h>		/* calloc */
# include <string.h>

# include <c_iraf.h>		/* for c_irafinit */

# include "../stis.h"
# include "../stiserr.h"
# include "../cs2.h"

int main (int argc, char **argv) {

	int status = 0;			/* zero is OK */

	char 	input[STIS_LINE], output[STIS_LINE];	/* file names */
	clpar 	par;			/* parameters used */
	int 	newpar[MAX_PAR+1];	/* user specifiable parameters */

	int cs2_command (int, char **, char *, char *, clpar *, int []);

	c_irafinit (argc, argv);

	/* Get input and output file names and switches in the command line. */
	cs2_command (argc, argv, input, output, &par, newpar);

	/* Reject cosmic rays. */
	if (status = CalStis2 (input, output, &par, newpar)) {
	    printf ("Error processing %s.\n", input);
	}

	if (status)
	    exit (ERROR_RETURN);
	else
	    exit (0);
}
