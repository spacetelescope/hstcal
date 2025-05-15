/* cs2 -- Reject cosmic rays */

# include <stdio.h>
# include <stdlib.h>		/* calloc */
# include <string.h>

# include "c_iraf.h"		/* for c_irafinit */

# include "stis.h"
# include "hstcalerr.h"
# include "cs2.h"

#include "hstcal_memory.h"

static struct TrlBuf trlbuf = {0};

int main (int argc, char **argv) {

	int status = 0;			/* zero is OK */

	char	*input;			/* list of input file names */
	char 	output[STIS_LINE];	/* output file name */
	clpar 	par;			/* parameters used */
	int 	newpar[MAX_PAR+1];	/* user specifiable parameters */

	int cs2_command (int, char **, char *, char *, clpar *, int []);

	c_irafinit (argc, argv);

	if (argc < 2) {
		printf("ERROR    missing input file name\n");
		exit(ERROR_RETURN);
	}

	if ((input = calloc (strlen(argv[1]) + 1, sizeof(char))) == NULL) {
	    printf("ERROR    out of memory in cs2.c\n");
	    exit (ERROR_RETURN);
	}

	PtrRegister ptrReg;
	initPtrRegister(&ptrReg);

	/* Initialize the structure for managing trailer file comments */
	InitTrlBuf ();
	addPtr(&ptrReg, &trlbuf , &CloseTrlBuf);
	trlGitInfo();


	/* Get input and output file names and switches in the command line. */
	status = cs2_command (argc, argv, input, output, &par, newpar);
	if (status != 0) {
	    if (input)
	        free(input);
	    exit (status);
	}

	/* Reject cosmic rays. */
	if ((status = CalStis2 (input, output, &par, newpar))) {
	    printf("Error processing %s.\n", input);
	}

	if (input)
	    free (input);

	freeOnExit(&ptrReg);

	if (status)
	    exit (ERROR_RETURN);
	else
	    exit (0);
}
