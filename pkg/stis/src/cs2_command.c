# include <stdio.h>
# include <stdlib.h>		/* atoi, atof */
# include <string.h>

# include "c_iraf.h"

# include "stis.h"
# include "cs2.h"

void cs2_reset (clpar *, int []); 
static int syntax_error (char *);
static int getArgS (char **, int, int *, short *);
static int getArgR (char **, int, int *, float *);
static int getArgT (char **, int, int *, char *);

/*  
   Processes command-line parameters for calstis-2.

   Revision history:
   ----------------
   20 Oct 97  -  Adapted from cs6 commline.c (JC Hsu)
   06 Jul 11  -  Add command-line option --version (Phil Hodge)
   10 Feb 12  -  Add command-line option -r (Phil Hodge)
*/

int cs2_command (int argc, char **argv, char *input, char *output,
              	 clpar *par, int newpar[]) 
{

/* arguments
int argc;               i: input command-line parameters
char **argv;
char *input;		o: input file name
char *output;		o: output file name
clpar *par;		o: user specified parameters
int newpar;		o: parameter been set by the user?
*/

	int 	ctoken;	/* current command-line token being processed */
	char	dummy[5];

	/* reset the parameters */
	cs2_reset (par, newpar); 

	for (ctoken = 1;  ctoken < argc;  ctoken++) {
	    if (strcmp (argv[ctoken], "--version") == 0) {
		PrVersion();
		exit (0);
	    }
	    if (strcmp (argv[ctoken], "-r") == 0) {
		PrFullVersion();
		exit (0);
	    }
	}

	/* not enough arguments */
	if (argc < 3)
	    return (syntax_error ("cs2 input output"));

	/* Get names of input and output files. These are mandatory. */
	strcpy (input,  argv[1]);
	strcpy (output, argv[2]);


	/* Scan remaining parameters. */
	if (argc > 3) {
	    ctoken = 3;
	    while (ctoken < argc) {

	        /* switch must begin with a "-" */
	        if (argv[ctoken][0] != '-')
	            return (syntax_error (argv[ctoken]));

	        else {

	            /* These do not require additional arguments. */

	            if (strcmp("t", argv[ctoken]+1) == 0) { 
	                par->printtime = 1;
	                ctoken++;

	            } else if (strcmp("v", argv[ctoken]+1) == 0) { 
	                par->verbose = 1;
	                ctoken++;

	            /* These require one additional argument, which is
                       handled by the getArg functions.
                    */

	            } else if (strcmp("crmask", argv[ctoken]+1) == 0) { 
			newpar[TOTAL]++;
			newpar[CRMASK] = 1;
	                if (getArgT (argv, argc, &ctoken, dummy))
	                    return (1);
			if (strcmp ("yes", dummy) == 0)
			    par->mask = 1;
			else
			    par->mask = 0;

	            } else if (strcmp("table", argv[ctoken]+1) == 0) { 
	                if (getArgT (argv, argc, &ctoken, par->tbname))
	                    return (1);

	            } else if (strcmp("scale", argv[ctoken]+1) == 0) { 
			newpar[TOTAL]++;
			newpar[SCALENSE] = 1;
	                if (getArgR (argv, argc, &ctoken, &par->scalenoise))
	                    return (1);

	            } else if (strcmp("init", argv[ctoken]+1) == 0) { 
			newpar[TOTAL]++;
			newpar[INITGUES] = 1;
	                if (getArgT (argv, argc, &ctoken, par->initial))
	                    return (1);

	            } else if (strcmp("sky", argv[ctoken]+1) == 0) { 
			newpar[TOTAL]++;
			newpar[SKYSUB] = 1;
	                if (getArgT (argv, argc, &ctoken, par->sky))
	                    return (1);

	            } else if (strcmp("sigmas", argv[ctoken]+1) == 0) { 
			newpar[TOTAL]++;
			newpar[CRSIGMAS] = 1;
	                if (getArgT (argv, argc, &ctoken, par->sigmas))
	                    return (1);

	            } else if (strcmp("radius", argv[ctoken]+1) == 0) { 
			newpar[TOTAL]++;
			newpar[CRRADIUS] = 1;
	                if (getArgR (argv, argc, &ctoken, &par->rej))
	                    return (1);

	            } else if (strcmp("thresh", argv[ctoken]+1) == 0) { 
			newpar[TOTAL]++;
			newpar[CRTHRESH] = 1;
	                if (getArgR (argv, argc, &ctoken, &par->psigma))
	                    return (1);

	            } else if (strcmp("pdq", argv[ctoken]+1) == 0) { 
			newpar[TOTAL]++;
			newpar[BADINPDQ] = 1;
	                if (getArgS (argv, argc, &ctoken, &par->badbits))
	                    return (1);

	            /* No match. */
	            } else
	                return (syntax_error (argv[ctoken]));
	        } 
	    }
	}

	return (0);
}


/*  Get a short parameter from the argv array. */

static int getArgS (char **argv, int argc, int *ctoken, short *value) {
	if (*ctoken <= (argc-2)) {
	    *value = (short) atoi(argv[++(*ctoken)]);
	    (*ctoken)++;
	    return (0);
	} else
	    return (syntax_error (argv[*ctoken]));
}

/*  Get a string parameter from the argv array. */

static int getArgT (char **argv, int argc, int *ctoken, char *value) {
	if (*ctoken <= (argc-2)) {
	    strcpy (value, argv[++(*ctoken)]);
	    (*ctoken)++;
	    return (0);
	} else
	    return (syntax_error (argv[*ctoken]));
}

/*  Get a float parameter from the argv array. */

static int getArgR (char **argv, int argc, int *ctoken, float *value) {
	if (*ctoken <= (argc-2)) {
	    *value = atof (argv[++(*ctoken)]);
	    (*ctoken)++;
	    return (0);
	} else
	    return (syntax_error (argv[*ctoken]));
}

static int syntax_error (char *msg) {
	printf ("ERROR  Syntax error: %s\n", msg);
	return (1);
}

