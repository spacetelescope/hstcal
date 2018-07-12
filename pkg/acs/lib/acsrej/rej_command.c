# include <stdio.h>
# include <stdlib.h>        /* atoi, atof, strtod, strtol */
# include <string.h>

# include "xtables.h"

# include "acs.h"
# include "acsrej.h"
# include "hstcalerr.h"

void rej_reset (clpar *, int []);
static int syntax_error (char *);
static int getArgS (char **, int, int *, short *);
static int getArgR (char **, int, int *, float *);
static int getArgT (char **, int, int *, char *);

/*
   Processes command-line parameters for AcsRej.

   Revision history:
   ----------------
   27 Aug 98  -  Adapted from CALSTIS2 cs2_command.c (WJ Hack)
   27 Jul 15  -  Dynamic size for input filenames (PL Lim)
*/

int rej_command (int argc, char **argv, char **input, char *output,
                clpar *par, int newpar[])
{

/* arguments
int argc;           i: input command-line parameters
char **argv;
char **input;       o: input file name or file list
char *output;       o: output file name
clpar *par;         o: user specified parameters
int newpar[];       o: array of parameters set by the user
*/

    extern int status;

    int ctoken;     /* current command-line token being processed */

    /* reset the parameters */
    rej_reset (par, newpar);

    /* not enough arguments
        List all options for user
    */
    if (argc < 3) {
        syntax_error ("acsrej.e input output [-t] [-v]");
        printf("                    [-shadcorr] [-crmask] [-readnoise_only] \n");
        printf("                    [-table <filename>] [-scale #]\n");
        printf("                    [-init (med|min)] [-sky (none|mode)]\n");
        printf("                    [-sigmas #] [-radius #] [-thresh #]\n");
        printf("                    [-pdq #]\n\n");
		printf("input             comma-delimited string (e.g., filename or filename1,filename2,filename3)\n");
		printf("output            string\n");
		printf("-t                print timestamps\n");
		printf("-v                turn on verbose mode\n");
		printf("-shadcorr         turn on shutter shading correction (intended for CCD images only)\n");
		printf("-crmask           flag CR-rejected pixels in the input files\n");
		printf("-readnoise_only   use read noise and not Poission noise to create ERR array (intended for BIAS images)\n");
		printf("-table filename   cosmic ray rejection table to use for default parameter values\n");
		printf("-scale #          multiplicative factor (in percent) applied to noise\n");
		printf("-init med|min     scheme for computing initial guess image (median or minimum)\n");
		printf("-sky none|mode    scheme for computing sky levels to be subtracted (none:0.0 or mode:most frequently occurring value)\n");
		printf("-sigmas #[,#...]  cosmic ray rejection thresholds, no. of thresholds are the no. of iterations\n");
		printf("-radius #         radius (in pixels) to propagate the cosmic ray\n");
		printf("-thresh #         cosmic ray rejection propagation threshold\n");
		printf("-pdq #            data quality flag used for cosmic ray rejection\n\n");
        return(status);
    }

    /* Get names of input and output files. These are mandatory.
      Input list is also limited by tables/c_imt.c */
    if ((*input = calloc(strlen(argv[1]) + 1, sizeof(char))) == NULL) {
        printf ("Can't even begin; out of memory.\n");
        return (status = OUT_OF_MEMORY);
    }
    strcpy (*input, argv[1]);
    strcpy (output, argv[2]);

    /* Scan remaining parameters. */
    if (argc > 3) {
        ctoken = 3;
        while (ctoken < argc) {

            /* switch must begin with a "-" */
            if (argv[ctoken][0] != '-') {
                return(syntax_error (argv[ctoken]));
            } else {

            /* These do not require additional arguments. */

            if (strcmp("t", argv[ctoken]+1) == 0) {
                par->printtime = 1;
                ctoken++;

            } else if (strcmp("v", argv[ctoken]+1) == 0) {
                par->verbose = 1;
                ctoken++;

            } else if (strcmp("shadcorr", argv[ctoken]+1) == 0) {
                par->shadcorr = 1;
                ctoken++;

            // This option/variable was previously named "newbias", and it has been renamed
            // to "readnoise_only" which is more informative.  The variable controls
            // the type of noise used in the computation of the error (ERR) array.  In particular,
            // this variable is intended to control the noise used when processing BIAS data.  When
            // using CALACS, readnoise_only will only be set to "1" (true) for BIAS data.  When using
            // ACSREJ as a standalone utility, the readnoise_only can be set to "1" for any data.
            } else if (strcmp("newbias", argv[ctoken]+1) == 0) {
                printf ("\n***************************************************************************************\n");
				printf ("     As of HSTCAL Version 2.0.0 (CALACS Version 10.0.0 HSTDP 2018_1) the 'newbias'\n ");
                printf ("    option for ACSREJ has been renamed 'readnoise_only' to clarify functionality.\n ");
                printf ("    Please update all programs which invoke ACSREJ with this option.\n");
                printf ("***************************************************************************************\n\n");
                return (status = INVALID_VALUE);

            } else if (strcmp("readnoise_only", argv[ctoken]+1) == 0) {
                par->readnoise_only = 1;
                ctoken++;

            } else if (strcmp("crmask", argv[ctoken]+1) == 0) {
                par->mask = 1;
                newpar[CRMASK] = 1;
                newpar[TOTAL]++;
                ctoken++;

            /* These require one additional argument, which is
               handled by the getArg functions.
            */

            } else if (strcmp("table", argv[ctoken]+1) == 0) {
                if (getArgT (argv, argc, &ctoken, par->tbname))
                    return (status = INVALID_VALUE);

            } else if (strcmp("scale", argv[ctoken]+1) == 0) {
                newpar[TOTAL]++;
                newpar[SCALENSE] = 1;
                if (getArgR (argv, argc, &ctoken, &par->scalense))
                    return (status = INVALID_VALUE);

            } else if (strcmp("init", argv[ctoken]+1) == 0) {
                newpar[TOTAL]++;
                newpar[INITGUES] = 1;
                if (getArgT (argv, argc, &ctoken, par->initgues))
                    return (status = INVALID_VALUE);

            } else if (strcmp("sky", argv[ctoken]+1) == 0) {
                newpar[TOTAL]++;
                newpar[SKYSUB] = 1;
                if (getArgT (argv, argc, &ctoken, par->sky))
                    return (status = INVALID_VALUE);

            } else if (strcmp("sigmas", argv[ctoken]+1) == 0) {
                newpar[TOTAL]++;
                newpar[CRSIGMAS] = 1;
                if (getArgT (argv, argc, &ctoken, par->sigmas))
                    return (status = INVALID_VALUE);

            } else if (strcmp("radius", argv[ctoken]+1) == 0) {
                newpar[TOTAL]++;
                newpar[CRRADIUS] = 1;
                if (getArgR (argv, argc, &ctoken, &par->radius))
                    return (status = INVALID_VALUE);

            } else if (strcmp("thresh", argv[ctoken]+1) == 0) {
                newpar[TOTAL]++;
                newpar[CRTHRESH] = 1;
                if (getArgR (argv, argc, &ctoken, &par->thresh))
                    return (status = INVALID_VALUE);

            } else if (strcmp("pdq", argv[ctoken]+1) == 0) {
                newpar[TOTAL]++;
                newpar[BADINPDQ] = 1;
                if (getArgS (argv, argc, &ctoken, &par->badinpdq))
                    return (status = INVALID_VALUE);

            /* No match. */
            } else
                return (syntax_error (argv[ctoken]));
            }
        }
    }

    return (0);
}


/* ------------------------------------------------------------------*/
/*                          getArgS                                  */
/* ------------------------------------------------------------------*/


/*  Get a short parameter from the argv array. */

static int getArgS (char **argv, int argc, int *ctoken, short *value) {
    extern int status;

    if (*ctoken <= (argc-2)) {
        *value = (short) strtol(argv[++(*ctoken)], (char **)NULL, 10);
        (*ctoken)++;
        return (status = ACS_OK);
    } else
        return (syntax_error (argv[*ctoken]));
}

/* ------------------------------------------------------------------*/
/*                          getArgT                                  */
/* ------------------------------------------------------------------*/


/*  Get a string parameter from the argv array. */

static int getArgT (char **argv, int argc, int *ctoken, char *value) {
    extern int status;

    if (*ctoken <= (argc-2)) {
        strcpy (value, argv[++(*ctoken)]);
        (*ctoken)++;
        return (status = ACS_OK);
    } else
        return (syntax_error (argv[*ctoken]));
}

/* ------------------------------------------------------------------*/
/*                          getArgR                                  */
/* ------------------------------------------------------------------*/


/*  Get a float parameter from the argv array. */

static int getArgR (char **argv, int argc, int *ctoken, float *value) {
    extern int status;

    if (*ctoken <= (argc-2)) {
        *value = strtod (argv[++(*ctoken)], (char **)NULL);
        (*ctoken)++;
        return (status = ACS_OK);
    } else
        return (syntax_error (argv[*ctoken]));
}

/* ------------------------------------------------------------------*/
/*                          syntax_error                             */
/* ------------------------------------------------------------------*/


static int syntax_error (char *msg) {
    extern int status;

    printf ("Syntax error: %s\n", msg);
    return (status = ERROR_RETURN);
}

/* ------------------------------------------------------------------*/
/*                          rej_reset                                */
/* ------------------------------------------------------------------*/


void rej_reset (clpar *par, int newpar[])
{

/* arguments
clpar *par;       o: user specified parameters
int newpar[];     o: array of parameters set by the user
*/

    /* First, set parameters specifiable in the command line to be
        non-existent. */
    par->tbname[0] = '\0';
    par->verbose = 0;
    par->printtime = 0;
    par->shadcorr = 0;
    par->readnoise_only = 0;

    newpar[TOTAL] = 0;
    newpar[SCALENSE] = 0;
    newpar[INITGUES] = 0;
    newpar[SKYSUB] = 0;
    newpar[CRSIGMAS] = 0;
    newpar[CRRADIUS] = 0;
    newpar[CRTHRESH] = 0;
    newpar[BADINPDQ] = 0;
    newpar[CRMASK] = 0;
}
