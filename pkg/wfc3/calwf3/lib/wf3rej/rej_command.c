# include <stdio.h>
# include <stdlib.h>        /* atoi, atof, strtod, strtol */
# include <string.h>

# include "xtables.h"

# include "wf3.h"
# include "wf3rej.h"
# include "hstcalerr.h"
# include "wf3version.h"

void rej_reset (clpar *, int []); 
static int syntax_error (char *);
static int getArgS (char **, int, int *, short *);
static int getArgR (char **, int, int *, float *);
static int getArgT (char **, int, int *, char *);

/*  
   Processes command-line parameters for Wf3Rej.

   Revision history:
   ----------------
   27 Aug 98  -  Adapted from CALSTIS2 cs2_command.c (WJ Hack)
   
   31  Jul 2015 - propagated memory fix from calacs (MLS)
   
*/

int rej_command (int argc, char **argv, char **input, char *output,
                clpar *par, int newpar[]) {

/* arguments
int argc;           i: number of input command-line parameters
char **argv;	    i: input command-line parameters
char **input;        o: input file name or file list, not asn
char *output;       o: output file name
clpar *par;         o: user specified parameters
int newpar[];       o: array of parameters set by the user
*/

    extern int status;
    
    int ctoken=0;     /* current command-line token being processed */

    /* Reset the parameters */
    rej_reset (par, newpar); 

    
    /* Not enough arguments or just wants version number; List all options for user */
    if (argc < 3  ) {
        if (argc >= 2){
            if (strcmp("-r", argv[1]) == 0) { 
                printf ("Current version: %s\n", WF3_CAL_VER);
                exit(0);
            }
            if (argv[ctoken][2] == '-') {
                printf("You need to speficy an output image name");
                exit(0);
            }
        } else {
            printf(" wf3rej input output [-t] [-v]");
            printf("                    [-shadcorr] [-crmask]\n ");
            printf("                    [-table <filename>] [-scale #]\n");
            printf("                    [-init (med|min)] [-sky (none|mode)]\n");
            printf("                    [-sigmas #] [-radius #] [-thresh #]\n");
            printf("                    [-pdq #]\n"); 
            printf("                    [-r]\n");
            exit(ERROR_RETURN);
        }
    }       
    /* Get names of input and output files. These are mandatory. */
    if ((*input = calloc(strlen(argv[1]) + 1, sizeof(char)))== NULL){
        printf("Cannot allocate memory for input\n");
        return(status=OUT_OF_MEMORY);
    }
    
    strcpy (*input,  argv[1]);
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
        return (status = WF3_OK);
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
        return (status = WF3_OK);
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
        return (status = WF3_OK);
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


void rej_reset (clpar *par, int newpar[]) {

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
