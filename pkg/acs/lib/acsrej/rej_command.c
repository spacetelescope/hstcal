# include <stdio.h>
# include <stdlib.h>        /* strtof, strtol */
# include <string.h>
# include <errno.h>
# include <ctype.h>
# include <argp.h>

# include "acs.h"
# include "acsrej.h"
# include "acsversion.h"
# include "hstcalerr.h"
# include "hstcalversion.h"

void rej_reset (clpar *, int []);
static int is_valid_filename(char *token);
static int getArgS (const struct argp_state *state, char *arg, short *value);
static int getArgR (const struct argp_state *state, char *arg, float *value);
static int getArgT (const struct argp_state *state, char *arg, char *value);

/*
   Processes command-line parameters for AcsRej.

   Revision history:
   ----------------
   27 Aug 98  -  Adapted from CALSTIS2 cs2_command.c (WJ Hack)
   27 Jul 15  -  Dynamic size for input filenames (PL Lim)
   04 Apr 25  -  Replace argument parser with argp
*/

extern int status;

// Arguments that exit early
#define OPT_GITINFO 1000
#define OPT_VERSION 1001

// Arguments that are toggleable
#define OPT_SHADCORR 2000
#define OPT_CRMASK 2001
#define OPT_READNOISE_ONLY 2002
#define OPT_NEWBIAS 2003

// Arguments that require input
#define OPT_TABLE 3000
#define OPT_SCALE 3001
#define OPT_INIT 3002
#define OPT_SKY 3003
#define OPT_SIGMAS 3004
#define OPT_RADIUS 3005
#define OPT_THRESH 3006
#define OPT_PDQ 3007

// Adjust argument group precedence (0 is argp's default)
#define OPT_GROUP_INFO 0
#define OPT_GROUP_TOGGLE 0
#define OPT_GROUP_INPUT 0

static char doc[] = "ACS cosmic ray rejection";
static char args_doc[] = "INPUT OUTPUT";
static struct argp_option options[] = {
    {
        .name = "version",
        .key = OPT_VERSION,
        .arg = 0,
        .flags = 0,
        .doc = "print version",
        .group = OPT_GROUP_INFO
    },
    {
        .name = "gitinfo",
        .key = OPT_GITINFO,
        .arg = 0,
        .flags = 0,
        .doc = "print git version information",
        .group = OPT_GROUP_INFO
    },
    {
        .name = "verbose",
        .key = 'v',
        .arg = 0,
        .flags = 0,
        .doc = "turn on verbose mode",
        .group = OPT_GROUP_TOGGLE
    },
    {
        .name = "print-timestamps",
        .key = 't',
        .arg = 0,
        .flags = 0,
        .doc = "print timestamps",
        .group = OPT_GROUP_TOGGLE
    },
    {
        .name = "shadcorr",
        .key = OPT_SHADCORR,
        .arg = 0,
        .flags = 0,
        .doc = "turn on shutter shading correction (intended for CCD images only)",
        .group = OPT_GROUP_TOGGLE
    },
    {
        .name = "crmask",
        .key = OPT_CRMASK,
        .arg = 0,
        .flags = 0,
        .doc = "flag CR-rejected pixels in the input files",
        .group = OPT_GROUP_TOGGLE
    },
    {
        .name = "readnoise_only",
        .key = OPT_READNOISE_ONLY,
        .arg = 0,
        .flags = 0,
        .doc = "use read noise and not Poission noise to create ERR array (intended for BIAS images)",
        .group = OPT_GROUP_TOGGLE},
    {
        .name = "newbias",
        .key = OPT_NEWBIAS,
        .arg = 0,
        .flags = OPTION_HIDDEN,
        .doc = "DEPRECATED",
        .group = OPT_GROUP_TOGGLE
    },
    {
        .name = "table",
        .key = OPT_TABLE,
        .arg = "FILE",
        .flags = 0,
        .doc = "cosmic ray rejection table to use for default parameter values",
        .group = OPT_GROUP_INPUT
    },
    {
        .name = "scale",
        .key = OPT_SCALE,
        .arg = "NUM",
        .flags = 0,
        .doc = "multiplicative factor (in percent) applied to noise",
        .group = OPT_GROUP_INPUT
    },
    {
        .name = "init",
        .key = OPT_INIT,
        .arg = "MODE",
        .flags = 0,
        .doc = "scheme for computing initial guess image (median or minimum)",
        .group = OPT_GROUP_INPUT
    },
    {
        .name = "sky",
        .key = OPT_SKY,
        .arg = "MODE",
        .flags = 0,
        .doc = "scheme for computing sky levels to be subtracted (none:0.0 or mode:most frequently occurring value)",
        .group = OPT_GROUP_INPUT
    },
    {
        .name = "sigmas",
        .key = OPT_SIGMAS,
        .arg = "NUM[,...]",
        .flags = 0,
        .doc = "cosmic ray rejection thresholds, no. of thresholds are the no. of iterations",
        .group = OPT_GROUP_INPUT
    },
    {
        .name = "radius",
        .key = OPT_RADIUS,
        .arg = "NUM",
        .flags = 0,
        .doc = "radius (in pixels) to propagate the cosmic ray",
        .group = OPT_GROUP_INPUT
    },
    {
        .name = "thresh",
        .key = OPT_THRESH,
        .arg = "NUM",
        .flags = 0,
        .doc = "cosmic ray rejection propagation threshold",
        .group = OPT_GROUP_INPUT
    },
    {
        .name = "pdq",
        .key = OPT_PDQ,
        .arg = "NUM",
        .flags = 0,
        .doc = "data quality flag used for cosmic ray rejection",
        .group = OPT_GROUP_INPUT
    },
    {.name = 0, .key = 0, .arg = 0, .flags = 0, .doc = 0, .group = 0} // END OF ARGUMENTS
};

struct arguments {
    char *args[2];
    clpar *par;
    int *newpar;
};

static error_t parse_opt(int key, char *arg, struct argp_state *state) {
    struct arguments *arguments = state->input;
    switch (key) {
        case 'v':
            arguments->par->verbose = 1;
            break;
        case 't':
            arguments->par->printtime = 1;
            break;
        case OPT_VERSION:
            printf("%s\n", ACS_CAL_VER);
            exit(0);
        case OPT_GITINFO:
            printGitInfo();
            exit(0);
        case OPT_SHADCORR:
            arguments->par->shadcorr = 1;
            break;
        case OPT_CRMASK:
            arguments->par->mask = 1;
            arguments->newpar[CRMASK] = 1;
            arguments->newpar[TOTAL]++;
            break;
        case OPT_READNOISE_ONLY:
            arguments->par->readnoise_only = 1;
            break;
        case OPT_NEWBIAS:
            fprintf(stderr, "\n***************************************************************************************\n");
            fprintf(stderr, "    As of HSTCAL Version 2.0.0 (CALACS Version 10.0.0 HSTDP 2018_1) the 'newbias'\n ");
            fprintf(stderr, "    option for ACSREJ has been renamed 'readnoise_only' to clarify functionality.\n ");
            fprintf(stderr, "    Please update all programs which invoke ACSREJ with this option.\n");
            fprintf(stderr, "***************************************************************************************\n\n");
            return (status = INVALID_VALUE);
        case OPT_TABLE:
            if (!is_valid_filename(arg)) {
                status = INVALID_FILENAME;
                argp_failure(state, status, errno, "'%s' is not a valid file name\n", arg);
            }
            getArgT(state, arg, arguments->par->tbname);
            break;
        case OPT_SCALE:
            arguments->newpar[TOTAL]++;
            arguments->newpar[SCALENSE] = 1;
            getArgR(state, arg, &arguments->par->scalense);
            break;
        case OPT_INIT:
            arguments->newpar[TOTAL]++;
            arguments->newpar[INITGUES] = 1;
            getArgT(state, arg, arguments->par->initgues);
            break;
        case OPT_SKY:
            arguments->newpar[TOTAL]++;
            arguments->newpar[INITGUES] = 1;
            getArgT(state, arg, arguments->par->sky);
            break;
        case OPT_SIGMAS:
            arguments->newpar[TOTAL]++;
            arguments->newpar[CRSIGMAS] = 1;
            getArgT(state, arg, arguments->par->sigmas);
            break;
        case OPT_RADIUS:
            arguments->newpar[TOTAL]++;
            arguments->newpar[CRRADIUS] = 1;
            getArgR(state, arg, &arguments->par->radius);
            break;
        case OPT_THRESH:
            arguments->newpar[TOTAL]++;
            arguments->newpar[CRTHRESH] = 1;
            arguments->par->thresh = strtof(arg, NULL);
            getArgR(state, arg, &arguments->par->thresh);
            break;
        case OPT_PDQ:
            arguments->newpar[TOTAL]++;
            arguments->newpar[BADINPDQ] = 1;
            getArgS(state, arg, &arguments->par->badinpdq);
            break;
        case ARGP_KEY_ARG:
            if (state->arg_num >= 2) {
                argp_usage(state);
            }
            if (!is_valid_filename(arg)) {
                status = INVALID_FILENAME;
                argp_failure(state, status, errno, "'%s' is not a valid file name\n", arg);
            }
            arguments->args[state->arg_num] = arg;
            break;
        case ARGP_KEY_END:
            if (state->arg_num < 2) {
                argp_usage(state);
            }
            break;
        default:
            return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

int rej_command (int argc, char **argv, char **input, char *output, clpar *par, int newpar[]) {
    /* arguments
    int argc;           i: input command-line parameters
    char **argv;
    char **input;       o: input file name or file list
    char *output;       o: output file name
    clpar *par;         o: user specified parameters
    int newpar[];       o: array of parameters set by the user
    */

    extern int status;

    /* reset the parameters */
    rej_reset (par, newpar);

    struct argp argp = {options, parse_opt, args_doc, doc, 0, 0, 0};
    struct arguments arguments = {
        .par = par,
        .newpar = newpar,
    };

    argp_parse(&argp, argc, argv, 0, 0, &arguments);
    *input = strdup(arguments.args[0]);
    if (!*input) {
        perror("unable to allocate bytes for input file name");
        return (status = OUT_OF_MEMORY);
    }
    output = strcpy(output, arguments.args[1]);

    return (0);
}

/* ------------------------------------------------------------------*/
/*                          is_valid_filename                        */
/* ------------------------------------------------------------------*/

/*  Determine the worthiness of a file name */

static int is_valid_filename(char *token) {
    char *token_p = token;
    if (!strlen(token_p)) {
        // empty string
        return 0;
    }

    // Check for invalid name
    int blanks = 0;
    int alnums = 0;
    while (*token_p != '\0') {
        if (isblank(*token_p)) {
            blanks++;
        } else if (isalnum(*token_p)) {
            alnums++;
        }
        token_p++;
    }

    if (blanks && !alnums) {
        // string consists of blank characters
        return 0;
    }

    return 1;
}

/* ------------------------------------------------------------------*/
/*                          getArgS                                  */
/* ------------------------------------------------------------------*/

/*  Get a short parameter from the argv array. */

static int getArgS (const struct argp_state *state, char *arg, short *value) {
    errno = 0;
    char *end = NULL;

    *value = (short) strtol(arg, &end, 10);
    if (((errno == EINVAL && *value == 0) || (end == arg && *value == 0)) || errno == ERANGE) {
        status = INVALID_VALUE;
        argp_failure(state, status, errno, "a valid integer is required: '%s'", arg);
    }
    return (status = ACS_OK);
}

/* ------------------------------------------------------------------*/
/*                          getArgT                                  */
/* ------------------------------------------------------------------*/

/*  Get a string parameter from the argv array. */

static int getArgT (const struct argp_state *state, char *arg, char *value) {

    if (!arg) {
        status = INVALID_VALUE;
        argp_failure(state, status, errno, "missing string argument");
    }
    strcpy(value, arg);
    return (status = ACS_OK);
}

/* ------------------------------------------------------------------*/
/*                          getArgR                                  */
/* ------------------------------------------------------------------*/


/*  Get a float parameter from the argv array. */

static int getArgR (const struct argp_state *state, char *arg, float *value) {
    char *end = NULL;
    errno = 0;
    *value = strtof(arg, &end);
    if ((end == arg && *value == 0.0f) || errno == ERANGE) {
        status = INVALID_VALUE;
        argp_failure(state, status, errno, "requires a valid integer argument, got '%s'", arg);
    }
    return (status = ACS_OK);
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
    par->initgues[0] = '\0';
    par->sigmas[0] = '\0';
    par->sky[0] = '\0';
    par->verbose = 0;
    par->printtime = 0;
    par->shadcorr = 0;
    par->readnoise_only = 0;
    par->badinpdq = 0;

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
