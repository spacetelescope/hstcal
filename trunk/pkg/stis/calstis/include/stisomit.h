/* These are the calibration steps that by default are not performed.
   This file is included only by lib/defswitch.c.
*/

# define STIS_N_OMIT   3

static char omitsw[STIS_N_OMIT][9] = {"atodcorr", "shadcorr", "sgeocorr"};
