#ifndef INCL_WF3OMIT_H
#define INCL_WF3OMIT_H

/* These are the calibration steps that by default are not performed.
   This file is included only by lib/defswitch.c.
*/

# define WF3_N_OMIT   2

static char omitsw[WF3_N_OMIT][9] = {"atodcorr", "shadcorr"};


#endif /* INCL_WF3OMIT_H */
