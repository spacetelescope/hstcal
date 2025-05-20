# include <stdio.h>
# include <string.h>

# include "hstcal.h"
# include "hstio.h"
# include "wf3.h"
# include "wf3dq.h"
# include "wf3info.h"
# include "hstcalerr.h"

/*
  Michele De La Pena, 2023 October 16

*/

/* Because the pixels representing the serial virtual overscan are present
   in the saturation image, it is necessary to skip over these pixels when
   determining the flagging for full-well saturation.

   This routine is based upon code in doblev.c.
*/
void
ComputeLimits(WF3Info *wf3, int xdim, int ydim, int begx[2], int begy[2], int endx[2], int endy[2]) {
    int amp;			/* Counter for amps used */
    int numamps;		/* Number of AMPS used to readout chip */
    char *ccdamp;		/* Amps which are used for this chip */
    char *amploc;
    int trimx1, trimx2;	/* width to trim off ends of each line */
    int trimx3, trimx4;	/* width to trim off middle of each line */
    int trimy1, trimy2;	/* amount to trim off ends of each column */
    int bias_loc;
    int bias_ampx, bias_ampy;
    int bias_orderx[4] = {0,1,0,1};
    int bias_ordery[4] = {0,0,1,1};

    void parseWFCamps (char *, int, char *);

    /* Copy out overscan size information for ease of reference in this function */
    trimx1 = wf3->trimx[0];	/* Left serial physical */
    trimx2 = wf3->trimx[1];	/* Right serial physical */
    trimx3 = wf3->trimx[2]; /* Left serial virtual */
    trimx4 = wf3->trimx[3]; /* Right serial virtual */
    trimy1 = wf3->trimy[0]; /* Bottom (Chip 1) parallel virtual */
    trimy2 = wf3->trimy[1]; /* Top (Chip 2) parallel virtual */

    /* Establish which amps from ccdamp string are appropriate
    ** for this chip */
    ccdamp = (char *) calloc (NAMPS+1, sizeof(char));
    ccdamp[0] = '\0';

    /* Set up the 2 AMPS used per line */
    if (wf3->detector == CCD_DETECTOR) {
        parseWFCamps (wf3->ccdamp, wf3->chip, ccdamp);
    }
    sprintf(MsgText,"Amp: %s chip: %d ccdamp: %s ", wf3->ccdamp, wf3->chip, ccdamp);
    trlmessage(MsgText);

    /* How many amps are used for this chip */
    numamps = strlen(ccdamp);
    sprintf(MsgText,"Number of amps %d", numamps);
    trlmessage(MsgText);

    for (amp = 0; amp < numamps; amp++) {
        bias_loc = 0;

        /* determine which section goes with the amp */
        /* bias_amp = 0 for BIASSECTA, = 1 for BIASSECTB */
        amploc = strchr (AMPSORDER, ccdamp[amp]);
        bias_loc = *amploc - *ccdamp;
        bias_ampx = bias_orderx[bias_loc];
        bias_ampy = bias_ordery[bias_loc];

        /* Compute range of pixels affected by each amp */
        begx[amp] = trimx1 + (wf3->ampx + trimx3 + trimx4) * bias_ampx;
        endx[amp] = (bias_ampx == 0 && wf3->ampx != 0) ? wf3->ampx + trimx1 :
                    xdim - trimx2;
        begy[amp] = trimy1 + wf3->ampy* bias_ampy;
        endy[amp] = (bias_ampy == 0 && wf3->ampy != 0) ? wf3->ampy +trimy1 :
                    ydim - trimy2;

        /* Make sure that endx and endy do not extend beyond the bounds of
        ** the image */
        if (endx[amp] > xdim) endx[amp] = xdim;
        if (endy[amp] > ydim) endy[amp] = ydim;
    }

    free(ccdamp);
}
