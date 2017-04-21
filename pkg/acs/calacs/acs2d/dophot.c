# include <ctype.h>
# include <stdio.h>
# include <stdlib.h>	/* malloc */
# include <string.h>    /* strchr */

# include "hstio.h"
# include "acs.h"
# include "acsinfo.h"
# include "err.h"

# define NPHOT      4		/* size of phot returned by c_phopar */
# define ARR_SIZE    10000    /* size of arrays used for throughputs */


/*static void MultFilter (PhotInfo *); */
static void Phot2Obs (char *, char *);

/* This routine gets the absolute flux conversion from PHOTTAB and the
 aperture throughput (filter throughput, really) from APERTAB.  These
 are multiplied together, and then the synphot routine phopar is called
 to determine the inverse sensitivity, reference magnitude (actually a
 constant), pivot wavelength, and RMS bandwidth.  These are written to
 keywords in the primary header.

 Warren Hack, 1998 June 12:
  Initial ACS version.
  Changes were made only to GetPhot, except for include file and
    variable names (sts -> acs2d, StisInfo1 -> ACSInfo).
 Warren Hack, 2002 Jan 30:
  Removed gain normalization of PHOTFLAM keyword value to be consistent
    with image datatype after flat-field calibration. Also, increased
    ARR_SIZE to match the default value used by SYNPHOT.
 Matt Davis, 2011 Apr 5:
  Update for removing synphot dependence, to be replaced by precomputed
    values in the IMPHTTAB table.
  Removed obsolete variables and function prototypes.
  Switched to use of getphttab functions.
 */

int doPhot (ACSInfo *acs2d, SingleGroup *x) {

  /* arguments:
   ACSInfo *acs2d     i: calibration switches, etc
   SingleGroup *x    io: image to be calibrated; primary header is modified
   */

  extern int status;

  PhotPar obs;

  char photmode[ACS_LINE],obsmode[ACS_LINE];

  /* function prototypes from lib */
  int GetKeyStr (Hdr *, char *, int, char *, char *, int);
  int PutKeyFlt (Hdr *, char *, float, char *);

  if (acs2d->photcorr == DUMMY)
    return (status);

  /* Extract photmode from sci extension header*/
  if (GetKeyStr(&x->sci.hdr,"PHOTMODE",USE_DEFAULT,"",photmode, ACS_LINE))
    return (status);

  /* Add commas to PHOTMODE string and change all strings to
   lower case for use in synphot. */
  Phot2Obs(photmode,obsmode);
  if (acs2d->verbose){
    sprintf(MsgText,"Created SYNPHOT obsmode of: %s",obsmode);
    trlmessage(MsgText);
  }

  /* Initialize PhotPar struct */
  InitPhotPar(&obs, acs2d->phot.name, acs2d->phot.pedigree);

  /* get phot values */
  if (GetPhotTab(&obs,obsmode)) {
    trlerror("Error return from GetPhotTab.");
    return status;
  }

  if (acs2d->verbose){
    sprintf(MsgText,"Computed PHOTFLAM value of %g",obs.photflam);
    trlmessage(MsgText);
  }

  /* Update the photometry keyword values in the SCI header.
   Revised 10 March 1999 WJH
   */

  if (PutKeyFlt (&x->sci.hdr, "PHOTFLAM", obs.photflam, "inverse sensitivity"))
    return (status);

  if (PutKeyFlt (&x->sci.hdr, "PHOTZPT", obs.photzpt, "zero point"))
    return (status);

  if (PutKeyFlt (&x->sci.hdr, "PHOTPLAM", obs.photplam, "pivot wavelength"))
    return (status);

  if (PutKeyFlt (&x->sci.hdr, "PHOTBW", obs.photbw, "RMS bandwidth"))
    return (status);

  FreePhotPar(&obs);

  /* check whether GetPhotTab returned -9999, which means it was asked
   to do extrapolation */
  if (obs.photflam == -9999) {
    trlwarn("IMPHTTAB extrapolation not supported, PHOTCORR skipped.");
    return (status = CAL_STEP_NOT_DONE);
  }

  return (status);
}

/*
 This function converts the PHOTMODE string into an OBSMODE
 string suitable for use with synphot functions.
 PHOTMODE - all upper case component names separated by blanks
 OBSMODE - all lower case names separated by commas
 NOTE: The operation occurs **in-place** on the photmode string.
 */
static void Phot2Obs(char *photmode, char *obsmode)
{
  char blank = 32, comma = 44;
  int i, len, c1;

  len = strlen(photmode);

  for (i=0; i < len; i++) {
    c1 = photmode[i];
    if (c1 == blank) {
      obsmode[i] = comma;
    } else
      obsmode[i] = tolower(c1);
  }
  obsmode[len] = '\0';
}
