# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include "hstio.h"
# include <time.h>
# include "acs.h"
# include "acscorr.h"
# include "calacs.h"
# include "hstcalerr.h"

/* This routine gets info from the primary header of the science file,
   and it calls routines to get calibration switches and reference file
   names.

   Warren Hack, 1998 May 26:
    Revised to work with ACS data.
*/

int GetSciInfo (ACSInfo *acs, CalSwitch *sci_sw, RefFileInfo *sciref) {

/* arguments:
ACSInfo *acs         i: calibration flags and other info
CalSwitch *sci_sw     o: all calibration switches (0 or 1) for science file
RefFileInfo *sciref  io: list of keyword,filename pairs for science file
*/

	extern int status;

	IODescPtr im;		/* descriptor for an image */
	Hdr phdr;		/* primary header */
	int nextend;		/* number of FITS extensions in rawfile */

    time_t date,date_limit;
    char dateobs[ACS_CBUF],targname[ACS_LINE];

    int GetKeyStr (Hdr *, char *, int, char *, char *, int);
	int GetKeyInt (Hdr *, char *, int, int, int *);
	int GetFlags (CalSwitch *, Hdr *);
	int SciFlags (ACSInfo *, CalSwitch *, Hdr *, RefFileInfo *);
        int parseObsDateVal (char *dateobs, time_t *date);

	/* Read primary header of rawfile into phdr. */
	initHdr (&phdr);
	im = openInputImage (acs->rawfile, "", 0);
	if (hstio_err()) {
		sprintf (MsgText, "Member \"%s\" is not present", acs->rawfile);
		trlerror (MsgText);
	    return (status = OPEN_FAILED);
	}
	getHeader (im, &phdr);		/* get primary header */
	if (hstio_err()) {
		sprintf (MsgText, "Could not open PRIMARY header for \"%s\" ", acs->rawfile);
		trlerror (MsgText);
        return (status = OPEN_FAILED);
    }
	closeImage (im);

	/* Get generic parameters. */

	/* Find out how many extensions there are in this file. */
	if (GetKeyInt (&phdr, "NEXTEND", USE_DEFAULT, EXT_PER_GROUP, &nextend))
	    return (status);
	acs->nchips = nextend / EXT_PER_GROUP;

	/* Get binning and gain info.  We really only need this for the CCD. */
	if (GetKeyInt (&phdr, "BINAXIS1", USE_DEFAULT, 1, &acs->scibin[0]))
	    return (status);
	if (GetKeyInt (&phdr, "BINAXIS2", USE_DEFAULT, 1, &acs->scibin[1]))
	    return (status);
    /*
	if (GetKeyInt (&phdr, "CCDGAIN", USE_DEFAULT, 1, &acs->scigain))
	    return (status);
    */
	acs->samebin = 1;	/* default */

    /* Get and parse DATE-OBS into a floating point value*/
    if (GetKeyStr (&phdr, "DATE-OBS", USE_DEFAULT, "", dateobs, ACS_CBUF))
        	    return (status);
    parseObsDateVal(dateobs, &date);
    /* Parse date used to check whether ACS HRC/WFC data was
       pre-SM4 or post-SM4 and turn into a float for comparison
       with exposure's date-obs
    */
    parseObsDateVal("2009-01-01", &date_limit);
    if (GetKeyStr (&phdr, "TARGNAME", USE_DEFAULT, "", targname, ACS_LINE))
        return (status);
    if (strncmp(targname,"BIAS",4) == 0 && date > date_limit){
        acs->newbias = 1;
    }

	/* Get calibration switches, and check that reference files exist. */
	if (GetFlags (sci_sw, &phdr))
	    return (status);
	if (SciFlags (acs, sci_sw, &phdr, sciref))
	    return (status);

	freeHdr (&phdr);
	return (status);
}
