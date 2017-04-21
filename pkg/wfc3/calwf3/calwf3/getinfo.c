# include <stdio.h>
# include <stdlib.h>
# include <string.h>

# include "hstio.h"

# include "wf3.h"
# include "wf3corr.h"
# include "calwf3.h"
# include "hstcalerr.h"

/* This routine gets info from the primary header of the science file,
   and it calls routines to get calibration switches and reference file
   names.

   Warren Hack, 1998 May 26:
    Revised to work with ACS data.

   Howard Bushouse, 2000 Aug 23:
    Revised to work with WFC3 data.

   H. Bushouse, 2003 Oct 16:
    Changed gain values from int to float for WFC3.

   H. Bushouse, 2007 Feb 15:
    Changed default gain for IR channel from 2.0 to 2.5 in GetIRInfo.
    
   M. Sosey, Aug 2015: 
    Some memory cleanup
*/

int GetCCDInfo (WF3Info *wf3, CCD_Switch *sci_sw, RefFileInfo *sciref) {

/* arguments:
WF3Info *wf3          i: calibration flags and other info
CCD_Switch *sci_sw    o: all calibration switches (0 or 1) for science file
RefFileInfo *sciref  io: list of keyword,filename pairs for science file
*/

	extern int status;

	IODescPtr im;		/* descriptor for an image */
	Hdr phdr;		/* primary header */
	int nextend;		/* number of FITS extensions in rawfile */

	int GetKeyInt (Hdr *, char *, int, int, int *);
	int GetKeyFlt (Hdr *, char *, int, float, float *);
	int GetCCDSws (CCD_Switch *, Hdr *);
	int GetCCDRef (WF3Info *, CCD_Switch *, Hdr *, RefFileInfo *);
	
	/* Open input raw data file. */
	initHdr (&phdr);
	im = openInputImage (wf3->rawfile, "", 0);
	if (hstio_err()) {
	    sprintf (MsgText, "Member \"%s\" is not present", wf3->rawfile);
	    trlerror (MsgText);
        freeHdr(&phdr);
	    return (status = OPEN_FAILED);
	}

	/* Read primary header into pdhr. */
	getHeader (im, &phdr);
	if (hstio_err()) {
	    sprintf (MsgText, "Could not open PRIMARY header for \"%s\" ",
		     wf3->rawfile);
	    trlmessage (MsgText);
        closeImage(im);
        freeHdr(&phdr);
	    return (status = OPEN_FAILED);
	}
	closeImage (im);
	
	/* Get generic parameters: */

	/* Find out how many extensions there are in this file. */
	if (GetKeyInt (&phdr, "NEXTEND", USE_DEFAULT, EXT_PER_GROUP, &nextend)){
        freeHdr(&phdr);
        closeImage(im);
	    return (status = KEYWORD_MISSING);
    }
        
	wf3->nchips = nextend / EXT_PER_GROUP;

	/* Get binning and gain info.  We really only need this for the CCD. */
	if (GetKeyInt (&phdr, "BINAXIS1", USE_DEFAULT, 1, &wf3->scibin[0])) {
        closeImage(im);
        freeHdr(&phdr);
	    return (status = KEYWORD_MISSING);
    }
	if (GetKeyInt (&phdr, "BINAXIS2", USE_DEFAULT, 1, &wf3->scibin[1])){
        closeImage(im);
        freeHdr(&phdr);
	    return (status = KEYWORD_MISSING);
    }
	if (GetKeyFlt (&phdr, "CCDGAIN",  USE_DEFAULT, 1.5, &wf3->scigain)){
        closeImage(im);
        freeHdr(&phdr);
	    return (status = KEYWORD_MISSING);
    }
    
	wf3->samebin = 1;	/* default */

	/* Get calibration switches, and check that reference files exist. */
	if (GetCCDSws (sci_sw, &phdr))
	    return (status = KEYWORD_MISSING);
	if (GetCCDRef (wf3, sci_sw, &phdr, sciref))
	    return (status = CAL_FILE_MISSING);

	freeHdr (&phdr);
	return (status);
}

int GetIRInfo (WF3Info *wf3, IR_Switch *sci_sw, RefFileInfo *sciref) {

/* arguments:
WF3Info *wf3          i: calibration flags and other info
IR_Switch *sci_sw     o: all calibration switches (0 or 1) for science file
RefFileInfo *sciref  io: list of keyword,filename pairs for science file
*/

	extern int status;

	IODescPtr im;		/* descriptor for an image */
	Hdr phdr;		/* primary header */
	int nextend;		/* number of FITS extensions in rawfile */

	int GetKeyInt (Hdr *, char *, int, int, int *);
	int GetKeyFlt (Hdr *, char *, int, float, float *);
	int GetIRSws (IR_Switch *, Hdr *);
	int GetIRRef (WF3Info *, IR_Switch *, Hdr *, RefFileInfo *);
	
	/* Open input raw data file. */
	initHdr (&phdr);
	im = openInputImage (wf3->rawfile, "", 0);
	if (hstio_err()) {
	    sprintf (MsgText, "Member \"%s\" is not present", wf3->rawfile);
	    trlerror (MsgText);
	    return (status = OPEN_FAILED);
	}
    
	/* Read primary header into pdhr. */
	getHeader (im, &phdr);
	if (hstio_err())
	    return (status = OPEN_FAILED);
	
	/* Get generic parameters: */

	/* Find out how many extensions there are in this file. */
	if (GetKeyInt (&phdr, "NEXTEND", USE_DEFAULT, EXT_PER_GROUP, &nextend)){
	    return (status = KEYWORD_MISSING);
    }
	
    wf3->nchips = nextend / EXT_PER_GROUP;

	/* Get binning and gain info.  We really only need this for the CCD. */
	if (GetKeyInt (&phdr, "BINAXIS1", USE_DEFAULT, 1, &wf3->scibin[0])){
	    return (status = KEYWORD_MISSING);
    }
	
    if (GetKeyInt (&phdr, "BINAXIS2", USE_DEFAULT, 1, &wf3->scibin[1])){
	    return (status = KEYWORD_MISSING);
    }
    
	if (GetKeyFlt (&phdr, "CCDGAIN",  USE_DEFAULT, 2.5, &wf3->scigain)){
	    return (status = KEYWORD_MISSING);
    }
    
	wf3->samebin = 1;	/* default */

	/* Get calibration switches, and check that reference files exist. */
	if (GetIRSws (sci_sw, &phdr)) {
	    return (status = KEYWORD_MISSING);
    }
    
	if (GetIRRef (wf3, sci_sw, &phdr, sciref)){
	    return (status = CAL_FILE_MISSING);
    }

    closeImage(im);
	freeHdr (&phdr);
	return (status);
}

