# include <stdio.h>
#include "hstcal.h"
# include "hstio.h"
# include "wf3.h"
# include "wf3corr.h"
# include "calwf3.h"

static void CCDSanity (int, char *);
static void IRSanity  (int, char *);

/* This routine gets the names of all the reference files (images and
   tables) that are needed for calibrating the science file, based on
   the switches.  The files are not checked at this time to see whether
   they actually exist.

   Some sanity checking is done, and warning messages are printed if
   switch values don't make sense for the current detector.

   Various flags are set (e.g. wf3->basic_2d), depending on which calwf3
   module will perform the calibration.  Basic 2-D reduction, derived from
   'calstis1', is divided into two parts: wf3ccd and wf32d.
   The steps dqicorr, atodcorr, flashcorr, and blevcorr constitute the
   wf3ccd part, and these are included with the sci_basic_ccd flag
   rather than with sci_basic_2d which is used for the wf32d steps.

	
   Warren Hack, 1998 May 26:
    Initial version for use with ACS; adapted from CALSTIS1

   Howard Bushouse, 2000 Aug 24:
    Adapted for use with WFC3
  
   H.Bushouse, 2001 May 8:
    Revised to add Post-Flash processing step.
 
   H.Bushouse, 2001 Nov 16:
    Replaced APERTAB and PHOTTAB with GRAPHTAB and COMPTAB.
 
   H.Bushouse, 2002 May 21:
    Removed call to CCDSanity for BLEVCORR step in GetCCDRef because that
    is a valid step for IR processing. Added call to IRSanity for NLINCORR
    step in GetIRRef because that is not a valid step for UVIS processing.
 
   H.Bushouse, 2003 Apr 25:
    Modified GetCCDRef so that CRREJTAB reference table is loaded if
    either CRCORR or RPTCORR is enabled, because it will be used by both
    processes for WFC3. Also modified GetIRRef to load CRREJTAB if
    RPTCORR is enabled.
  
   H.Bushouse, 2010 Sep 30:
    Modified GetIRRef to always check CCDTAB and OSCNTAB, because their
    contents are used in so many steps, and also add checks for zoffcorr,
    unitcorr, and crcorr switches. (PR #66081; Trac #608)

   H.Bushouse, 2011 Sep 7:
    Updated GetCCDRef and GetIRRef to use new IMPHTTAB instead of synphot 
    GRAPHTAB and CoMPTAB.
    
   M. sosey August 2014  PCETAB,  DRKCFILE, BIACFILE, SNKCFILE added for PCTECORR step
   and DQICORR step

*/

int GetCCDRef (WF3Info *wf3, CCD_Switch *sci_sw, Hdr *phdr,
	       RefFileInfo *sciref) {

/* arguments:
WF3Info *wf3          i: calibration flags and other info
CCD_Switch *sci_sw    i: all calibration switches for science file
Hdr *phdr             i: primary header of science file
RefFileInfo *sciref  io: list of keyword,filename pairs
*/

	extern int status;

	int refimage_used = 0;	/* = 1 if do bias, dark, flat, or shadcorr */

	int GetNewRef (Hdr *, char *, RefFileInfo *);

	/* */
	if (GetNewRef (phdr, "CCDTAB", sciref))
	    return (status);

    if (sci_sw->pctecorr == PERFORM) {
        
        CCDSanity (wf3->detector, "PCTECORR");
        wf3->sci_basic_cte = PERFORM;
        
        if (GetNewRef (phdr, "PCTETAB", sciref))
            return (status);
            
        if (GetNewRef (phdr, "DRKCFILE",  sciref))
            return (status);
        
        if (GetNewRef (phdr, "BIACFILE", sciref))
            return (status);
        
    }
        

	if (sci_sw->dqicorr == PERFORM) {
	    wf3->sci_basic_ccd = PERFORM;
        if (GetNewRef (phdr, "SNKCFILE", sciref))
            return (status);
	    if (GetNewRef (phdr, "BPIXTAB", sciref))
		    return (status);
	}

	if (sci_sw->crcorr == PERFORM) {
	    wf3->sci_crcorr = PERFORM;
	    if (wf3->nimages < 2) {
		trlwarn ("CRCORR will be omitted because there's only one image.");
		wf3->sci_crcorr = OMIT;
	    }
	    /* reference table is checked below, after getting rptcorr */
	}

	if (sci_sw->rptcorr == PERFORM) {
	    if (wf3->nimages < 2) {
		trlwarn 
		   ("RPTCORR will be omitted because there's only one image.");
		wf3->sci_rptcorr = OMIT;
	    } else {
		wf3->sci_rptcorr = PERFORM;
	    }
	}

	if (wf3->sci_crcorr == PERFORM || wf3->sci_rptcorr == PERFORM) {
	    if (GetNewRef (phdr, "CRREJTAB", sciref))
		    return (status);
	}


	/* Note that we set sci_basic_ccd rather than sci_basic_2d. */
	if (sci_sw->atodcorr == PERFORM) {
	    wf3->sci_basic_ccd = PERFORM;
	    CCDSanity (wf3->detector, "ATODCORR");
	    if (GetNewRef (phdr, "ATODTAB", sciref))
		    return (status);
	}

	if (sci_sw->biascorr == PERFORM) {
	    wf3->sci_basic_ccd = PERFORM;   refimage_used = 1;
	    CCDSanity (wf3->detector, "BIASCORR");
	    if (GetNewRef (phdr, "BIASFILE", sciref))
		    return (status);
	}

	if (sci_sw->flashcorr == PERFORM) {
	    wf3->sci_basic_ccd = PERFORM;   refimage_used = 1;
	    CCDSanity (wf3->detector, "FLSHCORR");
	    if (GetNewRef (phdr, "FLSHFILE", sciref))
		    return (status);
	}

	if (sci_sw->darkcorr == PERFORM) {
	    wf3->sci_basic_2d = PERFORM;   refimage_used = 1;
	    if (GetNewRef (phdr, "DARKFILE", sciref))
		    return (status);
	}

	if (sci_sw->flatcorr == PERFORM) {
	    wf3->sci_basic_2d = PERFORM;   refimage_used = 1;
	    if (GetNewRef (phdr, "PFLTFILE", sciref))
		    return (status);
	    if (GetNewRef (phdr, "DFLTFILE", sciref))
		    return (status);
	    if (GetNewRef (phdr, "LFLTFILE", sciref))
		    return (status);
	}

	if (sci_sw->shadcorr == PERFORM) {
	    wf3->sci_basic_2d = PERFORM;   refimage_used = 1;
	    CCDSanity (wf3->detector, "SHADCORR");
	    if (GetNewRef (phdr, "SHADFILE", sciref))
		    return (status);
	}

	/* Note that we set sci_basic_ccd rather than sci_basic_2d. */
	if (sci_sw->blevcorr == PERFORM) {
	    wf3->sci_basic_ccd = PERFORM;
	    if (GetNewRef (phdr, "OSCNTAB", sciref))
		    return (status);
	} else if (sci_sw->blevcorr == OMIT && refimage_used) {
	    /* Dark, flat, etc., assume the overscan has been subtracted. */
	    trlwarn 
		("For science file, should do BLEVCORR to remove overscan ");
	    /* Use trlmessage for second line of warning/error messages */
	    trlmessage ("before doing other steps that use reference images.");
	}

	if (sci_sw->photcorr == PERFORM) {
	    wf3->sci_basic_2d = PERFORM;
	    if (GetNewRef (phdr, "IMPHTTAB", sciref))
		    return (status);
	}

	if (sci_sw->fluxcorr == PERFORM) {
	    CCDSanity (wf3->detector, "FLUXCORR");
		return (status);
	}

	return (status);
}

int GetIRRef (WF3Info *wf3, IR_Switch *sci_sw, Hdr *phdr,
	      RefFileInfo *sciref) {

/* arguments:
WF3Info *wf3          i: calibration flags and other info
IR_Switch *sci_sw     i: all calibration switches for science file
Hdr *phdr             i: primary header of science file
RefFileInfo *sciref  io: list of keyword,filename pairs
*/

	extern int status;

	int GetNewRef (Hdr *, char *, RefFileInfo *);

	if (GetNewRef (phdr, "CCDTAB", sciref))
	    return (status);

	if (GetNewRef (phdr, "OSCNTAB", sciref))
	    return (status);

	if (sci_sw->dqicorr == PERFORM) {
	    wf3->sci_basic_ir = PERFORM;
	    if (GetNewRef (phdr, "BPIXTAB", sciref))
		return (status);
	}

	if (sci_sw->blevcorr == PERFORM)
	    wf3->sci_basic_ir = PERFORM;

	if (sci_sw->zoffcorr == PERFORM)
	    wf3->sci_basic_ir = PERFORM;

	if (sci_sw->noiscorr == PERFORM)
	    wf3->sci_basic_ir = PERFORM;

	if (sci_sw->darkcorr == PERFORM || sci_sw->zsigcorr == PERFORM) {
	    wf3->sci_basic_ir = PERFORM;
	    if (GetNewRef (phdr, "DARKFILE", sciref))
		    return (status);
	}

	if (sci_sw->nlincorr == PERFORM || sci_sw->zsigcorr == PERFORM) {
	    wf3->sci_basic_ir = PERFORM;
	    IRSanity (wf3->detector, "NLINCORR");
	    if (GetNewRef (phdr, "NLINFILE", sciref))
		return (status);
	}

	if (sci_sw->flatcorr == PERFORM) {
	    wf3->sci_basic_ir = PERFORM;
	    if (GetNewRef (phdr, "PFLTFILE", sciref))
		return (status);
	    if (GetNewRef (phdr, "DFLTFILE", sciref))
		return (status);
	    if (GetNewRef (phdr, "LFLTFILE", sciref))
		return (status);
	}

	if (sci_sw->photcorr == PERFORM) {
	    wf3->sci_basic_ir = PERFORM;
	    if (GetNewRef (phdr, "IMPHTTAB", sciref))
		return (status);
	}

	if (sci_sw->unitcorr == PERFORM)
	    wf3->sci_basic_ir = PERFORM;

	if (sci_sw->crcorr == PERFORM)
	    wf3->sci_basic_ir = PERFORM;

	if (sci_sw->rptcorr == PERFORM) {
	    if (wf3->nimages < 2) {
		trlwarn 
		   ("RPTCORR will be omitted because there's only one image.");
		wf3->sci_rptcorr = OMIT;
	    } else {
		wf3->sci_basic_ir = PERFORM;
	    }
	}

	if (sci_sw->crcorr == PERFORM || sci_sw->rptcorr == PERFORM) {
	    if (GetNewRef (phdr, "CRREJTAB", sciref))
		return (status);
	}

	return (status);
}

static void IRSanity (int detector, char *calswitch) {

	if (detector != IR_DETECTOR) {
	    sprintf (MsgText, "%s = PERFORM, but detector is UVIS.", calswitch);
	    trlwarn (MsgText);
	}
}

static void CCDSanity (int detector, char *calswitch) {

	if (detector == IR_DETECTOR) {
	    sprintf (MsgText, "%s = PERFORM, but detector is IR.", calswitch);
	    trlwarn (MsgText);
	}
}

