# include <stdio.h>
# include "hstio.h"
# include "acs.h"
# include "acscorr.h"
# include "calacs.h"


static void CCDSanity (int, char *);
static void MAMASanity (int, char *);


/* This routine gets the names of all the reference files (images and
   tables) that are needed for calibrating the science file, based on
   the switches.  The files are not checked at this time to see whether
   they actually exist.

   Some sanity checking is done, and warning messages are printed if
   switch values don't make sense for the current detector.

   Various flags are set (e.g. acs->basic_2d), depending on which calacs
   module will perform the calibration.  Basic 2-D reduction,
   derived from 'calstis1', is divided into 3 parts; acsccd, acscte, and acs2d.
   The steps dqicorr, atodcorr, blevcorr, and biascorr constitute the acsccd
   part, and these are included with the sci_basic_ccd flag rather than with
   sci_basic_2d which is used for the acs2d steps. The step pctecorr
   constitutes the acscte part.

   Warren Hack, 1998 May 26:
     Initial version for use with ACS
   Warren Hack, 2000 Sept 11:
     Revised to add Post-Flash processing step
   Warren Hack, 2001 Oct 17:
     Replaced APERTAB and PHOTTAB with GRAPHTAB and COMPTAB
   Pey Lian Lim, 2012 Dec 12:
     Moved FLSHCORR from ACSCCD to ACS2D
   Pey Lian Lim, 2013 Aug 9:
     Separated PCTECORR from ACSCCD.
*/
int SciFlags (ACSInfo *acs, CalSwitch *sci_sw, Hdr *phdr,
              RefFileInfo *sciref) {

    /* arguments:
       ACSInfo *acs          i: calibration flags and other info
       CalSwitch *sci_sw     i: all calibration switches for science file
       Hdr *phdr             i: primary header of science file
       RefFileInfo *sciref  io: list of keyword,filename pairs
    */

    extern int status;

    int refimage_used = 0;  /* = 1 if do bias, dark, flat, or shadcorr */

    int GetNewRef (Hdr *, char *, RefFileInfo *);

    if (acs->detector == WFC_CCD_DETECTOR || acs->detector == HRC_CCD_DETECTOR) {
        if (GetNewRef (phdr, "CCDTAB", sciref))
            return (status);
    }

    /* Note that we set sci_basic_2d for MAMA data rather than sci_basic_ccd. */
    if (sci_sw->dqicorr == PERFORM) {
        if (acs->detector == MAMA_DETECTOR) {
            acs->sci_basic_2d = PERFORM;
        } else {
            acs->sci_basic_ccd = PERFORM;
        }
        if (GetNewRef (phdr, "BPIXTAB", sciref))
            return (status);
    }

    if (sci_sw->glincorr == PERFORM) {
        acs->sci_basic_2d = PERFORM;
        MAMASanity (acs->detector, "GLINCORR");
    }
    if (sci_sw->lflgcorr == PERFORM) {
        acs->sci_basic_2d = PERFORM;
        MAMASanity (acs->detector, "LFLGCORR");
    }
    if (sci_sw->glincorr == PERFORM || sci_sw->lflgcorr == PERFORM) {
        if (GetNewRef (phdr, "MLINTAB", sciref))
            return (status);
    }

    if (sci_sw->crcorr == PERFORM) {
        acs->sci_crcorr = PERFORM;
        if (acs->detector == MAMA_DETECTOR) {
            trlwarn ("CRCORR will be omitted because detector is MAMA.");
            acs->sci_crcorr = OMIT;
        }
        if (acs->nimages < 2) {
            trlwarn ("CRCORR will be omitted because there's only one image.");
            acs->sci_crcorr = OMIT;
        }
        /* reference table is checked below, after getting rptcorr */
    }

    if (sci_sw->rptcorr == PERFORM) {
        if (acs->nimages < 2) {
            trlwarn ("RPTCORR will be omitted because there's only one image.");
            acs->sci_rptcorr = OMIT;
        } else {
            acs->sci_rptcorr = PERFORM;
        }
    }

    if (acs->sci_crcorr == PERFORM) {
        if (GetNewRef (phdr, "CRREJTAB", sciref))
            return (status);
    }

    /* Note that we set sci_basic_ccd rather than sci_basic_2d. */
    if (sci_sw->atodcorr == PERFORM) {
        acs->sci_basic_ccd = PERFORM;
        CCDSanity (acs->detector, "ATODCORR");
        if (GetNewRef (phdr, "ATODTAB", sciref))
            return (status);
    }

    if (sci_sw->biascorr == PERFORM) {
        acs->sci_basic_ccd = PERFORM;   refimage_used = 1;
        CCDSanity (acs->detector, "BIASCORR");
        if (GetNewRef (phdr, "BIASFILE", sciref))
            return (status);
    }

    if (sci_sw->pctecorr == PERFORM) {
        acs->sci_basic_cte = PERFORM;   refimage_used = 1;
        CCDSanity (acs->detector, "PCTECORR");
        if (GetNewRef (phdr, "PCTETAB", sciref))
            return (status);
    }

    if (sci_sw->darkcorr == PERFORM) {
        acs->sci_basic_2d = PERFORM;   refimage_used = 1;
        if (GetNewRef (phdr, "DARKFILE", sciref))
            return (status);
    }

    if (sci_sw->darkcorr == PERFORM && sci_sw->pctecorr == PERFORM) {
        acs->sci_basic_2d = PERFORM;   refimage_used = 1;
        if (GetNewRef (phdr, "DRKCFILE", sciref))
            return (status);
    }

    if (sci_sw->flashcorr == PERFORM) {
        acs->sci_basic_2d = PERFORM;
        refimage_used = 1;
        CCDSanity (acs->detector, "FLSHCORR");
        if (GetNewRef (phdr, "FLSHFILE", sciref))
            return (status);
    }

    if (sci_sw->flatcorr == PERFORM) {
        acs->sci_basic_2d = PERFORM;
        refimage_used = 1;
        if (GetNewRef (phdr, "PFLTFILE", sciref))
            return (status);
        if (GetNewRef (phdr, "DFLTFILE", sciref))
            return (status);
        if (GetNewRef (phdr, "LFLTFILE", sciref))
            return (status);
    }

    if (sci_sw->shadcorr == PERFORM) {
        acs->sci_basic_2d = PERFORM;
        refimage_used = 1;
        CCDSanity (acs->detector, "SHADCORR");
        if (GetNewRef (phdr, "SHADFILE", sciref))
            return (status);
    }

    /* Note that we set sci_basic_ccd rather than sci_basic_2d. */
    if (sci_sw->blevcorr == PERFORM) {
        acs->sci_basic_ccd = PERFORM;
        CCDSanity (acs->detector, "BLEVCORR");
        if (GetNewRef (phdr, "OSCNTAB", sciref))
            return (status);
    } else if (sci_sw->blevcorr == COMPLETE && acs->detector != MAMA_DETECTOR) {
        /* If it BLEVCORR has already been done, continue processing */
        CCDSanity(acs->detector, "BLEVCORR");
    } else if (sci_sw->blevcorr == OMIT && refimage_used &&
               acs->detector != MAMA_DETECTOR) {
        /* Dark, flat, etc., assume the overscan has been subtracted. */
        trlwarn ("For science file, should do BLEVCORR to remove overscan ");
        /* Use trlmessage for second line of warning/error messages */
        trlmessage ("before doing other steps that use reference images.");
    }

    if (sci_sw->photcorr == PERFORM) {
        acs->sci_basic_2d = PERFORM;
        if (GetNewRef (phdr, "IMPHTTAB", sciref))
            return (status);
    }

    return (status);
}


static void MAMASanity (int detector, char *calswitch) {
    if (detector != MAMA_DETECTOR) {
        sprintf (MsgText, "%s = PERFORM, but detector is CCD.", calswitch);
        trlwarn (MsgText);
    }
}


static void CCDSanity (int detector, char *calswitch) {
    if (detector == MAMA_DETECTOR) {
        sprintf (MsgText, "%s = PERFORM, but detector is MAMA.", calswitch);
        trlwarn (MsgText);
    }
}
