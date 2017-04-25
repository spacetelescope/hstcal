# include <stdio.h>

# include "wf3.h"
# include "calwf3.h"
# include "hstcalerr.h"
# include "wf3corr.h"

/* Initialize the CCD_Switch structure.  This includes information about the
   input and output images, calibration files, and flags to specify which
   calibration steps are to be done.  Get switches and reference file
   names from the primary headers, and check that the files exist.

   Howard Bushouse, 2000 August 23:
   	Initial WFC3 version
	[Adapted from acsinit.c by W. Hack]

   H.Bushouse, 2000-Sep-28: Added sci_basic_ir switch for IR chip
   H.Bushouse, 2002-Nov-26: Removed use of sflfile (following CALACS changes).
   H.Bushouse, 2008-Aug-20: Changed dithfile suffix from "dth" to "drz".
   M. Sosey, 2015 June: Updates for CTE implementation with UVIS 2.0
   M. Sosey, 2016 Jan: changed name of rac file to rac_tmp
*/

int CCDRefInit (WF3Info *wf3, CCD_Switch *sci_sw, RefFileInfo *sciref) {

/* arguments:
WF3Info *wf3        i: calibration flags and other info
CCD_Switch *sci_sw  o: all calibration switches (0 or 1) for science file
RefFileInfo *sciref o: reference file name and info
*/

	extern int status;

	int missing;		/* number of missing reference files */

	int  GetCCDInfo (WF3Info *, CCD_Switch *, RefFileInfo *);
	void RefExist (RefFileInfo *, int, int *);
	int  InsertWF3Suffix (WF3Info *);

	/* Construct output and temporary file names. */
	if (InsertWF3Suffix (wf3))
	    return (status);

	/* Get switches and reference file names for science file. */
	if (GetCCDInfo (wf3, sci_sw, sciref))
	    return (status);

	/* Verify that the reference files that we need do exist. */
	missing = 0;			/* initial value */
    
	RefExist (sciref, wf3->detector, &missing);	/* for science file */
    
	if (missing) {
        sprintf (MsgText, "%d reference file(s) missing.", missing);
        trlerror (MsgText);
    	return (status = CAL_FILE_MISSING);
    }

	return (status);
}

int IRRefInit (WF3Info *wf3, IR_Switch *sci_sw, RefFileInfo *sciref) {

/* arguments:
WF3Info *wf3        i: calibration flags and other info
IR_Switch *sci_sw   o: all calibration switches (0 or 1) for science file
RefFileInfo *sciref o: reference file name and info
*/

	extern int status;

	int missing;		/* number of missing reference files */

	int  GetIRInfo (WF3Info *, IR_Switch *, RefFileInfo *);
	void RefExist (RefFileInfo *, int, int *);
	int  InsertWF3Suffix (WF3Info *);

	/* Construct output and temporary file names. */
	if (InsertWF3Suffix (wf3))
	    return (status);

	/* Get switches and reference file names for science file. */
	if (GetIRInfo (wf3, sci_sw, sciref))
	    return (status);

	/* Verify that the reference files that we need do exist. */
	missing = 0;			/* initial value */
	RefExist (sciref, wf3->detector, &missing);	/* for science file */

	if (missing > 0) {
	    if (missing == 1) {
		    trlerror ("One reference file was missing.");
	    } else {
 		    sprintf (MsgText, "%d reference files were missing.", missing);
		    trlerror (MsgText);
	    }
	    return (status = CAL_FILE_MISSING);
	}

	return (status);
}

void WF3Defaults (WF3Info *wf3) {

	/* RAWFILE HAS BEEN ASSIGNED ALREADY */
	wf3->crj_root[0] = '\0';
    wf3->crc_root[0] = '\0';
	wf3->outroot[0]  = '\0';
	wf3->crjfile[0]  = '\0';
    wf3->crcfile[0]  = '\0';
	wf3->fltfile[0]  = '\0';
	wf3->imafile[0]  = '\0';
	wf3->blv_tmp[0]  = '\0';
    wf3->blc_tmp[0]  = '\0';
    wf3->crc_tmp[0]  = '\0';
    wf3->rac_tmp[0]  = '\0';
	wf3->crj_tmp[0]  = '\0';
	wf3->dthfile[0]  = '\0';
	wf3->mtype[0]    = '\0';
	wf3->rootname[0] = '\0';
	wf3->detector  = UNKNOWN_DETECTOR;
	wf3->nchips    = 1;
	wf3->nimages   = 1;
	wf3->scibin[0] = 0;
	wf3->scibin[1] = 0;
	wf3->scigain   = 0;
	wf3->samebin   = 0;

	/* INITIALIZE FLAGS TO NOT PERFORM THE STEP. */
	wf3->sci_basic_ccd = OMIT;
	wf3->sci_basic_2d  = OMIT;
	wf3->sci_basic_ir  = OMIT;
    wf3->sci_basic_cte = OMIT;
	wf3->sci_crcorr    = OMIT;
	wf3->sci_rptcorr   = OMIT;
	wf3->sci_dthcorr   = OMIT;
}

/* CONSTRUCT OUTPUT AND TEMPORARY FILE NAMES FROM OUTROOT. */

int InsertWF3Suffix (WF3Info *wf3) {

	extern int status;

	int MkName (char *, char *, char *, char *, char *, int);

	if (MkName (wf3->crj_root, "_crj", "_crj", "", wf3->crjfile, SZ_LINE))
	    return (status);

	if (MkName (wf3->crj_root, "_crj", "_crj_tmp", "", wf3->crj_tmp, SZ_LINE))
	    return (status);
        
    if (MkName (wf3->crc_root, "_crc", "_crc", "", wf3->crcfile, SZ_LINE))
        return(status);
        
    if (MkName (wf3->crc_root, "_crc", "_crc_tmp", "", wf3->crc_tmp, SZ_LINE))
        return (status);
            
    if (MkName (wf3->rootname, "_raw","_rac_tmp","", wf3->rac_tmp, SZ_LINE))
        return(status);
        		
	if (MkName (wf3->rootname, "_rac_tmp", "_blc_tmp", "", wf3->blc_tmp, SZ_LINE))
	    return (status);

	if (MkName (wf3->rootname, "_raw", "_ima", "", wf3->imafile, SZ_LINE))
	    return (status);

	if (MkName (wf3->rootname, "_raw", "_flt", "", wf3->fltfile, SZ_LINE))
	    return (status);

	if (MkName (wf3->rootname, "_raw", "_flc", "", wf3->flcfile, SZ_LINE))
	    return (status);

	if (MkName (wf3->asn_table, "_raw", "_drz", "", wf3->dthfile, SZ_LINE))
	    return (status);
        
	if (MkName (wf3->rootname, "_raw", "_blv_tmp", "", wf3->blv_tmp, SZ_LINE))
	    return (status);
            
            		
	return (status);
}

