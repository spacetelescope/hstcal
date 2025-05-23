# include <stdio.h>

#include "hstcal.h"
# include "acs.h"
# include "calacs.h"
# include "hstcalerr.h"
# include "acscorr.h"


/* Initialize the CalSwitch structure.  This includes information about the
   input and output images, calibration files, and flags to specify which
   calibration steps are to be done.  Get switches and reference file
   names from the primary headers, and check that the files exist.

   Warren Hack, 1998 May 26:
   Initial ACS version
*/
int ACSRefInit (CALACSInfo *acs, CalSwitch *sci_sw, RefFileInfo *sciref) {

    /* arguments:
       acs      i: calibration flags and other info
       sci_sw   o: all calibration switches (0 or 1) for science file
       sciref   o: reference file name and info
    */

    extern int status;

    int missing;  /* number of missing reference files */

    int GetSciInfo (CALACSInfo *, CalSwitch *, RefFileInfo *);
    void RefExist (RefFileInfo *, int, int *);
    int InsertACSSuffix (CALACSInfo *);

    /* Construct output and temporary file names. */
    if (InsertACSSuffix (acs))
        return (status);

    /* Get switches and reference file names for science file. */
    if (GetSciInfo (acs, sci_sw, sciref))
        return (status);

    /* Verify that the reference files that we need do exist. */
    missing = 0;  /* initial value */
    RefExist (sciref, acs->detector, &missing);  /* for science file */

    if (missing > 0) {
        if (missing == 1) {
            trlerror ("One reference file was missing.");
        } else {
            trlerror("%d reference files were missing.", missing);
        }
        return (status = CAL_FILE_MISSING);
    }

    return (status);
}


void ACSDefaults (CALACSInfo *acs) {

    /* rawfile has been assigned already */
    acs->crj_root[0] = '\0';
    acs->outroot[0] = '\0';
    acs->crjfile[0] = '\0';
    acs->crcfile[0] = '\0';
    acs->fltfile[0] = '\0';
    acs->flcfile[0] = '\0';
    acs->blv_tmp[0] = '\0';
    acs->blc_tmp[0] = '\0';
    acs->crj_tmp[0] = '\0';
    acs->crc_tmp[0] = '\0';
    acs->dthfile[0] = '\0';
    acs->mtype[0] = '\0';

    acs->rootname[0] = '\0';

    acs->detector = UNKNOWN_DETECTOR;
    acs->nchips = 1;
    acs->nimages = 1;

    acs->scibin[0] = 0;
    acs->scibin[1] = 0;
    acs->scigain = 0;
    acs->samebin = 0;
    acs->readnoise_only = 0;

    /* Initialize flags to not perform the step. */
    acs->sci_basic_ccd = OMIT;
    acs->sci_basic_cte = OMIT;
    acs->sci_basic_2d = OMIT;
    acs->sci_crcorr = OMIT;
    acs->sci_rptcorr = OMIT;
    acs->sci_dthcorr = OMIT;
}


/* Construct output and temporary file names from outroot.

   Note that we do not construct the wavecal (_wav) file name.  This is
   because it is either given on the command line or is gotten as the
   value of the WAVECAL keyword.
*/
int InsertACSSuffix (CALACSInfo *acs) {

    extern int status;

    int MkName (char *, char *, char *, char *, char *, int);

    if (MkName (acs->crj_root, "_crj", "_crj", "", acs->crjfile, CHAR_LINE_LENGTH))
        return (status);

    if (MkName (acs->crj_root, "_crj", "_crj_tmp", "", acs->crj_tmp, CHAR_LINE_LENGTH))
        return (status);

    if (MkName (acs->rootname, "_raw", "_flt", "", acs->fltfile, CHAR_LINE_LENGTH))
        return (status);

    if (MkName (acs->crj_root, "_crj", "_crc", "", acs->crcfile, CHAR_LINE_LENGTH))
        return (status);

    if (MkName (acs->crj_root, "_crj", "_crc_tmp", "", acs->crc_tmp, CHAR_LINE_LENGTH))
        return (status);

    if (MkName (acs->rootname, "_raw", "_flc", "", acs->flcfile, CHAR_LINE_LENGTH))
        return (status);

    if (MkName (acs->asn_table, "_raw", "_dth", "", acs->dthfile, CHAR_LINE_LENGTH))
        return (status);

    if (MkName (acs->rootname, "_raw", "_blv_tmp", "", acs->blv_tmp, CHAR_LINE_LENGTH))
        return (status);

    if (MkName (acs->rootname, "_raw", "_blc_tmp", "", acs->blc_tmp, CHAR_LINE_LENGTH))
        return (status);

    /*trlmessage("AcsInit: blv_tmp = %s ",acs->blv_tmp);
    */
    return (status);
}
