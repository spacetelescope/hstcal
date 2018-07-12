#ifndef CALACS_INCL
#define CALACS_INCL

#include "trlbuf.h"
#include "hstcal.h"

/* calacs -- integrated calacs reduction
Warren Hack, 1998 May 12: Initial version.
Pey Lian Lim, 2013 Aug 9: Separated PCTECORR from ACSCCD.
*/

typedef struct {
    /* name of association table exposure comes from */
    char asn_table[CHAR_LINE_LENGTH];
    char crj_root[CHAR_LINE_LENGTH];

    /* input, outroot, and temporary file names */
    char rawfile[CHAR_LINE_LENGTH];  /* uncalibrated science data */
    char outroot[CHAR_LINE_LENGTH];  /* file name _raw for output product */
    char crjfile[CHAR_LINE_LENGTH];  /* CR rejected, flat fielded science */
    char crcfile[CHAR_LINE_LENGTH];  /* crjfile + CTE correction */
    char fltfile[CHAR_LINE_LENGTH];  /* flat fielded science */
    char flcfile[CHAR_LINE_LENGTH];  /* fltfile + CTE correction */
    char blv_tmp[CHAR_LINE_LENGTH];  /* blevcorr, then CR flagged */
    char blc_tmp[CHAR_LINE_LENGTH];  /* blv_tmp + CTE correction */
    char crj_tmp[CHAR_LINE_LENGTH];  /* CR rejected, summed */
    char crc_tmp[CHAR_LINE_LENGTH];  /* crj_tmp + CTE correction */
    char dthfile[CHAR_LINE_LENGTH];  /* dither combined science data */
    char mtype[SZ_STRKWVAL+1];  /* Role of exposure in association */

    char rootname[CHAR_LINE_LENGTH];  /* root name for set of obs */

    /* info about input science file */
    int detector;  /* integer code for detector */
    int nchips;    /* number of IMSETs in file */
    int nimages;   /* number of images in this set */

    /* Info on the binning and gain of the science file. */
    int scibin[2];  /* binning factors */
    int scigain;    /* ccdgain values */
    int samebin;
    int readnoise_only;

    /* calibration switches */
    int sci_basic_ccd;  /* do acsccd? (dqicorr or blevcorr) */
    int sci_basic_cte;  /* do acscte? (PCTECORR) */
    int sci_basic_2d;   /* do acs2d for science file? */
    int sci_crcorr;     /* do cosmic-ray rejection for science file? */
    int sci_rptcorr;    /* combine repeatobs science data? */
    int sci_dthcorr;    /* dither combine science data? */

} ACSInfo;

#endif
