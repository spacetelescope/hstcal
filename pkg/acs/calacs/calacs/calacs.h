#ifndef CALACS_INCL
#define CALACS_INCL

/* calacs -- integrated calacs reduction
Warren Hack, 1998 May 12: Initial version.
Pey Lian Lim, 2013 Aug 9: Separated PCTECORR from ACSCCD.
*/

typedef struct {
    /* name of association table exposure comes from */
    char asn_table[ACS_LINE];
    char crj_root[ACS_LINE];

    /* input, outroot, and temporary file names */
    char rawfile[ACS_LINE];  /* uncalibrated science data */
    char outroot[ACS_LINE];  /* file name _raw for output product */
    char crjfile[ACS_LINE];  /* CR rejected, flat fielded science */
    char crcfile[ACS_LINE];  /* crjfile + CTE correction */
    char fltfile[ACS_LINE];  /* flat fielded science */
    char flcfile[ACS_LINE];  /* fltfile + CTE correction */
    char blv_tmp[ACS_LINE];  /* blevcorr, then CR flagged */
    char blc_tmp[ACS_LINE];  /* blv_tmp + CTE correction */
    char crj_tmp[ACS_LINE];  /* CR rejected, summed */
    char crc_tmp[ACS_LINE];  /* crj_tmp + CTE correction */
    char dthfile[ACS_LINE];  /* dither combined science data */
    char mtype[SZ_STRKWVAL+1];  /* Role of exposure in association */

    char rootname[ACS_LINE];  /* root name for set of obs */

    /* info about input science file */
    int detector;  /* integer code for detector */
    int nchips;    /* number of IMSETs in file */
    int nimages;   /* number of images in this set */

    /* Info on the binning and gain of the science file. */
    int scibin[2];  /* binning factors */
    int scigain;    /* ccdgain values */
    int samebin;
    int newbias;

    /* calibration switches */
    int sci_basic_ccd;  /* do acsccd? (dqicorr or blevcorr) */
    int sci_basic_cte;  /* do acscte? (PCTECORR) */
    int sci_basic_2d;   /* do acs2d for science file? */
    int sci_crcorr;     /* do cosmic-ray rejection for science file? */
    int sci_rptcorr;    /* combine repeatobs science data? */
    int sci_dthcorr;    /* dither combine science data? */

} ACSInfo;

#endif
