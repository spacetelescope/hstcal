#ifndef INCL_CALWF3_H
#define INCL_CALWF3_H

/* calwf3 -- integrated calwf3 reduction

	Howard Bushouse, 2000 August 22: Initial version (adapted from
                         calacs.h by W. Hack)

	H.Bushouse, 2000-Sep-28: Added "sci_basic_ir" calib switches for IR chip
	H.Bushouse, 2002-Nov-26: Removed "sflfile" char string.
	H.Bushouse, 2003-Oct-16: Changed scigain from int to float.
    M.Sosey, August 2014: added PCTECORR step names and switches
*/

typedef struct {
	/* name of association table exposure comes from */
	char asn_table[CHAR_LINE_LENGTH+1];
	char crj_root[CHAR_LINE_LENGTH+1];
    char crc_root[CHAR_LINE_LENGTH+1];

	/* input, outroot, and temporary file names */
	char rawfile[CHAR_LINE_LENGTH+1];	/* uncalibrated science data */
	char outroot[CHAR_LINE_LENGTH+1];	/* file name _raw for output product */
    char rac_tmp[CHAR_LINE_LENGTH+1];    /* PCTECORR corrected raw file in orig format */
	char crjfile[CHAR_LINE_LENGTH+1];	/* CR rejected, flat fielded science */
    char crcfile[CHAR_LINE_LENGTH+1];    /* CR reject, CTE fixed, flat fielded science */
	char imafile[CHAR_LINE_LENGTH+1];	/* IR intermediate file */
	char fltfile[CHAR_LINE_LENGTH+1];	/* flat fielded science, no CTE correction */
	char flcfile[CHAR_LINE_LENGTH+1];	/* flat fielded science with CTE correction */
	char blv_tmp[CHAR_LINE_LENGTH+1];	/* blevcorr,no CTE, then CR flagged */
    char blc_tmp[CHAR_LINE_LENGTH+1];    /* blevcorr, with CTE bias*/
	char crj_tmp[CHAR_LINE_LENGTH+1];	/* CR rejected, no CTE, summed */
    char crc_tmp[CHAR_LINE_LENGTH+1];    /* CR combined with CTE done*/
	char dthfile[CHAR_LINE_LENGTH+1];	/* dither combined science data */
	char mtype[SZ_FITS_VAL+1];	/* Role of exposure in association */

	char rootname[CHAR_LINE_LENGTH+1];	/* root name for set of obs */


	/* info about input science file */
	int detector;			/* integer code for detector */
	int nchips;			/* number of IMSETs in file */
	int nimages;			/* number of images in this set */

	/* Info on the binning and gain of the science file.	*/
	int scibin[2];			/* binning factors */
	float scigain;			/* ccdgain values */
	int samebin;

	/* calibration switches */
	int sci_basic_ccd;	/* do wf3ccd? (dqicorr or blevcorr) */
	int sci_basic_2d;	/* do wf32d for science file? */
	int sci_basic_ir;	/* do wf3ir for science file? */
    int sci_basic_cte;  /* CTE correction for uvis */
	int sci_crcorr;		/* do cosmic-ray rejection for science file? */
	int sci_rptcorr;	/* combine repeatobs science data? */
	int sci_dthcorr;	/* dither combine science data?    */

} WF3Info;


#endif /* INCL_CALWF3_H */
