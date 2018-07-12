#ifndef INCL_CALSTIS0_H
#define INCL_CALSTIS0_H

/* calstis0 -- integrated calstis reduction

   Phil Hodge, 1997 Dec 10:
	Include binning and gain info, and samebin flag.

   Phil Hodge, 1998 Aug 4:
	Add one to the length of each string buffer.

   Phil Hodge, 2000 Jan 5:
	Add echelle flag to StisInfo, for calstis4.

   Phil Hodge, 2000 Oct 5:
	Add sc2dcorr to CalSwitch struct.

   Phil Hodge, 2001 Feb 22:
	Add prism flag to StisInfo, for calstis4.

   Phil Hodge, 2004 Mar 1:
	Add ctecorr variable to CalSwitch struct.
*/

/* values based on OBSTYPE keyword */
# define UNKNOWN_TYPE      (-1)
# define SPECTROSCOPIC_TYPE  1
# define IMAGING_TYPE        2

typedef struct {

	/* input, outroot, and temporary file names */
	char rawfile[STIS_LINE+1];	/* uncalibrated science data */
	char outroot[STIS_LINE+1];	/* root name for output */
	char wavfile[STIS_LINE+1];	/* input wavecal file */
	char crjfile[STIS_LINE+1];	/* CR rejected, flat fielded science */
	char fltfile[STIS_LINE+1];	/* flat fielded science */
	char x1dfile[STIS_LINE+1];	/* extracted 1-D spectrum */
	char x2dfile[STIS_LINE+1];	/* rectified 2-D spectrum */
	char sx2file[STIS_LINE+1];	/* summed, rectified 2-D */
	char sx1file[STIS_LINE+1];	/* summed, extracted 1-D */
	char sflfile[STIS_LINE+1];	/* summed, not rectified 2-D */
	char blv_tmp[STIS_LINE+1];	/* blevcorr, then CR flagged */
	char crj_tmp[STIS_LINE+1];	/* CR rejected, summed */
	char fwv_tmp[STIS_LINE+1];	/* flat fielded wavecal */
	char cwv_tmp[STIS_LINE+1];	/* FF, source-subtracted wavecal */
	char w2d_tmp[STIS_LINE+1];	/* 2-D extracted wavecal */

	char rootname[STIS_CBUF+1];	/* root name for set of obs */

	/* info about input science or wavecal file */
	int detector;			/* integer code for detector */
	int obstype;			/* spectroscopic or imaging */
	int nimages;			/* number of "groups" in file */
	int echelle;			/* true if echelle data */
	int prism;			/* true if prism data */

	/* Info for comparing the binning and gain of the science file
	   and wavecal.  If they're the same (only relevant for the CCD),
	   the samebin flag should be set to one.
	*/
	int scibin[2], wavbin[2];	/* binning factors */
	int scigain, wavgain;		/* ccdgain values */
	int samebin;

	/* calibration switches */
	int sci_basic_2d_a;	/* do calstis1a? (dqicorr or blevcorr) */
	int sci_basic_2d;	/* do calstis1 for science file? */
	int sci_expscorr;	/* do calstis1 on CR flagged but non-summed? */
	int sci_crcorr;		/* do cosmic-ray rejection for science file? */
	int sci_rptcorr;	/* combine repeatobs science data? */
	int sci_wavecorr;	/* do calstis11, calstis4, and calstis12? */
	int sci_2d_rect;	/* do calstis7 for science file? */
	int sci_1d_extract;	/* do calstis6? */
	int sci_geocorr;	/* do 2-D image rectification? */
	int wav_basic_2d;	/* do calstis1 for wavecal? */
	int wav_subscicorr;	/* do calstis11?  (not a header switch) */

} StisInfo;

/* All the calibration switches. */

typedef struct {

	int atodcorr;
	int backcorr;
	int biascorr;
	int blevcorr;
	int crcorr;
	int darkcorr;
	int dispcorr;
	int doppcorr;
	int dqicorr;
	int expscorr;
	int flatcorr;
	int fluxcorr;
	int ctecorr;
	int geocorr;
	int glincorr;
	int helcorr;
	int lflgcorr;
	int lorscorr;
	int photcorr;
	int rptcorr;
	int sc2dcorr;
	int sgeocorr;
	int shadcorr;
	int wavecorr;
	int x1dcorr;
	int x2dcorr;
	int statcorr;	/* statflag */

} CalSwitch;

#endif /* INCL_CALSTIS0_H */
