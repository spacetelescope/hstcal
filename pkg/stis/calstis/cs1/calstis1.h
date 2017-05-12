#ifndef INCL_CALSTIS1_H
#define INCL_CALSTIS1_H

/* calstis1.h

  Phil Hodge, 1997 Sept 12:
	Change sdqflags from unsigned short to short.

  Phil Hodge, 1997 Dec 11:
	Include ref and apertab in StisInfo1, and include PhotInfo.

  Phil Hodge, 1998 June 8:
	Include obstype, aper_fov, and cdelt in StisInfo1.  Delete subarray.

   Phil Hodge, 1998 July 23:
	Include sts->bias_or_dark in StisInfo1.

   Phil Hodge, 1998 July 30:
	Add one to the length of each string buffer.
	Replace assoc_typ with wavecal.

   Phil Hodge, 1998 Sept 24:
	Remove bias_rej.

   Phil Hodge, 1998 Oct 6:
	Include temperature in StisInfo1.

   Phil Hodge, 1999 Nov 2:
	Include err_init_bias in StisInfo1.  This will either be zero or
	ccdbias, depending on whether blevcorr has been done or not.

   Ivo Busko, 2001 Oct 26:
	Include CCD housing average temperature in StisInfo1.

   Ivo Busko, 2002 Mar 19:
	Array stores darkscale factors for each IMSET.

   Phil Hodge, 2003 Jan 17:
	Add expend, ra_targ, dec_targ.

   Paul Barrett, 2003 Sep 18:
        Add epc_rows, epc_mjd, epc_temp.

   Paul Barrett, 2003 Sep 25:
        Add tdscorr and tdstab for TDS correction.

   Phil Hodge, 2004 July 28:
	Add blev_clip for rejecting outliers in virtual overscan region.
	Use "# define" for DEFAULT_TDC.

   Phil Hodge, 2004 Dec 27:
	Replace temperature and ccd_temperature with detector_temp.

   Phil Hodge, 2007 May 9:
	Add bias_exposure, to be set to 1 (true) if the observation is a bias.

   Phil Hodge, 2011 May 9:
	Remove filtcorr, tdscorr, apertab, and tdstab from StisInfo1.  The
        PhotInfo struct was deleted.

   Phil Hodge, 2011 Nov 17:
	Put tdscorr and tdstab back into StisInfo1.

   Phil Hodge, 2012 Oct 15:
	Remove DEFAULT_TDC.
*/

# define	MAX_IMSET	100

typedef struct {
	/* input and output image names */
	char input[STIS_LINE+1];	/* input image to be calibrated */
	char output[STIS_LINE+1];	/* output calibrated image */
	char outblev[STIS_LINE+1];	/* output text file for bias levels */

	char rootname[STIS_CBUF+1];	/* root name for set of obs */

	/* command-line flags */
	int printtime;			/* print time after each step? */
	int verbose;			/* print additional info? */

	/* keywords and file names for reference files */
	RefFileInfo *refnames;

	/* info about input image */
	char det[STIS_CBUF+1];		/* name of detector */
	char aperture[STIS_CBUF+1];	/* aperture name */
	char opt_elem[STIS_CBUF+1];	/* name of grating or mirror */
	char obsmode[STIS_CBUF+1];	/* e.g. ACCUM or TIMETAG */
	char obstype[STIS_CBUF+1];	/* SPECTROSCOPIC or IMAGING */
	char aper_fov[STIS_CBUF+1];	/* aperture size */
	int detector;			/* integer code for detector */
	int wavecal;			/* true if input file is a wavecal */
	int ncombine;			/* number previously summed together */
	int nimages;			/* number of "groups" in file */
	int bias_or_dark;		/* true if image is a BIAS or DARK */
	int bias_exposure;		/* true if image is a BIAS */
	double exptime;			/* exposure time */
	double expstart;		/* exposure start time (MJD) */
	double expend;			/* exposure end time (MJD) */
	double ra_targ;			/* right ascension of target */
	double dec_targ;		/* declination of target */

	short sdqflags;		/* serious data quality values */

	/* coordinate info */
	int bin[2];		/* size of pixel in detector coordinates */
	int dispaxis;		/* dispersion axis, 1 ==> X, 2 ==> Y */
	int dispsign;		/* sign of CD1_1 (dispersion), +1 or -1 */
	double cdelt[2];	/* pixel size, for imaging type */

	/* Detector temperature (or housing temperature for side-2 CCD)
	   is used for scaling the dark reference image and for temperature
	   dependence of the sensitivity.
	*/
	double detector_temp;	/* degrees Celsius */

	/* MAMA-specific info (regarding linearity) */
	double global_limit;	/* count rate for 10% global nonlinearity */
	double tau;		/* time constant for global nonlinearity */
	double local_limit;	/* count rate for 10% local nonlinearity */
	float expand;		/* radius for flagging local nonlinearity */

	/* MAMA-specific info (regarding Doppler shift) */
	double globrate;	/* global count rate */
	double doppzero;	/* Doppler shift zero phase time (MJD) */
	double doppmag;		/* Doppler shift magnitude (high-res pixels) */
	double orbitper;	/* Assumed HST orbital period (seconds) */

	/* CCD housing temperature is used for scaling CCD dark (side 2). */
        int    epc_rows;        /* Number of rows in EPC table */
        double *epc_mjd;        /* MJD from EPC table */
        double *epc_temp;       /* Temp from EPC table */

	/* Dark scaling factors for each IMSET */
	float darkscale[MAX_IMSET]; /* read from command line */
	int ndarkscale;             /* # of dark scale entries in array */

	/* CCD-specific info */
	char ccdamp[STIS_CBUF+1];	/* CCD amplifier read out (A,B,C,D) */
	int ccdgain;		/* commanded gain of CCD */
	int ccdoffset;		/* commanded offset */
	int binaxis[2];		/* BINAXIS1, BINAXIS2 from header */
	float atodgain;		/* actual gain of CCD */
	float ccdbias;		/* CCD bias offset number */
	float err_init_bias;	/* bias to subtract for error initialization */
	float readnoise;	/* readout noise */
	float saturate;		/* CCD saturation level */
	float blev_clip;	/* N * sigma rejection criterion */

	/* calibration flags (switches) */
	int doppcorr;		/* Doppler convolution needed (flats, etc)? */
	int lorscorr;		/* convert MAMA data to low-res */
	int dqicorr;		/* data quality initialization */
	int atodcorr;		/* analog to digital correction */
	int blevcorr;		/* subtract bias from overscan */
	int biascorr;		/* subtract bias image */
	int glincorr;		/* global nonlinearity correction */
	int lflgcorr;		/* flag local nonlinearity */
	int darkcorr;		/* subtract dark image */
	int flatcorr;		/* apply flat field */
	int pfltcorr;			/* apply pixel-to-pixel flat */
	int dfltcorr;			/* apply delta flat */
	int lfltcorr;			/* apply low-order flat */
	int shadcorr;		/* correct short exposures for shutter time */
	int photcorr;		/* add photometry header keyword values */
	int noisecorr;		/* initialize error array? */
	int statcorr;		/* compute statistics? */
	/* tdscorr is not a switch, it indicates whether the TDSTAB exists */
	int tdscorr;		/* temperature- and time-dependent corr */
        /* crcorr is not a switch, it indicates if the CRCORR is complete */
        int crcorr;             /* cosmic-ray rejection */

	/* calibration images and tables */
	RefImage bias;		/* bias image */
	RefImage dark;		/* dark image */
	RefImage pflt;		/* pixel-to-pixel flat field */
	RefImage dflt;		/* delta flat */
	RefImage lflt;		/* low-order flat */
	RefImage shad;		/* shutter shading correction image */
	RefTab bpix;		/* bad pixel table */
	RefTab ccdpar;		/* CCD parameters table */
        RefTab epctab;          /* Engineering parameters calibration table */
	RefTab mlin;		/* MAMA nonlinearity info table */
	RefTab atod;		/* analog to digital correction table */
	RefTab phot;		/* photometric throughput table */
	RefTab tdctab;		/* NUV dark temp- and time-dep corr table */
	RefTab tdstab;		/* Time-dependent corr. table (for temp.) */
} StisInfo1;

#endif /* INCL_CALSTIS1_H */
