#ifndef INCL_CALSTIS7_H
#define INCL_CALSTIS7_H

/* calstis7.h 2-D spectral extraction */

/* requires "../stis.h" */

/* values based on OBSTYPE keyword */
# define UNKNOWN_TYPE      (-1)
# define SPECTROSCOPIC_TYPE  1
# define IMAGING_TYPE        2


/* Image description and reference files for calstis7.

   Phil Hodge, 1997 Dec 11:
	Add pctcorr and pctab to the StisInfo7 struct;
	add ref to the StisInfo7 struct;
	add pcorr to the PhotInfo struct.

   Phil Hodge, 1998 Mar 17:
	Remove ref_aper from the ApInfo struct.

   Phil Hodge, 1998 Aug 4:
	Add one to the length of each string buffer (except ref_aper).

   Phil Hodge, 2000 Jan 13:
	Add one to ref_aper buffer size (in DispRelation).

   Phil Hodge, 2000 Apr 5:
	Add center_target to the StisInfo7 struct.

   Phil Hodge, 2000 Aug 9:
	Rename pc1 & pc2 --> xcoeff & ycoeff in DistInfo
	Remove all references to MAMA offset table or coefficients.

   Phil Hodge, 2000 Dec 5:
	Change the contents of the DistInfo struct.

   Phil Hodge, 2001 Apr 30:
	Add a4corr to StisInfo7.

   Phil Hodge, 2001 Aug 8:
	Add filter to StisInfo7.

   Ivo Busko, 2002 Jan 10:
	Add time-dependent sensitivity stuff: tdscorr and tdstab.

   Ivo Busko, 2002 Feb 04:
	Move PhotInfo structure definition to ../stispht.h.

   Ivo Busko, 2002 Apr 24:
        Add blazeshift parameter to StisInfo7 structure.

   Paul Barrett, 2004 Feb 11:
        Add slit angle to ApInfo structure.

   Phil Hodge, 2004 Dec 27:
	Add detector_temp to the StisInfo7 structure.

   Phil Hodge, 2006 Feb 7:
	Add err_algorithm to the StisInfo7 structure.

   Phil Hodge, 2006 Sept 21:
	Add wx2dcorr to the StisInfo7 structure (for output from wx2d).
*/

typedef struct {

	/* input and output image names */
	char input[STIS_LINE+1];	/* input image to be calibrated */
	char output[STIS_LINE+1];	/* output calibrated image */

	char rootname[STIS_CBUF+1];	/* root name for set of obs */

	int printtime;			/* print time after each step? */
	int verbose;			/* print additional info? */

	/* keywords and file names for reference files */
	RefFileInfo *refnames;

	/* info about input image */
	char obsmode[STIS_CBUF+1];
	char aperture[STIS_CBUF+1];
	char filter[STIS_CBUF+1];
	char opt_elem[STIS_CBUF+1];
	char det[STIS_CBUF+1];	/* name of detector */
	int detector;		/* NUV-MAMA, FUV-MAMA, or CCD */
	int obstype;		/* spectroscopic or imaging */
	int wavecal;		/* true if ASN_MTYP is "WAVECAL" */
	int nimages;		/* number of "groups" in file */
	int first_order;	/* true if first order grating */
	int center_target;	/* center target in output?  normally true */

	short sdqflags;		/* serious data quality values */

	/* CCD-specific info */
	double atodgain;	/* actual gain of CCD */

	/* values specific to each group */
	double exptime;		/* exposure time (sec) */
	double expstart;	/* exposure start time (MJD) */
	double expend;		/* exposure stop time (MJD) */
	double hfactor;		/* heliocentric correction factor */

	/* coordinate info */
	double ltm[2];		/* matrix part of transformation to det. */
	double ltv[2];		/* offset to detector coordinates */
	int cenwave;		/* central wavelength */
	int dispaxis;		/* dispersion axis, 1 ==> X, 2 ==> Y, or 0 */
	double ra_targ;		/* right ascension of target, degrees */
	double dec_targ;	/* declination of target, degrees */

	/* these crpix are from the SCI extension, not the reference table */
	double crpix[2];	/* reference pixel, zero-indexed */

	/* This is based on output cdelt, but scaled to be similar to input.
	   Note that this is only used for spectroscopic mode; for imaging
	   mode, the scale is given in the DistInfo struct.
	*/
	double plate_scale[2];	/* arcseconds per pixel, for spectroscopic */

	/* various offsets */
	double msm_slop[2];	/* SHIFTAi (measured by cs4) */
	double ap_offset[2];	/* aperture offset in pixels */
	double total_offset[2];	/* sum of the above offsets */

	/* correction to echelle dispersion solution for MSM shift in Y */
	double a4corr;

	/* blaze shift in pixels from command line */
	double blazeshift;

	/* algorithm for interpolating the error array (from command line) */
	int err_algorithm;

	/* Detector temperature (or housing temperature for side-2 CCD)
	   is used for the temperature dependence of the sensitivity.
	*/
	double detector_temp;	/* degrees Celsius */

	/* calibration switches */
	int wavecorr;		/* has wavecal processing been completed? */
	/* if wx2d was run, the trace should be set to zero */
	int wx2dcorr;
	int x2dcorr;		/* 2-D spectral extraction (or geocorr)? */
	int x2dcorr_o;		/* local to a given order (for pedigree) */
	int sgeocorr;		/* correct for small-scale distortions? */
	int heliocorr;		/* convert to heliocentric wavelengths? */
	int fluxcorr;		/* convert to absolute flux units? */
	int fc_corr;		/* conserve flux? */
	int statcorr;		/* compute statistics? */
	/* pctcorr and tdscorr are not header keywords */
	int pctcorr;		/* include PCT correction with fluxcorr? */
	int tdscorr;		/* include TDS correction with fluxcorr? */

	/* calibration images and tables */
	RefTab distntab;	/* SDC or IDC, coordinate & dist info table */
	RefImage sdstfile;	/* SSD, MAMA small-scale distortion file */
	RefTab apdestab;	/* APD, aperture description table */
	RefTab apertab;		/* APT, relative aperture throughput table */
	RefTab phottab;		/* PHT, photometric throughput table */
	RefTab tdstab;		/* TDS, time-dependent sensitivity table */
	RefTab disptab;		/* DSP, dispersion coefficients table */
	RefTab inangtab;	/* IAC, incidence-angle correction table */
	RefTab sptrctab;	/* 1DT, 1-D spectrum trace table */
	RefTab pctab;		/* PCT, photometric correction table */

	double trace_rotation;	/* trace rotation angle */

} StisInfo7;

/* CoordInfo describes the coordinates for the output image.  Note that
   the coordinate parameters (WCS info) describe the output rectified
   image, not the input or any reference image.
   These values are read from the _sdc table in spectroscopic mode.
   For imaging mode, npix[0] and npix[1] are read from the _idc table
   (and scaled depending on binning or subarray), but the other structure
   members are not used.
*/

typedef struct coo *CoordPtr;

typedef struct coo {
	int sporder;		/* order number */
	double a2center;	/* Y detector location corresp. to crpix[1] */
	int npix[2];		/* size of output image */
	double crpix[2];	/* reference pixel, zero-indexed */
	double crval[2];	/* coordinates at reference pixel */
	double cdelt[2];	/* increment per pixel */
	CoordPtr next;		/* pointer to next struct in list */
} CoordInfo;			/* for info from SDC */

/* Distortion information, used for imaging mode only, read from _idc.
   Unlike the CoordInfo struct, the members here are all in reference
   pixels; none of them will be scaled to image pixels.
*/

# define MAX_ORDER   5			/* maximum degree of polynomial */
# define DCM_SIZE  (MAX_ORDER + 1)	/* distortion coefficient matrix size */
# define MAX_NCOEFF (DCM_SIZE * DCM_SIZE)	/* number of coefficients */
# define WHICH_COEFF(i,j) (i*DCM_SIZE + j)

typedef struct {
	int allocated;		/* true if memory has been allocated */
	int npix[2];		/* size of output image */
	int norder;		/* order of the polynomial fit */
	double scale;		/* arcseconds per reference pixel */
	double xref, yref;	/* zero points in input (zero indexed) */
	double cxref, cyref;	/* zero points in output (zero indexed) */
	double *xcoeff, *ycoeff;	/* coefficients for X and Y */
} DistInfo;

/* Dispersion relation at A2CENTER, a Y location on the detector.  These
   values are read from the _dsp table.
*/
# define MAX_DISP_COEFF  10	/* length of array of dispersion coeff */

typedef struct dsp *DispPtr;

typedef struct dsp {
	double a2center;	/* Y location on detector */
	int ncoeff;		/* number of coefficients */
	double coeff[MAX_DISP_COEFF];	/* array of coefficients */
	char ref_aper[STIS_CBUF+1];
		/* reference aperture for dispersion relation */
	DispPtr next;		/* pointer to next struct in list */
} DispRelation;

/* Spectrum trace info, the displacement in the cross-dispersion
   direction at each detector pixel, read from the _1dt table.
*/
# define MAX_SP_TRACE  1024	/* length of array for spectrum trace */

typedef struct sptrc *TracePtr;

typedef struct sptrc {
	double a2center;	/* Y location on detector */
	double a1center;	/* X location on detector */
	int nelem;		/* actual size of array */
	double a2displ[MAX_SP_TRACE];	/* spectrum trace */
	TracePtr next;		/* pointer to next struct in list */
} SpTrace;

/* This contains the incidence-angle correction coefficients, which are
   read from the _iac table.
*/

typedef struct {
	int allocated;		/* true if memory has been allocated */
	int ncoeff1;		/* sizes of coefficient arrays */
	int ncoeff2;
	double *coeff1;		/* arrays of coefficients */
	double *coeff2;
} InangInfo;


/* This contains a description of the slit, including throughput.  Values
   are gotten from the aperture description table (_apd) and aperture
   throughput table (_apt).
*/

# define REF_ANGLE 0.315
# define MAX_BARS  3		/* maximum number of occulting bars */

typedef struct {
	/* aperture throughput _apt */
	int allocated;		/* true if memory has been allocated */
	int nelem;		/* size of wavelength and throughput arrays */
	double *wl;		/* array of wavelengths */
	double *thr;		/* array of fraction of light passed by slit */

	/* geometric description _apd */
        double angle;           /* angle of slit title, degrees */
	double width[2];	/* width of slit in X and Y directions */
	double ap_offset[2];	/* offset from nominal location in X and Y */
	double barlocn[MAX_BARS];	/* distance of bar from slit center */
	double barwidth[MAX_BARS];	/* width of bar */
	int nbars;			/* number of occulting bars */
} ApInfo;

#endif /* INCL_CALSTIS7_H */
