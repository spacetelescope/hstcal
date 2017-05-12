#ifndef INCL_CALSTIS4_H
#define INCL_CALSTIS4_H

/* requires "../stis.h"
   requires <stdio.h> for FILE
*/

/* Image description and reference files for calstis4.

   Phil Hodge, 1997 Dec 11:
	Add ref to the StisInfo4 struct.

   Phil Hodge, 1998 Aug 4:
	Add one to the length of each string buffer.

   Phil Hodge, 1998 Dec 11:
	Add parameters from the WCP table; add dbgfile and dbg pointer;
	add wl_sect1 & 2 and sp_sect1 & 2.

   Phil Hodge, 1999 Sept 23:
	Add exptime to StisInfo4 struct.

   Phil Hodge, 2000 Jan 12:
	Add variables for global offset for echelle.
	Move UNDEFINED_SHIFT to ../stis.h.

   Phil Hodge, 2001 Mar 7:
	Add disp_type to distinguish prism from grating.  Remove echelle.
	Add default values PRISM_CRPIX2 and PRISM_CDELT2.
	Add sdctab for prism.

   Phil Hodge, 2001 May 4:
	Add moffset2.

   Phil Hodge, 2004 July 23:
	Add slit_angle for long-slit echelle data.

   Nadia Dencheva, 2006 April 10:
	Add expstart for trace rotation.

   Phil Hodge, 2008 Nov 3:
	Add imset_ok.

   Phil Hodge, 2011 Jan 5:
	Add coeff_save, mref, yref, a4corr to DispRelation.  Delete moffset2.
*/

/* Codes for slit types. */
# define SHORT_ECHELLE_SLIT   1
# define MEDIUM_ECHELLE_SLIT  2
# define LONG_SLIT            3
# define UNKNOWN_SLIT       (-1)

# define DEGREES_TO_RADIANS (3.14159265358979323846 / 180.)

/* This info is for PRISM data prior to 2-D rectification.  The coordinate
   parameters in the header are for RA & Dec, but we need cross-dispersion
   spatial coordinates.  We get these from the SDC table, but if the SDC
   table doesn't have a row for prism, we can use these default values.
*/
# define PRISM_CRPIX2  (512.)			/* one indexed */
# define PRISM_CDELT2  (0.02915 / 3600.)	/* degrees per pixel */

typedef struct {

	/* input image name */
	char input[STIS_LINE+1];	/* input image */

	char rootname[STIS_CBUF+1];	/* root name for set of obs */

	char dbgfile[STIS_LINE+1];	/* file for debug output */
	FILE *dbg;			/* file handle for debug output */
	int printtime;			/* print time after each step? */
	int verbose;			/* print additional info? */
	double slit_angle;	/* (radians) for long slit with echelle */

	/* keywords and file names for reference files */
	RefFileInfo *refnames;

	/* info about input image */
	char obsmode[STIS_CBUF+1];	/* e.g. ACCUM or TIMETAG */
	char aperture[STIS_CBUF+1];	/* aperture name */
	char opt_elem[STIS_CBUF+1];	/* name of grating or mirror */
	char det[STIS_CBUF+1];		/* name of detector */
	char aper_fov[STIS_CBUF+1];	/* needed for echelle only */
	int detector;		/* NUV-MAMA, FUV-MAMA, or CCD */
	int nimages;		/* number of image sets in file */
	int disp_type;		/* disperser is grating, echelle, or prism */

	short sdqflags;		/* serious data quality values */

	/* wavecal-specific info */
	char sclamp[STIS_CBUF+1];	/* spectral calibration lamp in use */
	char lampset[STIS_CBUF+1];	/* current (mA) through cal. lamp */

	/* values specific to each group (except nimsets) */

	/* coordinate info */
	int cenwave;		/* central wavelength */
	int dispaxis;		/* dispersion axis, 1 ==> X, 2 ==> Y */
	double crpix[2];	/* reference pixel in input image, 0-indexed */
	double crval[2];	/* coordinates at reference pixel */
	double cdelt[2];	/* pixel spacing in physical units */
	double scale[2];	/* multiply to convert to reference pixels */

	/* these are for echelle only */
	/* Transformation from reference coordinates to image coordinates. */
	double ltm[2];		/* matrix part of transformation */
	double ltv[2];		/* linear part */
	int nimsets;		/* number of image sets in file */
	int nx, ny;		/* size of input image */

	/* we get the exposure time just to verify that it's > 0 */
        /* exposure start time is needed for trace rotation */
	double exptime;
        double expstart;

	int imset_ok;		/* true if exptime > 0 or max > min */

	/* calibration switches */
	int wavecorr;		/* process wavecal? */

	/* calibration tables */
	RefTab wcptab;		/* wavecal-processing parameters */
	RefTab lamptab;		/* LMP, template lamp table */
	RefTab apdestab;	/* APD, aperture description table */

	/* the following reference tables are needed only for
	   echelle or prism
	*/
	RefTab disptab;		/* DSP, dispersion coefficients table */
	RefTab inangtab;	/* IAC, incidence-angle correction table */
	RefTab sptrctab;	/* 1DT, 1-D spectrum trace table */

	/* the sdctab is only needed for prism */
	RefTab sdctab;		/* SDC, for cross-dispersion scale info */

	/* parameters from the WCP table, used to control processing */
	int wl_trim1;
	int wl_trim2;
	int sp_trim1;
	int sp_trim2;
	int wl_range;
	int sp_range;
	double nsigma_cr;
	double nsigma_illum;
	double mad_reject;
	double min_mad;

	/* image section to use when finding shift in wavelength (wl_sect)
	   or spatial (sp_sect) direction, derived from WCP parameters;
	   these are zero-indexed image pixel numbers
	*/
	int wl_sect1[2];	/* first axis, first and last pixel */
	int wl_sect2[2];	/* second axis, first and last pixel */
	int sp_sect1[2];	/* first axis */
	int sp_sect2[2];	/* second axis */

        double trace_rotation;  /* trace rotation angle */       
        

} StisInfo4;

/* This contains a description of the slit.  This differs from the
   similar struct in calstis7 in that this does not include the
   throughput info.
*/

# define MAX_BARS  3		/* maximum number of occulting bars */

typedef struct {

	/* The aperture throughput information in calstis7 is not
	   included here, and we also ignore the slit offsets.
	*/
	/* geometric description */
	double width[2];	/* width of slit in X and Y directions */
	int nbars;		/* number of occulting bars */
	double barlocn[MAX_BARS];	/* distance of bar from slit center */
	double barwidth[MAX_BARS];	/* width of bar */

} ApInfo;

/* This contains a spectrum of the calibration lamp. */

typedef struct {

	int allocated;		/* true if memory has been allocated */
	int nelem;		/* size of wavelength and flux arrays */
	double *wl;		/* array of wavelengths */
	double *flux;		/* array of fluxes */

} LampInfo;

/* Dispersion relation.  These values are read from the _dsp table. */
# define MAX_DISP_COEFF  10	/* length of array of dispersion coeff */

typedef struct {
	int ncoeff;		/* number of coefficients */
	double coeff[MAX_DISP_COEFF];	/* array of coefficients */
	double coeff_save[MAX_DISP_COEFF];	/* coeff. without a4corr */
	/* the following three are for the a4corr correction */
	int mref;
	double yref;
	double a4corr;
} DispRelation;

/* Spectrum trace info, the displacement in the cross-dispersion
   direction at each detector pixel, read from the _1dt table.
*/
# define MAX_SP_TRACE  1024	/* length of array for spectrum trace */

typedef struct sptrc *TracePtr;

typedef struct sptrc {
	double a2center;	/* Y location on detector */
	double a1center;	/* X location on detector */
	int sporder;		/* spectral order number */
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

#endif /* INCL_CALSTIS4_H */
