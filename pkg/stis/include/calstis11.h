#ifndef INCL_CALSTIS11_H
#define INCL_CALSTIS11_H

/* calstis11.h subtract science data */

/* include "../stis.h" */

/* Options for selecting the most appropriate group of the science
   file for subtracting from a given group of the wavecal file.
*/
# define STIS_PRECEDING  1
# define STIS_NEAREST    2
# define STIS_FOLLOWING  3

/* Image description and reference files for calstis11.

   Phil Hodge, 1998 Aug 13:
	Add one to the length of each string buffer.
	Add texpstrt to StisInfo11.

   Phil Hodge, 1999 Oct 1:
	Add gain to StisInfo11.
*/

typedef struct {

	/* image names */
	char input[STIS_LINE+1];
	char output[STIS_LINE+1];	/* output wavecal file or null */

	char rootname[STIS_CBUF+1];	/* root name for set of obs */

	/* info about input image */
	char sclamp[STIS_CBUF+1];	/* spectral calibration lamp in use */
	char obsmode[STIS_CBUF+1];	/* e.g. ACCUM or TIMETAG */
	char aperture[STIS_CBUF+1];	/* aperture name */
	char opt_elem[STIS_CBUF+1];	/* name of grating */
	char det[STIS_CBUF+1];		/* name of detector */
	int detector;			/* integer code for detector */
	int nimages;			/* number of image sets */
	int cenwave;			/* central wavelength */
	double texpstrt;		/* exposure start time */
	double gain;			/* CCD gain */

	int printtime;			/* print time after each step? */
	int verbose;			/* print additional info? */

	/* values specific to each group */
	double *exptime;		/* exposure time (sec) */
	double *midpt;			/* midpoint of exposure (MJD) */

} StisInfo11;

#endif /* INCL_CALSTIS11_H */
