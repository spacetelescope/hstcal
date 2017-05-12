#ifndef INCL_CALSTIS12_H
#define INCL_CALSTIS12_H

/* calstis12.h update SHIFTAi keywords in science file extensions */

/* include "../stis.h" */

/* Image description and reference files for calstis12.

   Phil Hodge, 1998 Mar 18:
	Remove offset, pedigree, and descrip.

   Phil Hodge, 1998 Aug 4:
	Add one to the length of each string buffer.

   Phil Hodge, 2000 Jan 12:
	Move interpolation options to ../cs12.h.
	Move UNDEFINED_SHIFT to ../stis.h.

   Phil Hodge, 2000 June 14:
	Update a comment.
*/

typedef struct {

	/* image name */
	char input[STIS_LINE+1];

	char rootname[STIS_CBUF+1];	/* root name for set of obs */

	int printtime;			/* print time after each step? */
	int verbose;			/* print additional info? */

	/* info about this image */
	char sclamp[STIS_CBUF+1];	/* spectral calibration lamp in use */
	char obsmode[STIS_CBUF+1];	/* e.g. ACCUM or TIMETAG */
	char aperture[STIS_CBUF+1];
	char opt_elem[STIS_CBUF+1];
	char det[STIS_CBUF+1];		/* name of detector */
	int nimages;			/* number of "groups" */
	int detector;			/* integer code for detector */
	int cenwave;			/* central wavelength */

	/* values specific to each group */
	double *midpt;			/* midpoint of exposure (MJD) */
	double *shift1, *shift2;	/* shift read from wavecal */

} StisInfo12;

#endif /* INCL_CALSTIS12_H */
