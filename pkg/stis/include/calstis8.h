#ifndef INCL_CALSTIS8_H
#define INCL_CALSTIS8_H

/* calstis8.h 2-D spectral extraction */

/* include "../stis.h" */

/* Image description for calstis8.

   Phil Hodge, 1998 Jan 22:
	Add echelle and nrptexp to the StisInfo8 struct.

   Phil Hodge, 1998 Aug 4:
	Add one to the length of each string buffer.
*/

typedef struct {

	char input[STIS_LINE+1];	/* input multi-group image */
	char output[STIS_LINE+1];	/* output sum */

	char rootname[STIS_CBUF+1];	/* root name for set of obs */

	/* info about input image */
	char obsmode[STIS_CBUF+1];	/* e.g. ACCUM or TIMETAG */
	char aperture[STIS_CBUF+1];	/* aperture name */
	char opt_elem[STIS_CBUF+1];	/* name of grating or mirror */
	char det[STIS_CBUF+1];		/* name of detector */

	int printtime;		/* print time after each step? */
	int verbose;		/* print additional info? */

	int echelle;		/* true if echelle data */

	int nimages;		/* number of imsets in file */
	int nrptexp;		/* number of imsets to add together */
	short sdqflags;		/* serious data quality flags */
	double exptime;		/* exposure time */

	/* calibration switches */
	int fluxcorr;		/* values have been converted to flux? */
	int statcorr;		/* compute statistics? */

} StisInfo8;

#endif /* INCL_CALSTIS8_H */
