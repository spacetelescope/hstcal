#ifndef INCL_ACSSUM_H
#define INCL_ACSSUM_H

/* acssum.h Sum Repeatobs data */

/* include "../acs.h" */
# include <c_iraf.h> 	/* Needed for definition of IRAFPointer */
# include <ximio.h>		/* Needed for c_imt* functions */

/* Image description for acssum.

*/

typedef struct {

	char **input;				/* input images */
	char output[CHAR_LINE_LENGTH];		/* output sum */

	char rootname[ACS_CBUF];	/* root name for set of obs */

	/* info about input image */
	char obsmode[ACS_CBUF];	/* e.g. ACCUM or TIMETAG */
	char aperture[ACS_CBUF];	/* aperture name */
	char filter1[ACS_CBUF];	/* name of filter or grating */
	char filter2[ACS_CBUF];	/* name of filter or grating */
	char det[ACS_CBUF];		/* name of detector */
	int detector;			/* detector variable */

	int printtime;		/* print time after each step? */
	int verbose;		/* print additional info? */

	int nimages;		/* number of images to be summed */
	int nimsets;		/* number of imsets per image */
	short sdqflags;		/* serious data quality flags */
	double exptime;		/* exposure time */
	double expend;		/* end of exposure time */

} AcsSumInfo;

#endif /* INCL_ACSSUM_H */
