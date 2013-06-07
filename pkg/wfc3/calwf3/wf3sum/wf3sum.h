/* wf3sum.h Sum Repeatobs data */

/* include "../wf3.h" */
# include "c_iraf.h" 	/* Needed for definition of IRAFPointer */
# include "ximio.h"	/* Needed for c_imt* functions */

/* Image description for wf3sum.  */

typedef struct {

	char **input;			/* input images */
	char output[SZ_LINE+1];		/* output sum */

	char rootname[SZ_CBUF+1];	/* root name for set of obs */

	/* info about input image */
	char obsmode[SZ_CBUF+1];	/* e.g. ACCUM */
	char aperture[SZ_CBUF+1];	/* aperture name */
	char filter[SZ_CBUF+1];		/* name of filter or grating */
	char det[SZ_CBUF+1];		/* name of detector */
	int detector;			/* detector variable */

	int printtime;		/* print time after each step? */
	int verbose;		/* print additional info? */

	int nimages;		/* number of images to be summed */
	int nimsets;		/* number of imsets per image */
	short sdqflags;		/* serious data quality flags */
	double exptime;		/* exposure time */
	double expend;		/* end of exposure time */

} Wf3SumInfo;

