/* This file contains:
	doBlev
	BlevOpen
	FitToOverscan
*/

# include <stdio.h>
# include <math.h>		/* sqrt */

# include "c_iraf.h"
# include "hstio.h"
# include "stis.h"
# include "calstis1.h"
# include "hstcalerr.h"
# include "stisdq.h"		/* for CALIBDEFECT */
# include "stisdef.h"

static FILE *BlevOpen (char *, int, int *);
static void FitToOverscan (SingleGroup *, short, int, int, int *, double);

/* This routine subtracts the bias determined from the overscan region
   on each line.  Only a portion of the trailing overscan region will be
   used to find the bias level.

   If an output text file for bias levels (third command-line argument)
   was specified, this routine also creates this file, writes the bias
   level for each image line, and then closes the file.

   The CRPIX and LTV keyword values in the output extension headers will be
   updated.

   NOTE:  If the input SingleGroup x is large enough to have overscan
   regions, a new SingleGroup out will be allocated that is smaller
   (because it will not have the overscan regions), a subset of x will be
   copied to out, x will be freed, and out (just the pointer) will be
   copied to x; done will be set to one.  If x is too small, however,
   the default ccdbias will be subtracted from x in-place, out will not be
   allocated, and x will not be freed; done will be set to zero.

   Phil Hodge, 1997 Nov 4:
	Change calling sequence of FindOverscan; move a section to
	FitToOverscan; call BlevDrift and DriftEval; don't write npix
	to the outblev file since npix and biaslevel are not both known
	at the same time.

   Phil Hodge, 1998 June 10:
	Move the declaration of FitToOverscan outside doBlev, since
	it's static.

   Phil Hodge, 1998 Oct 16:
	Change status value 1001 to GENERIC_ERROR_CODE.
	Include sdqflags in the calling sequences of FitToOverscan,
	FindBlev, and BlevDrift.

   Phil Hodge, 1999 Nov 2:
	Second argument is now SingleGroup **x instead of SingleGroup *x.
	Remove SingleGroup *y from calling sequence.

   Phil Hodge, 2004 July 28:
	Add blev_clip to calling sequence of BlevDrift.
*/

int doBlev (StisInfo1 *sts, SingleGroup **x, int extver,
		float *meanblev, int *done, int *driftcorr) {

/* arguments:
StisInfo1 *sts    i: calibration switches, etc
SingleGroup **x   io: input image to be calibrated; output stripped image
int extver        i: imset number
float *meanblev   o: mean value of bias levels that were subtracted
int *done         o: true if x has been copied to y, and x freed
int *driftcorr    o: true means correction for drift along lines was applied
*/

	SingleGroup *in;	/* local variable, for cleaner notation */
	SingleGroup *out;	/* copy of x with overscan stripped off */
	int status;

	double biaslevel;	/* value to subtract from a line of image */
	double averagedrift;	/* drift averaged over all (output) columns */
	double sumbias = 0.;	/* sum of bias levels (for getting mean) */
	int binx, biny;		/* bin factors */
	int trimx1, trimx2;	/* width to trim off ends of each line */
	int trimy1, trimy2;	/* amount to trim off ends of each column */
	int biassect[2];	/* section to use for finding bias level */
	int vx[2], vy[2];	/* location of virtual overscan region */
	int i, j;
	int avg = 0;		/* used by bin2d, but not signficant */
	FILE *ofp;		/* output file fp for bias levels */

	double BlevEval (int);
	int FindOverscan (char *, int, int, int *,
		int *, int *, int *, int *, int *,
		int *, int *);
	int BlevDrift (SingleGroup *, short, int *, int *, int, int *,
		float, int *);
	double DriftEval (int);
	double DriftMean (int);

	/* initial values */
	*done = 0;
	*driftcorr = 0;

	in = *x;		/* just for convenience of notation */

	binx = sts->bin[0];
	biny = sts->bin[1];

	if ((binx != 1 && binx != 2 && binx != 4) ||
	    (biny != 1 && biny != 2 && biny != 4)) {
	    printf ("ERROR    (doBlev) bin size must be 1, 2, or 4.\n");
	    return (GENERIC_ERROR_CODE);
	}

	/* Get biassect, the location of the region to use for determining
	   the bias level, and vx & vy, the location in the virtual
	   overscan to use for determining the drift with column number.
	   Also get the widths of the trim regions, the amount to remove
	   from each axis of the input when copying to output.
	*/
	if (FindOverscan (sts->ccdamp,
		in->sci.data.nx, in->sci.data.ny, sts->bin,
		&trimx1, &trimx2, &trimy1, &trimy2, biassect, vx, vy)) {

	    printf ("Warning  Image size is too small to do BLEVCORR; \\\n");
	    printf ("Warning  bias from CCDTAB will be subtracted.\n");
	    biaslevel = sts->ccdbias;
	    for (j = 0;  j < in->sci.data.ny;  j++) {
		for (i = 0;  i < in->sci.data.nx;  i++)
		    Pix (in->sci.data,i,j) = Pix (in->sci.data,i,j) - biaslevel;
	    }
	    *meanblev = biaslevel;
	    *done = 0;			/* x not freed; y not allocated */
	    return (0);
	}

	/* Create a new data set of the appropriate (smaller) size. */
	if ((out = calloc (1, sizeof (SingleGroup))) == NULL)
	    return (OUT_OF_MEMORY);
	initSingleGroup (out);
	allocSingleGroup (out,
		in->sci.data.nx - (trimx1 + trimx2),
		in->sci.data.ny - (trimy1 + trimy2), True);
	if (hstio_err())
	    return (OUT_OF_MEMORY);

	/* Copy out the non-overscan portion of x, putting the result in out.
	*/
	if ((status = bin2d (in, trimx1, trimy1, 1, 1, avg, out)))
	    return (status);

	/* Create an output file for individual bias levels, if requested. */
	ofp = BlevOpen (sts->outblev, extver, &status);
	if (status)
	    return (status);
	if (ofp != NULL)
	    fprintf (ofp, "# line  bias   extver=%d\n", extver);

	/* For each image line, determine the bias level from the overscan
	   in x, and fit to these values as a function of line number in
	   the output image.
	*/
	FitToOverscan (in, sts->sdqflags,
		out->sci.data.ny, trimy1, biassect, sts->ccdbias);

	/* Fit a line to the virtual overscan region as a function of
	   column number.
	*/
	if ((status = BlevDrift (in, sts->sdqflags,
                vx, vy, trimx1, biassect, sts->blev_clip, driftcorr)))
	    return (status);

	/* Evaluate the fit for each line, and subtract from the data. */
	averagedrift = DriftMean (out->sci.data.nx);
	for (j = 0;  j < out->sci.data.ny;  j++) {
	    biaslevel = BlevEval (j);
	    /* bias for current line plus average drift (constant) */
	    sumbias += (biaslevel + averagedrift);
	    for (i = 0;  i < out->sci.data.nx;  i++)
		Pix (out->sci.data,i,j) =
			Pix (out->sci.data,i,j) - biaslevel - DriftEval (i);
	    if (ofp != NULL)
		fprintf (ofp, "%6d  %0.6g\n", j+1, biaslevel + averagedrift);
	}

	if (ofp != NULL)
	    fclose (ofp);			/* done with output file */

	/* This is the mean value of all the bias levels subtracted. */
	*meanblev = sumbias / out->sci.data.ny;

	/* Free x, and copy out to x. */
	freeSingleGroup (*x);
	free (*x);
	*x = out;
	*done = 1;			/* overscan was actually removed */

	return (0);
}

/* This routine opens an output text file if the name is non-null. */

static FILE *BlevOpen (char *fname, int extver, int *status) {

	FILE *fp;			/* temp for value to return */

	*status = 0;

	if (fname[0] == '\0') {

	    fp = NULL;				/* no output file specified */

	} else if (extver == 1) {

	    fp = fopen (fname, "r");	/* does the file already exist? */

	    if (fp == NULL) {

		if ((fp = fopen (fname, "w")) == NULL) {
		    printf ("ERROR    Could not open %s for writing.\n", fname);
		    *status = OPEN_FAILED;
		}

	    } else {

		fclose (fp);
		printf ("ERROR    File %s already exists.\n", fname);
		*status = OPEN_FAILED;
		fp = NULL;
	    }

	} else {

	    if ((fp = fopen (fname, "a")) == NULL) {
		printf ("ERROR    Could not append to %s.\n", fname);
		*status = OPEN_FAILED;
	    }
	}

	return (fp);
}

/* This function determines the bias level from the overscan in the input
   image x, and it fits a straight line to these values as a function of
   line number in the output image.
*/

static void FitToOverscan (SingleGroup *x, short sdqflags,
	int ny, int trimy1,
	int *biassect, double ccdbias) {

/* arguments:
SingleGroup *x    i: input image
short sdqflags    i: "serious" data quality flags
int ny            i: size of second axis after overscan has been removed
int trimy1        i: offset between input and output line numbers
int biassect[2]   i: section to use for bias level determination
double ccdbias    i: bias level to subtract if we can't get it from overscan
*/

	double biaslevel;	/* bias level in one line of input image */
	int too_few = 0;	/* number of lines with too few good pixels */
	int npix;		/* number of pixels used to compute bias */
	int j;
	void BlevInit (int);
	void BlevAccum (int, double);
	int BlevFit (void);
	void BlevSet (double);
	int FindBlev (SingleGroup *, short, int, int *, double *, int *);

	BlevInit (ny / 2);			/* initialize for fitting */

	/* For each line, determine the bias level from the overscan in x.
	   Note that we loop over the number of pixels in the output image.
	   The argument j+trimy1 to FindBlev is the line number in the
	   input image x, with trimy1 the offset to the illuminated region.
	*/
	for (j = 0;  j < ny;  j++) {
	    if (FindBlev (x, sdqflags, j+trimy1, biassect, &biaslevel, &npix)) {
		too_few++;			/* not fatal */
	    } else {
		/* Note:  j is line number in output image. */
		BlevAccum (j, biaslevel);	/* increment sums */
	    }
	}
	if (too_few > 0) {
	    printf ("Warning  (blevcorr) %d image line", too_few);
	    if (too_few == 1)
		printf (" has");
	    else
		printf ("s have");
	    printf (" too few overscan pixels.\n");
	}

	/* Fit a curve to the bias levels found. */
	if (BlevFit()) {
	    printf ("Warning  No bias level data, or singular fit; \\\n");
	    printf ("Warning  bias from CCDTAB will be subtracted.\n");
	    BlevSet (ccdbias);		/* assign the default value */
	}
}
