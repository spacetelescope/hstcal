# include <stdlib.h>
# include <string.h>
# include "err.h"

/* If the slit is narrower than this, the spectrum does not need to be
   convolved with the slit.
*/
# define NARROW_SLIT  3

/* This routine convolves the template spectrum with a rectangle of about
   the same size as the slit that was used for the observation.

   Phil Hodge, 2000 Jan 5:
	Extracted from xcwave.c.
*/

int ConvSlit (double slitwidth, double tspec[], int nwl) {

/* arguments:
double slitwidth    i: width of slit (dispersion direction) in pixels
double tspec[nwl]   io: template spectrum, convolved in place
int nwl             i: size of tspec array
*/

	double *buf;	/* scratch for template spectrum convolved with slit */
	double sum;	/* sum of template spectrum over slit at one pixel */
	int iwidth;	/* slit width, odd integer */
	int half;	/* half of iwidth, truncated down */
	int i;		/* loop index in buf */
	int j;		/* loop index in tspec */

	/* iwidth should be the nearest odd integer to slitwidth. */
	iwidth = (int) slitwidth;
	if (iwidth / 2 * 2 == iwidth)
	    iwidth++;

	if (iwidth < NARROW_SLIT)
	    return (0);

	half = iwidth / 2;		/* truncate */

	if ((buf = calloc (nwl + 2*half, sizeof (double))) == NULL)
	    return (OUT_OF_MEMORY);

	memcpy (&buf[half], tspec, nwl*sizeof(double));

	sum = 0.;
	for (i = 0;  i < iwidth;  i++)
	    sum += buf[i];
	tspec[0] = sum / (double)iwidth;

	/* Convolve with the slit.  The normalization results in the ends
	   being tapered off, but that should be OK.
	*/
	for (j = 1;  j < nwl;  j++) {
	    i = j + half;
	    sum += (buf[i+half] - buf[i-half-1]);
	    tspec[j] = sum / (double)iwidth;
	}

	free (buf);

	return (0);
}
