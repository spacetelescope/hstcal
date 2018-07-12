# include <stdio.h>
# include <stdlib.h>		/* calloc */
# include <string.h>		/* strncmp */
# include <math.h>		/* sqrt */
# include "c_iraf.h"
# include "hstio.h"

# include "stis.h"
# include "hstcalerr.h"	/* SIZE_MISMATCH */
# include "stisdef.h"

/* This routine takes an input data array, extracts a subset, bins it
   by averaging within rectangular bins, and assigns the values to an
   output data array.  The calling routine must allocate the output
   SingleGroup (setting its size) and free it when done.
   The coordinate keywords in the output extension headers will be updated.

   Phil Hodge, 1998 Oct 5:
	Change status value 1001 to HEADER_PROBLEM.
*/

int bin2d (SingleGroup *a, int xcorner, int ycorner, int binx, int biny,
	int avg, SingleGroup *b) {

/* arguments:
SingleGroup *a        i: input data
int xcorner, ycorner  i: starting location of subarray (zero indexed,
			input pixel coordinates)
int binx, biny        i: number of input pixels for one output pixel
int avg               i: == 0 means we should sum the values within a bin;
			> 0 means we should average the values
SingleGroup *b        o: output data
*/

	int status;

	double block[2];	/* number of input pixels for one output */
	double offset[2];	/* offset of binned image */
	float weight;		/* binx * biny (or sqrt for errors) */
	float sum;		/* for summing data values */
	float sum_err;		/* for summing error array */
	short sum_dq;		/* for ORing data quality array */
	int nx, ny;		/* size of output array */
	int m, n;		/* pixel index in output array */
	int i, j;		/* pixel index in input array */
	int i0, j0;		/* starting location of (m,n) in input */

	nx = b->sci.data.nx;
	ny = b->sci.data.ny;

	weight = binx * biny;

	block[0] = binx;
	block[1] = biny;
	offset[0] = xcorner;
	offset[1] = ycorner;

	/* Check the subset of the input corresponding to the output. */
	if (xcorner < 0 || ycorner < 0 ||
		xcorner + nx*binx > a->sci.data.nx ||
		ycorner + ny*biny > a->sci.data.ny) {
	    printf ("ERROR    (bin2d)  subset is out of bounds: \\\n");
	    printf ("         input is %d x %d, output is %d x %d \\\n",
		a->sci.data.nx, a->sci.data.ny,
		b->sci.data.nx, b->sci.data.ny);
	    printf ("         start = (%d,%d), binx = %d, biny = %d.\n",
		xcorner+1, ycorner+1, binx, biny);
	    return (SIZE_MISMATCH);
	}

	if (binx == 1 && biny == 1) {

	    for (n = 0, j = ycorner;  n < ny;  n++, j++)
		for (m = 0, i = xcorner;  m < nx;  m++, i++)
		    Pix (b->sci.data, m, n) = Pix (a->sci.data, i, j);

	    /* Extract the error array. */
	    for (n = 0, j = ycorner;  n < ny;  n++, j++)
		for (m = 0, i = xcorner;  m < nx;  m++, i++)
		    Pix (b->err.data, m, n) = Pix (a->err.data, i, j);

	    /* Extract the data quality array. */
	    for (n = 0, j = ycorner;  n < ny;  n++, j++)
		for (m = 0, i = xcorner;  m < nx;  m++, i++)
		    DQSetPix (b->dq.data, m, n, DQPix(a->dq.data,i,j));

	} else {

	    /* Average the science data array. */

	    j0 = ycorner;				/* zero indexed */
	    /* for each line of output ... */
	    for (n = 0;  n < ny;  n++) {
		i0 = xcorner;				/* zero indexed */
		/* for each sample in current output line ... */
		for (m = 0;  m < nx;  m++) {
		    /* loop over each pixel in input for current output pixel */
		    sum = 0.;
		    for (j = j0;  j < j0+biny;  j++)
			for (i = i0;  i < i0+binx;  i++)
			    sum += Pix (a->sci.data, i, j);
		    if (avg)
			Pix (b->sci.data, m, n) = sum / weight;
		    else
			Pix (b->sci.data, m, n) = sum;
		    i0 += binx;
		}
		j0 += biny;
	    }

	    /* Average the error array (square root of sum of squares). */

	    j0 = ycorner;
	    for (n = 0;  n < ny;  n++) {
		i0 = xcorner;
		for (m = 0;  m < nx;  m++) {
		    sum_err = 0.;
		    for (j = j0;  j < j0+biny;  j++)
			for (i = i0;  i < i0+binx;  i++)
			    sum_err += Pix (a->err.data, i, j) *
				       Pix (a->err.data, i, j);
		    if (avg)
			Pix (b->err.data, m, n) = sqrt (sum_err) / weight;
		    else
			Pix (b->err.data, m, n) = sqrt (sum_err);
		    i0 += binx;
		}
		j0 += biny;
	    }

	    /* Bitwise OR the data quality array. */
	    j0 = ycorner;
	    for (n = 0;  n < ny;  n++) {
		i0 = xcorner;
		for (m = 0;  m < nx;  m++) {
		    sum_dq = 0;				/* zero is OK */
		    for (j = j0;  j < j0+biny;  j++)
			for (i = i0;  i < i0+binx;  i++)
			    sum_dq |= DQPix (a->dq.data, i, j);
		    DQSetPix (b->dq.data, m, n, sum_dq);
		    i0 += binx;
		}
		j0 += biny;
	    }
	}

	/* Copy the headers. */
	copyHdr (b->globalhdr, a->globalhdr);
	if (hstio_err())
	    return (HEADER_PROBLEM);
	copyHdr (&b->sci.hdr, &a->sci.hdr);
	if (hstio_err())
	    return (HEADER_PROBLEM);
	copyHdr (&b->err.hdr, &a->err.hdr);
	if (hstio_err())
	    return (HEADER_PROBLEM);
	copyHdr (&b->dq.hdr, &a->dq.hdr);
	if (hstio_err())
	    return (HEADER_PROBLEM);

	/* Update the coordinate parameters that depend on the binning. */
	if ((status = BinCoords (&a->sci.hdr, block, offset,
                                 &b->sci.hdr, &b->err.hdr, &b->dq.hdr)))
	    return (status);

	return (0);
}
