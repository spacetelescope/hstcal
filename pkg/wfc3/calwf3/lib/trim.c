# include <stdlib.h>		/* calloc */
# include <string.h>		/* strncmp */
# include <math.h>		/* sqrt */

# include "hstio.h"

# include "wf3.h"
# include "hstcalerr.h"	/* SIZE_MISMATCH */

/* This routine takes an input data array, bins it in 1d, extracts a subset, 
   and assigns the values to an
   output data array.  The calling routine must allocate the output
   SingleGroupLine (setting its size) and free it when done.
   The coordinate keywords in the output extension headers will be updated.
   
** Need to add 'ystart' as input to keep header information correct in output
    array.
    
    This function does not treat any expansion or reduction in Y!
*/

int trim1d (SingleGroupLine *a, int xstart, int ystart, int binx, int avg,
	    int update, SingleGroupLine *b) {

/* arguments:
SingleGroupLine *a    i: input data
int xstart	      i: starting location of subarray (zero indexed,
			 input pixel coordinates)
int ystart	      i: starting location of subarray (zero indexed,
			 input pixel coordinates)
int binx	      i: number of input pixels for one output pixel
int avg               i: == 0 means we should sum the values within a bin;
			 > 0 means we should average the values
int updated	      i: == 1 means have already updated the header info
			 == 0 means we need to update the header info
SingleGroupLine *b    o: output data
*/

	extern int status;

	double block[2];	/* number of input pixels for one output */
	double offset[2];	/* offset of binned image */
	int nx;			/* size of output array */
	int m;			/* pixel index in output array */
	int i;			/* pixel index in input array */

	int BinCoords (Hdr *, double *, double *, Hdr *, Hdr *, Hdr *);

	/* Initialize internal variables... */
	nx = b->sci.tot_nx;
	block[0] = binx;
	block[1] = 1;
	offset[0] = xstart;
	offset[1] = ystart;  /* Need to modify to pass ystart */
	
	/* Check the subset of the input corresponding to the output. */
	if (xstart < 0 || xstart + nx*binx > a->sci.tot_nx) {
	    trlerror ("(trim1d)  subset is out of bounds:");
	    sprintf (MsgText,"         input is %d pixels, output is %d pixels",
		     a->sci.tot_nx, b->sci.tot_nx);
	    trlmessage (MsgText);
	    sprintf (MsgText, "         start = (%d), binx = %d.", xstart+1,
		     binx);
	    trlmessage (MsgText);
	    return (status = SIZE_MISMATCH);
	}

	for (m = 0, i = xstart;  m < nx;  m++, i++) {
	     b->sci.line[m] = a->sci.line[i];

	     /* Extract the error array. */
	     b->err.line[m] = a->err.line[i];

	     /* Extract the data quality array. */
             b->dq.line[m] = a->dq.line[i];
	}

	/* Copy header info only if requested... */
	if (update == YES) { 

	    /* Copy the headers. */
	    copyHdr (b->globalhdr, a->globalhdr);
	    if (hstio_err())
		return (status = 1001);
	    copyHdr (&b->sci.hdr, &a->sci.hdr);
	    if (hstio_err())
		return (status = 1001);
	    copyHdr (&b->err.hdr, &a->err.hdr);
	    if (hstio_err())
		return (status = 1001);
	    copyHdr (&b->dq.hdr, &a->dq.hdr);
	    if (hstio_err())
		return (status = 1001);

	    /* Update the coordinate parameters that depend on the binning. */
	    if (BinCoords (&a->sci.hdr, block, offset, &b->sci.hdr,
			   &b->err.hdr, &b->dq.hdr))
		return (status);
	}
	return (status);
}
