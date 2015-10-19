# include <stdio.h>
# include <string.h>

# include "hstio.h"
# include "wf3.h"
# include "wf3info.h"
# include "wf3err.h"
# include "cte.h"


/* This routine subtracts the special BIACFILE image from x (in-place).
   For WF3 CTE corrected science data.

   This is a special bias image created for the CTE correction and 
   follows the same format as other reference bias images. The bias subtraction
   is done before the RAZ image is made. Instead of calling the regular routine.
   The regular bias subtraction is still performed.
   
   Megan Sosey, 2015

 */

int doCteBias (WF3Info *wf3, SingleGroup *x) {

	/* arguments:
	   WF3Info *wf3     i: calibration switches, etc
	   SingleGroup *x	io: image to be calibrated; written to in-place
	 */

	extern int status;

	SingleGroupLine y, z;	/* y and z are scratch space */
	int rx, ry;		/* for binning biac image down to size of x */
	int x0, y0;		/* offsets of sci image */
	int same_size;		/* true if no binning of ref image required */
	int avg = 0;		/* bin2d should sum within each bin */
	int scilines; 		/* number of lines in science image */
	int i, j;
	int update;
    

	sprintf(MsgText,"CTE: Subtracting BIACFILE: %s for group %d",wf3->biac.name, x->group_num);
	trlmessage(MsgText);

	initSingleGroupLine (&y);
	scilines = x->sci.data.ny;


	/* Initialize local variables */
	rx = 1;
	ry = 1;
	x0 = 0;
	y0 = 0;
	same_size = 1;


	/* Get the first line of biac image data. */    
	openSingleGroupLine (wf3->biac.name, x->group_num, &y);
	if (hstio_err())
		return (status = OPEN_FAILED);

	/* 
	   Reference image should already be selected to have the
	   same binning factor as the science image.  All we need to
	   make sure of is whether the science array is a sub-array of
	   the biac image.  
	 */

	if (FindLine (x, &y, &same_size, &rx, &ry, &x0, &y0))
		return (status);

	/* Return with error if reference data not binned same as input */
	if (rx != 1 || ry != 1) {
		closeSingleGroupLine (&y);
		freeSingleGroupLine (&y);
		sprintf (MsgText,
				"BIAC image and input are not binned to the same pixel size!");
		trlerror (MsgText);
		return (status = SIZE_MISMATCH);
	}

	/* Subtract the biasc image from x. */

	/* If the science image is binned, it will be assumed to have
	   same size as the reference image, since reading subarrays
	   of a binned chip is not supported in current flight software.
	 */
	if (same_size) {

		/* Loop over all the lines in the science image */
		for (i=0; i < scilines; i++) {
			if ( getSingleGroupLine (wf3->biac.name, i, &y)) {
				sprintf(MsgText,"Could not read line %d from biac image.",
						i+1);
				trlerror(MsgText);
			}

			/* No trimming required. */
			if (sub1d(x, i, &y)) {
				trlerror ("(BIACFILE) size mismatch.");
				return (status);
			}
		}

	} else {

		/* Loop over all the lines in the science array, and
		 ** match them to the appropriate line in the reference image. */
		/* 
		   i - index for line in science image
		   j - index for line in reference image
		   y0 - line in reference image corresponding to
		   line in input image
		 */
		initSingleGroupLine (&z);
		allocSingleGroupLine (&z, x->sci.data.nx);
		for (i=0, j=y0; i < scilines; i++,j++) { 

			/* We are working with a sub-array and need to apply the
			 ** proper section from the reference image to the science
			 ** image.  */
			if (getSingleGroupLine (wf3->biac.name, j, &y)) {
				sprintf (MsgText,"Could not read line %d from biac image.",
						j+1);
				trlerror(MsgText);
			}			

			update = NO;

			if (trim1d (&y, x0, y0, rx, avg, update, &z)) {
				trlerror ("(BIACFILE) size mismatch.");
				return (status);
			}

			if (sub1d (x, i, &z))
				return (status);
		}
		freeSingleGroupLine (&z);			/* done with z */
	}

	closeSingleGroupLine (&y);
	freeSingleGroupLine (&y);

	if(wf3->printtime){
		TimeStamp("Finished subtracting BIAC file: ",wf3->rootname);
	}
    
	return (status);
}
