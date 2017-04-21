# include <stdio.h>
# include <string.h>

# include "hstio.h"
# include "wf3.h"
# include "wf3info.h"
# include "err.h"
# include "cte.h"


/* This routine subtracts the special BIACFILE image from x (in-place).
   For WF3 CTE corrected science data.

   This is a special bias image created for the CTE correction and
   follows the same format as other reference bias images. The bias subtraction
   is done before the RAZ image is made. Instead of calling the regular routine.
   The regular bias subtraction is still performed.

   Megan Sosey, 2015

 */
int sub1dreform (SingleGroup *, int, int, SingleGroupLine *);

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
	int dimx, dimy;
	int straddle=0;
	int overstart=0;
	int offsetx=0;
	int offsety=0;

	dimx = x->sci.data.nx;
	dimy = x->sci.data.ny;

	sprintf(MsgText,"CTE: Subtracting BIACFILE: %s for imset %d",wf3->biac.name, x->group_num);
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

	 if (wf3->verbose){
 		sprintf(MsgText,"Image has starting location of %d,%d in the reference image",x0,y0);
 		trlmessage(MsgText);
 	}

	 /* Subtract the bias image from x. */

	 /* If the science image is binned, it will be assumed to have
	    same size as the reference image, since reading subarrays
	    of a binned chip is not supported in current flight software.
	  */
	 if (same_size) {

	 	/* Loop over all the lines in the science image */
	 	for (i=0; i < scilines; i++) {
	 		status = getSingleGroupLine (wf3->biac.name, i, &y);
	 		if (status) {
	 			sprintf(MsgText,"Could not read line %d from bias image.",
	 					i+1);
	 			trlerror(MsgText);
	 		}

	 		/* No trimming required. */
	 		status = sub1d(x, i, &y);
	 		if (status) {
	 			trlerror ("(biascorr) size mismatch.");
	 			return (status);
	 		}
	 	}

	 } else {

		 if (x0 >= 2072){ /*image starts in B or D regions and we can just shift the starting pixel*/
		 	if (wf3->verbose){
		 		sprintf(MsgText,"Subarray starts in B or D region, moved from (%d,%d) to ",x0,y0);
		 		trlmessage(MsgText);
		 	}
		 		x0 += 60;
		 	if (wf3->verbose){
		 		sprintf(MsgText,"(%d,%d) to avoid virtual overscan in reference",x0,y0);
		 		trlmessage(MsgText);
		 	}
		 } else { /*the subarray starts somewhere in A or C and might straddle the virtual overscan region */

		 	if ( (x0 + dimx) >= 2072){
		 		straddle=1;
		 		overstart=2073-x0;
		 	}

		 }

		 if (wf3->verbose){
	 	    sprintf(MsgText,"ccdamp=%s, straddle=%d, offset=(%d,%d),ampx,ampy=(%d,%d),x0,y0=(%d,%d)",wf3->ccdamp,straddle,offsetx,offsety,wf3->ampx,wf3->ampy,x0,y0);
	         trlmessage(MsgText);
	     }

	 	/* Loop over all the lines in the science array, and
	 	 ** match them to the appropriate line in the reference image. */
	 	/*
	 	   i - index for line in science image
	 	   j - index for line in reference image
	 	   y0 - line in reference image corresponding to
	 	   line in input image
	 	 */
		initSingleGroupLine (&z);
 	    if (straddle){
 	    	allocSingleGroupLine (&z, x->sci.data.nx+60);
 	    } else {
 	        allocSingleGroupLine (&z, x->sci.data.nx);
 	    }

	 	for (i=0, j=y0; i < scilines; i++,j++) {

	 		/* We are working with a sub-array and need to apply the
	 		 ** proper section from the reference image to the science
	 		 ** image.  */
	 		status = getSingleGroupLine (wf3->biac.name, j, &y);
	 		if (status) {
	 			sprintf (MsgText,"Could not read line %d from biac image.",
	 					j+1);
	 			trlerror(MsgText);
	 		}

	 		update = NO;

		    if (trim1d (&y, x0, j, rx, avg, update, &z)) {
				trlerror ("(ctebiascorr) reference file size mismatch.");
				return (status);
		    }
			if (straddle) {
				status = sub1dreform (x, i, overstart,&z);
			} else {
				status = sub1d (x, i, &z);
			}
			if (status)
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
