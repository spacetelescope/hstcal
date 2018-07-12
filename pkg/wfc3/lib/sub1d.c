# include <math.h>
# include <stdio.h>
# include "hstio.h"
# include "hstcalerr.h"	/* SIZE_MISMATCH */
# include "wf3.h"


/* Subtract the second SingleGroupLine from the first, leaving the
   result in the first.

   (*a) -= (*b)

   The science data arrays are subtracted; the error arrays are combined;
   the data quality arrays are combined.
	
    
   M. Sosey, 2013 Sept 09
   Added new routine to deal with uvis subarrays  for postflash correction
   which knows what to do with different reference file / science sizes and the
   uvis overscan region 
*/

int sub1d (SingleGroup *a, int line, SingleGroupLine *b) {

/* arguments:
SingleGroup *a		io: input data; output difference
int line		 i: line of input data to subtract 1-d data from
SingleGroupLine *b	 i: second input data
*/

	extern int status;

	int i;
	float da, db;		/* errors for a and b */
	short dqa, dqb, dqab;	/* data quality for a, b, combined */
	int dimx;

	if (a->sci.data.nx != b->sci.tot_nx)
	    return (status = SIZE_MISMATCH);

	/* science, error, and DQ data */
	dimx = a->sci.data.nx;
	for (i = 0;  i < dimx;  i++) {

    	     /* science array */
	     Pix(a->sci.data, i, line) =
			Pix (a->sci.data, i, line) - b->sci.line[i];

    	     /* error array */
	     da = Pix (a->err.data, i, line);
	     db = b->err.line[i];
	     Pix (a->err.data, i, line) = sqrt (da * da + db * db);

    	     /* data quality */
	     dqa = DQPix (a->dq.data, i, line);
	     dqb = b->dq.line[i];
	     dqab = dqa | dqb;
	     DQSetPix (a->dq.data, i, line, dqab);
	}

	return (status);
}

int sub1dreform (SingleGroup *a, int line, int overstart, SingleGroupLine *b) {

/* arguments:
SingleGroup *a		io: input data; output difference
int line		 i: line of input data to subtract 1-d data from
SingleGroupLine *b	 i: second input data
int overstart : this is where the overscan starts in the cut image

This is the same as sub1d, except that it is meant to be used with the postflash
correction for subarrays only. This is because the reference file is always a full frame 4-amp readout.
GO users are restricted to subarrays which are entirely contained in one amp (C amp I think). Group
members can place subarrays anywhere for calibration programs and have cases where the subrarray bisects
the virtual overscan region in the postflash reference file. This means that just that striped section needs to
be removed from the reference image and a new stitched line created for subtraction

So in this case, the second single group line is larger than the first and the virtual overscan
strip is removed from the data so that they are the same size before the subtraction

In this case the size of *a does NOT match the size of *b
*/

	extern int status;

	int i,j;
	float da, db;		/* errors for a and b */
	short dqa, dqb, dqab;	/* data quality for a, b, combined */
    int sizea;
    int sizeb;
    
    sizea=a->sci.data.nx;
    sizeb=b->sci.tot_nx;    
    

	/* science, error, and DQ data */
	for (i=0,j=0;  i < sizea;  i++,j++) {
        if (i == overstart){
            j+=60;
         }
    	     /* science array */
	     Pix(a->sci.data, i, line) =
			Pix (a->sci.data, i, line) - b->sci.line[j];

    	     /* error array */
	     da = Pix (a->err.data, i, line);
	     db = b->err.line[j];
	     Pix (a->err.data, i, line) = sqrt (da * da + db * db);

    	     /* data quality */
	     dqa = DQPix (a->dq.data, i, line);
	     dqb = b->dq.line[j];
	     dqab = dqa | dqb;
	     DQSetPix (a->dq.data, i, line, dqab);
	}

	return (status);
}

