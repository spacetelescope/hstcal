# include <stdio.h>

# include "xtables.h"
# include "hstio.h"

# include "stis.h"
# include "hstcalerr.h"
# include "calstis6.h"



/*
   RebinData  -  This function handles the case of hires data.

   The IDT algorithm works with lowres data only. In case of hires
   input data, it is rebinned to a work array/table, the algorithm
   is applied to the work array/table, and at the end the image data
   is rebinned back to hires and then the final 1-D data is extracted.
   In case of lowres data, it is just copied to the output array/table.

   The 'rebinfactor' parameter can have three possible values:

   1 - No rebinning, just copy.
   2 - Squeeze by a factor 2.
  -2 - Expand by a factor 2.

   The SingleGroup output image and RowContents output array must be
   allocated by the caller. This function will take care of allocating
   memory for the internal arrays in each RowContents element.


   Revision history:
   ----------------
   23 Feb 00  -  Implemented (I.Busko)
*/

int RebinData (SingleGroup *iim, SingleGroup *oim, RowContents **ix1d,
               RowContents **ox1d , int rebinfactor, int nrows) {


/* arguments
SingleGroup *iim;		i: input image
SingleGroup *oim;		o: output image
RowContents **ix1d;		i: input array with table data.
RowContents **ox1d;		o: output array with table data.
int rebinfactor;		i: rebinning factor
int nrows;			i: number of rows
*/
	int i, j, i2, j2, rf;
	double hold;

	/* Alloc memory in output x1d array. */

	rf = (rebinfactor > 0) ? rebinfactor : -rebinfactor;
	for (j = 0; j < nrows; j++) {
	    ox1d[j]->sporder = ix1d[j]->sporder;
	    ox1d[j]->npts    = ix1d[j]->npts / rf;
	    ox1d[j]->wave     = (double *) calloc (ox1d[j]->npts,
                                                  sizeof(double));
	    ox1d[j]->gross    = (float *) calloc (ox1d[j]->npts,
                                                  sizeof(float));
	    ox1d[j]->net      = (float *) calloc (ox1d[j]->npts,
                                                  sizeof(float));
	    ox1d[j]->extrlocy = (float *) calloc (ox1d[j]->npts,
                                                 sizeof(float));
	    if (ox1d[j]->wave == NULL || ox1d[j]->gross    == NULL ||
                ox1d[j]->net  == NULL || ox1d[j]->extrlocy == NULL) {
                printf ("Not enough memory to allocate data arrays.\n");
	        return (OUT_OF_MEMORY);
	    }
	}

	switch (rebinfactor) {

	case 1:

	    /* Just copy. */

	    for (j = 0; j < iim->sci.data.ny; j++) {
	        for (i = 0; i < iim->sci.data.nx; i++) {
	            Pix (oim->sci.data, j, i)   = Pix (iim->sci.data, j, i);
	            Pix (oim->err.data, j, i)   = Pix (iim->err.data, j, i);
	            DQPix (oim->dq.data, j, i)  = DQPix (iim->dq.data, j, i);
	        }
	    }
	    for (j = 0; j < nrows; j++) {
	        for (i = 0; i < iim->sci.data.nx; i++) {
	            ox1d[j]->wave[i] = ix1d[j]->wave[i];
	            ox1d[j]->gross[i] = ix1d[j]->gross[i];
	            ox1d[j]->net[i] = ix1d[j]->net[i];
	            ox1d[j]->extrlocy[i] = ix1d[j]->extrlocy[i];
	        }
	    }
	    return (STIS_OK);

	case 2:

	    /* Squeeze by a factor 2. Only the SCI data is processed since
               the algorithm doesn't propagate ERR and DQ info anyway.
            */
	    for (j = 0; j < oim->sci.data.ny-1; j++) {
	        j2 = 2 * j;
	        for (i = 0; i < oim->sci.data.nx-1; i++) {
	            i2 = 2 * i;
	            Pix (oim->sci.data, j, i) = (Pix (iim->sci.data, j2, i2) +
                        Pix (iim->sci.data, j2, i2+1) +
                        Pix (iim->sci.data, j2+1, i2) +
                        Pix (iim->sci.data, j2+1, i2+1)) / 4.0;
	        }
	    }

	    for (j = 0; j < nrows; j++) {
	        for (i = 0; i < oim->sci.data.nx-1; i++) {
	            i2 = 2 * i;
	            ox1d[j]->wave[i]     = (ix1d[j]->wave[i2] +
                                            ix1d[j]->wave[i2+1]) / 2.0 ;
	            ox1d[j]->gross[i]    = (ix1d[j]->gross[i2] +
                                            ix1d[j]->gross[i2+1]) / 2.0 ;
	            ox1d[j]->net[i]      = (ix1d[j]->net[i2] +
                                            ix1d[j]->net[i2+1]) / 2.0 ;
	            ox1d[j]->extrlocy[i] = (ix1d[j]->extrlocy[i2] +
                                            ix1d[j]->extrlocy[i2+1]) / 2.0 ;
	            ox1d[j]->extrlocy[i] /= (double)rebinfactor;
	        }
	        ox1d[j]->wave[oim->sci.data.nx-1] =
                    ox1d[j]->wave[oim->sci.data.nx-2];
	        ox1d[j]->gross[oim->sci.data.nx-1] =
                    ox1d[j]->gross[oim->sci.data.nx-2];
	        ox1d[j]->net[oim->sci.data.nx-1] =
                    ox1d[j]->net[oim->sci.data.nx-2];
	        ox1d[j]->extrlocy[oim->sci.data.nx-1] =
                    ox1d[j]->extrlocy[oim->sci.data.nx-2];
	    }

	    return (STIS_OK);

	case -2:

	    /* Expand by a factor 2. Only the SCI data is processed since
               the algorithm doesn't propagate ERR and DQ info anyway. It
               is assumed that the output image already stores this info.
            */
	    for (j = 0; j < iim->sci.data.ny; j++) {
	        j2 = 2 * j;
	        for (i = 0; i < iim->sci.data.nx; i++) {
	            i2 = 2 * i;
	            hold = Pix (iim->sci.data,j,i) / 2.0;
	            Pix (oim->sci.data, j2,i2)     = hold;
	            Pix (oim->sci.data, j2+1,i2)   = hold;
	            Pix (oim->sci.data, j2,i2+1)   = hold;
	            Pix (oim->sci.data, j2+1,i2+1) = hold;
	        }
	    }
	    /* No need to expand x1d array. */

	    return (STIS_OK);

	default:
	    return (ERROR_RETURN);
	}
}
