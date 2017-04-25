# include <stdio.h>
# include <stdlib.h>

# include "xtables.h"

# include "stis.h"
# include "calstis6.h"
# include "hstcalerr.h"
# include "stispht.h"


/*
   Routines to manage memory.






   Revision history:
   ----------------
   20 Feb 97  -  Implemented, with some code adapted from similar routines
                 in calstis7 (I.Busko)
   13 Apr 98  -  Add extrlocy array (IB)
   24 Jun 98  -  Profile memory management (IB)
   29 Jun 98  -  Profile rejection ranges (IB)
   10 Jul 98  -  Fixed inconsistent extrsize definition (IB)
   18 Sep 98  -  Optimal extraction support (IB)
   06 Dec 99  -  Profile DQ array (IB)
   17 Dec 99  -  Flux normalization factors from profile (IB)
   25 Apr 00  -  Extra room for profile array (IB)
   16 Jun 00  -  2-D array with profile Y positions (IB)
   13 Jul 00  -  Redefined extrsz as double (IB)
   01 Nov 00  -  Profile offset (IB)
   28 Nov 00  -  Subsampling in profile builder (IB)
   30 Jan 01  -  Centroid output (IB)
   08 Apr 05  -  In FreeThroughput6, free memory for GAC info.
   21 Oct 11  -  Allocate/deallocate memory for new NET_ERROR column (PEH)
*/




/* Routines to handle memory related to the output arrays associated with
   the row structure. The 'npts' structure member must be set by the caller
   to the number of elemnents in each array prior to calling the allocation
   function.

   calloc is used instead of malloc so arrays that are left empty by ommited
   steps won't crash the output table routines with garbage.
*/

int AllocOutArrays (RowContents *row) {

	if ((row->wave = (double *)calloc(row->npts, sizeof (double))) == NULL)
	    return (OUT_OF_MEMORY);
	if ((row->gross = (float *)calloc(row->npts, sizeof (float))) == NULL)
	    return (OUT_OF_MEMORY);
	if ((row->back = (float *)calloc(row->npts, sizeof (float))) == NULL)
	    return (OUT_OF_MEMORY);
	if ((row->net = (float *)calloc(row->npts, sizeof (float))) == NULL)
	    return (OUT_OF_MEMORY);
	if ((row->flux = (float *)calloc(row->npts, sizeof (float))) == NULL)
	    return (OUT_OF_MEMORY);
	if ((row->error = (float *)calloc(row->npts, sizeof (float))) == NULL)
	    return (OUT_OF_MEMORY);
	if ((row->net_error = (float *)calloc(row->npts, sizeof (float))) ==
		NULL)
	    return (OUT_OF_MEMORY);
	if ((row->dq = (short *)calloc(row->npts, sizeof (short))) == NULL)
	    return (OUT_OF_MEMORY);
	if ((row->extrlocy = (float *)calloc(row->npts, sizeof (float))) ==
            NULL)
	    return (OUT_OF_MEMORY);
	return (STIS_OK);

}



/*  Alloc memory for profile array. */

int AllocProfile (StisInfo6 *sts, int imxsize, double extrsz) {

	int extrsize, nalloc, i;

	extrsize = (int)extrsz;

	/* Translate to image pixels. */
	extrsize *= sts->ltm[1];

	/* If building a profile file, increase extraction box size
           by one pixel. Because the box center is non-integer, the
           actual number of pixels covered by the box when laid out
           on top of the detector pixel grid is one more than its
           specified size. Once retrieved, the pixels will be sampled
           back into the proper extraction box size.

           This is not necessary anymore since we subsample.
        */
	/*
	if (sts->do_profile)
	    extrsize++;
	*/

	/* Profile-related arrays are allocated independent of the
           do_profile and optimal switches. This avoids excessive
           testing in AccumPix. The amount of allocated memory does
           depend on switch settings tough. This is necessary to avoid
           rui errors when passing the addresses as function parameters
           in X1DSpec.
        */
	if ((sts->profile = (double **) malloc (imxsize * sizeof (double *)))
             == NULL)
	    return (OUT_OF_MEMORY);

	if ((sts->profile_dq = (short **) malloc (imxsize * sizeof (short *)))
             == NULL)
	    return (OUT_OF_MEMORY);

	if ((sts->subprof = (double **) malloc (imxsize *  sizeof(double)))
             == NULL)
	    return (OUT_OF_MEMORY);

	if ((sts->profile_rejf = (double *) calloc (imxsize, sizeof(double)))
             == NULL)
	    return (OUT_OF_MEMORY);

	if ((sts->profile_offset = (double *) calloc (imxsize, sizeof(double)))
             == NULL)
	    return (OUT_OF_MEMORY);

	if ((sts->profile_centroid = (double *) calloc (imxsize,
                                     sizeof(double))) == NULL)
	    return (OUT_OF_MEMORY);

	if ((sts->profile_rej = (short *) calloc (imxsize, sizeof(short)))
             == NULL)
	    return (OUT_OF_MEMORY);

	/* The PROF_MARGIN additional pixels added to compute the subpixel
           array size are meant to provide a safety margin to let the
           trace wander around a little. An alternate way of fixing this
           would be to find out beforehand the maximum amount of
           trace wandering based on the actual trace file. This would
           require an extra loop over the input array columns, with all
           the trace handling machinery duplicated from somewhere else
           (actually triplicated by this point...). And it probably would
           still require the addition of a couple pixels anyway.

           Also, notice that the definition in here MUST be the same
           used in Do1Dx to define the table array size. This code
           duplication cannot be avoided due to the convoluted way the
           optimal extraction code was developed.
        */
	if (sts->do_profile || sts->optimal) {
	    nalloc  = extrsize + 1;
	    sts->subprof_size = (int)sts->subscale *
                                (extrsize * sts->ltm[1] + PROF_MARGIN);
	} else {
	    nalloc  = 1;
	    sts->subprof_size = 1;
	}

	for (i = 0; i < imxsize; i++) {
	    if ((sts->profile[i] = (double *) calloc (nalloc,
                                   sizeof (double))) == NULL)
	        return (OUT_OF_MEMORY);
	}
	for (i = 0; i < imxsize; i++) {
	    if ((sts->profile_dq[i] = (short *) calloc (nalloc,
                                      sizeof (short))) == NULL)
	        return (OUT_OF_MEMORY);
	}
	for (i = 0; i < imxsize; i++) {
	    if ((sts->subprof[i] = (double *) calloc (sts->subprof_size,
                                    sizeof (double))) == NULL)
	        return (OUT_OF_MEMORY);
	}

	sts->profile_x = imxsize;
	sts->profile_y = extrsize;

	return (0);
}



/* Free profile array memory. */

void FreeProfile (StisInfo6 *sts, int imxsize) {

	int i;

	for (i = 0; i < imxsize; free (sts->profile[i++]));
	for (i = 0; i < imxsize; free (sts->profile_dq[i++]));
	for (i = 0; i < imxsize; free (sts->subprof[i++]));
	free (sts->profile);
	free (sts->profile_dq);
	free (sts->subprof);

	free (sts->profile_offset);
	free (sts->profile_centroid);
	free (sts->profile_rejf);
	free (sts->profile_rej);
}



/* Free memory associated with the output arrays. */

void FreeOutArrays (RowContents *row) {

	if (row->wave      != NULL) free (row->wave);
	if (row->gross     != NULL) free (row->gross);
	if (row->back      != NULL) free (row->back);
	if (row->net       != NULL) free (row->net);
	if (row->flux      != NULL) free (row->flux);
	if (row->error     != NULL) free (row->error);
	if (row->net_error != NULL) free (row->net_error);
	if (row->dq        != NULL) free (row->dq);
	if (row->extrlocy  != NULL) free (row->extrlocy);
}



/* Free memory related to reference data, if it has been allocated,
   and resets the flag to indicate that memory is no longer allocated.
*/

void FreePhot6 (PhotInfo *phot) {

	if (phot->allocated) {
	    free (phot->wl);
	    free (phot->thru);
	    free (phot->error);
	    if (phot->pcorr != NULL) {
		free (phot->pcorr);
		phot->pcorr = NULL;
	    }
	    phot->allocated = 0;
	}
}

void FreeThroughput6 (ApInfo *slit) {

	if (slit->allocated) {
	    free (slit->wl);
	    free (slit->thr);
	    slit->allocated = 0;
	}
	if (slit->gac_allocated) {
	    free (slit->gac_wl);
	    free (slit->gac_thr);
	    slit->gac_nelem = 0;
	    slit->gac_allocated = 0;
	}
}

void FreeInang6 (InangInfo *iac) {

	if (iac->allocated) {
	    if (iac->ncoeff1 > 0)
		free (iac->coeff1);
	    if (iac->ncoeff2 > 0)
		free (iac->coeff2);
	    iac->allocated = 0;
	}
}

void FreeIntensity (IntensArray *inta) {

	if (inta->allocated) {
	    free (inta->wave);
	    free (inta->intens);
	    inta->allocated = 0;
	}
}


