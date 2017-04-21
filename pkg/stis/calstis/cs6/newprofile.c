
/* This file contains:
	NewProfile	  adds a ProfileArray record to the list
	CheckProfile      checks for invalid entries in list
	InterpProfile	  returns an interpolated ProfileArray
	FreeProfileArray  frees memory
        DebugProfile
internal:
	CopyProfile	copy out one ProfileArray




   Revision history:
   ----------------
   17 Sep 98  -  Implemented (I.Busko)
   06 Dec 00  -  Subsampled profile (IB)
*/

# include <stdio.h>
# include <stdlib.h>
# include <string.h>

# include "xtables.h"

# include "stis.h"
# include "calstis6.h"
# include "err.h"

static void CopyProfile (ProfileArray *, ProfileArray *);

/* This routine adds a ProfileArray record to the list.  The list is
   sorted by minp.  Before calling NewProfile for the first record
   to be inserted in the list, *profa should have been initialized to NULL.
*/

int NewProfile (ProfileArray **profa, ProfileArray *input) {

/* arguments:
ProfileArray **profa  io: identifies the first record in the list
ProfileArray *input   i: a new record to be inserted into the list
*/
	ProfileArray *current, *next, *new;
	int done = 0;

	/* Allocate space for the new record, and copy input to new. */
	if ((new = (ProfileArray *) malloc (sizeof (ProfileArray))) == NULL) {
	    printf ("ERROR    Can't allocate memory.\n");
	    return (OUT_OF_MEMORY);
	}
	new->profoff = (double *) malloc (input->nptsoff * sizeof (double));
	new->prof    = (double *) malloc (input->npts    * sizeof (double));
	if (new->profoff == NULL || new->prof == NULL) {
	    printf ("ERROR    Can't allocate memory.\n");
	    return (OUT_OF_MEMORY);
	}
	CopyProfile (input, new);

	/* Compute and store mid-point. */
	new->xpix = (new->maxp - new->minp) / 2 + new->minp;

	/* Find where this record belongs in the list, sorting by minp,
	   and insert the record.
	*/
	current = *profa;

	if ((*profa) == NULL) {
	    /* The list is empty; new becomes the first record. */
	    *profa = new;
	    done = 1;
	} else if (new->minp <= (*profa)->minp) {
	    /* new is before the beginning of the list; new becomes first. */
	    new->next = *profa;
	    *profa = new;
	    done = 1;
	}

	while (!done) {

	    if (current->next == NULL) {
		current->next = new;		/* add to end of list */
		done = 1;

	    } else {

		next = current->next;

		if (new->minp <= next->minp) {
		    /* insert here */
		    new->next = current->next;
		    current->next = new;
		    done = 1;
		}
	    }
	    if (!done) {
		current = next;
	    }
	}

	return (0);
}




/* This routine checks that there are no duplicate values of minp
   in the ProfileArray list. This should be called after reading all
   rows of the profile table and adding the relevant ones to the list.
   If the list is empty, this routine returns -1.
*/

int CheckProfile (ProfileArray **profa) {

	ProfileArray *current, *previous;

	current = *profa;
	if (current == NULL)
	    return (-1);		/* empty list */

	previous = current;
	current = current->next;

	while (current != NULL) {
	    if (current->minp == previous->minp) {
		printf (
	"ERROR    Duplicate values of MIN_PIX = %d in OPROFTAB.\n",
			current->minp);
		return (INVALID);
	    }
	    previous = current;
	    current = current->next;
	}

	return (0);
}



/* This routine frees all memory for the ProfileArray list, including
   profa itself.
*/

void FreeProfileArray (ProfileArray **profa) {

	ProfileArray *current, *next;

	current = *profa;
	while (current != NULL) {
	    next = current->next;
	    free (current->profoff);
	    free (current->prof);
	    free (current);
	    current = next;
	}

	*profa = NULL;
}




/* This routine interpolates the profile elements on the X coordinate.
   If X is outside the range of pixel values in the list, the first (or
   last) record will be copied to output without interpolation; otherwise,
   the two records that have their X coordinate bracketing the input
   X will be used, and the profile values will be linearly interpolated.

   The X coordinate of a given record is definde by

   X = (maxp - minp) / 2 + minp

   where maxp and minp are the values actually stored in the record.
   Notice that the X value is not stored in the output record. It is
   the caller's responsibility to keep track of X pixel coordinates.

   If *output is NULL, memory will be allocated for that record; otherwise,
   it is assumed that *output is actually a valid pointer to a ProfileArray,
   i.e. that the memory has been allocated.

   The offsets array is not copied into the output profile structure.
   Instead, the function sets the 'offset' variable to the actual offset
   corresponding to the X column in the image.
*/

int InterpProfile (ProfileArray **profa, int xpix, ProfileArray **output,
                   double *offset) {

/* arguments:
ProfileArray **profa  i: identifies the first record in the list
int xpix              i: the X pixel coordinate where to compute output
ProfileArray **output o: the profile interpolated to xpix
double *offset;       o: offset corresponding to position X
*/

	ProfileArray *current, *next;
	float p;
	int i, done;

	current = *profa;

	/* Allocate space for the output record. */
	if (*output == NULL) {
	    if ((*output = (ProfileArray *) malloc (sizeof (ProfileArray))) == 
                            NULL) {
		printf ("ERROR    Can't allocate memory.\n");
		return (OUT_OF_MEMORY);
	    }
	    (*output)->profoff = (double *) malloc (current->nptsoff * 
                                  sizeof (double));
	    (*output)->prof    = (double *) malloc (current->npts * 
                                 sizeof (double));
	    if ((*output)->profoff == NULL || (*output)->prof == NULL) {
	        printf ("ERROR    Can't allocate memory.\n");
	        return (OUT_OF_MEMORY);
	    }
	}

	done = (current == NULL);

	while (!done) {

	    next = current->next;

            /* offset corresponding to X position. Note that in here
               "current" doesn't mean exactly that. This code was
               adapted from other similar routines in calstis that
               handle linked lists, but here we have the twist that
               effective ranges for interpolation purposes begin and
               end at the *mid point* of the physical ranges, thus the
               additional and confusing set of tests. Note also that
               xpix is zero-indexed and range endpoints are 1-indexed.
            */
	    if (xpix < current->maxp)
	        *offset = current->profoff[xpix + 1 - current->minp];
	    else
	        *offset = next->profoff[xpix + 1 - next->minp];

	    if (xpix <= current->xpix) {

		/* before beginning of list, or exact match */

		CopyProfile (current, *output);	/* no interpolation needed */
		done = 1;

	    } else if (next == NULL) {

		/* end of list */

		CopyProfile (current, *output);
		done = 1;

	    } else if (current->xpix < xpix && xpix < next->xpix) {

		/* Linear interpolation of profile values. */

		p = (float)(xpix - current->xpix) / 
                    (next->xpix - current->xpix);

		(*output)->npts    = current->npts;
		(*output)->nptsoff = current->nptsoff;

		for (i = 0;  i < (*output)->npts;  i++)
		    (*output)->prof[i] =
			(1. - p) * current->prof[i] + p * next->prof[i];

		(*output)->next = NULL;

		done = 1;

	    } else {

		current = next;
		next = current->next;
	    }

	    if (current == NULL)
		done = 1;
	}

	return (0);
}




/* This routine copies the current ProfileArray to output, which must
   actually point to allocated space.  Note that in output, the pointer
   to the next record is set to NULL.
*/

static void CopyProfile (ProfileArray *current, ProfileArray *output) {

/* arguments:
ProfileArray *current   i: the current ProfileArray record
ProfileArray *output    o: receives a copy of current
*/

	int i;
	output->npts    = current->npts;
	output->nptsoff = current->nptsoff;
	output->minw    = current->minw;
	output->maxw    = current->maxw;
	output->minp    = current->minp;
	output->maxp    = current->maxp;
	output->sn      = current->sn;
	for (i = 0;  i < current->npts;  i++)
	    output->prof[i]  = current->prof[i];
	for (i = 0;  i < current->nptsoff;  i++)
	    output->profoff[i] = current->profoff[i];
	output->next = NULL;
}




/* This routine prints the minp coordinate from the ProfileArray list. */

void DebugProfile (ProfileArray **profa) {

/* argument:
ProfileArray **profa  i: identifies the first record in the list
*/

	ProfileArray *current;
	int i = 1;

	current = *profa;

	while (current != NULL) {
	    printf ("%3d:  %d\n", i, current->minp);
	    current = current->next;
	    i++;
	}
}



