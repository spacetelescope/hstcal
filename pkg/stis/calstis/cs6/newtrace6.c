
/* This file contains:
	NewTrace6	add an SpTrace record to the list
	CheckTrace6	check for invalid entries in list
	InterpTrace	return an SpTrace, interpolated on a2center
	FreeTrace6	free memory
internal:
	CopyTrace	copy out one SpTrace
	PrintTrace	print info from full SpTrace linked list (debugging)




   Revision history:
   ----------------
   20 Feb 97  -  Borrowed from calstis7 (I.Busko)
   20 Feb 97  -  CheckTrace modified to allow skipping of orders (IB)
   20 Feb 97  -  Removed DebugTrace (IB)
   24 Feb 97  -  Rename routines to avoid conflicts with cs7 (IB)
   11 Apr 97  -  Changes after code review (IB):
                 - explicit cast for malloc-returned pointer.
   08 May 97  -  Conform to new _trl standard (IB)
   02 Jul 97  -  Fixed logic that zeros coefficients (IB)
   28 Jul 97  -  Added InterpTrace routine (IB)
   25 Sep 97  -  Added PrintTrace function (IB)
   30 Jul 98  -  Add 1 to duplicate A2CENTER error message (IB)
   02 Sep 98  -  Removed fflush(stdout) calls (IB)
   27 Dec 04  -  Remove the "!= 0" from
		"if (current->a1center != previous->a1center != 0) {" (PEH)
*/

# include <stdio.h>
# include <stdlib.h>
# include <string.h>

# include "xtables.h"
# include "stis.h"
# include "calstis6.h"
# include "err.h"

# define MAX(x,y)  ((x) >= (y) ? (x) : (y))

static void CopyTrace (SpTrace *, SpTrace *);

/* This routine adds an SpTrace record to the list.  The list is sorted
   by a2center.  Before calling NewTrace6 for the first record to be
   inserted in the list, *trace should have been initialized to NULL.
*/

int NewTrace6 (SpTrace **trace, SpTrace *input) {

/* arguments:
SpTrace **trace  io: identifies the first record in the list
SpTrace *input   i: a new record to be inserted into the list
*/

	SpTrace *current, *next, *new;
	int done = 0;

	if (input->nelem > MAX_SP_TRACE) {
	    printf ("ERROR    %d trace coefficients.\n",
		input->nelem);
	    return (TABLE_ERROR);
	}

	/* Allocate space for the new record, and copy input to new. */
	if ((new = (SpTrace *) malloc (sizeof (SpTrace))) == NULL) {
	    printf ("ERROR    Can't allocate memory.\n");
	    return (OUT_OF_MEMORY);
	}
	CopyTrace (input, new);

	/* Find where this record belongs in the list, sorting by a2center,
	   and insert the record.
	*/
	current = *trace;

	if (*trace == NULL) {
	    /* The list is empty; new becomes the first record. */
	    *trace = new;
	    done = 1;
	} else if (new->a2center <= (*trace)->a2center) {
	    /* new is before the beginning of the list; new becomes first. */
	    new->next = *trace;
	    *trace = new;
	    done = 1;
	}

	while (!done) {

	    if (current->next == NULL) {
		current->next = new;		/* add to end of list */
		done = 1;

	    } else {

		next = current->next;

		if (new->a2center <= next->a2center) {
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


/* This routine checks that there are no duplicate valus of a2center
   in the SpTrace list, and it checks that all values of a1center
   are the same.  This should be called after reading all rows of the
   table of spectrum traces and adding the relevant ones to
   the list.  If the list is empty, this routine returns NO_TRACE.
*/

int CheckTrace6 (SpTrace **trace) {

/* arguments:
SpTrace **trace  i: identifies the first record in the list
*/

	SpTrace *current, *previous;

	current = *trace;
	if (current == NULL)
	    return (NO_TRACE);		/* empty list */

	previous = current;
	current = current->next;
	while (current != NULL) {
	    if (current->a2center == previous->a2center) {
	        /* a2center is stored internally in zero-indexed format */
		printf (
	"ERROR    Duplicate values of A2CENTER=%.8g in SPTRCTAB.\n",
			current->a2center + 1.0);
		return (ERROR_TRACE);
	    }
	    if (current->a1center != previous->a1center) {
		printf (
	"ERROR    Different A1CENTER values in different rows in SPTRCTAB.\n");
		return (ERROR_TRACE);
	    }
	    previous = current;
	    current = current->next;
	}

	return (0);
}


/* This routine interpolates the spectrum traces on a2center.
   If a2center is outside the range of values in the list, the first (or
   last) record will be copied to output without interpolation; otherwise,
   the two records that have a2center bracketing the input a2center will
   be used, and the spectrum traces will be linearly interpolated.

   If *output is NULL, memory will be allocated for that record; otherwise,
   it is assumed that *output is actually a valid pointer to an SpTrace,
   i.e. that the memory has been allocated.
*/

int InterpTrace6 (SpTrace **trace, double a2center,
		 SpTrace **output) {

/* arguments:
SpTrace **trace   i: identifies the first record in the list
double a2center   i: line number in input image
SpTrace **output  o: the spectrum trace interpolated to a2center
*/

	SpTrace *current, *next;
	double p;
	int i;
	int done;

	current = *trace;

	/* Allocate space for the output. */
	if (*output == NULL) {
	    if ((*output = malloc (sizeof (SpTrace))) == NULL) {
		printf ("ERROR    Can't allocate memory in InterpTrace.\n");
		return (OUT_OF_MEMORY);
	    }
	}

	done = (current == NULL);

	while (!done) {

	    next = current->next;

	    if (a2center <= current->a2center) {

		/* before beginning of list, or exact match */
		CopyTrace (current, *output);	/* no interpolation needed */
		(*output)->a2center = a2center;
		done = 1;

	    } else if (next == NULL) {

		/* end of list */
		CopyTrace (current, *output);
		(*output)->a2center = a2center;
		done = 1;

	    } else if (current->a2center < a2center &&
				a2center < next->a2center) {

		/* linear interpolation of spectrum traces */
		p = (a2center - current->a2center) /
			(next->a2center - current->a2center);

		(*output)->a2center = a2center;
		(*output)->a1center = current->a1center;
		(*output)->nelem = MAX (current->nelem, next->nelem);
		for (i = 0;  i < (*output)->nelem;  i++) {
		    (*output)->a2displ[i] =
			(1. - p) * current->a2displ[i] + p * next->a2displ[i];
		}
		for (i = (*output)->nelem;  i < MAX_SP_TRACE;  i++)
		    (*output)->a2displ[i] = 0.;
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


/* This routine frees all memory for the SpTrace list.  *trace itself
   will be freed as well, and *trace will be set to NULL.
*/

void FreeTrace6 (SpTrace **trace) {

/* argument:
SpTrace **trace  io: identifies the first record in the list; set to NULL
*/

	SpTrace *current, *next;

	current = *trace;

	while (current != NULL) {
	    next = current->next;
	    free (current);
	    current = next;
	}

	*trace = NULL;
}


/* This routine copies the current SpTrace to output, which must
   actually point to allocated space.  Note that in output, the pointer
   to the next record is set to NULL.
*/

static void CopyTrace (SpTrace *current, SpTrace *output) {

/* arguments:
SpTrace *current   i: the current SpTrace record
SpTrace *output    o: receives a copy of current
*/

	int i;

	output->a2center = current->a2center;
	output->a1center = current->a1center;
	output->nelem = current->nelem;
	for (i = 0;  i < current->nelem;  i++)
	    output->a2displ[i] = current->a2displ[i];
	for (i = current->nelem;  i < MAX_SP_TRACE;  i++)
	    output->a2displ[i] = 0.;
	output->next = NULL;
}




/* This routine transverses a linked list of SpTrace structures,
   printing at stdout the a2center element of each entry. It is
   mean for debugging purposes only.
*/

void PrintTrace6 (SpTrace **trace) {

/* arguments:
SpTrace **trace   i: identifies the first record in the list
*/
	SpTrace *current;
	int done;

	current = *trace;
	done = (current == NULL);
	printf ("Beginning of trace list.\n");

	while (!done) {

	    printf ("%g  %g\n", current->a2center, current->a2displ[0]);
	    current = current->next;
	    if (current == NULL)
		done = 1;
	}

	printf ("End of trace list.\n");
}

