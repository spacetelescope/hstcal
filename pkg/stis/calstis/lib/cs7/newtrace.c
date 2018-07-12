/* This file contains:
	NewTrace	add an SpTrace record to the list
	CheckTrace	check for invalid entries in list
	InterpTrace	return an SpTrace, interpolated on a2center
	FreeTrace	free memory
internal:
	CopyTrace	copy out one SpTrace
*/

# include <stdio.h>
# include <stdlib.h>
# include <string.h>

# include "stis.h"
# include "calstis7.h"
# include "hstcalerr.h"

# define MAX(x,y)  ((x) >= (y) ? (x) : (y))

static void CopyTrace (SpTrace *, SpTrace *);

/* This routine adds an SpTrace record to the list.  The list is sorted
   by a2center.  Before calling NewTrace for the first record to be
   inserted in the list, *trace should have been initialized to NULL.

   Phil Hodge, 1998 Oct 5:
	Change status values 1212 and 1214 to GENERIC_ERROR_CODE.

   Phil Hodge, 2000 Jan 13:
	NewTrace now adds newrec directly to the list, rather than
	making a copy.

   Phil Hodge, 2001 Mar 7:
	Remove DebugTrace.

   Phil Hodge, 2004 Dec 27:
	In CheckTrace, modify
	"if (current->a1center != previous->a1center != 0) {"
	to remove the "!= 0".
*/

int NewTrace (SpTrace **trace, SpTrace *newrec) {

/* arguments:
SpTrace **trace  io: identifies the first record in the list
SpTrace *newrec   i: a new record to be inserted into the list
*/

	SpTrace *current, *next;
	int done = 0;

	if (newrec->nelem > MAX_SP_TRACE) {
	    printf ("ERROR    (NewTrace) %d elements in array.\n",
		newrec->nelem);
	    return (TABLE_ERROR);
	}

	/* Find where this record belongs in the list, sorting by a2center,
	   and insert the record.
	*/
	current = *trace;

	if (*trace == NULL) {
	    /* The list is empty; newrec becomes the first record. */
	    *trace = newrec;
	    done = 1;
	} else if (newrec->a2center <= (*trace)->a2center) {
	    /* newrec is before the beginning of the list;
		newrec becomes the first.
	    */
	    newrec->next = *trace;
	    *trace = newrec;
	    done = 1;
	}

	while (!done) {

	    if (current->next == NULL) {
		current->next = newrec;		/* add to end of list */
		done = 1;

	    } else {

		next = current->next;

		if (newrec->a2center <= next->a2center) {
		    /* insert here */
		    newrec->next = current->next;
		    current->next = newrec;
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
   the list.  If the list is empty, this routine returns -1.
*/

int CheckTrace (SpTrace **trace) {

/* arguments:
SpTrace **trace  i: identifies the first record in the list
*/

	SpTrace *current, *previous;

	current = *trace;
	if (current == NULL)
	    return (-1);		/* empty list */

	previous = current;
	current = current->next;

	while (current != NULL) {
	    if (current->a2center == previous->a2center) {
		printf (
	"ERROR    Duplicate values of A2CENTER=%.8g in SPTRCTAB.\n",
			current->a2center);
		return (GENERIC_ERROR_CODE);
	    }
	    if (current->a1center != previous->a1center) {
		printf (
	"ERROR    Different A1CENTER values in different rows in SPTRCTAB.\n");
		return (GENERIC_ERROR_CODE);
	    }
	    previous = current;
	    current = current->next;
	}

	return (0);
}

/* This routine frees all memory for the SpTrace list.  *trace itself
   will be freed as well, and *trace will be set to NULL.
*/

void FreeTrace (SpTrace **trace) {

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

/* This routine interpolates the spectrum traces on a2center.
   If a2center is outside the range of values in the list, the first (or
   last) record will be copied to output without interpolation; otherwise,
   the two records that have a2center bracketing the input a2center will
   be used, and the spectrum traces will be linearly interpolated.

   If *output is NULL, memory will be allocated for that record; otherwise,
   it is assumed that *output is actually a valid pointer to an SpTrace,
   i.e. that the memory has been allocated.
*/

int InterpTrace (SpTrace **trace, double a2center,
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
