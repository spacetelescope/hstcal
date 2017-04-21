/* This file contains:
	NewCoord	add a CoordInfo record to the list
	ReturnCoord	extract a CoordInfo record from the list
	RangeCoord	find the range of sporders in the list
	FreeCoord	free memory
internal:
	CopyCoord	copy out one CoordInfo
*/

# include <stdio.h>
# include <stdlib.h>

# include "stis.h"
# include "calstis7.h"
# include "err.h"

static void CopyCoord (CoordInfo *, CoordInfo *);

/* This routine adds a CoordInfo record to the list.  The list is sorted
   by sporder.  Before calling NewCoord for the first record to be inserted
   in the list, *coords should have been initialized to NULL.

   Phil Hodge, 1998 Oct 5:
	Change status values 1212, 1213, 1214 to GENERIC_ERROR_CODE.

   Phil Hodge, 2001 Mar 7:
	Remove DebugCoord.
*/

int NewCoord (CoordInfo **coords, CoordInfo *input) {

/* arguments:
CoordInfo **coords  io: identifies the first record in the list
CoordInfo *input    i: a new record to be inserted into the list
*/

	CoordInfo *current, *next, *newrec;
	int done = 0;

	/* Allocate space for the new record, and copy to memory. */
	if ((newrec = malloc (sizeof (CoordInfo))) == NULL) {
	    printf ("ERROR    Can't allocate memory in NewCoord.\n");
	    return (OUT_OF_MEMORY);
	}

	CopyCoord (input, newrec);

	/* Find where this record belongs in the list, sorting by sporder,
	   and insert the record.
	*/
	current = *coords;

	if ((*coords) == NULL) {
	    /* The list is empty; newrec becomes the first record. */
	    *coords = newrec;
	    done = 1;
	} else if (newrec->sporder <= (*coords)->sporder) {
	    /* newrec is before the beginning of the list;
		newrec becomes the first.
	    */
	    newrec->next = *coords;
	    *coords = newrec;
	    done = 1;
	}

	while (!done) {

	    if (current->next == NULL) {
		current->next = newrec;		/* add to end of list */
		done = 1;

	    } else {

		next = current->next;

		if (newrec->sporder <= next->sporder) {
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

/* This routine returns a CoordInfo record from the list.  The record is
   selected by sporder.

   If *output is NULL, memory will be allocated for that record; otherwise,
   it is assumed that *output is actually a valid pointer to a CoordInfo,
   i.e. that the memory has been allocated.
*/

int ReturnCoord (CoordInfo **coords, int sporder, CoordInfo **output) {

/* arguments:
CoordInfo **coords  i: identifies the first record in the list
CoordInfo **output   o: the record that matches sporder
*/

	CoordInfo *current;
	int foundit = 0;

	current = *coords;

	/* Allocate space for the output. */
	if (*output == NULL) {
	    if ((*output = malloc (sizeof (CoordInfo))) == NULL) {
		printf ("ERROR    Can't allocate memory in ReturnCoord.\n");
		return (OUT_OF_MEMORY);
	    }
	}

	while (current != NULL) {
	    if (current->sporder == sporder) {
		foundit = 1;
		break;
	    }
	    current = current->next;
	}

	if (foundit)
	    CopyCoord (current, *output);
	else
	    /* sporder not found, or empty list */
	    return (GENERIC_ERROR_CODE);

	return (0);
}


/* This routine finds the minimum and maximum values of sporder in the
   CoordInfo list.  The list is also checked to ensure that the sporder
   numbers are consecutive and there are no duplicates.
*/

int RangeCoord (CoordInfo **coords, int *minorder, int *maxorder) {

/* arguments:
CoordInfo **coords        i: identifies the first record in the list
int *minorder, *maxorder  o: minimum and maximum values of sporder
*/

	CoordInfo *current, *previous;

	current = *coords;
	if (current == NULL)
	    return (-1);		/* empty list */

	*minorder = current->sporder;		/* initial values */
	*maxorder = current->sporder;

	previous = current;
	current = current->next;

	while (current != NULL) {

	    if (current->sporder != previous->sporder + 1) {
		if (current->sporder == previous->sporder) {
		    printf (
		"ERROR    Duplicate order number %d in SDC table.\n",
			current->sporder);
		    return (GENERIC_ERROR_CODE);
		} else {
		    printf (
		"ERROR    Order numbers in SDC table jump from %d to %d.\n",
		    previous->sporder, current->sporder);
		    return (GENERIC_ERROR_CODE);
		}
	    }
	    *maxorder = current->sporder;
	    previous = current;
	    current = current->next;
	}

	return (0);
}

/* This routine frees all memory for the CoordInfo list.  The value of
   coords is not changed.
*/

void FreeCoord (CoordInfo **coords) {

/* argument:
CoordInfo **coords  io: identifies the first record in the list
*/

	CoordInfo *current, *next;

	current = *coords;

	while (current != NULL) {
	    next = current->next;
	    free (current);
	    current = next;
	}
}

/* This routine copies the current CoordInfo to output, which must
   actually point to allocated space.  Note that in output, the pointer
   to the next record is set to NULL.
*/

static void CopyCoord (CoordInfo *current, CoordInfo *output) {

/* argument:
CoordInfo *current  i: the current CoordInfo record
CoordInfo *output   o: receives a copy of current
*/

	output->sporder  = current->sporder;
	output->a2center = current->a2center;
	output->npix[0]  = current->npix[0];
	output->npix[1]  = current->npix[1];
	output->crpix[0] = current->crpix[0];
	output->crpix[1] = current->crpix[1];
	output->crval[0] = current->crval[0];
	output->crval[1] = current->crval[1];
	output->cdelt[0] = current->cdelt[0];
	output->cdelt[1] = current->cdelt[1];
	output->next     = NULL;
}
