/* This file contains:
	NewCoord6	add a CoordInfo record to the list
	ReturnCoord6	extract a CoordInfo record from the list
	RangeCoord6	find the range of sporders in the list
	FreeCoord6	free memory
	DebugCoord6	print debugging info
internal:
	CopyCoord	copy out one CoordInfo




   Revision history:
   ----------------
   14 Nov 97  -  Borrowed from calstis7 (I.Busko)
   14 Nov 97  -  Rename routines to avoid conflicts with cs7 (IB)
   07 Jan 00  -  Add npix handling (IB)
*/

# include <stdio.h>
# include <stdlib.h>
# include <string.h>

# include "xtables.h"
# include "stis.h"
# include "calstis6.h"
# include "err.h"

static void CopyCoord (CoordInfo *, CoordInfo *);

/* This routine adds a CoordInfo record to the list.  The list is sorted
   by sporder.  Before calling NewCoord6 for the first record to be inserted
   in the list, *coords should have been initialized to NULL.
*/

int NewCoord6 (CoordInfo **coords, CoordInfo *input) {

/* arguments:
CoordInfo **coords  io: identifies the first record in the list
CoordInfo *input    i: a new record to be inserted into the list
*/

	CoordInfo *current, *next, *newrec;
	int done;

	done = 0;

	/* Allocate space for the new record, and copy to memory. */
	if ((newrec = (CoordInfo *) malloc (sizeof (CoordInfo))) == NULL) {
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

int ReturnCoord6 (CoordInfo **coords, int sporder, CoordInfo **output) {

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
	    (*output)->next = NULL;
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
	    return (INVALID);	/* sporder not found, or empty list */

	return (0);
}


/* This routine finds the minimum and maximum values of sporder in the
   CoordInfo list.  The list is also checked to ensure that the sporder
   numbers are consecutive and there are no duplicates.
*/

int RangeCoord6 (CoordInfo **coords, int *minorder, int *maxorder) {

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
		    return (1212);
		} else {
		    printf (
		"ERROR    Order numbers in SDC table jump from %d to %d.\n",
		    previous->sporder, current->sporder);
		    return (1214);
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

void FreeCoord6 (CoordInfo **coords) {

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

	*coords = NULL;
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
	output->cdelt2   = current->cdelt2;
	output->npix     = current->npix;
	output->next     = NULL;
}

/* This routine prints the spectral order number from the CoordInfo list. */

void DebugCoord6 (CoordInfo **coords) {

/* argument:
CoordInfo **coords  i: identifies the first record in the list
*/

	CoordInfo *current;
	int i = 1;

	current = *coords;

	while (current != NULL) {
	    printf ("%3d:  %3d\n", i, current->sporder);
	    current = current->next;
	    i++;
	}
}
