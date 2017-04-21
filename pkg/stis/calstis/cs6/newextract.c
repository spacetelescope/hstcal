
/* This file contains:
	NewXtract	add a XtractInfo record to the list
	ReturnXtract	extract a XtractInfo record from the list
	RangeXtract	find the range of sporders in the list
	FreeXtract	free memory
	DebugXtract	print debugging info
internal:
	CopyXtract	copy out one XtractInfo




   Revision history:
   ----------------
   20 Feb 97  -  Implemented, with code adapted from similar routines
                 in calstis7 (I.Busko)
   11 Apr 97  -  Changes after code review (IB):
                 - replaced literals by INVALID constant.
                 - explicit cast for malloc-returned pointer.
                 - set last pointer to NULL when freeing XtractInfo list.
   16 Apr 97  -  Replaced scalar bktilt by a polynomial description. Also
                 implemented tilted spectrum extraction box (what for ?) (IB)
   08 May 97  -  Conform to new _trl standard (IB)
*/


# include <stdio.h>
# include <stdlib.h>
# include <string.h>

# include "xtables.h"
# include "stis.h"
# include "calstis6.h"
# include "err.h"


static void CopyXtract (XtractInfo *, XtractInfo *);


/* This routine adds a XtractInfo record to the list.  The list is sorted
   by sporder.  Before calling NewXtract for the first record to be inserted
   in the list, *extract should have been initialized to NULL.
*/

int NewXtract (XtractInfo **extract, XtractInfo *input) {

/* arguments:
XtractInfo **extract  io: identifies the first record in the list
XtractInfo *input     i:  a new record to be inserted into the list
*/

	XtractInfo *current, *next, *new;
	int done;

	done = 0;

	if (input->ncoeffsl > MAX_SLIT_COEFF) {
	    printf ("ERROR    %d slit coefficients.\n",
		input->ncoeffsl);
	    return (TABLE_ERROR);
	}
	if (input->ncoeffbk > MAX_BACK_COEFF) {
	    printf ("ERROR    %d background coefficients.\n",
		input->ncoeffbk);
	    return (TABLE_ERROR);
	}

	/* Allocate space for the new record, and copy to memory. */
	if ((new = (XtractInfo *) malloc (sizeof (XtractInfo))) == NULL) {
	    printf ("ERROR    Can't allocate memory.\n");
	    return (OUT_OF_MEMORY);
	}

	CopyXtract (input, new);

	/* Find where this record belongs in the list, sorting by sporder,
	   and insert the record.
	*/
	current = *extract;

	if ((*extract) == NULL) {
	    /* The list is empty; new becomes the first record. */
	    *extract = new;
	    done = 1;
	} else if (new->sporder <= (*extract)->sporder) {
	    /* new is before the beginning of the list; new becomes first. */
	    new->next = *extract;
	    *extract = new;
	    done = 1;
	}

	while (!done) {

	    if (current->next == NULL) {
		current->next = new;		/* add to end of list */
		done = 1;

	    } else {

		next = current->next;

		if (new->sporder <= next->sporder) {
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



/* This routine returns a XtractInfo record from the list.  The record is
   selected by sporder.

   If *output is NULL, memory will be allocated for that record; otherwise,
   it is assumed that *output is actually a valid pointer to a XtractInfo,
   i.e. that the memory has been allocated.
*/

int ReturnXtract (XtractInfo **extract, int sporder, XtractInfo **output) {

/* arguments:
XtractInfo **extract  i: identifies the first record in the list
XtractInfo **output   o: the record that matches sporder
*/

	XtractInfo *current;
	int foundit = 0;

	current = *extract;

	/* Allocate space for the output. */
	if (*output == NULL) {
	    if ((*output = (XtractInfo *) malloc (sizeof (XtractInfo))) == 
                           NULL) {
		printf ("ERROR    Can't allocate memory.\n");
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
	    CopyXtract (current, *output);
	else
	    return (INVALID);	/* sporder not found, or empty list */

	return (0);
}



/* This routine finds the minimum and maximum values of sporder in the
   XtractInfo list.  The list is also checked to ensure that the sporder
   numbers are consecutive and there are no duplicates.
*/

int RangeXtract (XtractInfo **extract, int *minorder, int *maxorder) {

/* arguments:
XtractInfo **extract      i: identifies the first record in the list
int *minorder, *maxorder  o: minimum and maximum values of sporder
*/

	XtractInfo *current, *previous;

	current = *extract;
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
		"ERROR    Duplicate order number %d in XTRACTAB \n",
			current->sporder);
		    return (INVALID);
		} else {
		    printf (
		"ERROR    Order numbers in XTRACTAB jump from %d to %d\n",
		    previous->sporder, current->sporder);
		    return (INVALID);
		}
	    }
	    *maxorder = current->sporder;
	    previous = current;
	    current = current->next;
	}

	return (0);
}



/* This routine frees all memory for the XtractInfo list.  The value of
   extract is not changed.
*/

void FreeXtract (XtractInfo **extract) {

/* argument:
XtractInfo **extract  io: identifies the first record in the list
*/

	XtractInfo *current, *next;

	current = *extract;

	while (current != NULL) {
	    next = current->next;
	    free (current);
	    current = next;
	}

	*extract = NULL;
}



/* This routine copies the current XtractInfo to output, which must
   actually point to allocated space.  Note that in output, the pointer
   to the next record is set to NULL.
*/

static void CopyXtract (XtractInfo *current, XtractInfo *output) {

/* argument:
XtractInfo *current  i: the current XtractInfo record
XtractInfo *output   o: receives a copy of current
*/
	int i;

	output->sporder     = current->sporder;
	output->extrsize    = current->extrsize;
	output->bksize[0]   = current->bksize[0];
	output->bksize[1]   = current->bksize[1];
	output->bkoffset[0] = current->bkoffset[0];
	output->bkoffset[1] = current->bkoffset[1];
	output->ncoeffsl    = current->ncoeffsl;
	output->ncoeffbk    = current->ncoeffbk;
	for (i = 0;  i < current->ncoeffsl;  i++)
	    output->sltcoeff[i] = current->sltcoeff[i];
	for (i = 0;  i < current->ncoeffbk;  i++)
	    output->bktcoeff[i] = current->bktcoeff[i];
	output->backord     = current->backord;
	output->maxsearch   = current->maxsearch;
	strcpy (output->xtracalg, current->xtracalg);
	output->next = NULL;
}



/* This routine prints the spectral order number from the XtractInfo list. */

void DebugXtract (XtractInfo **extract) {

/* argument:
XtractInfo **extract  i: identifies the first record in the list
*/

	XtractInfo *current;
	int i = 1;

	current = *extract;

	while (current != NULL) {
	    printf ("%3d:  %3d\n", i, current->sporder);
	    current = current->next;
	    i++;
	}
}
