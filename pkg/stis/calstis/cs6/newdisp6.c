/* This file contains:
	NewDisp6	add a DispRelation record to the list
	CheckDisp6	check for invalid entries in list
	InterpDisp6	return a DispRelation, interpolated on a2center
	FreeDisp6	free memory
	DebugDisp6	print debugging info
internal:
	CopyDisp	copy out one DispRelation




   Revision history:
   ----------------
   20 Feb 97  -  Borrowed from calstis7 (I.Busko)
   24 Feb 97  -  Rename routines to avoid conflicts with cs7 (IB)
   11 Apr 97  -  Changes after code review (IB):
                 - replaced literals by INVALID constant.
                 - explicit cast for malloc-returned pointer.
   08 May 97  -  Conform to new _trl standard (IB)
   02 Jul 97  -  Fixed logic that zeros coefficients (IB)
*/

# include <stdio.h>
# include <stdlib.h>
# include <string.h>

# include "xtables.h"

# include "stis.h"
# include "calstis6.h"
# include "hstcalerr.h"

# define MAX(x,y)  ((x) >= (y) ? (x) : (y))

static void CopyDisp (DispRelation *, DispRelation *);

/* This routine adds a DispRelation record to the list.  The list is
   sorted by a2center.  Before calling NewDisp6 for the first record
   to be inserted in the list, *disp should have been initialized to NULL.
*/

int NewDisp6 (DispRelation **disp, DispRelation *input) {

/* arguments:
DispRelation **disp  io: identifies the first record in the list
DispRelation *input  i: a new record to be inserted into the list
*/

	DispRelation *current, *next, *new;
	int done = 0;

	if (input->ncoeff > MAX_DISP_COEFF) {
	    printf ("ERROR    %d dispersion coefficients.\n",
		input->ncoeff);
	    return (TABLE_ERROR);
	}

	/* Allocate space for the new record, and copy input to new. */
	if ((new = (DispRelation *) malloc (sizeof (DispRelation))) == NULL) {
	    printf ("ERROR    Can't allocate memory.\n");
	    return (OUT_OF_MEMORY);
	}
	CopyDisp (input, new);

	/* Find where this record belongs in the list, sorting by a2center,
	   and insert the record.
	*/
	current = *disp;

	if ((*disp) == NULL) {
	    /* The list is empty; new becomes the first record. */
	    *disp = new;
	    done = 1;
	} else if (new->a2center <= (*disp)->a2center) {
	    /* new is before the beginning of the list; new becomes first. */
	    new->next = *disp;
	    *disp = new;
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
   in the DispRelation list, and it checks that all values of ref_aper
   are the same.  This should be called after reading all rows of the
   table of dispersion coefficients and adding the relevant ones to
   the list.  If the list is empty, this routine returns -1.
*/

int CheckDisp6 (DispRelation **disp) {

/* arguments:
DispRelation **disp  i: identifies the first record in the list
*/

	DispRelation *current, *previous;

	current = *disp;
	if (current == NULL)
	    return (-1);		/* empty list */

	previous = current;
	current = current->next;

	while (current != NULL) {
	    if (current->a2center == previous->a2center) {
		printf (
	"ERROR    Duplicate values of A2CENTER=%.8g in DISPTAB.\n",
			current->a2center);
		return (INVALID);
	    }
	    if (strcmp (current->ref_aper, previous->ref_aper) != 0) {
		printf (
	"ERROR    Different REF_APER values in different rows in DISPTAB.\n");
		return (INVALID);
	    }
	    previous = current;
	    current = current->next;
	}

	return (0);
}

/* This routine frees all memory for the DispRelation list, including
   disp itself.
*/

void FreeDisp6 (DispRelation **disp) {

/* argument:
DispRelation **disp  io: identifies the first record in the list; set to NULL
*/

	DispRelation *current, *next;

	current = *disp;

	while (current != NULL) {
	    next = current->next;
	    free (current);
	    current = next;
	}

	*disp = NULL;
}

/* This routine interpolates the dispersion coefficients on a2center.
   If a2center is outside the range of values in the list, the first (or
   last) record will be copied to output without interpolation; otherwise,
   the two records that have a2center bracketing the input a2center will
   be used, and the coefficients will be linearly interpolated.

   If *output is NULL, memory will be allocated for that record; otherwise,
   it is assumed that *output is actually a valid pointer to a DispRelation,
   i.e. that the memory has been allocated.
*/

int InterpDisp6 (DispRelation **disp, double a2center,
		DispRelation **output) {

/* arguments:
DispRelation **disp   i: identifies the first record in the list
double a2center       i: line number in input image
DispRelation **output o: the dispersion relation interpolated to a2center
*/

	DispRelation *current, *next;
	double p;
	int i;
	int done;

	current = *disp;

	/* Allocate space for the output. */
	if (*output == NULL) {
	    if ((*output = (DispRelation *) malloc (sizeof (DispRelation))) == 
                            NULL) {
		printf ("ERROR    Can't allocate memory.\n");
		return (OUT_OF_MEMORY);
	    }
	}

	done = (current == NULL);

	while (!done) {

	    next = current->next;

	    if (a2center <= current->a2center) {
		/* before beginning of list, or exact match */
		CopyDisp (current, *output);	/* no interpolation needed */
		(*output)->a2center = a2center;
		done = 1;

	    } else if (next == NULL) {

		/* end of list */
		CopyDisp (current, *output);
		(*output)->a2center = a2center;
		done = 1;

	    } else if (current->a2center < a2center &&
				a2center < next->a2center) {

		/* linear interpolation of coefficients */
		p = (a2center - current->a2center) /
			(next->a2center - current->a2center);

		(*output)->a2center = a2center;
		(*output)->ncoeff = MAX (current->ncoeff, next->ncoeff);
		for (i = 0;  i < (*output)->ncoeff;  i++) {
		    (*output)->coeff[i] =
			(1. - p) * current->coeff[i] + p * next->coeff[i];
		}
		for (i = (*output)->ncoeff;  i < MAX_DISP_COEFF;  i++)
		    (*output)->coeff[i] = 0.;
		strcpy ((*output)->ref_aper, current->ref_aper);
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

/* This routine copies the current DispRelation to output, which must
   actually point to allocated space.  Note that in output, the pointer
   to the next record is set to NULL.
*/

static void CopyDisp (DispRelation *current, DispRelation *output) {

/* arguments:
DispRelation *current   i: the current DispRelation record
DispRelation *output    o: receives a copy of current
*/

	int i;

	output->a2center = current->a2center;
	output->ncoeff = current->ncoeff;
	for (i = 0;  i < current->ncoeff;  i++)
	    output->coeff[i] = current->coeff[i];
	for (i = current->ncoeff;  i < MAX_DISP_COEFF;  i++)
	    output->coeff[i] = 0.;
	strcpy (output->ref_aper, current->ref_aper);
	output->next = NULL;
}

/* This routine prints the Y coordinate from the DispRelation list. */

void DebugDisp6 (DispRelation **disp) {

/* argument:
DispRelation **disp  i: identifies the first record in the list
*/

	DispRelation *current;
	int i = 1;

	current = *disp;

	while (current != NULL) {
	    printf ("%3d:  %.8g\n", i, current->a2center);
	    current = current->next;
	    i++;
	}
}
