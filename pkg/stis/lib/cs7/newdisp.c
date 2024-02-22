/* This file contains:
	NewDisp		add a DispRelation record to the list
	CheckDisp	check for invalid entries in list
	InterpDisp	return a DispRelation, interpolated on a2center
	FreeDisp	free memory
internal:
	CopyDisp	copy out one DispRelation
*/

# include <stdio.h>
# include <stdlib.h>
# include <string.h>

# include "stis.h"
# include "calstis7.h"
# include "hstcalerr.h"

# define MAX(x,y)  ((x) >= (y) ? (x) : (y))

static void CopyDisp (DispRelation *, DispRelation *);

/* This routine adds a DispRelation record to the list.  The list is
   sorted by a2center.  Before calling NewDisp for the first record
   to be inserted in the list, *disp should have been initialized to NULL.

   Phil Hodge, 1998 Oct 5:
	Change status values 1212 and 1214 to GENERIC_ERROR_CODE.

   Phil Hodge, 2001 Mar 7:
	Remove DebugDisp.
*/

int NewDisp (DispRelation **disp, DispRelation *input) {

/* arguments:
DispRelation **disp  io: identifies the first record in the list
DispRelation *input  i: a new record to be inserted into the list
*/

	DispRelation *current, *next, *newrec;
	int done = 0;

	if (input->ncoeff > MAX_DISP_COEFF) {
	    printf ("ERROR    (NewDisp) %d dispersion coefficients.\n",
		input->ncoeff);
	    return (TABLE_ERROR);
	}

	/* Allocate space for the new record, and copy input to newrec. */
	if ((newrec = malloc (sizeof (DispRelation))) == NULL) {
	    printf ("ERROR    Can't allocate memory in NewDisp.\n");
	    return (OUT_OF_MEMORY);
	}
	CopyDisp (input, newrec);

	/* Find where this record belongs in the list, sorting by a2center,
	   and insert the record.
	*/
	current = *disp;

	if ((*disp) == NULL) {
	    /* The list is empty; newrec becomes the first record. */
	    *disp = newrec;
	    done = 1;
	} else if (newrec->a2center <= (*disp)->a2center) {
	    /* newrec is before the beginning of the list;
		newrec becomes the first.
	    */
	    newrec->next = *disp;
	    *disp = newrec;
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
   in the DispRelation list, and it checks that all values of ref_aper
   are the same.  This should be called after reading all rows of the
   table of dispersion coefficients and adding the relevant ones to
   the list.  If the list is empty, this routine returns -1.
*/

int CheckDisp (DispRelation **disp) {

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
		return (GENERIC_ERROR_CODE);
	    }
	    if (strcmp (current->ref_aper, previous->ref_aper) != 0) {
		printf (
	"ERROR    Different REF_APER values in different rows in DISPTAB.\n");
		return (GENERIC_ERROR_CODE);
	    }
	    previous = current;
	    current = current->next;
	}

	return (0);
}

/* This routine frees all memory for the DispRelation list, including
   disp itself.
*/

void FreeDisp (DispRelation **disp) {

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

int InterpDisp (DispRelation **disp, double a2center,
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
	    if ((*output = malloc (sizeof (DispRelation))) == NULL) {
		printf ("ERROR    Can't allocate memory in InterpDisp.\n");
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
