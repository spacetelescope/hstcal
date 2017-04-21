/* This file contains:
	InitRefFile
	NewRefFile
	FindRefFile
	FreeRefFile
*/

# include <stdlib.h>
# include <string.h>

# include "acs.h"
# include "err.h"

/* These routines are for managing a list of keyword & value pairs
   for reference file names.

   Warren Hack, 1998 May 26:
	Revised for ACS, although still basically the STIS code.
*/

/* This initialization routine must be called before any others in this
   file.  The first record is a dummy that just points to the first
   real record.
*/

void InitRefFile (RefFileInfo *ref) {

	strcpy (ref->keyword, "dummy");
	strcpy (ref->filename, "N/A");
	ref->next = NULL;
}

/* Add a new keyword,filename pair to the list.  If the keyword is already
   in the list, then filename will replace the current value.
*/

int NewRefFile (RefFileInfo *ref, char *keyword, char *filename) {

/* arguments:
RefFileInfo *ref   io: the list of reference info
char *keyword      i: keyword name of reference file
char *filename     i: name of reference file
*/

	extern int status;

	RefFileInfo *previous, *current, *newrec;
	int foundit = 0;

	/* Search for the keyword in the list. */
	previous = ref;
	current = ref->next;		/* skip over dummy first record */
	while (!foundit) {
	    if (current == NULL) {
		current = previous;	/* back up to last record */
		break;
	    } else if (strcmp (keyword, current->keyword) == 0) {
		foundit = 1;
	    } else {
		previous = current;
		current = current->next;
	    }
	}

	if (foundit) {

	    /* Replace current filename. */
	    strcpy (current->filename, filename);

	} else {

	    /* Allocate space for a new record, and copy info into it. */
	    if ((newrec = malloc (sizeof (RefFileInfo))) == NULL)
		return (status = OUT_OF_MEMORY);

	    strcpy (newrec->keyword, keyword);
	    strcpy (newrec->filename, filename);
	    newrec->next = NULL;

	    /* Add the new record to the end of the list. */
	    current->next = newrec;
	}

	return (status);
}

/* Find a reference file keyword in the list, and return the name.
   If the keyword is found in the list, foundit will be set to 1;
   otherwise, foundit will be set to 0.
*/

void FindRefFile (RefFileInfo *ref, char *keyword,
	char *filename, int *foundit) {

/* arguments:
RefFileInfo *ref   i: the list of reference info
char *keyword      i: keyword name of reference file
char *filename     o: name of reference file
int *foundit       o: true if keyword was found in list
*/

	RefFileInfo *current;
	int done = 0;

	filename[0] = '\0';
	*foundit = 0;
	current = ref;
	while (!done) {
	    if (current == NULL) {
		done = 1;
	    } else if (strcmp (keyword, current->keyword) == 0) {
		strcpy (filename, current->filename);
		*foundit = 1;
		done = 1;
	    } else {
		current = current->next;
	    }
	}
}

/* Free memory allocated for this list. */

void FreeRefFile (RefFileInfo *ref) {

/* arguments:
RefFileInfo *ref   io: the list of reference info
*/

	RefFileInfo *current, *next;

	current = ref->next;		/* don't free the first record */
	while (current != NULL) {
	    next = current->next;
	    free (current);
	    current = next;
	}
	ref->next = NULL;
}
