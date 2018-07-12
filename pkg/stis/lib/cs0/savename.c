# include <stdio.h>
# include <stdlib.h>		/* malloc */
# include <string.h>
# include "hstcalerr.h"

# define MAX_NAMES_INCR   20

/* These routines are for saving a string (a file name) in a list of
   names.  If the file name is already present in the list of names,
   SaveName sets oldname to 1; if the name is new, oldname will be 0.

	int InitNames (void);
	int SaveName (char *filename, int *oldname);
	void FreeNames (void);
*/

static char **fnames;		/* array of pointers to file names */
static int maxnames = 0;	/* allocated size of array */
static int nnames = 0;		/* current number of names */
static int ReallocNames (void);

int InitNames (void) {

	int i;

	nnames = 0;
	fnames = malloc (MAX_NAMES_INCR * sizeof (char *));
	if (fnames == NULL) {
	    printf ("ERROR    (InitNames) can't allocate memory.\n");
	    return (OUT_OF_MEMORY);
	}
	maxnames = MAX_NAMES_INCR;

	for (i = 0;  i < maxnames;  i++)
	    fnames[i] = NULL;

	return (0);
}

int SaveName (char *filename, int *oldname) {

/* arguments:
char *filename   i: name to be saved in list
int *oldname     o: true if name already present in list
*/

	int status;
	int i;

	/* Check whether the current filename is already in the list. */
	for (i = 0;  i < nnames;  i++) {
	    if (strcmp (filename, fnames[i]) == 0) {
		*oldname = 1;		/* name was found */
		return (0);
	    }
	}
	*oldname = 0;			/* name not found */

	/* Need more space for names? */
	if (nnames >= maxnames) {
	    if ((status = ReallocNames()))
		return (status);
	}

	/* Add the current filename to the list of names. */
	i = nnames;
	fnames[i] = malloc ((strlen (filename) + 1) * sizeof (char));
	if (fnames[i] == NULL) {
	    printf ("ERROR    (SaveName) can't allocate memory.\n");
	    return (OUT_OF_MEMORY);
	}
	strcpy (fnames[i], filename);
	nnames++;

	return (0);
}

static int ReallocNames (void) {

	int i;
	int newmax;		/* new max number of names */

	newmax = maxnames + MAX_NAMES_INCR;

	fnames = realloc (fnames, newmax * sizeof (char *));
	if (fnames == NULL) {
	    printf ("ERROR    (SaveName) can't allocate memory.\n");
	    return (OUT_OF_MEMORY);
	}

	for (i = maxnames;  i < newmax;  i++)
	    fnames[i] = NULL;

	maxnames = newmax;

	return (0);
}

void FreeNames (void) {

	int i;

	for (i = 0;  i < nnames;  i++)
	    free (fnames[i]);
	free (fnames);
	nnames = 0;
	maxnames = 0;
}
