# include <stdio.h>
# include <stdlib.h>
# include <string.h>

# include "c_iraf.h"
# include "xtables.h"

# include "stis.h"
# include "calstis0.h"
# include "hstcalerr.h"

static int SpecialCheck (char *, char *, int, int *, int *);
static void MissingFile (char *, char *, int *);
static int DoesFileExist (char *, int *);

/* This function checks for the existence of the reference files that
   we need.  There are a couple of cases where a reference file need not
   be specified (dqitab, if CCD, and flat fields), but if a file name
   was given, we require that the file exist.

   This function may be called more than once, e.g. for science and
   wavecal RefFileInfo structs.  Note that InitNames must have been
   called prior to calling this routine, as SaveName is called.

   The keywords that are given special treatment are BPIXTAB,
   PFLTFILE, DFLTFILE, LFLTFILE, PCTAB, APERTAB, and WCPTAB.

   Phil Hodge, 1997 Dec 10:
	File created.

   Phil Hodge, 1998 Jan 13:
	Include PCTAB and APERTAB as "special" files.

   Phil Hodge, 1998 July 14:
	Include MOFFTAB as a "special" file.

   Phil Hodge, 1998 Nov 20:
	Include WCPTAB as a "special" file.

   Phil Hodge, 1999 Dec 6:
	Use tbtopn instead of openInputImage to check for the existence
	of both tables and (by appending "[0]") images.

   Phil Hodge, 2001 May 4:
	Remove MOFFTAB as a "special" file.  We no longer use this table.

   Phil Hodge, 2003 Jun 20:
	Remove the check on wbiafile.
*/

void RefExist (RefFileInfo *ref, int detector, int *missing) {

/* arguments:
RefFileInfo *ref   i: the list of reference info
int detector       i: detector code, to check for CCD
int *missing       o: a count of missing (or blank) reference files
*/

	RefFileInfo *current;
	int name_specified;		/* true if name is not blank or N/A */
	int exist;			/* EXISTS_YES or EXISTS_NO */
	int flat_not_specified = 0;

	current = ref->next;		/* skip dummy first record */
	while (current != NULL) {

	    /* Check if file exists. */
	    exist = DoesFileExist (current->filename, &name_specified);
	    if (exist != EXISTS_YES) {

		/* If a name was specified, the file must exist.
		   If this is a "special" keyword, treat it individually;
		   otherwise, log it as missing.
		*/
		if (name_specified ||
		    !SpecialCheck (current->keyword, current->filename,
				detector,
				&flat_not_specified, missing))
		    MissingFile (current->keyword, current->filename, missing);
	    }

	    current = current->next;
	}

	/* All three flats not specified?  That's not OK. */
	if (flat_not_specified == 3) {
	    (*missing)++;
	    printf ("Warning  PFLTFILE, DFLTFILE, LFLTFILE are all blank\n");
	}
}

static void MissingFile (char *keyword, char *filename, int *missing) {

	printf ("Warning  %s `%s' not found or can't open.\n",
				keyword, filename);
	(*missing)++;
}

/* This routine checks whether a file exists by trying to open it
   read-only.  If it does exist, the macro EXISTS_YES will be returned
   as the function value; if not, EXISTS_NO will be returned.  In addition,
   the argument name_specified will be set to true unless the filename
   is blank or "N/A", as determined by calling GotFileName.
*/

static int DoesFileExist (char *filename, int *name_specified) {

	char *pname;		/* filename with [0] appended */
	int oldname;		/* already included in our list of names? */
	int flag;
	IRAFPointer tp;		/* for checking if file exists */
	int SaveName (char *, int *);

	if (!GotFileName (filename)) {

	    *name_specified = 0;
	    flag = EXISTS_NO;

	} else {

	    *name_specified = 1;

	    /* Save this name in the list, so we don't need to
		check it more than once.
	    */
	    if (SaveName (filename, &oldname))
		return (EXISTS_NO);

	    if (oldname) {
		/* If it really doesn't exist, we already counted it. */
		flag = EXISTS_YES;
	    } else {
		tp = c_tbtopn (filename, IRAF_READ_ONLY, 0);
		if (c_iraferr()) {
		    clear_cvoserr();
		    /* +10 is more than we need; +4 would do for "[0]" */
		    pname = calloc (strlen (filename) + 10, sizeof (char));
		    if (pname == NULL)
			return (OUT_OF_MEMORY);
		    strcpy (pname, filename);
		    strcat (pname, "[0]");	/* open primary header */
		    tp = c_tbtopn (pname, IRAF_READ_ONLY, 0);
		    free (pname);
		    if (c_iraferr()) {
			flag = EXISTS_NO;
			clear_cvoserr();
		    } else {
			flag = EXISTS_YES;
			c_tbtclo (tp);
		    }
		} else {
		    flag = EXISTS_YES;
		    c_tbtclo (tp);
		}
	    }
	}

	return (flag);
}

/* This routine checks whether the current keyword is "special".  If it is,
   it is handled individually.

   The function value will be one if it is a special keyword, and it will
   be zero otherwise.
*/

static int SpecialCheck (char *keyword, char *filename,
		int detector,
		int *flat_not_specified, int *missing) {

	int flag = 0;		/* to be returned as the function value */

	if (streq_ic (keyword, "BPIXTAB")) {

	    /* We could do DQICORR just to check saturation,
		if the detector is the CCD.
	    */
	    if (detector != CCD_DETECTOR)
		MissingFile (keyword, filename, missing);
	    flag = 1;

	} else if (streq_ic (keyword, "PFLTFILE") ||
		   streq_ic (keyword, "DFLTFILE") ||
		   streq_ic (keyword, "LFLTFILE")) {

	    /* It's OK for some flats to not be specified, but keep count. */
	    (*flat_not_specified)++;
	    flag = 1;

	} else if (streq_ic (keyword, "PCTAB")) {

	    flag = 1;		/* this keyword might be absent */

	} else if (streq_ic (keyword, "WCPTAB")) {

	    flag = 1;		/* this keyword might be absent */

	} else if (streq_ic (keyword, "APERTAB")) {

	    /* It's OK for this to not be specified only for PHOTCORR. */
	    flag = 1;
	}

	return (flag);
}
