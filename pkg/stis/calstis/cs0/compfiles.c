# include <stdio.h>
# include <string.h>

# include "stis.h"
# include "calstis0.h"

static int SpecialKey (char *);
static void SpecialComp (RefFileInfo *, RefFileInfo *, int);
static void NoComparison (char *, char *, char *);

/* This routine compares reference file names read from the science file
   header and the wavecal header.  For each entry in the wavecal list,
   the keyword is searched for in the science file list.  If the keyword
   is found, the corresponding file names in the wavecal and science file
   lists are compared; if they differ, a warning is printed, and the name
   from the science file list replaces the name in the wavecal list.

   The samebin argument is used when comparing BIASFILE values.

   The keywords that are given special treatment are BIASFILE and
   LAMPTAB.

   Phil Hodge, 1997 Dec 10:
	Function created.

   Phil Hodge, 2003 Jun 20:
	Remove the test on wbiafile.
*/

void CompFiles (RefFileInfo *sciref, RefFileInfo *wavref, int samebin) {

/* arguments:
RefFileInfo *sciref  i: the list of reference info from the science header
RefFileInfo *wavref  i: the list of reference info from the wavecal header
int samebin          i: true if same binning and gain in science & wavcal
*/

	RefFileInfo *wcurrent;		/* current record in wavecal list */
	char filename[STIS_FITS_REC+1];
	int foundit;

	wcurrent = wavref;
	while (wcurrent != NULL) {

	    if (SpecialKey (wcurrent->keyword)) {

		SpecialComp (sciref, wcurrent, samebin);

	    } else {

		FindRefFile (sciref, wcurrent->keyword, filename, &foundit);

		if (foundit &&
		    !streq_ic (filename, wcurrent->filename)) {

		    NoComparison (wcurrent->keyword,
				filename, wcurrent->filename);
		    strcpy (wcurrent->filename, filename);
		}
	    }

	    wcurrent = wcurrent->next;
	}
}

static int SpecialKey (char *keyword) {

	if (streq_ic (keyword, "dummy"))
	    return (1);

	if (streq_ic (keyword, "LAMPTAB"))
	    return (1);

	if (streq_ic (keyword, "BIASFILE"))
	    return (1);

	return (0);
}

static void SpecialComp (RefFileInfo *sciref, RefFileInfo *wcurrent,
		int samebin) {

	char filename[STIS_FITS_REC+1];
	int foundit;

	if (streq_ic (wcurrent->keyword, "dummy"))
	    return;

	if (streq_ic (wcurrent->keyword, "LAMPTAB")) {

	    /* If LAMPTAB is blank or N/A in science file header, that's
		OK, just replace it in the science file list with the value
		from the wavecal header.
	    */

	    FindRefFile (sciref, wcurrent->keyword, filename, &foundit);

	    if (foundit) {
		if (GotFileName (filename)) {
		    if (!streq_ic (filename, wcurrent->filename)) {
			NoComparison (wcurrent->keyword,
				filename, wcurrent->filename);
			strcpy (wcurrent->filename, filename);
		    }
		} else {
		    /* copy lamptab name from wavecal to science list */
		    if (NewRefFile (sciref,
				wcurrent->keyword, wcurrent->filename))
			return;
		}
	    }
	}

	if (streq_ic (wcurrent->keyword, "BIASFILE")) {

	    if (samebin) {
		/* Compare BIASFILE in science and wavecal. */
		FindRefFile (sciref, "BIASFILE", filename, &foundit);
		if (foundit) {
		    if (!streq_ic (filename, wcurrent->filename)) {
			NoComparison (wcurrent->keyword,
				filename, wcurrent->filename);
			strcpy (wcurrent->filename, filename);
		    }
		}
	    }
	}
}


/* This routine prints warning messages. */

static void NoComparison (char *keyword, char *sciname, char *wavname) {

	printf (
	"Warning  %s is not the same in the science and wavecal headers; \\\n",
		keyword);
	printf (
	"Warning  values in the science and wavecal are respectively: \\\n");
	printf ("Warning  `%s' \\\n", sciname);
	printf ("Warning  `%s' \\\n", wavname);
	printf (
"Warning  The value from the science header will be used for the wavecal.\n");
}
