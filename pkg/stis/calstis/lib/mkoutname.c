/* This file contains:
	MkOutName	constructs an output name from the input
	DefaultExtn	appends the default filename extension
*/

# include <stdio.h>
# include <stdlib.h>
# include <string.h>

# include "stis.h"
# include "hstcalerr.h"

# define FITS_EXTN  ".fits"	/* default extension */

static int FindExtn (char *);
static int strcatN (char *, char *, int);

/* This routine constructs the output file name based on the input name
   if the output name is null.  If output is not null but lacks a filename
   extension, the default extension will appended.  If output already has
   an extension, then no change will be made to output.

   If output is null, the input name will first be copied to output.
   The argument isuffix contains an array of recognized suffixes for
   the input file name, and osuffix contains an array of corresponding
   suffixes for the output file name.  Both arrays should be the same
   length, nsuffix.  For each suffix in isuffix, the input name will be
   checked to see if The end of the input name (before the filename
   extension) will be compared with each suffix in isuffix until a
   match is found, in which case that isuffix will be replaced by the
   corresponding osuffix in the output name.  If no match is found
   with any of the isuffix strings, the first string in osuffix will be
   appended to the root name in the output, rather than replacing the
   current suffix.  If input does not contain a filename extension
   (but it should), the default extension will be appended to output.

   For example, given these suffixes:

	char *isuffix[] = {"_flt", "_crj"};
	char *osuffix[] = {"_x1d", "_sx1"};
	int nsuffix = 2;

   here are some examples of input and resulting output names:

           input            output
           -----            ------
          a                a_x1d.fits
          a_raw.fits       a_raw_x1d.fits
          a_flt.fits       a_x1d.fits
          a_crj.fits       a_sx1.fits

   Phil Hodge, 1998 Feb 5:
	File created.
*/

int MkOutName (char *input, char **isuffix, char **osuffix, int nsuffix,
		char *output, int maxch) {

/* arguments:
char *input        i: name of input FITS file
char *isuffix[]    i: suffixes expected for input, e.g. "_raw"
char *osuffix[]    i: suffixes to append for output, e.g. "_flt"
int nsuffix        i: length of isuffix and osuffix arrays
char *output       io: name of output FITS file
int maxch          i: maximum size of output
*/

	int status;

	char *extn;		/* extension on input (or default) */

	int is_len;		/* length of current isuffix */
	int tr_len;		/* length of truncated input name */
	int i;			/* loop index */
	int dotlocn;		/* location of '.' in input name */

	if (output[0] == '\0') {

	    extn = calloc (1 + strlen(input) + strlen(FITS_EXTN),
			sizeof (char));
	    if (extn == NULL)
		return (OUT_OF_MEMORY);

	    /* Find the extension (if any) on the input name. */
	    dotlocn = FindExtn (input);

	    if ((status = strcatN (output, input, maxch)))
		return (status);

	    if (dotlocn >= 0) {
		strcpy (extn, &input[dotlocn]);	/* save extension from input */
		output[dotlocn] = '\0';		/* truncate at '.' */
	    } else {
		strcpy (extn, FITS_EXTN);	/* default extension */
	    }

	    /* If the input name ends in one of the expected suffixes,
		chop it off of the output name before appending the output
		suffix and extension.
	    */
	    for (i = 0;  i < nsuffix;  i++) {
		is_len = strlen (isuffix[i]);
		tr_len = strlen (output);	/* length of truncated name */
		if (tr_len >= is_len) {
		    if (strcmp (output+tr_len-is_len, isuffix[i]) == 0) {
			output[tr_len-is_len] = '\0';
			if ((status = strcatN (output, osuffix[i], maxch)))
			    return (status);
			break;
		    }
		}
	    }
	    /* Append first output suffix if none of the input suffixes
		was found.
	    */
	    if (i >= nsuffix) {
		if ((status = strcatN (output, osuffix[0], maxch)))
		    return (status);
	    }
	    if ((status = strcatN (output, extn, maxch)))
		return (status);
	    free (extn);

	} else {

	    /* An output name was specified.  Check that it includes an
		extension, and if not, append the default extension.
	    */
	    dotlocn = FindExtn (output);
	    if (dotlocn < 0) {
		if ((status = strcatN (output, FITS_EXTN, maxch)))
		    return (status);
	    }		/* else the output name is OK as is */
	}

	return (0);
}

/* This routine appends the default extension if no extension is found
   in the input name.  The input name is modified in-place.  The only
   error condition is when the sum of the string lengths exceeds maxch.
*/

int DefaultExtn (char *input, int maxch) {

	int status;
	int dotlocn;		/* location of '.' in input name */

	status = 0;
	dotlocn = FindExtn (input);
	if (dotlocn < 0)
	    status = strcatN (input, FITS_EXTN, maxch);

	return (status);
}

/* This function returns the index of '.' in fname, or -1 if it wasn't
  found.  The searching starts at the end of fname and works backward.
  Note that searching for '.' will stop if '$', '/', or ']' is
  encountered, as these could be used as parts of a directory name.
  Note also that we do not require the extension to have any particular
  value.  fname may begin with '.'.
*/

static int FindExtn (char *fname) {

	int dotlocn = -1;
	int i;

	for (i = strlen(fname)-1;  i >= 0;  i--) {
	    if (fname[i] == '.') {
		dotlocn = i;
		break;
	    }
	    /* Have we reached a directory prefix? */
	    if (fname[i] == '$' || fname[i] == '/' || fname[i] == ']')
		break;
	}

	return (dotlocn);
}

/* This routine concatenates instr to outstr, but only if the sum of
   their lengths is no longer than maxch.
*/

static int strcatN (char *outstr, char *instr, int maxch) {

	int status;

	if (strlen (instr) + strlen (outstr) > maxch) {
	    printf ("ERROR    (MkOutName) strings are too long: \\\n");
	    printf ("ERROR    `%s' + `%s'\n", outstr, instr);
	    status = 2011;
	} else {
	    strcat (outstr, instr);
	    status = 0;
	}

	return (status);
}
