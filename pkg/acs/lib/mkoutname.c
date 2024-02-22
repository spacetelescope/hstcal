/* This file contains:
	MkOutName	constructs an output name from the input
	DefaultExtn	appends the default filename extension
	MkNewExtn	replaces current extension with new extension
*/

# include <stdio.h>
# include <stdlib.h>
# include <string.h>
#include "hstcal.h"
# include "hstcalerr.h"
# include "acs.h"	/* for message output */

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
	
   Warren Hack, 1998 Nov 17:
   	Added MkNewExtn to support trailer file name conventions
   Warren Hack, 2004 Mar 19:
    Added calls to free() to clean up memory upon exiting for errors.	
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

	extern int status;

	char *extn;		/* extension on input (or default) */

	int is_len;		/* length of current isuffix */
	int tr_len;		/* length of truncated input name */
	int i;			/* loop index */
	int dotlocn;		/* location of '.' in input name */
    
	if (output[0] == '\0') {

	    extn = calloc (1 + strlen(input) + strlen(FITS_EXTN),
			sizeof (char));
	    if (extn == NULL)
		return (status = OUT_OF_MEMORY);

	    /* Find the extension (if any) on the input name. */
	    dotlocn = FindExtn (input);

	    if (strcatN (output, input, maxch)){
            free(extn);
		    return (status);
        }

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
					if (strcatN (output, osuffix[i], maxch)) {
                        free(extn);
			    		return (status);
                    }
					break;
		    	}
			}
	    }
	    /* Append first output suffix if none of the input suffixes
		was found.
	    */
	    if (i >= nsuffix) {
			if (strcatN (output, osuffix[0], maxch)) {
                free(extn);
		    	return (status);
            }
	    }
	    if (strcatN (output, extn, maxch)) {
            free(extn);
			return (status);
        }
	    free (extn);

	} else {

	    /* An output name was specified.  Check that it includes an
		extension, and if not, append the default extension.
	    */
	    dotlocn = FindExtn (output);
	    if (dotlocn < 0) {
			if (strcatN (output, FITS_EXTN, maxch))
		    	return (status);
	    }		/* else the output name is OK as is */
	}

	return (status);
}

/* This routine appends the default extension if no extension is found
   in the input name.  The input name is modified in-place.  The only
   error condition is when the sum of the string lengths exceeds maxch.
*/

int DefaultExtn (char *input, int maxch) {

	extern int status;
	int dotlocn;		/* location of '.' in input name */

	dotlocn = FindExtn (input);
	if (dotlocn < 0)
	    strcatN (input, FITS_EXTN, maxch);

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

	extern int status;

	if (strlen (instr) + strlen (outstr) > maxch) {
	    trlerror ("(MkOutName) strings are too long:");
	    sprintf (MsgText, "`%s' + `%s'", outstr, instr);
	    trlerror (MsgText);
		status = INVALID_FILENAME;
	} else {
	    strcat (outstr, instr);
	}

	return (status);
}

/* This function replaces the extension found in 'input' with
	the extension specified by 'newext'.
*/
int MkNewExtn (char *input, char *newext) {
	
	extern int status;
	
	int dotlocn;
	char *output;
	
	/* Set up internal buffer for working with filename */
	output = (char *) calloc(strlen(input)+strlen(newext)+1, sizeof(char) );
	if (output == NULL) return (status = OUT_OF_MEMORY);
	
	output[0] = '\0';
	
	/* Find the extension (if any) on the input name. */
	dotlocn = FindExtn (input);
        
	strcpy (output, input);		/* copy input to output */
        
	if (dotlocn >= 0) {
		output[dotlocn] = '\0';		/* truncate at '.' */
		strcat (&output[dotlocn], newext);	/* append new extension */
	} else {
		/* If no extension was found in 'input',
			simply append given extension onto 'output'
		*/
		strcat (output, newext);
	}
    
	/* Copy new filename back into 'input' */
	strcpy (input, output);
    
	free (output);
	return (status);
}
