# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include "hstcalerr.h"
# include "acs.h"	/* for message output */

int MkName (char *input, char *isuffix, char *osuffix, char *outextn,
		char *output, int maxch) {

/* arguments:
char *input        i: name of input FITS file
char *output       o: name of output FITS file
char *isuffix      i: suffix expected for input, e.g. "_raw"
char *osuffix      i: suffix to append for output, e.g. "_flt"
char *outextn	   i: extension to append for output, e.g. ".trl" (optional)
int maxch          i: maximum size of output
*/

	extern int status;

	char fitsext[] = ".fits";	/* default extension */
	char *extn;		/* extension on input (or default) */

	size_t i_len;		/* length of input name */
	size_t is_len;		/* length of isuffix */
	size_t fits_len;	/* length of fitsext */
	size_t tr_len;		/* length of truncated input name */
	size_t dotlocn;		/* location of '.' in input name */

	i_len = strlen (input);
	is_len = strlen (isuffix);
	fits_len = strlen (fitsext);

	// TODO: maxch should be an unsigned type
	if (i_len > (size_t) maxch) {
	    trlerror ("(MkName) Input name is too long.");
	    return (status = INVALID_FILENAME);
	}

	extn = calloc (1 + is_len + fits_len, sizeof (char));
	if (extn == NULL)
	    return (status = OUT_OF_MEMORY);

	/* Find the extension (if any) on the input name. */
	dotlocn = 0;
	for (size_t i = i_len ? i_len - 1 : 0;  i != 0;  i--) {
	    if (input[i] == '.') {
		dotlocn = i;
		break;
	    }
	    /* Check for various special characters. */
	    if (input[i] == '$')
		break;
	    else if (input[i] == '/')
		break;
	    else if (input[i] == ']')
		break;
	}

	strcpy (output, input);
	
	if (dotlocn > 0) {
	    strcpy (extn, &input[dotlocn]);	/* save extension from input */
	    output[dotlocn] = '\0';		/* truncate at '.' */
	} else {
	    strcpy (extn, fitsext);		/* default extension */
	}
	
	/* If an extension was specified for the output, use it */
	if (strcmp(outextn,"") != 0) {
		strcpy (extn, outextn);		/* use specified extension */
	}
	
	/* If the input name ends in the expected suffix, chop it off of
	   the output name before appending the output suffix and extension.
	*/
	tr_len = (int)strlen (output);		/* length of truncated name */
	if (tr_len >= is_len) {
	    if (strcmp (output+tr_len-is_len, isuffix) == 0)
		output[tr_len-is_len] = '\0';
	}

	strcat (output, osuffix);
	strcat (output, extn);
	
	free (extn);
	return (status);
}
