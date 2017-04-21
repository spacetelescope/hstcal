# include <stdio.h>
# include <stdlib.h>
# include <string.h>

# include "stis.h"
# include "err.h"

int MkName (char *input, char *isuffix, char *osuffix,
		char *output, int maxch) {

/* arguments:
char *input        i: name of input FITS file
char *output       o: name of output FITS file
char *isuffix      i: suffix expected for input, e.g. "_raw"
char *osuffix      i: suffix to append for output, e.g. "_flt"
int maxch          i: maximum size of output
*/

	char fitsext[] = ".fits";	/* default extension */
	char *extn;		/* extension on input (or default) */

	int i_len;		/* length of input name */
	int is_len;		/* length of isuffix */
	int fits_len;		/* length of fitsext */
	int tr_len;		/* length of truncated input name */
	int i;			/* loop index */
	int dotlocn;		/* location of '.' in input name */

	i_len = strlen (input);
	is_len = strlen (isuffix);
	fits_len = strlen (fitsext);

	if (i_len > maxch) {
	    printf ("ERROR    (MkName) Input name is too long.\n");
	    return (2011);
	}

	extn = calloc (1 + is_len + fits_len, sizeof (char));
	if (extn == NULL)
	    return (OUT_OF_MEMORY);

	/* Find the extension (if any) on the input name. */
	dotlocn = 0;
	for (i = i_len-1;  i >= 0;  i--) {
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

	/* If the input name ends in the expected suffix, chop it off of
	   the output name before appending the output suffix and extension.
	*/
	tr_len = strlen (output);		/* length of truncated name */
	if (tr_len >= is_len) {
	    if (strcmp (output+tr_len-is_len, isuffix) == 0)
		output[tr_len-is_len] = '\0';
	}

	strcat (output, osuffix);
	strcat (output, extn);

	free (extn);
	return (0);
}
