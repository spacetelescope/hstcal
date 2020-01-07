# include <stdio.h>
# include <string.h>
# include "stis.h"
# include "stisaper.h"	/* pseudo_ap, for pseudo-aperture suffixes */

/* This routine checks whether propaper ends in one of the "pseudo-aperture"
   suffixes, and if it does, it appends that string to aperture (unless
   aperture already includes that suffix).

   Phil Hodge, 2003 Mar 21  Function created.
*/

void pseudoap (char *propaper, char *aperture, int verbose) {

/* arguments:
char *propaper   i: proposed aperture name, which may include "E1", etc.
char *aperture  io: aperture name; may be modified in-place
int verbose      i: true --> print message if aperture is modified
*/

	int i;

	i = 0;
	while (strcmp (pseudo_ap[i], "NULL") != 0) {
	    if (strcmptail (propaper, pseudo_ap[i]) == 0 &&
		strcmptail (aperture, pseudo_ap[i]) != 0) {
		strcat (aperture, pseudo_ap[i]);
		if (verbose) {
		    printf ("INFO     '%s' has been appended to aperture\n",
			pseudo_ap[i]);
		}
		break;
	    }
	    i++;
	}
}
