# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include "stis.h"
# include "calstis4.h"

/* These two bracket 6, which is the length of the medium echelle slit. */
# define MEDIUM_SHORT  5.	/* boundary between short and medium */
# define MEDIUM_LONG   7.	/* boundary between medium and long */

/* This routine reads the slit length from the aperture name and assigns
   an integer code for the type of slit.  The length is read from the
   beginning of the aperture name if dispaxis is one, but it is read from
   the portion of the name after the 'X' if dispaxis is two.
*/

void WhichSlit (char *aperture, int dispaxis, int *slittype) {

/* arguments:
char *aperture   i: aperture name, e.g. 52X0.1
int dispaxis     i: dispersion axis, 1 or 2
int *slittype    o: code for slit type:  long, medium, short
*/

	double length;	/* slit length in arcsec, from aperture name */
	char *startchar;

	if (dispaxis == 1) {
	    length = atof (aperture);

	} else {
	    /* Skip past the 'X' in aperture name before reading length. */
	    startchar = strchr (aperture, 'X');
	    if (startchar == NULL)
		startchar = strchr (aperture, 'x');
	    if (startchar == NULL) {
		*slittype = UNKNOWN_SLIT;
		return;
	    }
	    startchar++;		/* increment past the 'X' */
	    length = atof (startchar);
	}

	if (length == 0.)	/* atof conversion error -- bad aperture name */
	    *slittype = UNKNOWN_SLIT;

	else if (length < MEDIUM_SHORT)
	    *slittype = SHORT_ECHELLE_SLIT;

	else if (length < MEDIUM_LONG)
	    *slittype = MEDIUM_ECHELLE_SLIT;

	else
	    *slittype = LONG_SLIT;
}
