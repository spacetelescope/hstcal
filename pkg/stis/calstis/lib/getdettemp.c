# include <stdio.h>
# include <string.h>

# include "c_iraf.h"
# include "hstio.h"

# include "stis.h"
# include "hstcalerr.h"
# include "stisdef.h"
# include "stistemperature.h"	/* defines BEGIN_SIDE2 */

/* This routine gets the detector temperature from a header keyword.
   The keyword name depends on the detector.  If the detector is
   invalid or the keyword isn't found, or for side-1 CCD data, the
   temperature will be set to -1.  (For the CCD, the keyword gives the
   housing temperature, not the detector temperature.)

   Phil Hodge, 2004 Dec 27:
        Initial version.
*/

int GetDetTemp (Hdr *hdr, int detector, double *temperature) {

/* arguments:
Hdr *hdr             i: header of current extension
int detector         i: integer code for detector
double *temperature  o: temperature (degrees C), or -1 if not found
*/

	/* name of keyword (depending on detector) for temperature */
	char keyword[STIS_CBUF+1];
	int use_def = 1;	/* use default if missing keyword */
	int status;
	FitsKw key;		/* for testing whether the keyword exists */

	*temperature = -1.;		/* default */
	if (detector == CCD_DETECTOR) {
	    double expstart;	/* MJD of exposure start time */
	    if ((status = Get_KeyD (hdr, "EXPSTART", use_def, -1., &expstart)))
		return (status);
	    if (expstart < BEGIN_SIDE2)
		return (0);
	    strcpy (keyword, "OCCDHTAV");
	} else if (detector == FUV_MAMA_DETECTOR) {
	    strcpy (keyword, "OM1CAT");
	} else if (detector == NUV_MAMA_DETECTOR) {
	    strcpy (keyword, "OM2CAT");
	}
	key = findKw (hdr, keyword);
	if (key == NotFound) {
	    printf ("Warning  keyword %s not found;" \
		" no temperature correction applied to sensitivity.\n",
		keyword);
	} else {
	    if ((status = Get_KeyD (hdr, keyword, use_def, -1., temperature)))
		return (status);
	    if (*temperature == -1.) {
		printf (
	"Warning  temperature (keyword %s) has default value;\n", keyword);
		printf (
	"         no temperature correction applied to sensitivity.\n");
	    }
	}

	return (0);
}
