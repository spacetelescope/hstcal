# include "c_iraf.h"
# include "hstio.h"

# include "stis.h"
# include "stisdef.h"
# include "stisversion.h"	/* defines CAL_VER */

/* This routine updates the CAL_VER primary header keyword, or adds it
   to the header if it's not already present.  This keyword gives the
   current version number and date that was used.
*/

void UCalVer (Hdr *phdr) {

/* argument:
Hdr *phdr       io: pointer to header; CAL_VER will be updated
*/

	Put_KeyS (phdr, "CAL_VER", STIS_CAL_VER, "version number and date");
}
