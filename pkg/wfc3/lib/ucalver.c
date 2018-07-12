# include "hstio.h"
# include "wf3version.h"	/* defines WF3_CAL_VER */

/* This routine updates the CAL_VER primary header keyword, or adds it
   to the header if it's not already present.  This keyword gives the
   current version number and date that was used.
   This routine is used to limit the number of files which need to 
    be recompiled when the version number is changed.  WJH 27 July 1999 
*/

void UCalVer (Hdr *phdr) {

/* argument:
Hdr *phdr       io: pointer to header; CAL_VER will be updated
*/

	int PutKeyStr (Hdr *, char *, char *, char *);

	PutKeyStr (phdr, "CAL_VER", WF3_CAL_VER, "version number and date");
}
