# include "hstio.h"
# include "acsversion.h"	/* defines ACS_CAL_VER */

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

	PutKeyStr (phdr, "CAL_VER", ACS_CAL_VER, "version number and date");
}


/* This adds the version number of the CTE algorithm to the image header.
 * MRD 29 Sept 2011 */

void UCteVer(Hdr *phdr, char * cteName, char * cteVersion) {
  int PutKeyStr (Hdr *, char *, char *, char *);
  
  PutKeyStr(phdr, "CTE_NAME", cteName, "name of CTE algorithm");
  PutKeyStr(phdr, "CTE_VER", cteVersion, "version of CTE algorithm");
}
