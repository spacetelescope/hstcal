# include <stdio.h>
# include <string.h>
# include "hstio.h"
# include "wf3.h"
# include "err.h"

/* This routine determines the chip associated with the
	EXTVER.

   Warren Hack, 18 June 1998

    Revised to loop over all reference file extensions to search
    for the appropriate CCDCHIP. This allows subarrays to be supported
    more transparently, by not relying on them having the same number
    of extensions as the reference file.
     WJH     12 April 2000 
*/

int DetCCDChip (char *fname, int chip, int nimsets, int *extver) {

/*
char *fname		i: name of file
int chip		i: CHIP ID to be found
int nimsets;		i: number of IMSETS in image
int *extver		o: extension (IMSET) from file corresponding
			   to input image chip ID
*/

	Hdr scihdr, prihdr;
	IODescPtr ip;
	
	extern int status;
	int ccdchip;		/* CHIP id from reference header */
	int n, foundit;
	int nextend, ngrps;
	
	int GetKeyInt (Hdr *, char *, int, int, int *);
	
	initHdr (&scihdr);
	initHdr (&prihdr);
	*extver = 0;
	ip = NULL;
	foundit = NO;
	ip = openInputImage (fname, "", 0);
	getHeader (ip, &prihdr);
	closeImage (ip);

	/* Find out how many extensions there are in this file. */
	if (GetKeyInt(&prihdr, "NEXTEND", USE_DEFAULT, EXT_PER_GROUP, &nextend))
	    return (status);
	ngrps = nextend / EXT_PER_GROUP;

	/* Loop over all the extensions in the reference file
	** to search for the extension which corresponds to the desired 
	** CCDCHIP id of the exposure. */

	for (n=1; n <= ngrps; n++) {
	     ip = openInputImage (fname, "SCI", n);

	     getHeader (ip, &scihdr);

	     if (ip != NULL)
		 closeImage (ip);

	     /* Get CCD-specific parameters. */
	     if (GetKeyInt (&scihdr, "CCDCHIP", USE_DEFAULT, 1, &ccdchip))
		 return (status);

	     if (ccdchip == chip) {
		 /* Found it! */		
		 *extver = n;
        	 foundit = YES;
		 break;
	     } else {
		 /* Check next extension for CHIP id */
		 ccdchip = 0;
	     }
	}
        
	freeHdr(&scihdr);
	freeHdr(&prihdr);

	if (foundit == NO){
	    sprintf (MsgText, "No Reference Data found for chip %d",chip);
	    trlerror (MsgText);
	    return (status = NO_CHIP_FOUND);
	}

	return (status);
}

