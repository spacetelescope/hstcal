# include <stdio.h>

# include "hstio.h"
# include "wf3.h"
# include "hstcalerr.h"

/*  Load primary header from input image

   Warren Hack, 1998 July 15:	
   	Initial ACS Version.
	
	Warren Hack, 1998 Nov 17:
		Revised to support trailer files.

*/

int LoadHdr (char *input, Hdr *phdr) {

	extern int status;
	
	IODescPtr im;		/* descriptor for input image */
   	
	sprintf(MsgText, "Trying to open %s...",input);
	trlmessage (MsgText);

	/* Open input image in order to read its primary header. */
	im = openInputImage (input, "", 0);				

	if (hstio_err()) {
	    trlopenerr (input);
	    return (status = OPEN_FAILED);
	}
		
	initHdr (phdr);	

	/* get primary header */
	if (getHeader (im, phdr) )
	    status = HEADER_PROBLEM;	
	if (hstio_err() || status) {
	    trlreaderr (input);
	    closeImage (im);
	    freeHdr (phdr);
	    return (status = OPEN_FAILED);
	}
	
	closeImage (im);
    
	sprintf(MsgText, "Read in Primary header from %s...",input);
	trlmessage (MsgText);

	return (status);
}	
