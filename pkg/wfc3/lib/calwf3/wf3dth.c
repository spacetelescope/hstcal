# include <stdio.h>
# include <stdlib.h>	/* calloc */
# include <string.h>

/* This section for EmptyGroup functions */
# if defined(VMS)
# include <stat.h>
# else
# include <sys/types.h>
# include <sys/stat.h>
# endif

#include "hstcal.h"
# include "hstio.h"
# include "ximio.h"

# include "wf3.h"
# include "calwf3.h"
# include "hstcalerr.h"
# include "trlbuf.h"

/* wf3dth -- produce DTH product (empty)

   Howard A. Bushouse, 2000 Aug 23:
	Initial version.
	[Adapted from acsdth.c by W. Hack]

   H.Bushouse, 2001 May 8:
	Updated to sync with changes to acsdth.c: Modified trailer file
	creation to include ALL input files; Also modified SPT file creation
	to include all extensions from input SPT files.

   H.Bushouse, 2002 June 19:
	Updated to track CALACS acsdth changes: Moved trailer file
	initialization outside this routine so that a trailer file can be
	created/appended for all ASN tables regardless of whether they contain
	actual dither members or not. This initialization is now handled by
	CalWf3Run in "calwf3.c". Also changed output product suffix from
	"dth" to "drz". Eliminated creation of dummy "_drz.fits" products.

   H.Bushouse, 2006 June 20:
	Updated to track CALACS acsdth changes: Fixed bug InitDthTrl in 
	reallocating memory for list of trailer file names (trl_in).

   H.Bushouse, 2008 June 11:
	Restored old acsdth code for creating dummy dth product, until we get
	MultiDrizzle implemented in WFC3 pipeline.

   H.Bushouse, 2009 Aug 13:
	Updated to always set NEXTEND=3 in dummy drz file header.

   H.Bushouse, 2009 Dec 22:
	Eliminated creation of dummy '_drz.fits' product, now that MultiDrizzle
	is being used to produce them. (calwf3 v2.0)

   H.Bushouse, 2010 Oct 20:
	Updated InitDthTrl to return if the input list is empty, in order to
	handle associations with missing members. (PR 66366)

  M Sosey, 2012 December 27:
      Updated to account for a memory leak on linux machines during BuildDth 
      when RPTCORR is off and a new spt is being constructed (#967)       

  M Sosey, 2015 May: 
      Updated for UVIS2
      
*/


int Wf3Dth ( char *in_list, char *output, int dthcorr, int printtime,
	    int verbose){
	
	extern int status;
	
	char mtype[SZ_CBUF+1];		/* role of exposure in Association */
	IRAFPointer tpin;

	void PrBegin (char *);
	void PrEnd (char *);
	void TimeStamp (char *, char *);	

	int  mkNewSpt (char *, char *, char *);
	
/* ----------------------- Start Code --------------------------------*/

	/* Start the task... */
	PrBegin ("WF3DTH");
		
	if (printtime)
	    TimeStamp ("WF3DTH started", "");

	sprintf (MsgText, "Astrodrizzle needs to be run in order to generate");
	trlmessage (MsgText);
	sprintf (MsgText, "a geometrically corrected, drizzle-combined product.");
	trlmessage (MsgText);

	/* Open the input file list */
	tpin = c_imtopen (in_list);

	/* create new SPT file for output product */
	sprintf(mtype,"PROD-DTH");
	if (mkNewSpt (in_list, mtype, output)) {
	    return(status);
	}

	c_imtclose (tpin);

	PrEnd ("WF3DTH");

	/* Write out temp trailer file to final file */
	WriteTrlFile ();

	return (status);	
}


void InitDthTrl (char *inlist, char *output) {
	
	extern int status;
	
	IRAFPointer tpin;
	int n, nfiles;

	char *trl_in;			/* trailer filename for input */
	int  trl_len;
	char trl_out[CHAR_LINE_LENGTH+1]; 	/* output trailer filename */
	char input[CHAR_FNAME_LENGTH+1];		/* Name of image in list */
	char out_name[CHAR_FNAME_LENGTH+1];
	
	char *isuffix[]={"_sfl", "_crj", "_flt", "_flc", "_crc", "_sfl"};
	char *osuffix[]={"_drz", "_drz", "_drz", "_drc", "_drc", "_drc"};
	char *trlsuffix[]={"", "", "", "", "", ""};
	int nsuffix = 6;
	
	int MkOutName (char *, char **, char **, int, char *, int);
	int MkNewExtn ( char *, char *);
	void WhichError (int);

	/* Allocate space for trailer file input list */
	trl_in = (char *) calloc((CHAR_LINE_LENGTH+1), sizeof(char));
	trl_len = CHAR_LINE_LENGTH+1;

	/* Make TRL filenames */
	trl_in[0]  = '\0';
	trl_out[0] = '\0';
	out_name[0] = '\0';
	tpin = c_imtopen(inlist);
	nfiles = c_imtlen(tpin);

	/* If the input list is null, then just return */
	if (nfiles == 0) {
	    c_imtclose (tpin);
	    free(trl_in);
        printf("\nNothing in DTH input list\n");
	    return;
	}

	for (n = 0; n < nfiles; ++n) {
	     c_imtgetim (tpin, input, CHAR_FNAME_LENGTH);
	     /* Start by stripping off suffix from input/output filenames */
	     if (MkOutName (input, isuffix, trlsuffix, nsuffix, out_name,
			    CHAR_LINE_LENGTH)) {
         
		    WhichError (status);
		    sprintf (MsgText, "Couldn't determine trailer filename for %s",input);
		    trlmessage (MsgText);
	     }
         
	     /* Now convert trailer filename extension from '.fits' to '.trl' */
	     if (MkNewExtn (out_name, TRL_EXTN)) {
		    sprintf (MsgText, "Error creating input trailer filename %s",out_name);
		    trlerror (MsgText);
		    WhichError (status);
	     }
         
         
	     if ( (strlen(out_name) + strlen(trl_in) + 1) >= trl_len) {
		    trl_len += CHAR_LINE_LENGTH;
		    trl_in = realloc (trl_in, trl_len);
	     }

	     /* Append each filename to create list of input trailer files */
	     strcat(trl_in, out_name);
         
	     /* But don't put a comma after the last filename */
	     if (n < (nfiles-1)) strcat (trl_in, ",");
	     /* Reset value for the output filename for the next image */
	     out_name[0] = '\0';
	}

    
	if (MkOutName (output, osuffix, trlsuffix, nsuffix, trl_out, CHAR_LINE_LENGTH)) {
	    WhichError (status);
	    sprintf (MsgText, "Couldn't create trailer filename for %s\n",output);
	    trlmessage (MsgText);
	}

	if (MkNewExtn (trl_out, TRL_EXTN) ) {
	    sprintf (MsgText, "Error creating output trailer filename %s\n",trl_out);
	    trlerror (MsgText);
	    WhichError (status);
	}

    
	/* Sets up temp trailer file for output and copies input
	** trailer file into it. */
	InitTrlFile (trl_in, trl_out);

	/* Deallocate memory */
	free(trl_in);
	c_imtclose (tpin);
}

