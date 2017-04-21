/* ACS2d -- basic 2-D image reduction

   This file contains:
	ACS2d
*/

# include <stdio.h>
# include <time.h>
# include <string.h>

# include "hstio.h"

# include "acs.h"
# include "acsinfo.h"
# include "err.h"
# include "acscorr.h"		/* calibration switch names for calacs */

void Init2DTrl (char *, char *);

/* Do basic 2-D calibration.

   	Warren Hack, 1998 June 9:
		Initial ACS version.
	
	Warren Hack, 1998 Nov 17:
		Revised to support trailer file comments
*/

int ACS2d (char *input, char *output, CalSwitch *acs2d_sw, RefFileInfo *refnames,
		int printtime, int verbose) {

	extern int status;

	ACSInfo acs2d;	/* calibration switches, reference files, etc */
	int extver;

	Hdr phdr;		/* primary header for input image */

	int Do2D (ACSInfo *, int);
	int FileExists (char *);
	int Get2dFlags (ACSInfo *, Hdr *);
	int GetACSKeys (ACSInfo *, Hdr *);
	int GetLinTab (ACSInfo *);
	void Sanity2d (ACSInfo *);
	void TimeStamp (char *, char *);
	int LoadHdr (char *, Hdr *);
	void ACSInit (ACSInfo *);
	void PrBegin (char *label);
	void PrEnd (char *label);
	void PrFileName (char *label, char *filename);
	void PrHdrInfo (char *aperture, char *filter1, char *filter2, char *detector);
	void PrGrpBegin (char *label, int n);
	void PrGrpEnd (char *label, int n);
	void Init2DTrl (char *, char *);
	
/* ----------------------- Start Code --------------------------------*/

	/* Determine the names of the trailer files based on the input
		and output file names, then initialize the trailer file buffer
		with those names.
	*/
	Init2DTrl (input, output);
	/* If we had a problem initializing the trailer files, quit... */
	if (status != ACS_OK) 
		return (status);
	
	PrBegin ("ACS2D");

	if (printtime)
	    TimeStamp ("ACS2D started", "");

	/* Initialize structure containing calacs information. */
	ACSInit (&acs2d);

	/* Copy command-line arguments into acs2d. */
	strcpy (acs2d.input, input);
	strcpy (acs2d.output, output);
	acs2d.dqicorr  = acs2d_sw->dqicorr;
	acs2d.glincorr = acs2d_sw->glincorr;
	acs2d.lflgcorr = acs2d_sw->lflgcorr;
	acs2d.darkcorr = acs2d_sw->darkcorr;
  acs2d.flashcorr = acs2d_sw->flashcorr;
	acs2d.flatcorr = acs2d_sw->flatcorr;
	acs2d.shadcorr = acs2d_sw->shadcorr;
	acs2d.photcorr = acs2d_sw->photcorr;
	acs2d.expscorr = acs2d_sw->expscorr;
  acs2d.pctecorr = acs2d_sw->pctecorr;
	acs2d.noisecorr = PERFORM;
	acs2d.printtime = printtime;
	acs2d.verbose = verbose;
    
	acs2d.refnames = refnames;

	PrFileName ("input", acs2d.input);
	PrFileName ("output", acs2d.output);
	
	/* Check whether the output file already exists. */
	if (FileExists (acs2d.output))
	    return (status);

	/* Open input image in order to read its primary header. */
	if (LoadHdr (acs2d.input, &phdr) )
		return (status);

	/* Get keyword values from primary header using same function 
		used in ACSCCD
	*/
	if (GetACSKeys (&acs2d, &phdr)) {
		freeHdr (&phdr);
	    return (status);
	}

	/* Print information about this image. */
	PrHdrInfo (acs2d.aperture, acs2d.filter1, acs2d.filter2, acs2d.det);

	/* Reset switches that are inappropriate for current detector. */
	Sanity2d (&acs2d);

	/* Get reference file names from input image header.  Pedigree is
	   checked, and the calibration switch (an internal flag, not the
	   header keyword) will be reset from PERFORM to DUMMY if the
	   reference file has pedigree = "DUMMY".  Switches that are
	   currently set to PERFORM will be reset to OMIT if the value
	   in the header is COMPLETE.
	*/
	if (Get2dFlags (&acs2d, &phdr)) {
		freeHdr (&phdr);
	    return (status);
	}

	freeHdr (&phdr);

	if (acs2d.glincorr == PERFORM || acs2d.lflgcorr == PERFORM) {
	    if (GetLinTab (&acs2d))
		return (status);
	}

	/* Do basic 2-D image reduction. */

	if (acs2d.printtime)
	    TimeStamp ("Begin processing", acs2d.rootname);
	
	if (verbose) {
		sprintf(MsgText,"Processing %d IMSETs... ",acs2d.nimsets);
		trlmessage(MsgText);
	}
	
	for (extver = 1;  extver <= acs2d.nimsets;  extver++) {
	    trlmessage ("\n");
	    PrGrpBegin ("imset", extver);
	    if (Do2D (&acs2d, extver))
		return (status);
	    PrGrpEnd ("imset", extver);
	}

	trlmessage ("\n");
	PrEnd ("ACS2D");

	if (acs2d.printtime)
	    TimeStamp ("ACS2D completed", acs2d.rootname);

	/* Write out temp trailer file to final file */
	WriteTrlFile ();

	return (status);
}



void Init2DTrl (char *input, char *output) {
	
	extern int status;
	int exist;

	char trl_in[ACS_LINE+1]; 	/* trailer filename for input */
	char trl_out[ACS_LINE+1]; 	/* output trailer filename */
	
	int MkOutName (char *, char **, char **, int, char *, int);
	int MkNewExtn (char *, char *);
	void WhichError (int);
	int TrlExists (char *);
	void SetTrlOverwriteMode (int);

	/* Input and output suffixes. */
	char *isuffix[] = {"_raw", "_blv_tmp", "_blc_tmp", "_crj_tmp", "_crc_tmp"};
	char *osuffix[] = {"_flt", "_flt",     "_flc",     "_crj",     "_crc"};
	char *trlsuffix[] = {"", "", "", "", ""};

	int nsuffix = 5;
	
	/* Initialize internal variables */
	trl_in[0] = '\0';
	trl_out[0] = '\0';
	exist = EXISTS_UNKNOWN;

	/* Start by stripping off suffix from input/output filenames */
	if (MkOutName (input, isuffix, trlsuffix, nsuffix, trl_in, ACS_LINE)) {
		WhichError (status);
		sprintf (MsgText, "Couldn't determine trailer filename for %s", input);
		trlmessage (MsgText);
	}
	if (MkOutName (output, osuffix, trlsuffix, nsuffix, trl_out, ACS_LINE)) {
		WhichError (status);
		sprintf (MsgText, "Couldn't create trailer filename for %s", output);
		trlmessage (MsgText);
	}
    
	/* Now, convert trailer filename extensions from '.fits' to '.trl' */
	if (MkNewExtn (trl_in, TRL_EXTN) ) {
		sprintf(MsgText, "Error with input trailer filename %s", trl_in);
		trlerror (MsgText);
		WhichError (status);
	}
	if (MkNewExtn (trl_out, TRL_EXTN) ) {
		sprintf(MsgText, "Error with output trailer filename %s", trl_out);
		trlerror (MsgText);
		WhichError (status);
	}

	/* If we are working with a RAW file, then see if a TRL file 
		needs to be overwritten after the generic conversion comments.
	*/
	if (strstr(input, isuffix[0]) != NULL) {
		/* Test whether the output file already exists */
		exist = TrlExists(trl_out);            
		if (exist == EXISTS_YES) {
			/* The output file exists, so we want to overwrite them with
				the new trailer comments.
			*/
			SetTrlOverwriteMode (YES);	
		}
	}

	/* Sets up temp trailer file for output and copies input
		trailer file into it.
	*/
	InitTrlFile (trl_in, trl_out);
	/* Check value of STATUS for errors in calling routine */
}
