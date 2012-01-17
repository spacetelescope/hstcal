/* WF32d -- basic 2-D image reduction

   This file contains:
	WF32d
*/

# include <stdio.h>
# include <time.h>
# include <string.h>

# include "hstio.h"

# include "wf3.h"
# include "wf3info.h"
# include "wf3err.h"
# include "wf3corr.h"		/* calibration switch names for calwf3 */


/* Do basic 2-D calibration.

   	Warren Hack, 1998 June 9:
		Initial ACS version.
	
	Warren Hack, 1998 Nov 17:
		Revised to support trailer file comments

	Howard Bushouse, 2000 Aug 28:
		Initial WFC3 version.

	Howard Bushouse, 2001 Nov 15:
		Removed use of filtcorr for new PHOTCORR methods
		(in accordance with CALACS updates).

	Howard Bushouse, 2002 June 17:
	      Removed use of statcorr - now done by default
	      (in accordance with CALACS updates).

	Howard Bushouse, 2003 Oct 24:
	      Added expscorr switch as command argument (CALACS changes).
*/

int WF32d (char *input, char *output, CCD_Switch *wf32d_sw,
	   RefFileInfo *refnames, int printtime, int verbose) {

	extern int status;

	WF3Info wf32d;	/* calibration switches, reference files, etc */
	int extver;

	Hdr phdr;		/* primary header for input image */

	int Do2D (WF3Info *, int);
	int FileExists (char *);
	int Get2dFlags (WF3Info *, Hdr *);
	int GetKeys (WF3Info *, Hdr *);
	void Sanity2d (WF3Info *);
	void TimeStamp (char *, char *);
	int LoadHdr (char *, Hdr *);
	void WF3Init (WF3Info *);
	void PrBegin (char *label);
	void PrEnd (char *label);
	void PrFileName (char *label, char *filename);
	void PrHdrInfo (char *aperture, char *filter, char *detector);
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
	if (status != WF3_OK) 
	    return (status);
	
	PrBegin ("WF32D");

	if (printtime)
	    TimeStamp ("WF32D started", "");

	/* Initialize structure containing calwf3 information. */
	WF3Init (&wf32d);

	/* Copy command-line arguments into wf32d. */
	strcpy (wf32d.input, input);
	strcpy (wf32d.output, output);
	wf32d.dqicorr  = wf32d_sw->dqicorr;
	wf32d.darkcorr = wf32d_sw->darkcorr;
	wf32d.flatcorr = wf32d_sw->flatcorr;
	wf32d.shadcorr = wf32d_sw->shadcorr;
	wf32d.photcorr = wf32d_sw->photcorr;
	wf32d.expscorr = wf32d_sw->expscorr;
	wf32d.noiscorr = PERFORM;
	wf32d.printtime = printtime;
	wf32d.verbose = verbose;

	wf32d.refnames = refnames;

	PrFileName ("input", wf32d.input);
	PrFileName ("output", wf32d.output);
	
	/* Check whether the output file already exists. */
	if (FileExists (wf32d.output))
	    return (status);

	/* Open input image in order to read its primary header. */
	if (LoadHdr (wf32d.input, &phdr) )
		return (status);

	/* Get keyword values from primary header using same function 
		used in WF3CCD
	*/
	if (GetKeys (&wf32d, &phdr)) {
	    freeHdr (&phdr);
	    return (status);
	}

	/* Print information about this image. */
	PrHdrInfo (wf32d.aperture, wf32d.filter, wf32d.det);

	/* Reset switches that are inappropriate for current detector. */
	Sanity2d (&wf32d);

	/* Get reference file names from input image header.  Pedigree is
	   checked, and the calibration switch (an internal flag, not the
	   header keyword) will be reset from PERFORM to DUMMY if the
	   reference file has pedigree = "DUMMY".  Switches that are
	   currently set to PERFORM will be reset to OMIT if the value
	   in the header is COMPLETE.
	*/
	if (Get2dFlags (&wf32d, &phdr)) {
	    freeHdr (&phdr);
	    return (status);
	}

	freeHdr (&phdr);

	/* Do basic 2-D image reduction. */

	if (wf32d.printtime)
	    TimeStamp ("Begin processing", wf32d.rootname);
	
	if (verbose) {
		sprintf(MsgText,"Processing %d IMSETs... ",wf32d.nimsets);
		trlmessage(MsgText);
	}
	
	for (extver = 1;  extver <= wf32d.nimsets;  extver++) {
	    trlmessage ("\n");
	    PrGrpBegin ("imset", extver);
	    if (Do2D (&wf32d, extver))
		return (status);
	    PrGrpEnd ("imset", extver);
	}

	trlmessage ("\n");
	PrEnd ("WF32D");

	if (wf32d.printtime)
	    TimeStamp ("WF32D completed", wf32d.rootname);

	/* Write out temp trailer file to final file */
	WriteTrlFile ();

	return (status);
}


void Init2DTrl (char *input, char *output) {
	
	extern int status;
	int exist;

	char trl_in[SZ_LINE+1]; 	/* trailer filename for input */
	char trl_out[SZ_LINE+1]; 	/* output trailer filename */
	
	int MkOutName (char *, char **, char **, int, char *, int);
	int MkNewExtn (char *, char *);
	void WhichError (int);
	int TrlExists (char *);
	void SetTrlOverwriteMode (int);

	/* Input and output suffixes. */
	char *isuffix[] = {"_raw", "_blv_tmp", "_crj_tmp"};
	char *osuffix[] = {"_flt", "_flt",     "_crj"};
	char *trlsuffix[] = {"", "", ""};

	int nsuffix = 3;
	
	/* Initialize internal variables */
	trl_in[0] = '\0';
	trl_out[0] = '\0';
	exist = EXISTS_UNKNOWN;

	/* Start by stripping off suffix from input/output filenames */
	if (MkOutName (input, isuffix, trlsuffix, nsuffix, trl_in, SZ_LINE)) {
	    WhichError (status);
	    sprintf (MsgText, "Couldn't determine trailer filename for %s",
		     input);
	    trlmessage (MsgText);
	}
	if (MkOutName (output, osuffix, trlsuffix, nsuffix, trl_out, SZ_LINE)) {
	    WhichError (status);
	    sprintf (MsgText, "Couldn't create trailer filename for %s",
		     output);
	    trlmessage (MsgText);
	}

	/* Now, convert trailer filename extensions from '.fits' to '.trl' */
	if (MkNewExtn (trl_in, TRL_EXTN) ) {
	    sprintf (MsgText, "Error with input trailer filename %s", trl_in);
	    trlerror (MsgText);
	    WhichError (status);
	}
	if (MkNewExtn (trl_out, TRL_EXTN) ) {
	    sprintf (MsgText, "Error with output trailer filename %s", trl_out);
	    trlerror (MsgText);
	    WhichError (status);
	}

	/* If we are working with a RAW file, then see if a TRL file 
	   needs to be overwritten after the generic conversion comments.  */
	if (strstr(input, isuffix[0]) != NULL) {
	    /* Test whether the output file already exists */
	    exist = TrlExists(trl_out);            
	    if (exist == EXISTS_YES) {
		/* The output file exists, so we want to overwrite them with
		** the new trailer comments.  */
		SetTrlOverwriteMode (YES);	
	    }
	}

	/* Sets up temp trailer file for output and copies input
		trailer file into it.
	*/
	InitTrlFile (trl_in, trl_out);
	/* Check value of STATUS for errors in calling routine */
}
