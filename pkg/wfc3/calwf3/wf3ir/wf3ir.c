/* wf3ir -- basic IR image reduction

   This file contains:
	WF3ir
*/

# include <stdio.h>
# include <stdlib.h>
# include <string.h>

# include "hstio.h"	/* defines HST I/O functions */

# include "wf3.h"
# include "wf3info.h"
# include "hstcalerr.h"     /* error values */
# include "wf3corr.h"    /* calibration switch names */
# include "trl.h"        /* trailer functions */

extern int status;

static void mkNames (char *, char *, char *);

/* Do basic IR calibration.
**
** Howard Bushouse, 2000 Oct 11:
**	Initial WFC3 version, created from modified calnica.c in
**	NICMOS pipeline.
** H.Bushouse	May 2001	Migrated from CALNIC to CALWF3 style
**				handling of reference data, use of
**				trailer file, etc.
** H.Bushouse	21 June 2002	Removed use of statcorr (always perform doStat).
** H.Bushouse	21 Feb 2007	Added use of new putCalDataSect routine to
**				write out calibrated images that have the
**				ref pixels trimmed off.
*/

int WF3ir (char *raw_file, char *flt_file, IR_Switch *ir_sw,
	   RefFileInfo *refnames, int printtime, int verbose,
	   float cr_thresh, float zs_thresh, int samp_rej) {

/* Arguments:
**	raw_file	i: input raw file name
**	flt_file	i: output flt file name
**	ir_sw		i: ir processing switches
**	refnames	i: reference file name list
**	printtime	i: print time stamps?
**	verbose		i: verbose messages?
**	cr_thresh	i: CR rejection threshold
**	zs_thresh	i: Zero-read signal threshold
**	samp_rej	i: Number of initial samples to reject
*/

	/* Local variables */
	Hdr phdr;			/* Primary header */
	MultiNicmosGroup  allinput;	/* Input image data */
	SingleNicmosGroup crimage;	/* CR reject output image */
	char ima_file[SZ_FNAME+1];	/* Intermediate MultiAccum file name */
	int sizex, sizey;		/* Output trimmed image sizes */

	WF3Info wf3;	/* calibration switches, reference files, etc. */

	/* Function definitions */
	int  DoIR (WF3Info *, MultiNicmosGroup *, SingleNicmosGroup *);
	int  FileExists (char *);
	int  GetIRFlags (WF3Info *, Hdr *);
	int  GetKeys (WF3Info *, Hdr *);
	void TimeStamp (char *, char *);
	void PrBegin (char *);
	void PrEnd   (char *);
	void PrFileName (char *, char *);
	void PrHdrInfo (char *, char *, char *);
	int  LoadHdr (char *, Hdr *);
	void InitIRTrl (char *, char *);
	void WF3Init   (WF3Info *);

	int  getRawData (WF3Info *, MultiNicmosGroup *);
	int  putMultiCalData (MultiNicmosGroup *, char *);
	int  putCalDataSect (SingleNicmosGroup *, char *, int, int, int, int);

/* ----------------------------------------------------------------------- */

	/* Determine the name of the trailer file based on the input and
	** output file names, then initialize the trailer file buffer
	** with those names. */
	InitIRTrl (raw_file, flt_file);

	/* If there was a problem initializing the trailer file, quit. */
	if (status != WF3_OK)
	    return (status);

	/* Print startup message */
	PrBegin ("WF3IR");

	if (printtime)
	    TimeStamp ("WF3IR started", "");

	/* Initialize structure containing calwf3 information. */
	WF3Init (&wf3);

	/* Build the complete input and output file names */
	mkNames (raw_file, ima_file, flt_file);

	/* Copy command-line arguments into wf3 structure. */
	strcpy (wf3.input,  raw_file);
	strcpy (wf3.output, flt_file);

	wf3.dqicorr   = ir_sw->dqicorr;
	wf3.blevcorr  = ir_sw->blevcorr;
	wf3.zsigcorr  = ir_sw->zsigcorr;
	wf3.zoffcorr  = ir_sw->zoffcorr;
	wf3.noiscorr  = ir_sw->noiscorr;
	wf3.darkcorr  = ir_sw->darkcorr;
	wf3.nlincorr  = ir_sw->nlincorr;
	wf3.flatcorr  = ir_sw->flatcorr;
	wf3.photcorr  = ir_sw->photcorr;
	wf3.unitcorr  = ir_sw->unitcorr;
	wf3.crcorr    = ir_sw->crcorr;
	wf3.printtime = printtime;
	wf3.verbose   = verbose;
	wf3.crthresh  = cr_thresh;
	wf3.zsthresh  = zs_thresh;
	wf3.samp_rej  = samp_rej;

	wf3.refnames  = refnames;

	PrFileName ("input", raw_file);
	PrFileName ("output", flt_file);

	/* Check whether the output files already exist */
	if (FileExists (flt_file))
	    return (status);
	if (FileExists (ima_file))
	    return (status);

	/* Open input image and read its primary header. */
	if (LoadHdr (wf3.input, &phdr))
	    return (status);

	/* Get keyword values from primary header */
	if (GetKeys (&wf3, &phdr)) {
	    freeHdr (&phdr);
	    return (status);
	}

	/* If we have UVIS data, do not even proceed here... */
	if (wf3.detector == CCD_DETECTOR) {
	    trlerror ("Can NOT process UVIS data with WF3IR...");
	    freeHdr (&phdr);
	    return (status = ERROR_RETURN);
	}

        /* Print information about this image */
	PrHdrInfo (wf3.aperture, wf3.filter, wf3.det);

	/* Get reference file names from input image header. */
	if (GetIRFlags (&wf3, &phdr)) {
	    freeHdr (&phdr);
	    return (status);
	}
	
	/* Free the primary header. */
	freeHdr (&phdr);

	/* Load all input data groups */
	sprintf (MsgText, "Reading data from %s ...", wf3.input);
	trlmessage (MsgText);
	if (getRawData (&wf3, &allinput))
	    return (status);

	/* Do IR image reduction. */

	if (wf3.printtime)
	    TimeStamp ("Begin processing", wf3.rootname);

	if (verbose) {
	    sprintf (MsgText, "Processing %d imsets... ",wf3.ngroups);
	    trlmessage (MsgText);
	}

	/* Apply the calibration steps */
	if (DoIR (&wf3, &allinput, &crimage)) {
	    freeMultiNicmosGroup (&allinput);
	    return (status);
	}

	/* Write the intermediate calibrated data to the IMA file */
	sprintf (MsgText, "Writing calibrated readouts to %s", ima_file);
	trlmessage (MsgText);
	if (putMultiCalData (&allinput, ima_file)) {
	    freeMultiNicmosGroup (&allinput);
	    freeSingleNicmosGroup (&crimage);
	    return (status);
	}

	/* Compute size of trimmed output image */
	sizex = allinput.group[0].sci.data.nx - (wf3.trimx[0]+wf3.trimx[1]);
	sizey = allinput.group[0].sci.data.ny - (wf3.trimy[0]+wf3.trimy[1]);

	/* Write the final trimmed image data to FLT file */
	sprintf (MsgText, "Writing final image to %s", flt_file);
	trlmessage (MsgText);
	sprintf (MsgText, " with trimx = %d,%d, trimy = %d,%d", wf3.trimx[0],
		 wf3.trimx[1], wf3.trimy[0], wf3.trimy[1]);
	trlmessage (MsgText);

	if (wf3.crcorr == PERFORM) {
	    if (putCalDataSect (&crimage, flt_file, wf3.trimx[0],
				wf3.trimy[0], sizex, sizey)) {
		freeSingleNicmosGroup (&crimage);
		return (status);
	    }
	    freeSingleNicmosGroup (&crimage);
	} else {
	    if (putCalDataSect (&(allinput.group[0]), flt_file, wf3.trimx[0],
				wf3.trimy[0], sizex, sizey)) {
		freeMultiNicmosGroup (&allinput);
		return (status);
	    }
	}
	if (wf3.printtime)
	    TimeStamp ("Output written to disk", wf3.rootname);

	/* Free memory for the input data groups */
	freeMultiNicmosGroup (&allinput);

	trlmessage ("\n");
	PrEnd ("WF3IR");

	if (wf3.printtime)
	    TimeStamp ("WF3IR completed", raw_file);

	/* Write out temp trailer file to final file */
	WriteTrlFile ();

	/* Successful return */
	return (status = 0);
}

/* MKNAMES: Build the complete input and output file names from
** whatever the user supplied as arguments */

static void mkNames (char *rawfile, char *imafile, char *fltfile) {

/* Arguments:
**	rawfile	io: input raw file name
**	imafile	 o: output intermediate file name (used for MULTIACCUM only)
**	fltfile	io: output flat file name
*/

	/* Local variables */
	Bool ext1;			/* Was a "fits" extension specified? */
	Bool ext2;			/* Was a "fit"  extension specified? */
	Bool israw;			/* Was a raw file suffix specified? */
	char inroot[SZ_FNAME+1];	/* root name of input file */
	char outroot[SZ_FNAME+1];	/* root name of output file */
	int  rlen;			/* length of root name */
	char fitsext1[] = ".fits";	/* FITS file extension (default) */
	char fitsext2[] = ".fit";	/* FITS file extension (alternate) */
	char rawsuf[] = "_raw";		/* raw file suffix */
	char imasuf[] = "_ima";		/* intermediate file suffix */
	char fltsuf[] = "_flt";		/* flt file suffix */

	/* Initialize */
	ext1   = False;
	ext2   = False;
	israw = False;
	strcpy (inroot,  rawfile);
	strcpy (outroot, fltfile);

	/* Search for "fits" extension on input file name */
	rlen = strlen(inroot);
	if (strncmp (inroot+rlen-5, fitsext1, 5) == 0) {
	    ext1 = True;
	    inroot[rlen-5] = '\0';
	    rlen = strlen (inroot);

	/* Search for "fit" extension on input file name */
	} else if (strncmp (inroot+rlen-4, fitsext2, 4) == 0) {
	    ext2 = True;
	    inroot[rlen-4] = '\0';
	    rlen = strlen (inroot);
	}

	/* Search for "raw" suffix on input file name */
	if (strncmp (inroot+rlen-4, rawsuf, 4) == 0) {
	    israw = True;
	    inroot[rlen-4] = '\0';
	    rlen = strlen (inroot);
	}
	
	/* Add any missing pieces to the input raw file name */
	if (israw && !ext1 && !ext2) {
	    strcat (rawfile, fitsext1);
	    ext1 = True;
	} else if (!israw && !ext1 && !ext2) {
	    strcat (rawfile, rawsuf);
	    strcat (rawfile, fitsext1);
	    ext1 = True;
	}

	/* If no output name was specified, build the flt and ima file
	** file names from the input root name */
	if (fltfile[0] == '\0') {

	    strcpy (fltfile, inroot);
	    strcpy (imafile, inroot);

	} else {

	    /* Search for "fits" extension on output file name */
	    rlen = strlen (outroot);
	    if (strncmp (outroot+rlen-5, fitsext1, 5) == 0) {
		ext1 = True;
		ext2 = False;
		outroot[rlen-5] = '\0';
		rlen = strlen (outroot);

	    /* Search for "fit" extension on output file name */
	    } else if (strncmp (outroot+rlen-4, fitsext2, 4) == 0) {
		ext1 = False;
		ext2 = True;
		outroot[rlen-4] = '\0';
		rlen = strlen (outroot);
	    }

	    /* Search for "flt" suffix on output file name */
	    if (strncmp (outroot+rlen-4, fltsuf, 4) == 0) {
		outroot[rlen-4] = '\0';
		rlen = strlen (outroot);
	    }

	    strcpy (fltfile, outroot);
	    strcpy (imafile, outroot);
	}

	/* Build the flt and ima file names from their roots */
	strcat (imafile, imasuf);
	strcat (fltfile, fltsuf);

	if (ext1) {
	    strcat (fltfile, fitsext1);
	    strcat (imafile, fitsext1);
	} else if (ext2) {
	    strcat (fltfile, fitsext2);
	    strcat (imafile, fitsext2);
	}

}

void InitIRTrl (char *input, char *output) {

        extern int status;

        char trl_in[SZ_LINE+1];         /* trailer filename for input */
        char trl_out[SZ_LINE+1];        /* output trailer filename */
        int exist;

        char isuffix[] = "_raw";
        char osuffix[] = "_flt";

        int MkName (char *, char *, char *, char *, char *, int);
        void WhichError (int);
        int TrlExists (char *);
        void SetTrlOverwriteMode (int);

        /* Initialize internal variables */
        trl_in[0] = '\0';
        trl_out[0] = '\0';
        exist = EXISTS_UNKNOWN;

        /* Start by stripping off suffix from input/output filenames */
        if (MkName (input, isuffix, "", TRL_EXTN, trl_in, SZ_LINE)) {
            WhichError (status);
            sprintf (MsgText, "Couldn't determine trailer filename for %s",
                     input);
            trlmessage (MsgText);
        }

        if (MkName (output, osuffix, "", TRL_EXTN, trl_out, SZ_LINE)) {
            WhichError (status);
            sprintf (MsgText, "Couldn't create trailer filename for %s",
                     output);
            trlmessage (MsgText);
        }

        /* Test whether the output file already exists */
        exist = TrlExists(trl_out);
        if (exist == EXISTS_YES) {
            /* The output file exists, so we want to overwrite them with
            ** the new trailer comments.  */
            SetTrlOverwriteMode (YES);
        }

        /* Sets up temp trailer file for output and copies input
        ** trailer file into it.  */
        InitTrlFile (trl_in, trl_out);

}

