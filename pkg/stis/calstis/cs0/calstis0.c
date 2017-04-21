# include <stdio.h>
# include <stdlib.h>	/* calloc */
# include <string.h>

# include "hstio.h"

# include "stis.h"
# include "calstis0.h"
# include "hstcalerr.h"

/* calstis0 -- integrated calstis processing

   Phil Hodge, 1997 Oct 30:
	Add interfaces for calling CalStis2, CalStis4, CalStis6, CalStis7;
	these are used in order to set up the calling sequences for the
	default parameters.  Add ResetCs1Sw to reset dqicorr & blevcorr.

   Phil Hodge, 1997 Nov 13:
	Replace sts.sx1file with sts.x1dfile in the call to CalStis6_0
	for the case of non-CR-split data.  ("sx1file" was just a typo)

   Phil Hodge, 1997 Dec 11:
	Change calling sequence of StisInit;
	initialize and free reference file lists;
	pass sciref or wavref to CalStis1, etc.

   Phil Hodge, 1998 Mar 18:
	Remove RefFileInfo sciref from calling sequence for CalStis12.

   Phil Hodge, 1998 Apr 10:
	Remove any directory prefix from input file name before appending
	to outroot.  This is done by AppendRawfile.

   Phil Hodge, 1998 Apr 13:
	Change the calling sequence for CalStis6:  remove debugFile,
	replace debug by extrloc and set it to 1.

   Phil Hodge, 1998 June 3:
	Include BIASCORR among the steps that should be done prior to
	cosmic ray rejection.  This affects SetCs1Sw and ResetCs1Sw.

   Phil Hodge, 1998 July 13:
	Add arguments to the calling sequence of CalStis6.

   Phil Hodge, 1998 Sept 24:
	Add three char * arguments to the calling sequence of CalStis6.

   Phil Hodge, 1998 Oct 5:
	Change status values 1021, 1023, 1025 to GENERIC_ERROR_CODE.

   Phil Hodge, 1998 Dec 11:
	Add dbgfile to the calling sequence of CalStis4.
	Add backval, backerr, variance to the calling sequence of CalStis6.

   Phil Hodge, 1999 May 25:
	Add fflux to the calling sequence of CalStis6.

   Phil Hodge, 1999 Dec 13:
	Add psclip and sclip to the calling sequence of CalStis6.

   Phil Hodge, 2000 Jan 19:
	If echelle data, pass fwv_tmp instead of w2d_tmp to calstis4 and
	calstis12.  If basic 2-D is not to be done for the wavecal,
	copy the wavecal file to fwv_tmp, so that calstis4 can modify it
	without affecting the input wavecal itself.  Add CalStis12_0,
	an interface for calstis12, include the cs12 include file for
	the interpolation option macros, and change the default from
	nearest neighbor to linear interpolation.

   Phil Hodge, 2000 Feb 8:
	Add lfilter to the calling sequence of CalStis6.

   Phil Hodge, 2000 Apr 5:
	Add center_target to the calling sequence of CalStis7.

   Phil Hodge, 2000 July 26:
	Add idtfile to the calling sequence of CalStis6.
	In the call to CalStis6, variance and fflux were out of order
	with respect to calstis6.c in ../cs6/; they have been interchanged
	so that the call and the function are now consistent.

   Phil Hodge, 2000 Oct 5:
	Add sc2dcorr to the calling sequence of CalStis6.

   Phil Hodge, 2001 Jan 4:
	Add subscale to the calling sequence of CalStis6.

   Phil Hodge, 2001 Feb 22:
	Pass fwv_tmp to calstis4 for prism data.

   Phil Hodge, 2001 Aug 27:
	If basic 2-D is not being done but cosmic-ray rejection is,
	set fltexists to false.

   Phil Hodge, 2002 Mar 1:
	Remove the test that prevents calling CalStis8 on the flt file to
	produce sfl if x2dcorr or geocorr is perform (i.e. if CalStis8 has
	already been called to produce the sx2 file).

   Ivo Busko, 2002 April 16:
	Add argument to the calling sequence of CalStis6.

   Ivo Busko, 2002 Aug 01:
	Initialize darkscale_string.

   Ivo Busko, 2002 Oct 29:
	Turn on background smoothing switch.

   Ivo Busko, 2002 Dec 16:
	Add ref. star position parameter to CalStis6 calling sequence.

   Paul Barrett, 2003 Jun 23:
        Changed the calling sequence for CalStis6.  Added PERFORM constant
        ctecorr.

   Phil Hodge, 2004 July 23:
	Add ctecorr variable, set from primary header keyword (default
	value is PERFORM).
	Add slit_angle = 0. to the calling sequence for CalStis4.

   Phil Hodge, 2005 April 21:
	In CalStis6_0, replace variable stpos with xoffset, and initialize
	the value to 0.  stpos was a position, and xoffset is an offset.

   Phil Hodge, 2006 Jan 11:
	Include wecfcorr = SKIPPED in the calling sequences for CalStis4 and
	CalStis6.

   Phil Hodge, 2006 Feb 7:
	Include err_algorithm in the calling sequence for CalStis7.

   Phil Hodge, 2007 May 9:
	Ignore status = NOTHING_TO_DO from CalStis8.
*/

static int CalStis2_0 (char *, char *, int, int);
static int CalStis4_0 (char *, RefFileInfo *, int, int);
static int CalStis6_0 (char *, char *, CalSwitch *, int, int);
static int CalStis7_0 (char *, char *, CalSwitch *, RefFileInfo *, int, int);
static int CalStis12_0 (char *, char *, int, int);
static void AppendRawfile (StisInfo *, char *, char *);
static int CopyFFile (char *, char *);
static void SetCs1Sw (CalSwitch *, CalSwitch *,
		cs1_switch *, cs1_switch *, cs1_switch *);
static void ResetCs1Sw (cs1_switch *, cs1_switch *);

int CalStis0 (char *rawfile, char *wavfile, char *outroot,
		int printtime, int save_tmp, int verbose) {

/* arguments:
char *rawfile    i: name of the FITS file to be calibrated
char *wavfile    i: name of the associated wavecal, or null
char *outroot    i: root for output file names, or null
int printtime    i: true --> print time stamps at intermediate steps
int save_tmp     i: true --> save temporary files
int verbose      i: true --> print info in individual csN
*/

	int status;

	StisInfo sts;		/* calibration switches, etc */
	RefFileInfo sciref, wavref;	/* ref file keywords and names */
	CalSwitch sci_sw;	/* all cal switches for science file */
	CalSwitch wav_sw;	/* all cal switches for wavecal */
	int save_crj;		/* don't delete crj_tmp file? */
	int save_fwv;		/* don't delete fwv_tmp file? */
	int fltexists;		/* fltfile will be created? (yes, normally) */
	cs1_switch cs1a_sci_sw;	/* calstis1a switches for science file */
	cs1_switch cs1_sci_sw;	/* calstis1 switches for science file */
	cs1_switch cs1_wav_sw;	/* calstis1 switches for wavecal */

	/* This will be either the fwv_tmp or x2d_tmp file name, depending
	   on whether the wavecal observation was echelle or first order.
	*/
	char *wavecal_file;

	int StisInit (StisInfo *, CalSwitch *, CalSwitch *,
		RefFileInfo *, RefFileInfo *);

	PrBegin (0);	/* *** CALSTIS-0 -- Version ... *** */

	if (printtime)
	    TimeStamp ("CALSTIS started", "");

	strcpy (sts.rawfile, rawfile);
	strcpy (sts.wavfile, wavfile);
	AppendRawfile (&sts, rawfile, outroot);	/* append rawfile to outroot */

	save_crj = save_tmp;		/* initial values; */
	save_fwv = save_tmp;		/*   we may need to save these */

	/* Print image name.  We may not have the wavecal name yet. */
	printf ("\n");
	PrFileName ("input", sts.rawfile);
	if (outroot[0] != '\0')
	    PrFileName ("outroot", sts.outroot);

	/* Normally, the fltfile will be either be created or the rawfile
	   name will be copied to fltfile so that subsequent processing will
	   use rawfile as input.  The only case where fltfile really doesn't
	   exist is when crcorr = perform and expscorr = omit.
	*/
	fltexists = 1;

	/* Initialize the lists of reference file keywords and names. */
	InitRefFile (&sciref);
	InitRefFile (&wavref);

	/* Initialize variables, and get info (calibration switches and
	   reference file names) from primary headers.
	   If wavfile was not specified on the command line, then the
	   default name will be gotten from the science file header;
	   the wavecal name will be printed.
	*/
	if ((status = StisInit (&sts, &sci_sw, &wav_sw, &sciref, &wavref)))
	    return (status);

	/* Copy switch values for calstis1 argument list. */
	SetCs1Sw (&sci_sw, &wav_sw, &cs1a_sci_sw, &cs1_sci_sw, &cs1_wav_sw);

	if (sts.sci_basic_2d != PERFORM &&
	    sts.sci_basic_2d_a != PERFORM &&
	    sts.sci_crcorr != PERFORM &&
	    sts.sci_rptcorr != PERFORM &&
	    sts.sci_wavecorr != PERFORM &&
	    sts.sci_2d_rect != PERFORM &&
	    sts.sci_1d_extract != PERFORM &&
	    sts.sci_geocorr != PERFORM &&
	    sts.wav_basic_2d != PERFORM) {
	    printf ("Warning  No calibration switch was set to PERFORM.\n");
	    return (NOTHING_TO_DO);
	}

	/* Either do cosmic ray rejection followed by calstis1, or just
	   run calstis1.  If the flags are such that calstis1 does not
	   need to be run, however, we will copy the names of permanent
	   files to the variables used for temporary file names.  This is
	   so we can use those temp names in subsequent steps.
	*/
	if (sts.sci_crcorr == PERFORM) {

	    /* Cosmic ray rejection, followed by basic 2-D processing. */

	    /* First do atodcorr, dqicorr, and blevcorr (or copy the name). */
	    if (sts.sci_basic_2d_a == PERFORM) {
		/* This is calstis1a; note that we use cs1a_sci_sw. */
		if ((status = CalStis1 (sts.rawfile, sts.blv_tmp, "",
                         &cs1a_sci_sw, &sciref, printtime, verbose)))
		    return (status);
		ResetCs1Sw (&cs1a_sci_sw, &cs1_sci_sw);	/* reset switches */
	    } else {
		/* Copy the input file to blv_tmp.  We need to do this
		   because calstis2 can modify the DQ flags in its input file.
		*/
		if ((status = CopyFFile (sts.rawfile, sts.blv_tmp)))
		    return (status);
	    }

	    if (sts.sci_basic_2d != PERFORM) {
		/* Copy the name so that calstis2 will write crjfile
		   directly, because we won't be calling calstis1.
		*/
		strcpy (sts.crj_tmp, sts.crjfile);
		save_crj = 1;
	    }

	    /* Reject cosmic rays. */
	    if ((status = CalStis2_0 (sts.blv_tmp, sts.crj_tmp,
                                      printtime, verbose)))
		return (status);

	    if (sts.sci_basic_2d == PERFORM) {

		/* Flat field the summed, cosmic ray rejected image. */
		if ((status = CalStis1 (sts.crj_tmp, sts.crjfile, "",
                         &cs1_sci_sw, &sciref, printtime, verbose)))
		    return (status);

		if (sts.sci_expscorr == PERFORM) {
		    /* Flat field the cosmic ray flagged but non-summed data. */
		    if ((status = CalStis1 (sts.blv_tmp, sts.fltfile, "",
                             &cs1_sci_sw, &sciref, printtime, verbose)))
			return (status);
		} else {
		    fltexists = 0;	/* fltfile not created */
		}
	    } else {
		fltexists = 0;
	    }
	    if (!save_crj)
		remove (sts.crj_tmp);
	    if (!save_tmp)
		remove (sts.blv_tmp);

	    /* blv_tmp and crj_tmp are not used beyond this point. */

	} else if (sts.sci_basic_2d_a == PERFORM ||
		   sts.sci_basic_2d == PERFORM) {

	    /* Basic 2-D processing (flat field, etc). */

	    if ((status = CalStis1 (sts.rawfile, sts.fltfile, "",
                     &cs1_sci_sw, &sciref, printtime, verbose)))
		return (status);

	} else {

	    strcpy (sts.fltfile, sts.rawfile);
	}

	if (sts.obstype == SPECTROSCOPIC_TYPE) {

	    /* Process wavecal to get MSM shift. */

	    if (sts.sci_wavecorr == PERFORM) {

		/* wavfile --> fwv_tmp */

		if (sts.wav_basic_2d == PERFORM) {
		    if ((status = CalStis1 (sts.wavfile, sts.fwv_tmp, "",
                            &cs1_wav_sw, &wavref, printtime, verbose)))
			return (status);
		} else {
		    /* Copy the wavecal (or just its name) to fwv_tmp. */
		    if (sts.echelle || sts.prism) {
			if ((status = CopyFFile (sts.wavfile, sts.fwv_tmp)))
			    return (status);
		    } else {
			strcpy (sts.fwv_tmp, sts.wavfile);
			save_fwv = 1;		/* no longer a temp file */
		    }
		}

		/* Now do the wavecal processing, fwv_tmp --> w2d_tmp. */

		if (sts.wav_subscicorr == PERFORM) {

		    /* Subtract science from wavecal, then rectify. */
		    if (sts.sci_crcorr == PERFORM)
			status = CalStis11 (sts.fwv_tmp, sts.crjfile,
				sts.cwv_tmp, printtime, verbose);
		    else
			status = CalStis11 (sts.fwv_tmp, sts.fltfile,
				sts.cwv_tmp, printtime, verbose);

		    if (status == 0) {
			if ((status = CalStis7_0 (sts.cwv_tmp, sts.w2d_tmp,
                                  &wav_sw, &wavref, printtime, verbose)))
			    return (status);
			if (!save_tmp)
			    remove (sts.cwv_tmp);
		    } else if (status == NOTHING_TO_DO) {
			/* The science file was not subtracted; rectify
				the wavecal itself.
			*/
			status = 0;
			if ((status = CalStis7_0 (sts.fwv_tmp, sts.w2d_tmp,
                                  &wav_sw, &wavref, printtime, verbose)))
			    return (status);
		    } else {
			return (status);
		    }
		} else {
		    /* Rectify the wavecal, except for echelle or prism. */
		    if (!sts.echelle && !sts.prism) {
			if ((status = CalStis7_0 (sts.fwv_tmp, sts.w2d_tmp,
                                &wav_sw, &wavref, printtime, verbose)))
			    return (status);
		    }
		}

		/* Now compute the MSM shift and update the fltfile header. */

		/* select the appropriate file */
		if (sts.echelle || sts.prism)
		    wavecal_file = sts.fwv_tmp;
		else
		    wavecal_file = sts.w2d_tmp;

		if ((status = CalStis4_0 (wavecal_file, &sciref,
                                          printtime, verbose)))
		    return (status);
		if (sts.sci_crcorr == PERFORM) {
		    /* Update header of crjfile. */
		    if ((status = CalStis12_0 (wavecal_file, sts.crjfile,
                                               printtime, verbose)))
			return (status);
		}
		if (fltexists) {
		    if ((status = CalStis12_0 (wavecal_file, sts.fltfile,
                                               printtime, verbose)))
			return (status);
		}
		if (!save_fwv)
		    remove (sts.fwv_tmp);
		if (!save_tmp && !sts.echelle && !sts.prism) {
		    /* wasn't created for echelle or prism data */
		    remove (sts.w2d_tmp);
		}
	    }

	    if (sts.sci_2d_rect == PERFORM) {
		if (sts.sci_crcorr == PERFORM) {
		    if ((status = CalStis7_0 (sts.crjfile, sts.sx2file,
                              &sci_sw, &sciref, printtime, verbose)))
			return (status);
		} else {
		    if ((status = CalStis7_0 (sts.fltfile, sts.x2dfile,
                              &sci_sw, &sciref, printtime, verbose)))
			return (status);
		    if (sts.sci_rptcorr == PERFORM) {
			status = CalStis8 (sts.x2dfile, sts.sx2file,
					printtime, verbose);
			if (status != 0) {
			    if (status == NOTHING_TO_DO) {
				printf (
		"Warning  Output file %s was not written.\n", sts.sx2file);
				status = 0;
			    } else {
				return (status);
			    }
			}
		    }
		}
	    }

	    if (sts.sci_1d_extract == PERFORM) {

		if (sts.sci_crcorr == PERFORM) {
		    if ((status = CalStis6_0 (sts.crjfile, sts.sx1file,
                                              &sci_sw, printtime, verbose)))
			return (status);
		} else {
		    if ((status = CalStis6_0 (sts.fltfile, sts.x1dfile,
                                              &sci_sw, printtime, verbose)))
			return (status);
/*
### calstis13 is not written yet ...
		    if (sts.sci_rptcorr == PERFORM) {
			if ((status = CalStis13 (sts.x1dfile, sts.sx1file,
					printtime, verbose)))
			    return (status);
		    }
*/
		}
	    }

	} else if (sts.obstype == IMAGING_TYPE && sts.sci_geocorr == PERFORM) {

	    /* 2-D rectification for image. */
	    if (sts.sci_crcorr == PERFORM) {
		if ((status = CalStis7_0 (sts.crjfile, sts.sx2file,
                        &sci_sw, &sciref, printtime, verbose)))
		    return (status);
	    } else {
		if ((status = CalStis7_0 (sts.fltfile, sts.x2dfile,
                        &sci_sw, &sciref, printtime, verbose)))
		    return (status);
		if (sts.sci_rptcorr == PERFORM) {
		    status = CalStis8 (sts.x2dfile, sts.sx2file,
				    printtime, verbose);
			if (status != 0) {
			    if (status == NOTHING_TO_DO) {
				printf (
		"Warning  Output file %s was not written.\n", sts.sx2file);
				status = 0;
			    } else {
				return (status);
			    }
			}
		}
	    }
	}

	if (sts.sci_rptcorr == PERFORM && fltexists) {
	    status = CalStis8 (sts.fltfile, sts.sflfile, printtime, verbose);
	    if (status != 0) {
		if (status == NOTHING_TO_DO) {
		    printf (
		"Warning  Output file %s was not written.\n", sts.sflfile);
		    status = 0;
		} else {
		    return (status);
		}
	    }
	}

	/* Done with lists of reference file keywords and names. */
	FreeRefFile (&sciref);
	FreeRefFile (&wavref);

	printf ("\n");
	PrEnd (0);		/* *** CALSTIS-0 complete *** */

	if (printtime)
	    TimeStamp ("CALSTIS completed", sts.rootname);

	return (0);
}

/* This routine assigns a value to sts->outroot, the name which will be
   used to construct the output file names (by replacing _raw with _flt,
   etc).  If outroot is null, the name of the input file (rawfile) will be
   assigned directly to sts->outroot.  If outroot is a directory, then the
   file name in rawfile (skipping the directory, if any was specified)
   will be appended to outroot.

   The test on directory is to look for /, $, or ].  In order to be
   regarded as a directory, outroot must end in one of those characters
   (e.g. "./" for the default directory, "abc/" for subdirectory abc).
   When checking for a directory portion in rawfile, we first look for
   the last '/' in rawfile.  If none is found, we check for a VMS directory
   name; the test is that rawfile contains a ], but there's something
   following it (i.e. a name such as abc_raw.fits[noinh] would not qualify
   as a filename preceded by a VMS directory).  If none is found, we look
   for the last '$' in rawfile.  If none is found, we assume rawfile does
   not contain a directory prefix.
*/

static void AppendRawfile (StisInfo *sts, char *rawfile, char *outroot) {

/* arguments:
StisInfo *sts     i: outroot is assigned a value in this struct
char *rawfile     i: input file name
char *outroot     i: output root name, or null
*/

	int rlen, olen;	/* lengths of rawfile and outroot strings */
	char *infile;	/* points to beginning of file name in rawfile */

	olen = strlen (outroot);

	if (olen == 0) {

	    strcpy (sts->outroot, rawfile);	/* default value */

	} else {

	    strcpy (sts->outroot, outroot);

	    /* Is outroot a directory name? */
	    if (outroot[olen-1] == '/' ||	/* local subdirectory */
		outroot[olen-1] == '$' ||	/* iraf environment variable */
		outroot[olen-1] == ']') {	/* vms directory */

		rlen = strlen (rawfile);

		/* Append the input file name (rawfile) to sts->outroot,
		   but exclude the directory prefix, if any.
		*/
		infile = strrchr (rawfile, '/');
		if (infile != NULL) {
		    infile++;		/* directory; skip past the '/' */
		} else {
		    /* Check for a VMS directory name (there's a ], but it's
			not the last character in rawfile).
		    */
		    infile = strrchr (rawfile, ']');
		    if (infile != NULL && rawfile[rlen-1] != ']') {
			infile++;	/* yes, VMS; skip past the ']' */
		    } else {
			infile = strrchr (rawfile, '$');
			if (infile != NULL) {
			    infile++;		/* IRAF environment variable */
			} else {
			    infile = rawfile;	/* no directory found */
			}
		    }
		}

		strcat (sts->outroot, infile);
	    }
	}
}


/* This routine copies switch values from sci_sw and wav_sw to
   cs1a_sci_sw, cs1_sci_sw, and cs1_wav_sw.
*/

static void SetCs1Sw (CalSwitch *sci_sw, CalSwitch *wav_sw,
	cs1_switch *cs1a_sci_sw, cs1_switch *cs1_sci_sw,
	cs1_switch *cs1_wav_sw) {

/* arguments:
CalSwitch *sci_sw           i: all calibration switches for science file
CalSwitch *wav_sw           i: all calibration switches for wavecal
cs1_switch *cs1a_sci_sw     o: calstis1a switches for science file
cs1_switch *cs1_sci_sw      o: calstis1 switches for science file
cs1_switch *cs1_wav_sw      o: calstis1 switches for wavecal
*/

	/* These are the switches for calstis1 prior to running calstis2
	   (cosmic-ray rejection).
	*/
	cs1a_sci_sw->dqicorr = sci_sw->dqicorr;
	cs1a_sci_sw->atodcorr = sci_sw->atodcorr;
	cs1a_sci_sw->blevcorr = sci_sw->blevcorr;
	cs1a_sci_sw->biascorr = sci_sw->biascorr;
	cs1a_sci_sw->doppcorr = OMIT;
	cs1a_sci_sw->lorscorr = OMIT;
	cs1a_sci_sw->glincorr = OMIT;
	cs1a_sci_sw->lflgcorr = OMIT;
	cs1a_sci_sw->darkcorr = OMIT;
	cs1a_sci_sw->flatcorr = OMIT;
	cs1a_sci_sw->shadcorr = OMIT;
	cs1a_sci_sw->photcorr = OMIT;
	cs1a_sci_sw->statcorr = OMIT;
	strcpy (cs1a_sci_sw->darkscale_string, "");

	/* These are the switches for calstis1 for the science file. */
	cs1_sci_sw->dqicorr = sci_sw->dqicorr;
	cs1_sci_sw->atodcorr = sci_sw->atodcorr;
	cs1_sci_sw->blevcorr = sci_sw->blevcorr;
	cs1_sci_sw->doppcorr = sci_sw->doppcorr;
	cs1_sci_sw->lorscorr = sci_sw->lorscorr;
	cs1_sci_sw->glincorr = sci_sw->glincorr;
	cs1_sci_sw->lflgcorr = sci_sw->lflgcorr;
	cs1_sci_sw->biascorr = sci_sw->biascorr;
	cs1_sci_sw->darkcorr = sci_sw->darkcorr;
	cs1_sci_sw->flatcorr = sci_sw->flatcorr;
	cs1_sci_sw->shadcorr = sci_sw->shadcorr;
	cs1_sci_sw->photcorr = sci_sw->photcorr;
	cs1_sci_sw->statcorr = sci_sw->statcorr;
	strcpy (cs1_sci_sw->darkscale_string, "");

	/* These are the switches for calstis1 for the wavecal. */
	cs1_wav_sw->dqicorr = wav_sw->dqicorr;
	cs1_wav_sw->atodcorr = wav_sw->atodcorr;
	cs1_wav_sw->blevcorr = wav_sw->blevcorr;
	cs1_wav_sw->doppcorr = wav_sw->doppcorr;
	cs1_wav_sw->lorscorr = wav_sw->lorscorr;
	cs1_wav_sw->glincorr = wav_sw->glincorr;
	cs1_wav_sw->lflgcorr = wav_sw->lflgcorr;
	cs1_wav_sw->biascorr = wav_sw->biascorr;
	cs1_wav_sw->darkcorr = wav_sw->darkcorr;
	cs1_wav_sw->flatcorr = wav_sw->flatcorr;
	cs1_wav_sw->shadcorr = wav_sw->shadcorr;
	cs1_wav_sw->photcorr = wav_sw->photcorr;
	cs1_wav_sw->statcorr = wav_sw->statcorr;
	strcpy (cs1_wav_sw->darkscale_string, "");
}

/* This routine resets dqicorr, atodcorr, blevcorr, and biascorr for
   calstis1.  These steps are done prior to cosmic ray rejection, then
   the other steps for calstis1 are done.  It isn't essential to reset
   these switches, since calstis1 will recognize that blevcorr and
   biascorr have already been done and will turn off those switches
   locally, but dqicorr will be done a second time, so it's a good idea
   to explicitly turn that off.
*/

static void ResetCs1Sw (cs1_switch *cs1a_sci_sw, cs1_switch *cs1_sci_sw) {

/* arguments:
cs1_switch *cs1a_sci_sw     i: calstis1a switches for science file
cs1_switch *cs1_sci_sw      o: calstis1 switches for science file
*/

	/* These are the switches performed by calstis1 prior to running
	   calstis2 (cosmic-ray rejection).
	*/
	if (cs1a_sci_sw->dqicorr == PERFORM)
	    cs1_sci_sw->dqicorr = OMIT;

	if (cs1a_sci_sw->atodcorr == PERFORM)
	    cs1_sci_sw->atodcorr = OMIT;

	if (cs1a_sci_sw->blevcorr == PERFORM)
	    cs1_sci_sw->blevcorr = OMIT;

	if (cs1a_sci_sw->biascorr == PERFORM)
	    cs1_sci_sw->biascorr = OMIT;
}

# define FITS_BUFSIZE  2880	/* size of a FITS block */

/* This routine copies a FITS file. */

static int CopyFFile (char *infile, char *outfile) {

/* arguments:
char *infile    i: name of input file
char *outfile   i: name of output file
*/

	FILE *ifp, *ofp;	/* for input and output files */
	void *buf;		/* buffer for copying blocks */
	int nin, nout;		/* number read and written */
	int done;

	if ((buf = calloc (FITS_BUFSIZE, sizeof(char))) == NULL)
	    return (OUT_OF_MEMORY);

	if ((ofp = fopen (outfile, "wb")) == NULL) {
	    printf ("ERROR    Can't create temporary file %s.\n", outfile);
	    free (buf);
	    return (GENERIC_ERROR_CODE);
	}

	if ((ifp = fopen (infile, "rb")) == NULL) {
	    printf ("ERROR    Can't open %s.\n", infile);
	    fclose (ofp);
	    remove (outfile);
	    free (buf);
	    return (OPEN_FAILED);
	}

	done = 0;
	while (!done) {
	    nin = fread (buf, sizeof(char), FITS_BUFSIZE, ifp);
	    if (ferror (ifp)) {
		printf ("ERROR    Can't read from %s (copying to %s).\n",
				infile, outfile);
		fclose (ofp);
		fclose (ifp);
		free (buf);
		return (GENERIC_ERROR_CODE);
	    }
	    if (feof (ifp))
		done = 1;

	    nout = fwrite (buf, sizeof(char), nin, ofp);
	    if (nout < nin) {
		printf ("ERROR    Can't copy %s to %s.\n", infile, outfile);
		fclose (ofp);
		fclose (ifp);
		free (buf);
		return (GENERIC_ERROR_CODE);
	    }
	}

	fclose (ofp);
	fclose (ifp);
	free (buf);

	return (0);
}

/* This function is an interface for calling calstis2.  It sets up the
   default parameters, copies printtime and verbose into the par structure,
   and calls calstis2.
*/

# include "cs2.h"

static int CalStis2_0 (char *input, char *output, int printtime, int verbose) {

	int status;
	clpar 	par;			/* parameters used */
	int 	newpar[MAX_PAR+1];	/* user specifiable parameters */
	void cs2_reset (clpar *, int []);

	cs2_reset (&par, newpar);
	par.printtime = printtime;
	par.verbose = verbose;

	status = CalStis2 (input, output, &par, newpar);

	return (status);
}

/* This function is an interface for calling calstis4. */

static int CalStis4_0 (char *input, RefFileInfo *refnames,
	int printtime, int verbose) {

	int status;
	char dbgfile[] = "";		/* (not used by pipeline) */
	double slit_angle = 0.;		/* (not used by pipeline) */

	if ((status = CalStis4 (input, dbgfile,
                                refnames, printtime, verbose, slit_angle)))
	    return (status);

	return (0);
}

/* This function is an interface for calling calstis6.  It sets up the
   default parameters and calls calstis6.
*/

static int CalStis6_0 (char *input, char *output,
	CalSwitch *sci_sw, int printtime, int verbose) {

	int status;

	/* Definitions and defaults for calstis6. */
	double cl_a2center = -1.;
	int maxsearch = -1;
	double extrsize = -1., bk1size = -1., bk2size = -1.;
	double bk1offset = -1., bk2offset = -1., bktilt = -100.;
	int bkord = -1;
	int sporder = -1;
	char xtracalg[] = "";
	int extrloc = 1;	/* yes, add extra info to table */
	int ccglobal = 0;	/* use global crosscor everywhere ? */
	double ccthresh = 0.;	/* crosscor theshold */
	int do_profile = 0;	/* use calstis6 as profile generator ? */
	int pstep = 1;		/* pixel step used to compress profiles */
	double wstep = 1.;	/* wavelength step used to compress profiles */
	double minsn = 1.;	/* minimum acceptable S/N */
	char rejranges[] = "";	/* rejection ranges (not used by pipeline) */
	char profilefile[] = "";	/* (not used by pipeline) */
	char fluxfile[] = "";		/* (not used by pipeline) */
	char outw[] = "";		/* (not used by pipeline) */
	double backval = 0.;		/* (not used by pipeline) */
	double backerr = 0.;		/* (not used by pipeline) */
	int variance = 0;		/* (not used by pipeline) */
	int fflux = 1;			/* (not used by pipeline) */
	double psclip = 0.;		/* (not used by pipeline) */
	double sclip = 0.;		/* (not used by pipeline) */
	int lfilter = 17;		/* (not used by pipeline) */
	char idtfile[] = "";		/* (not used by pipeline) */
	double subscale = 10.;		/* (not used by pipeline) */
	double blazeshift = NO_VALUE;	/* (not used by pipeline) */
	int bks_mode = BKS_MEDIAN;	/* activates background smoothing */
	int bks_order = 3;		/* background smoothing poly. order */
	double xoffset = 0.;	/* offset for slitless data */
	int pipeline = 1;	/* calstis6 is being run from the pipeline */

	if ((status = CalStis6 (input, output,
		sci_sw->backcorr, sci_sw->dispcorr, sci_sw->fluxcorr,
		sci_sw->helcorr, sci_sw->sgeocorr, sci_sw->ctecorr,
                sci_sw->sc2dcorr,
		cl_a2center, maxsearch, extrsize,
		bk1size, bk2size, bk1offset, bk2offset, bktilt, bkord,
		sporder, xtracalg, printtime, verbose,
		extrloc,
		ccglobal, ccthresh,
		do_profile, pstep, wstep, minsn,
		rejranges,
		profilefile, fluxfile, outw,
		backval, backerr, variance, fflux,
		psclip, sclip,
		lfilter, idtfile,
		subscale, blazeshift,
		bks_mode, bks_order, xoffset,
                pipeline)))
	    return (status);

	return (0);
}

/* This function is an interface for calling calstis7. */

static int CalStis7_0 (char *input, char *output,
		CalSwitch *sci_sw, RefFileInfo *refnames,
		int printtime, int verbose) {

	int status;
	int center_target = 0;	/* center target in output image */
	int err_algorithm = WGT_VARIANCE;

	if ((status = CalStis7 (input, output,
			sci_sw->sgeocorr, sci_sw->helcorr,
			sci_sw->fluxcorr, sci_sw->statcorr,
			refnames, printtime, verbose, center_target,
                        NO_VALUE, err_algorithm)))
	    return (status);

	return (0);
}

/* This function is an interface for calling calstis12. */

# include "cs12.h"		/* for interpolation options */

static int CalStis12_0 (char *wavecal, char *science,
		int printtime, int verbose) {

	int status;
	/* linear interpolation, if there are multiple wavecal imsets */
	int w_option = STIS_LINEAR;

	if ((status = CalStis12 (wavecal, science,
                                 w_option, printtime, verbose)))
	    return (status);

	return (0);
}
