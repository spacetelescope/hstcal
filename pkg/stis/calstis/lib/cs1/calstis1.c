/* calstis1 -- basic 2-D image reduction

   This file contains:
	CalStis1
	StisInit1
	InitRefImage
*/

# include <stdio.h>
# include <time.h>
# include <string.h>

# include "hstio.h"

# include "stis.h"
# include "calstis1.h"
# include "hstcalerr.h"

static void StisInit1 (StisInfo1 *);
static void InitRefImage (RefImage *);
static int ParseDarkscaleString (cs1_switch *, StisInfo1 *);

/* Do basic 2-D calibration.

   Phil Hodge, 1997 Dec 10:
	Include refnames in calling sequence, and copy to sts;
	initialize filtcorr switch; initialize apertab info.

   Phil Hodge, 1998 June 8:
	In StisInit1, add obstype and aper_fov, and delete subarray.

   Phil Hodge, 1998 July 30:
	In StisInit1, delete assoc_typ and add wavecal.

   Phil Hodge, 1998 Sept 24:
	Remove bias_rej.

   Phil Hodge, 1998 Oct 5:
	Add InitRefImage and InitRefTab to this file, and call them.

   Phil Hodge, 1999 July 26:
	In StisInit1, initialize the temperature to -1.

   Ivo Busko, 2002 Mar 19:
	Function that parses darkscale command line string.

   Paul Barrett, 2003 Sep 18:
        Add GetEPCTab.

   Paul Barrett, 2003 Sep 25:
        Add tdscorr and tdstab for TDS correction.

   Paul Barrett, 2004 Jun 25:
        Moved GetEPCTab to GetGrpInfo1.

   Phil Hodge, 2004 Sept 3:
	In StisInit1, initialize ccd_temperature to -1.

   Phil Hodge, 2004 Dec 27:
	Remove InitRefTab from this file.
	Initialize detector_temp instead of temperature and ccd_temperature.

   Phil Hodge, 2007 May 9:
	Add output_extver to the calling sequence for Do2D; this will be
	different from extver only in the case when one or more imsets are
	blank (zero exposure time or constant value).
	Update NEXTEND if the final output_extver is less than nimages.
	If there are no output imsets, however, return NOTHING_TO_DO.

   Phil Hodge, 2008 Nov 3:
	Rename variable output_extver to ngood_extver.  Imsets with zero
	exposure time or constant value will be written to output, but with
	a new header keyword IMSET_OK set to F.  ngood_extver will be a
	counter of the number of good (i.e. not exptime = 0) imsets.
	NOTHING_TO_DO will be returned if ngood_extver is 0.

   Phil Hodge, 2011 May 9:
	CalStis1 was modified to not set sts.filtcorr or sts.tdscorr.
	StisInit1 was modified to not call InitRefTab for tdstab or for
	apertab.

   Phil Hodge, 2011 Nov 17:
	Include tdscorr and tdstab.

   Robert Jedrzejewski, 2017 Feb 13:
        Initialize tdctab
*/

int CalStis1 (char *input, char *output, char *outblev,
		cs1_switch *cs1_sw, RefFileInfo *refnames,
		int printtime, int verbose) {

	int status;

	StisInfo1 sts;	/* calibration switches, reference files, etc */
	int extver;
	int ngood_extver;	/* counter for "good" image sets (imsets) */

	IODescPtr im;		/* descriptor for input image */
	Hdr phdr;		/* primary header for input image */

	int Do2D (StisInfo1 *, int, int *);
	int GetCCDTab (StisInfo1 *);
	int GetFlags1 (StisInfo1 *, Hdr *);
	int GetKeyInfo1 (StisInfo1 *, Hdr *);
	int GetLinTab (StisInfo1 *);
	void Sanity1 (StisInfo1 *);

	PrBegin (1);

	if (printtime)
	    TimeStamp ("CALSTIS-1 started", "");

	/* Initialize structure containing calstis information. */
	StisInit1 (&sts);

	/* Copy command-line arguments into sts. */
	strcpy (sts.input, input);
	strcpy (sts.output, output);
	strcpy (sts.outblev, outblev);
	sts.dqicorr  = cs1_sw->dqicorr;
	sts.atodcorr = cs1_sw->atodcorr;
	sts.blevcorr = cs1_sw->blevcorr;
	sts.doppcorr = cs1_sw->doppcorr;
	sts.lorscorr = cs1_sw->lorscorr;
	sts.glincorr = cs1_sw->glincorr;
	sts.lflgcorr = cs1_sw->lflgcorr;
	sts.biascorr = cs1_sw->biascorr;
	sts.darkcorr = cs1_sw->darkcorr;
	sts.flatcorr = cs1_sw->flatcorr;
	sts.shadcorr = cs1_sw->shadcorr;
	sts.photcorr = cs1_sw->photcorr;
	sts.tdscorr  = sts.photcorr;		/* not a header switch */
	sts.statcorr = cs1_sw->statcorr;
	sts.noisecorr = PERFORM;
	sts.printtime = printtime;
	sts.verbose = verbose;

	sts.refnames = refnames;

	/* Parse darkscale string. */
        if ((status = ParseDarkscaleString (cs1_sw, &sts)))
            return (status);

	PrFileName ("input", sts.input);
	PrFileName ("output", sts.output);
	if (sts.outblev[0] != '\0')
	    PrFileName ("outblev", sts.outblev);

	/* Check whether the output file already exists. */
	if ((status = FileExists (sts.output)))
	    return (status);

	initHdr (&phdr);

	/* Open input image in order to read its primary header. */
	im = openInputImage (sts.input, "", 0);
	if (hstio_err())
	    return (OPEN_FAILED);
	getHeader (im, &phdr);		/* get primary header */
	if (hstio_err())
	    return (OPEN_FAILED);
	closeImage (im);

	/* Get keyword values from primary header. */
	if ((status = GetKeyInfo1 (&sts, &phdr)))
	    return (status);

	/* Check match between number of IMSETS and number of darkscale
	   entries in command line.
	*/
	if (sts.ndarkscale > 0) {
	    if (sts.ndarkscale != sts.nimages) {
	        printf (
"Warning: number of IMSETS = %d, number of darkscale entries = %d\n",
                        sts.nimages, sts.ndarkscale);
	        fflush (stdout);
	    }
	}

	/* Print information about this image. */
	PrHdrInfo (sts.obsmode, sts.aperture, sts.opt_elem, sts.det);

	/* Reset switches that are inappropriate for current detector. */
	Sanity1 (&sts);

	/* Get reference file names from input image header.  Pedigree is
	   checked, and the calibration switch (an internal flag, not the
	   header keyword) will be reset from PERFORM to DUMMY if the
	   reference file has pedigree = "DUMMY".  Switches that are
	   currently set to PERFORM will be reset to OMIT if the value
	   in the header is COMPLETE.
	*/
	if ((status = GetFlags1 (&sts, &phdr)))
	    return (status);

	freeHdr (&phdr);

	/* Get values from tables. */

	if (sts.detector == CCD_DETECTOR) {
	    if ((status = GetCCDTab (&sts)))
		return (status);
	}

	if (sts.glincorr == PERFORM || sts.lflgcorr == PERFORM) {
	    if ((status = GetLinTab (&sts)))
		return (status);
	}

	/* Do basic 2-D image reduction. */

	if (sts.printtime)
	    TimeStamp ("Begin processing", sts.rootname);

	ngood_extver = 0;
	for (extver = 1;  extver <= sts.nimages;  extver++) {
	    printf ("\n");
	    PrGrpBegin ("imset", extver);
	    if ((status = Do2D (&sts, extver, &ngood_extver)) != 0)
		return (status);
	    PrGrpEnd ("imset", extver);
	}
	if (ngood_extver <= 0) {
	    printf ("Warning  No good data were written to output.\n");
	    return (NOTHING_TO_DO);
	}

	printf ("\n");
	PrEnd (1);

	if (sts.printtime)
	    TimeStamp ("CALSTIS-1 completed", sts.rootname);

	return (0);
}

/* Assign initial values.  The default values may be set in GetKeyInfo1
   or GetGrpInfo1, and they may be different from these.
*/

static void StisInit1 (StisInfo1 *sts) {

	sts->input[0] = '\0';
	sts->output[0] = '\0';
	sts->outblev[0] = '\0';
	sts->rootname[0] = '\0';
	sts->det[0] = '\0';
	sts->aperture[0] = '\0';
	sts->opt_elem[0] = '\0';
	sts->obsmode[0] = '\0';
	sts->obstype[0] = '\0';
	sts->aper_fov[0] = '\0';
	sts->printtime = 0;
	sts->detector = UNKNOWN_DETECTOR;
	sts->wavecal = 0;
	sts->nimages = 1;
	sts->bin[0] = 1;
	sts->bin[1] = 1;
	sts->dispaxis = 2;
	sts->dispsign = 1;
	sts->exptime = 1.;

	sts->sdqflags = 32767;			/* 15 bits set */

	sts->detector_temp = -1.;

	/* MAMA-specific info */
	sts->globrate = 0.;
	sts->expstart = 0.;
	sts->doppzero = 0.;
	sts->doppmag = 0.;
	sts->orbitper = 96.;

	/* CCD-specific info */
	sts->ccdamp[0] = '\0';
	sts->ccdgain = 1;
	sts->atodgain = 1.;
	sts->ccdbias = 0.;
	sts->readnoise = 0.;

	/* No info yet about whether reference files exist, and
	   therefore no values for pedigree and descrip.
	*/

	/* images */
	InitRefImage (&sts->bias);
	InitRefImage (&sts->dark);
	InitRefImage (&sts->pflt);
	InitRefImage (&sts->dflt);
	InitRefImage (&sts->lflt);
	InitRefImage (&sts->shad);

	/* tables */
	InitRefTab (&sts->bpix);
	InitRefTab (&sts->ccdpar);
        InitRefTab (&sts->epctab);
	InitRefTab (&sts->mlin);
	InitRefTab (&sts->atod);
	InitRefTab (&sts->phot);
	InitRefTab (&sts->tdstab);
	InitRefTab (&sts->tdctab);
}

/* Initialize the elements of a RefImage structure. */

static void InitRefImage (RefImage *image) {

	image->name[0] = '\0';
	image->pedigree[0] = '\0';
	image->descrip[0] = '\0';
	image->exists = EXISTS_UNKNOWN;
	image->goodPedigree = PEDIGREE_UNKNOWN;
}

/* Parses the darkscale string passed via command line.
*/

static int ParseDarkscaleString (cs1_switch *cs1_sw, StisInfo1 *sts) {

	int i;

	int strtor (char *, float[]);

	for (i = 0; i < MAX_IMSET; sts->darkscale[i++] = 0.0);
	sts->ndarkscale = 0;

	sts->ndarkscale = strtor (cs1_sw->darkscale_string, sts->darkscale);

	if (sts->ndarkscale < 0)
	    return (ERROR_RETURN);
	else
	    return (STIS_OK);
}
