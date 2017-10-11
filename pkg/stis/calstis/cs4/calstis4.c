/* calstis4 -- determine offsets from wavecal

   This file contains:
	CalStis4
     internal:
	StisInit4
*/

# include <stdio.h>
# include <string.h>

# include "hstio.h"

# include "stis.h"
# include "calstis4.h"
# include "hstcalerr.h"

static void StisInit4 (StisInfo4 *);

/* This determines the offset in each axis due to MSM uncertainty.

   Phil Hodge, 1997 Dec 11:
	Include refnames in calling sequence, and copy to sts.

   Phil Hodge, 1998 Dec 11:
	Add dbgfile and dbg pointer.  Add wcptab in StisInit4, and
	assign default values for those parameters.

   Phil Hodge, 2000 Jan 14:
	Open the debug file here only if the input file does not contain
	echelle data.  (For echelle data, the debug file will be a FITS
	file, and we'll write it in MakeTemplate.)  Use InitRefTab.

   Phil Hodge, 2000 Mar 1:
	Write the input file name and rootname to the debug file.

   Phil Hodge, 2001 Feb 27:
	Set disp_type in StisInit4.

   Phil Hodge, 2004 July 23:
	Add slit_angle to the calling sequence for CalStis4.  Initialize it
	to zero in StisInit4.

   Phil Hodge, 2004 Dec 27:
	Remove InitRefTab from this file.
*/

int CalStis4 (char *input, char *dbgfile,
	RefFileInfo *refnames, int printtime, int verbose, double slit_angle) {

/* arguments:
char *input            i: name of input (modified in-place)
char *dbgfile          i: output text file for debug info
RefFileInfo *refnames  i: reference file keywords and names
int printtime          i: true means that TimeStamp should be called
int verbose            i: print extra info?
double slit_angle      i: angle of long slit used with echelle; this is
                          given in degrees but will be converted to radians
*/

	int status;

	StisInfo4 sts;	/* calibration switches, reference files, etc */

	IODescPtr im;		/* descriptor for input image */
	Hdr phdr;		/* primary header for input image */
	int i;

	int GetFlags4 (StisInfo4 *, Hdr *);
	int GetKeyInfo4 (StisInfo4 *, Hdr *);
	int WaveCal (StisInfo4 *);

	PrBegin (4);

	if (printtime)
	    TimeStamp ("CALSTIS-4 started", "");

	/* Initialize structure containing calstis information. */
	StisInit4 (&sts);

	/* Copy command-line arguments into sts. */
	strcpy (sts.input, input);
	sts.printtime = printtime;
	sts.verbose = verbose;
	sts.slit_angle = DEGREES_TO_RADIANS * slit_angle;

	sts.refnames = refnames;

	for (i = 0;  i < strlen (dbgfile);  i++) {
	    if (dbgfile[i] != ' ') {
		strcpy (sts.dbgfile, dbgfile+i);
		sts.verbose = 1;	/* turn on verbose if debug specified */
		break;
	    }
	}

	PrFileName ("input", sts.input);

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
	if ((status = GetKeyInfo4 (&sts, &phdr)))
	    return (status);

	/* Now we know whether we have echelle data or not.  If an echelle
	   was used, append ".fits" to the debug file name if it doesn't
	   already end in ".fit" or ".fits".  If we have first-order data,
	   open the debug file now.
	*/
	if (sts.dbgfile[0] != '\0') {

	    if (sts.disp_type == ECHELLE_DISP) {

		/* append ".fits" if necessary */
		int len;
		len = strlen (sts.dbgfile);
		if (len < 5 ||
		    (strcmp (&sts.dbgfile[len-4], ".fit") != 0 &&
		     strcmp (&sts.dbgfile[len-5], ".fits") != 0)) {

		    strcat (sts.dbgfile, ".fits");
		}

	    } else {

		sts.dbg = fopen (sts.dbgfile, "a");
		if (sts.dbg == NULL) {
		    printf ("Warning  Can't open debug file %s\n", sts.dbgfile);
		    sts.dbgfile[0] = '\0';
		}
	    }
	}

	if (sts.dbgfile[0] != '\0')
	    PrFileName ("debugfile", sts.dbgfile);
	if (sts.dbg != NULL) {
	     fprintf (sts.dbg, "# Input file %s, rootname %s ###\n",
			sts.input, sts.rootname);
	}

	/* Print information about this image. */
	PrHdrInfo (sts.obsmode, sts.aperture, sts.opt_elem, sts.det);

	/* Get file names from input image header. */
	if ((status = GetFlags4 (&sts, &phdr)))
	    return (status);

	freeHdr (&phdr);

	if (sts.printtime)
	    TimeStamp ("Begin processing", sts.rootname);

	/* Do wavecal processing */
	if ((status = WaveCal (&sts)))
	    return (status);

	printf ("\n");
	PrEnd (4);

	if ((status = fcloseWithStatus(&sts.dbg)))
	    return status;

	if (sts.printtime)
	    TimeStamp ("CALSTIS-4 completed", sts.rootname);

	return (0);
}

/* Initialize the calstis4 structure.  This includes information about the
   input image, calibration files, and flags to specify which calibration
   steps are to be done.
*/

static void StisInit4 (StisInfo4 *sts) {

	/* Assign default values. */

	sts->input[0] = '\0';
	sts->dbgfile[0] = '\0';
	sts->dbg = NULL;
	sts->slit_angle = 0.;
	sts->rootname[0] = '\0';
	sts->aperture[0] = '\0';
	sts->detector = UNKNOWN_DETECTOR;
	sts->opt_elem[0] = '\0';
	sts->nimages = 1;
	sts->sdqflags = 32767;			/* 15 bits set */
	sts->cenwave = 0;
	sts->dispaxis = 1;
	sts->disp_type = RECTIFIED;

	/* Default values for parameters from the WCP table.
	   These will be used if there is no WCP table.
	*/
	sts->wl_trim1 = 0;
	sts->wl_trim2 = 300;
	sts->sp_trim1 = 200;
	sts->sp_trim2 = 0;
	sts->wl_range = 63;
	sts->sp_range = 61;
	sts->nsigma_cr = 3.;
	sts->nsigma_illum = 2.;
	sts->mad_reject = 3.;
	sts->min_mad = 1.;

	/* Initialize flag to not perform the step. */
	sts->wavecorr = OMIT;

	/* tables */
	InitRefTab (&sts->wcptab);
	InitRefTab (&sts->lamptab);
	InitRefTab (&sts->apdestab);
	InitRefTab (&sts->disptab);
	InitRefTab (&sts->inangtab);
	InitRefTab (&sts->sptrctab);

	/* Initialize trace rotation angle */
	sts->trace_rotation = 0.;
}
