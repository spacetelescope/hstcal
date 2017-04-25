/* calstis7 -- 2-D spectral rectification

   This file contains:
	CalStis7
	StisInit7
*/

# include <stdio.h>
# include <string.h>

# include "hstio.h"

# include "stis.h"
# include "calstis7.h"
# include "hstcalerr.h"

static void StisInit7 (StisInfo7 *);

/* Do basic 2-D calibration.

   Phil Hodge, 1997 Dec 11:
	Include refnames in calling sequence, and copy to sts;
	initialize pctcorr and pctab.

   Phil Hodge, 1998 Aug 24:
	The arguments sgeocorr, helcorr, fluxcorr, and statcorr were
	being used as if they were boolean variables.  In fact, they
	take on the values PERFORM, OMIT, etc., so just assign them
	to the StisInfo7 structure members.

   Phil Hodge, 2000 Apr 5:
	Add center_target to the calling sequence, and copy to sts.

   Phil Hodge, 2000 Aug 9:
	Remove references to MAMA offset table or coefficients.

   Ivo Busko, 2002 Apr 24:
        Add blazeshift command line argument.

   Phil Hodge, 2004 Dec 27:
	Remove InitRefTab from this file.  Initialize detector_temp to -1.

   Phil Hodge, 2006 Jan 26:
	Add err_algorithm argument.
*/

int CalStis7 (char *input, char *output,
		int sgeocorr, int helcorr, int fluxcorr, int statcorr,
		RefFileInfo *refnames, int printtime, int verbose,
		int center_target, double blazeshift, int err_algorithm) {

/* arguments:
char *input            i: name of input file
char *output           i: name of output file
int sgeocorr, helcorr, fluxcorr, statcorr    i: calibration switches
RefFileInfo *refnames  i: reference file keywords and names
int printtime          i: true means that TimeStamp should be called
int verbose            i: print extra info?
int center_target      i: center the target in output image?
int blazeshift         i: blaze shift
int err_algorithm      i: specifies how interpolation of error estimates
				should be done
*/

	int status;

	StisInfo7 sts;	/* calibration switches, reference files, etc. */

	IODescPtr im;		/* descriptor for input image */
	Hdr phdr;		/* primary header for input image */

	int Do2Dx (StisInfo7 *);
	int FileExists (char *);
	int GetFlags7 (StisInfo7 *, Hdr *);
	int GetKeyInfo7 (StisInfo7 *, Hdr *);

	PrBegin (7);

	if (printtime)
	    TimeStamp ("CALSTIS-7 started", "");

	/* Initialize structure containing calstis information. */
	StisInit7 (&sts);

	/* Copy command-line arguments into sts. */
	strcpy (sts.input, input);
	strcpy (sts.output, output);
	sts.sgeocorr = sgeocorr;
	sts.heliocorr = helcorr;
	sts.fluxcorr = fluxcorr;
	sts.pctcorr = fluxcorr;		/* not a header switch */
	sts.statcorr = statcorr;
	sts.printtime = printtime;
	sts.verbose = verbose;
	sts.center_target = center_target;
	sts.blazeshift = blazeshift;
	sts.err_algorithm = err_algorithm;

	sts.refnames = refnames;

	PrFileName ("input", sts.input);
	PrFileName ("output", sts.output);

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
	if ((status = GetKeyInfo7 (&sts, &phdr)))
	    return (status);

	/* Print information about this image. */
	PrHdrInfo (sts.obsmode, sts.aperture, sts.opt_elem, sts.det);

	/* Get calibration file names from input image header. */
	if ((status = GetFlags7 (&sts, &phdr)))
	    return (status);

	freeHdr (&phdr);

	if (sts.printtime)
	    TimeStamp ("Begin processing", sts.rootname);

	/* Do 2-D spectral extraction. */
	if ((status = Do2Dx (&sts)))
	    return (status);

	printf ("\n");
	PrEnd (7);

	if (sts.printtime)
	    TimeStamp ("CALSTIS-7 completed", sts.rootname);

	return (0);
}

/* Initialize the calstis7 structure.  This includes information about the
   input and output images, calibration files, and flags to specify which
   calibration steps are to be done.
*/

static void StisInit7 (StisInfo7 *sts) {

	/* Assign default values. */

	sts->input[0] = '\0';
	sts->output[0] = '\0';
	sts->rootname[0] = '\0';
	sts->obsmode[0] = '\0';
	sts->aperture[0] = '\0';
	sts->opt_elem[0] = '\0';
	sts->det[0] = '\0';
	sts->detector = UNKNOWN_DETECTOR;
	sts->obstype = UNKNOWN_TYPE;
	sts->wavecal = 0;
	sts->nimages = 1;
	sts->sdqflags = 32767;			/* 15 bits set */
	sts->exptime = 1.;
	sts->hfactor = 1.;
	sts->ltm[0] = 1.;
	sts->ltm[1] = 1.;
	sts->ltv[0] = 0.;
	sts->ltv[1] = 0.;
	sts->cenwave = 0;
	sts->dispaxis = 1;
	sts->atodgain = 1.;
	sts->detector_temp = -1.;

	sts->msm_slop[0] = 0.;
	sts->msm_slop[1] = 0.;
	sts->ap_offset[0] = 0.;
	sts->ap_offset[1] = 0.;
	sts->total_offset[0] = 0.;
	sts->total_offset[1] = 0.;

	/* Initialize flags to not perform most steps. */
	sts->wavecorr = OMIT;
	sts->x2dcorr = PERFORM;		/* else why run cs7? */
	sts->sgeocorr = OMIT;
	sts->heliocorr = OMIT;
	sts->fluxcorr = OMIT;
	sts->pctcorr = OMIT;	/* not a header switch; goes with fluxcorr */
	sts->fc_corr = OMIT;		/* currently not a header switch */
	sts->statcorr = OMIT;		/* header flag is STATFLAG */
	sts->tdscorr = OMIT;

	/* No info yet about whether reference files exist, and
	   therefore no values for pedigree and descrip.
	*/

	sts->sdstfile.name[0] = '\0';			/* image */
	sts->sdstfile.pedigree[0] = '\0';
	sts->sdstfile.descrip[0] = '\0';
	sts->sdstfile.exists = EXISTS_UNKNOWN;
	sts->sdstfile.goodPedigree = PEDIGREE_UNKNOWN;

	/* tables */
	InitRefTab (&sts->distntab);
	InitRefTab (&sts->apdestab);
	InitRefTab (&sts->apertab);
	InitRefTab (&sts->phottab);
	InitRefTab (&sts->tdstab);
	InitRefTab (&sts->disptab);
	InitRefTab (&sts->inangtab);
	InitRefTab (&sts->sptrctab);
	InitRefTab (&sts->pctab);

	sts->trace_rotation = 0.;
}
