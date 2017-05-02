/* This file includes:
	Do2Dx
	OpenSGeo
	AddOffsets
	x2dMsg
	FreeDist
	FreePhot
	FreeThroughput
*/

# include <stdio.h>
# include <stdlib.h>
# include <math.h>		/* fabs */
# include <string.h>

# include "c_iraf.h"
# include "hstio.h"

# include "stis.h"
# include "calstis7.h"
# include "hstcalerr.h"
# include "stisdef.h"
# include "stispht.h"
# include "stistds.h"

static int MOCAdjustDisp (StisInfo7 *, DispRelation *);
static void AddOffsets (StisInfo7 *, ApInfo *);
static void FreeDist (DistInfo *);
static void FreePhot (PhotInfo *);
static void FreeThroughput (ApInfo *);
static int OpenSGeo (char *, FloatHdrData *, FloatHdrData *);
static void SetCrpix2 (CoordInfo *);
static void x2dMsg (StisInfo7 *, int);
static void FluxMsg (StisInfo7 *, int);

/* This routine extracts a 2-D spectrum from the input image.  Each
   extver will be extracted separately.  The primary header will be
   updated with history information, and the calibration switches in
   the header will be reset from "PERFORM" to "COMPLETE".

   Phil Hodge, 1997 Sept 12:
	Change calling sequence of ScaleWCS to include input image size.

   Phil Hodge, 1997 Dec 12:
	Remove dispaxis from calling sequence of AbsFlux;
	replace dispaxis with obstype in PutGrpInfo7.
	Call GetPCT; print PCTAB info; free phot->pcorr in FreePhot.
	Move call to FreePhot inside loop on sporder.
	Assign first_order, depending on maxorder.

   Phil Hodge, 1998 Jan 26:
	Change message about PCTAB not being included in FLUXCORR to
	not being included in DIFF2PT.

   Phil Hodge, 1998 Feb 4:
	Print messages about DISPCORR = PERFORM and COMPLETE.

   Phil Hodge, 1998 Mar 17:
	Add a comment to AddOffsets to say that total_offset[0] is not used
	for spectroscopic type data.

   Phil Hodge, 1998 Oct 5:
	Change status value 1001 to GENERIC_ERROR_CODE.

   Phil Hodge, 1999 Sept 23:
	If exptime = 0, reset fluxcorr (for current imset) to omit.
	Allocate memory using malloc for in and out SingleGroups,
	and add a comment where the input can be transposed if dispaxis = 2.

   Phil Hodge, 2000 Jan 12:
	Remove coord_o and slit from calling sequence of GeoCorr7.
	Free coords and coord_o for both spectroscopic and imaging data.

   Phil Hodge, 2000 Apr 5:
	Add SetCrpix2; call this for first-order data if center_target is true.

   Phil Hodge, 2000 Aug 9:
	Move the call to SetCrpix2 from just before to just after the call
	to ScaleWCS.
	Only call GetApDes7 for spectroscopic data.
	In FreeDist, the coefficients are now xcoeff & ycoeff.
	For imaging data, set total offset to zero.
	Remove ssgx and ssgy from the calling sequence of X2DCorr7.
	Remove all references to MAMA offset table or coefficients.

   Phil Hodge, 2001 May 4:
	Move the call to GetDisp to here from X2DCorr7.
	Call GetMOC, and include this correction for echelle data.

   Ivo Busko, 2002 Jan 10:
	Add TdsInfo structure to store time-dependent sensitivity data.
	Add call the GetTds and FreeTds functions.
	Modify calling sequence to AbsFlux to include tds structure.

   Ivo Busko, 2002 Feb 14:
	Add modified PhtInfo structure to store MSM/blaze shift data.
	Add code to support MSM/blaze shift algorithm.

   Phil Hodge, 2003 Jan 17:
	Change the calling sequence of HelioFactor; update V_HELIO keyword.

   Phil Hodge, 2003 Feb 24:
	Include units for v_helio in the comment field for the keyword.

   Phil Hodge, 2003 July 7:
	Include a check on status == ROW_NOT_FOUND after calling GetTds,
	so that tdscorr will be set to omit if no matching row is found
	in the TDS table.  In FluxMsg, change the warning about TDS
	correction to remove the reference to DIFF2PT.

   Phil Hodge, 2006 Sept 21:
	In x2dMsg, if the input file is output from the wx2d task, print
	a message saying that the trace was applied earlier by wx2d.
*/

int Do2Dx (StisInfo7 *sts) {

/* argument:
StisInfo7 *sts    i: calibration switches and info
*/

	int status;

	IODescPtr oim;		/* descriptor for output image */
	Hdr phdr;		/* primary header for output image */
	DistInfo dist;		/* distortion info */
	DispRelation *disp;	/* list of dispersion relations */
	DispRelation *disp_y;	/* dispersion relation interpolated at y */
	CoordInfo *coords;	/* list of size & coordinate info records */
	CoordInfo *coord_o;	/* one record extracted from coords */
	PhotInfo phot;		/* for conversion to absolute flux */
	PhotInfo photb;		/* for MSM/blaze correction */
	TdsInfo tds;		/* for time-dependent sensitivity correction */
	ApInfo slit;		/* description of slit */
	InangInfo iac;		/* incidence-angle correction coeff */
	SpTrace *trace;		/* spectral trace for spectral order mref */
	SingleGroup *in, *out;	/* input and output data */
	FloatHdrData ssgx;	/* small-scale distortion in X */
	FloatHdrData ssgy;	/* small-scale distortion in Y */
	int fluxcorr_extver;	/* local to each imset */
	int extver;		/* loop index for extension version number */
	int o_extver;		/* extension version number for output */
	int sporder;		/* current spectral order number */
	int minorder, maxorder;	/* min & max spectral order */
	int disp_type;		/* grating or prism ? */
	int mref;		/* MSM/blaze correction reference order */
	double ypos;		/* MSM/blaze correction reference position */
	double wpos;		/* MSM/blaze correction reference wavelength */
	double ddisp;		/* MSM/blaze correction reference dispersion*/
	int i;
	double blazeshift;	/* MSM/blaze shift from command line */
	double delta;
	double v_helio;		/* radial velocity in direction of target */
	int warn, warn1, warn2;

	int AbsFlux (StisInfo7 *, SingleGroup *, CoordInfo *, PhotInfo *,
		ApInfo *, TdsInfo *, double *, double, double, double, int,
		double *);
	void AdjustDisp (StisInfo7 *, DispRelation *,double, InangInfo *,
	                 int *, int *);
	int InterpDisp (DispRelation **, double, DispRelation **);
	int doStat (SingleGroup *, short);
	int GeoCorr7 (StisInfo7 *, DistInfo *,
		FloatHdrData *, FloatHdrData *, SingleGroup *, SingleGroup *);
	int GetAbsPhot (StisInfo7 *, int, PhotInfo *, int, int *);
	int GetApDes7 (StisInfo7 *, ApInfo *);
	int GetApThr (StisInfo7 *, ApInfo *);
	int GetApOffset (StisInfo7 *, ApInfo *, char *, double *);
	int GetDisp (StisInfo7 *, DispRelation **);
	double GetWavelength (DispRelation *, int, double, int, double, double);
	void FreeDisp (DispRelation **);
	void FreeTds (TdsInfo *);
	int GetGrpInfo7 (StisInfo7 *, Hdr *);
	int GetIDC (StisInfo7 *, CoordInfo **, DistInfo *);
	int GetInang (StisInfo7 *, RefTab *, int, InangInfo *);
	int GetSDC (StisInfo7 *, CoordInfo **, int *, int *);
	int GetPCT (StisInfo7 *, PhotInfo *);
	int GetTds (char *, char *, TdsInfo *);
	int GetTrace (StisInfo7 *, int, SpTrace **);
	double HelioFactor (StisInfo7 *, double *);
	int History7 (StisInfo7 *, Hdr *);
	int PutGrpInfo7 (SingleGroup *, SingleGroup *, CoordInfo *, int, int);
	int ReturnCoord (CoordInfo **, int, CoordInfo **);
	void FreeCoord (CoordInfo **);
	void FreeTrace (SpTrace **);
	int ScaleWCS (StisInfo7 *, int, int, CoordInfo *);
	int X2DCorr7 (StisInfo7 *, CoordInfo *, DispRelation *, ApInfo *,
		int, SingleGroup *, SingleGroup *);

	o_extver = 0;

	sts->x2dcorr_o = PERFORM;	/* before reading ref tables */

	/* Set flags to indicate that memory has not been allocated yet. */
	slit.allocated  = 0;
	dist.allocated  = 0;
	phot.allocated  = 0;
	photb.allocated = 0;
	tds.allocated   = 0;
	iac.allocated   = 0;
	phot.pcorr  = NULL;
	photb.pcorr = NULL;
	trace       = NULL;

	disp    = NULL;
	disp_y  = NULL;
	coords  = NULL;
	coord_o = NULL;

	/* Allocate memory for SingleGroup structures. */
	in = malloc (sizeof (SingleGroup));
	out = malloc (sizeof (SingleGroup));
	if (in == NULL || out == NULL)
	    return (OUT_OF_MEMORY);

	/* Get aperture description from APD. */
	if (sts->obstype == SPECTROSCOPIC_TYPE) {
	    if ((status = GetApDes7 (sts, &slit)))
		return (status);
	}

	/* Get aperture throughput info. */
	if (sts->fluxcorr == PERFORM) {
	    if ((status = GetApThr (sts, &slit)))
		return (status);
	}

	if (sts->tdscorr == PERFORM) {
	    /* Read time-dependent sensitivity info into memory. */
	    status = GetTds (sts->tdstab.name, sts->opt_elem, &tds);
            if (status == DUMMY || status == ROW_NOT_FOUND)
                sts->tdscorr = OMIT;
            else if (status)
		return (status);
	}

	if (sts->sgeocorr == PERFORM) {
	    /* Read the small-scale distortion data into memory. */
	    if ((status = OpenSGeo (sts->sdstfile.name, &ssgx, &ssgy)))
		return (status);
	}

	if (sts->obstype == SPECTROSCOPIC_TYPE) {

	    /* Get output coordinate parameters for all spectral orders. */
	    if ((status = GetSDC (sts, &coords, &minorder, &maxorder)))
		return (status);

	    /* Get dispersion coefficients from the DSPTAB table. */
	    if ((status = GetDisp (sts, &disp)))
		return (status);
	}

	if (sts->obstype == IMAGING_TYPE) {

	    /* Get output coordinate parameters and distortion; should be
		the same for all imsets.
	    */
	    if ((status = GetIDC (sts, &coords, &dist)))
		return (status);
	    minorder = 1;  maxorder = 1;	/* loop limits */
	}

	sts->first_order = (maxorder <= 1);	/* first order grating? */

	if (sts->x2dcorr_o != PERFORM) {
	    if (sts->x2dcorr_o == DUMMY) {
		printf (
		"ERROR    X2DCORR not performed due to dummy pedigree.\n");
	    } else {
		printf (
		"ERROR    X2DCORR not performed due to missing row.\n");
	    }
	    return (NOTHING_TO_DO);
	}

	for (extver = 1;  extver <= sts->nimages;  extver++) {

	    PrGrpBegin ("imset", extver);

	    initSingleGroup (in);
	    initSingleGroup (out);

	    /* The input data may contain more than one order. */
	    getSingleGroup (sts->input, extver, in);
	    if (hstio_err())
		return (OPEN_FAILED);

	    /* Get keyword values from extension header. */
	    if ((status = GetGrpInfo7 (sts, &in->sci.hdr)))
		return (status);

	    /* If exptime = 0, don't do fluxcorr. */
	    fluxcorr_extver = sts->fluxcorr;	/* local to current imset */
	    if (fluxcorr_extver == PERFORM) {
		if (sts->exptime <= 0.) {
		    printf (
	"Warning  FLUXCORR skipped for imset %d because exptime = %.6g\n",
			extver, sts->exptime);
		    fluxcorr_extver = SKIPPED;
		}
	    }

	    /* This has to be done inside the loop over imsets because it
		uses ltm2_2, which is an extension header keyword.  It must
		be done only once, because it modifies the disp coefficients.
	    */
	    if (extver == 1) {
		if (sts->obstype == SPECTROSCOPIC_TYPE && !sts->first_order) {
		    /* Adjust dispersion coefficients to account for offset. */
		    if ((status = MOCAdjustDisp (sts, disp)))
			return (status);
		}
	    }

	    /* Transpose the input data, if a cross disperser was used. */
/* ...
	    if (sts->dispaxis == 2) {
		if ((status = Transpose (&in)))
		    return (status);
	    }
... */

	    if (sts->obstype == SPECTROSCOPIC_TYPE && !sts->wavecal &&
		sts->wavecorr != COMPLETE && extver == 1)
		printf
		("Warning  Wavecal processing has not been performed.\n");

	    /* Get heliocentric correction factor. */
	    if (sts->heliocorr == PERFORM && !sts->wavecal) {
		printf ("\n");
		PrSwitch ("helcorr", PERFORM);
		sts->hfactor = HelioFactor (sts, &v_helio);
		/* Note:  written to input in memory, copied later to output */
		if ((status = Put_KeyD (&in->sci.hdr, "V_HELIO", v_helio,
                                        "heliocentric radial velocity (km/s)")))
		    return (status);
		PrSwitch ("helcorr", COMPLETE);
	    } else {
		sts->hfactor = 1.;
	    }

	    blazeshift = sts->blazeshift;

	    /* In order to apply the MSM/blaze correction, we must get
	       some information up front. The first item to get is the
	       reference spectral order number. It is gotten by scanning
	       the _pht table. Then we have to open the trace table (_1dt),
	       as well as poke into the dispersion solution structure (disp)
	       in order to retrieve the remaining information necessary
	       to run the algorithm.
	    */
	    mref = 0;
	    warn = 1;
	    photb.mref = mref;
	    if (sts->obstype  == SPECTROSCOPIC_TYPE &&
	        sts->fluxcorr == PERFORM) {

	        /* Get reference order number. */
	        for (i = minorder; i <= maxorder; i++) {
	            status = GetAbsPhot (sts, i, &photb, 0, &warn);
	            if (status == ROW_NOT_FOUND) {
	                FreePhot (&photb);
	                continue;
	            } else if (status) {
	                FreePhot (&photb);
	                return (status);
	            } else {
	                mref = photb.mref;
	                FreePhot (&photb);
	                break;
	            }
	        }
	        if (status)
	            return (status);

	        if (mref > 0) {

	            /* Get the trace for the reference order */
	            if ((status = GetTrace (sts, mref, &trace)))
	                return (status);
	            ypos = trace->a2center;
	            FreeTrace (&trace);

	            /* Add all pertinet offsets. */
	            ypos += sts->msm_slop[1] + sts->ap_offset[1];

	            /* Get center wavelength and dispersion for the reference
	               order. Note that if no dispersion solution is activated, no phot
	               correction will take place, so the wpos and ddisp values
	               will never get used.
	            */
	            if ((status = GetApOffset (sts, &slit, disp->ref_aper,
                                               &delta)))
	                return (status);
	            if ((status = GetInang (sts, &sts->inangtab, mref, &iac)))
	                return (status);
	            if ((status = InterpDisp (&disp, ypos, &disp_y)))
	                return (status);

	            warn1 = 0;  warn2 = 0;
	            AdjustDisp (sts, disp_y, delta, &iac, &warn1, &warn2);

	            if (strcmp (sts->opt_elem, "PRISM") == 0)
	                disp_type = PRISM_DISP;
	            else
	                disp_type = GRATING_DISP;

	            wpos = GetWavelength (disp_y, mref, 512.0, disp_type,
	                                  sts->cenwave - 100.,
	                                  sts->cenwave + 100.);
	            ddisp = wpos - GetWavelength (disp_y, mref, 511.0,
	                                  disp_type,
	                                  sts->cenwave - 100.,
	                                  sts->cenwave + 100.);
	            FreeDisp (&disp_y);
	        }
	    }

	    /* Do for each spectral order in input image. */
	    for (sporder = minorder;  sporder <= maxorder;  sporder++) {

		if (sts->obstype == SPECTROSCOPIC_TYPE) {
		    printf ("\n");
		    PrGrpBegin ("order", sporder);
		}

		/* Extract the CoordInfo record for this order. */
		if ((status = ReturnCoord (&coords, sporder, &coord_o)))
		    return (status);

		/* Scale some values (in-place in coord_o) according to
		   binning, and update plate_scale in sts.
		*/
		if ((status = ScaleWCS (sts,
                        in->sci.data.nx, in->sci.data.ny, coord_o)))
		    return (status);

		/* If center_target is true, reset CRPIX2 so the output
		   image will be shifted to center the target in Y.
		*/
		if (sts->first_order && sts->center_target)
		    SetCrpix2 (coord_o);

		/* Total offset = MSM offset + aperture. */
		AddOffsets (sts, &slit);

		/* Create output imset. */
		allocSingleGroup (out, coord_o->npix[0], coord_o->npix[1], True);
		if (hstio_err())
		    return (OUT_OF_MEMORY);

		/* Copy the headers from input to output.  For spectroscopic
		   mode, update order number and coordinates.  For imaging
		   mode, GeoCorr7 will update crpix.
		*/
		if ((status = PutGrpInfo7 (in, out,
                                           coord_o, extver, sts->obstype)))
		    return (status);

		/* Perform geometric correction.
		   (incl &ssgx, &ssgy in call for small-scale distortion)
		*/
		if (sts->obstype == SPECTROSCOPIC_TYPE) {
		    if ((status = X2DCorr7 (sts, coord_o, disp, &slit,
                                            o_extver, in, out))) {
			if (status > 0)
			    return (status);		/* a real error */
			status = 0;
			printf ("Warning  Skipping spectral order %d.\n",
					sporder);
			freeSingleGroup (out);
			continue;
		    }
		} else if (sts->obstype == IMAGING_TYPE) {
                        if ((status = GeoCorr7 (sts, &dist, &ssgx, &ssgy, in, out)))
                            return (status);
		}

		x2dMsg (sts, o_extver);

		if (sts->printtime)
		    TimeStamp ("2-D rectification complete", sts->rootname);

		/* Convert to absolute flux. */
		if (fluxcorr_extver == PERFORM) {
		    if ((status = GetAbsPhot (sts, sporder, &phot, 1, &warn)))
			return (status);
		    /* PCT is actually indpendent of sporder. */
		    if ((status = GetPCT (sts, &phot)))
			return (status);
		    FluxMsg (sts, o_extver);
		    if (sts->x2dcorr_o != PERFORM) {
			printf (
		"Warning  Skipping spectral order %d \\\n", sporder);
			if (sts->x2dcorr_o == DUMMY) {
			    printf (
			"Warning  due to dummy pedigree in PHOTTAB row.\n");
			} else {
			    printf (
			"Warning  due to missing row in PHOTTAB.\n");
			}
			freeSingleGroup (out);
			continue;
		    }

	            /* Store data required by MSM/blaze correction */
	            if (mref > 0) {
	                phot.wpos = wpos;
	                phot.ypos = ypos;
	                phot.disp = ddisp;
	            } else {
	                phot.wpos = 0.0;
	                phot.ypos = 0.0;
	                phot.disp = 0.0;
	            }

		    if ((status = AbsFlux (sts, out, coord_o, &phot, &slit, &tds,
				sts->plate_scale, sts->atodgain,
				sts->exptime, sts->hfactor, sporder,
                                &blazeshift)))
			return (status);
		    PrSwitch ("fluxcorr", COMPLETE);
		    if (sts->printtime)
			TimeStamp ("FLUXCORR complete", sts->rootname);
		}

		/* Compute statistics and add to output headers. */
		if (sts->statcorr == PERFORM) {
		    if ((status = doStat (out, sts->sdqflags)))
			return (status);
		    if (sts->printtime)
			TimeStamp ("STATFLAG complete", sts->rootname);
		}

		/* Write blazeshif info. */

                if (blazeshift != NO_VALUE) {
	            if ((status = Put_KeyF (&out->sci.hdr, "BLZSHIFT",
                                            (float) blazeshift,
                                            "Blaze shift (in pixels)"))) {
	                return (status);
	            }
	        }

		o_extver += 1;			/* output extver */

		/* Update primary header:  CAL_VER, FILENAME, and history. */
		if (o_extver == 1) {
		    UCalVer (out->globalhdr);
		    UFilename (sts->output, out->globalhdr);
		    if ((status = History7 (sts, out->globalhdr)))
			return (status);
		}

		if (sts->obstype == SPECTROSCOPIC_TYPE) {
		    PrSwitch ("x2dcorr", COMPLETE);
		    PrSwitch ("dispcorr", COMPLETE);
		    PrGrpEnd ("order", sporder);
		} else {
		    PrSwitch ("geocorr", COMPLETE);
		}

		putSingleGroup (sts->output, o_extver, out, 0);
		if (hstio_err()) {
		    printf ("ERROR    Couldn't write imset %d.\n", o_extver);
		    return (GENERIC_ERROR_CODE);
		}
		freeSingleGroup (out);

		FreePhot (&phot);	/* free memory for PhotInfo */

		if (sts->printtime)
		    TimeStamp ("Output written to disk", sts->rootname);
	    }

	    freeSingleGroup (in);
	    PrGrpEnd ("imset", extver);

	    /* Have we actually written any output for this imset? */
	    if (o_extver <= 0) {
		printf ("ERROR    No output written for any order.\n");
		return (NOTHING_TO_DO);
	    }
	}

	if (sts->verbose && sts->obstype == SPECTROSCOPIC_TYPE &&
                       sts->trace_rotation != 0.)
                   printf ("         trace was rotated by = %.6g degree.\n",
                       sts->trace_rotation);

	/* free memory */

	FreeThroughput (&slit);
	FreeTds (&tds);
	FreeCoord (&coords);
	FreeCoord (&coord_o);

	if (sts->obstype == IMAGING_TYPE)
	    FreeDist (&dist);		/* free distortion arrays */

	if (sts->sgeocorr == PERFORM) {
	    freeFloatHdrData (&ssgx);
	    freeFloatHdrData (&ssgy);
	}

	/* Update nextend in primary header, if the number has changed. */
	if (o_extver > extver-1) {	/* extver was incremented in for loop */
	    initHdr (&phdr);
	    oim = openUpdateImage (sts->output, "", 0, &phdr);
	    if (hstio_err())
		return (OPEN_FAILED);
	    if ((status = Put_KeyI (&phdr, "NEXTEND", o_extver * EXT_PER_GROUP,
                                    "number of extensions")))
		return (status);
	    putHeader (oim);
	    if (hstio_err())
		return (OPEN_FAILED);
	    closeImage (oim);
	    freeHdr (&phdr);
	    if (sts->printtime)
		TimeStamp ("NEXTEND updated", sts->rootname);
	}

	FreeDisp (&disp);

	free (in);
	free (out);

	/* Have we written any output for any imset? */
	if (o_extver <= 0) {
	    printf ("Warning  No output written for any imset.\n");
	    return (NOTHING_TO_DO);
	}

	return (0);
}

/* This routine modifies the dispersion coefficients in-place to account
   for the displacement of the image on the detector from the location
   that was used for measuring the coefficients.
   This only needs to be applied for echelle data.

   Note:  This uses extension header keywords ltm2_2 and shifta2, so it
   should be called within the loop over imsets.  However, it should only
   be called once, because it modifies the coefficients in-place.
*/

static int MOCAdjustDisp (StisInfo7 *sts, DispRelation *disp) {

	int mref;		/* spectral order number of reference order */
	double yref;		/* Y location of reference order */
	double a4corr;		/* correction factor */
	double ydiff;		/* Y position offset */
	double a2center;	/* from the trace for order mref */
	double r_shifta2;	/* shifta2 converted to reference pixel size */
	int status;

	SpTrace *trace;		/* spectral trace for spectral order mref */
	int GetTrace (StisInfo7 *, int, SpTrace **);
	void FreeTrace (SpTrace **);

	trace = NULL;

	/* Get MAMA offset parameters from the dispersion coeff. table. */
	if ((status = GetMOC (&sts->disptab, sts->opt_elem, sts->cenwave,
                              &mref, &yref, &a4corr)))
	    return (status);
	if (a4corr == 0.)
	    return (0);		/* nothing further to do */

	/* Get the trace for spectral order mref, and get its a2center. */
	if ((status = GetTrace (sts, mref, &trace)))
	    return (status);
	a2center = trace->a2center;

	/* Here's what we're doing:
		ydiff = ypos - yref
		shifta2 = ypos - a2center, and shifta2 is called msm_slop[1]
	   so:
		ydiff = shifta2 + a2center - yref
	*/
	r_shifta2 = sts->msm_slop[1] / sts->ltm[1];	/* reference pixels */
	ydiff = r_shifta2 + a2center - yref;

	/* Adjust dispersion coefficients using a4corr. */
	disp->coeff[0] -= ydiff * sts->cenwave * a4corr;
	disp->coeff[4] += ydiff * a4corr;

	FreeTrace (&trace);

	return (0);
}


/* This routine opens the two images in a small-scale geometric
   correction reference file.
*/

static int OpenSGeo (char *sdstfile, FloatHdrData *ssgx, FloatHdrData *ssgy) {

/* arguments:
char *sdstfile          i: name of file to open
FloatHdrData *ssgx      o: small-scale distortion in X
FloatHdrData *ssgy      o: small-scale distortion in Y
*/

	initFloatHdrData (ssgx);
	initFloatHdrData (ssgy);

	getFloatHD (sdstfile, "AXIS1", 1, ssgx);
	if (hstio_err())
	    return (OPEN_FAILED);
	getFloatHD (sdstfile, "AXIS2", 1, ssgy);
	if (hstio_err())
	    return (OPEN_FAILED);

	return (0);
}

static void x2dMsg (StisInfo7 *sts, int o_extver) {

	printf ("\n");
	if (sts->obstype == SPECTROSCOPIC_TYPE) {
	    PrSwitch ("x2dcorr", PERFORM);
	    PrSwitch ("dispcorr", PERFORM);
	} else {
	    PrSwitch ("geocorr", PERFORM);
	}

	if (o_extver == 0) {

	    if (sts->obstype == SPECTROSCOPIC_TYPE) {

		PrRefInfo ("apdestab", sts->apdestab.name,
			sts->apdestab.pedigree,
			sts->apdestab.descrip, sts->apdestab.descrip2);

		PrRefInfo ("sdctab", sts->distntab.name,
			sts->distntab.pedigree,
			sts->distntab.descrip, sts->distntab.descrip2);

		PrRefInfo ("disptab", sts->disptab.name,
			sts->disptab.pedigree,
			sts->disptab.descrip, sts->disptab.descrip2);

		PrRefInfo ("inangtab", sts->inangtab.name,
			sts->inangtab.pedigree,
			sts->inangtab.descrip, sts->inangtab.descrip2);

		PrRefInfo ("sptrctab", sts->sptrctab.name,
			sts->sptrctab.pedigree,
			sts->sptrctab.descrip, sts->sptrctab.descrip2);
		if (sts->wx2dcorr == COMPLETE)
		    printf ("SPTRCTAB trace was applied earlier, by wx2d\n");

	    } else {

		PrRefInfo ("idctab", sts->distntab.name,
			sts->distntab.pedigree,
			sts->distntab.descrip, sts->distntab.descrip2);
	    }
	}
}

static void FluxMsg (StisInfo7 *sts, int o_extver) {

	printf ("\n");
	PrSwitch ("fluxcorr", PERFORM);

	if (o_extver == 0) {

	    PrRefInfo ("phottab", sts->phottab.name,
			sts->phottab.pedigree,
			sts->phottab.descrip, sts->phottab.descrip2);

	    PrRefInfo ("apertab", sts->apertab.name,
			sts->apertab.pedigree,
			sts->apertab.descrip, sts->apertab.descrip2);

	    if (sts->pctcorr == PERFORM) {
		PrRefInfo ("pctab", sts->pctab.name,
			sts->pctab.pedigree,
			sts->pctab.descrip, sts->pctab.descrip2);
	    } else {
		printf (
		"Warning  PCTAB correction was not included in DIFF2PT.\n");
	    }
	    if (sts->tdscorr == PERFORM) {
		PrRefInfo ("tdstab", sts->tdstab.name,
			             sts->tdstab.pedigree,
			             sts->tdstab.descrip,
			             sts->tdstab.descrip2);
	    } else {
		printf ("Warning  TDS correction not performed.\n");
	    }
	}
}

/* This routine should only be called for first-order data, never for echelle
   data.

   For FUV G140L and G140M spectra, the target has at different times been
   offset to a location either above or below the middle of the detector,
   to avoid the repeller wire.  Normally, the A2CENTER and CRPIX2 values in
   the SDC table result in the input image being copied to approximately the
   center of the 2-D rectified output image.  This routine can be called to
   modify CRPIX2 so that the image will be shifted vertically as it is copied
   into the output, so that the target ends up in the middle of the output
   image.  This may be desirable because the target would then be at the same
   location in the rectified image regardless of when the data were taken.
   For an extended target or a wavecal, however, some of the data could be
   shifted off the rectified image and be lost.

   Regardless of whether this routine is called or not, the location in the
   cross-dispersion direction of the target in the output, 2-D rectified
   image is given by the CRPIX2 keyword in the extension headers.  The
   exception to this rule is that POSTARG is not taken into account, so if
   POSTARG2 is non-zero, the target will be offset from CRPIX2.
*/

static void SetCrpix2 (CoordInfo *coord_o) {

	coord_o->crpix[1] = (coord_o->npix[1] + 1.) / 2. - 1.;	/* 600 */
}

/* This routine adds together the offsets:  MAMA "dither", MSM offset,
   and aperture.  The aperture offset will be converted from arcseconds
   to pixels.  Both the total and aperture offsets will be saved in sts.
   Note that the "pixel" referred to here is reference pixel size.

   For imaging data, the total offset is set to zero.
*/

static void AddOffsets (StisInfo7 *sts, ApInfo *slit) {

/* arguments:
StisInfo7 *sts      i: calibration switches and info
ApInfo *slit        i: description of slit (aperture offset)
*/

	if (sts->obstype == SPECTROSCOPIC_TYPE) {

	    sts->ap_offset[0] = slit->ap_offset[0] / sts->plate_scale[0];
	    sts->ap_offset[1] = slit->ap_offset[1] / sts->plate_scale[1];

	    /* Actually, sts->total_offset[0] is never used. */
	    sts->total_offset[0] = sts->msm_slop[0] + sts->ap_offset[0];

	    sts->total_offset[1] = sts->msm_slop[1] + sts->ap_offset[1];

	} else {

	    sts->total_offset[0] = 0.;
	    sts->total_offset[1] = 0.;

	    sts->ap_offset[0] = 0.;
	    sts->ap_offset[1] = 0.;
	}
}

/* This routine frees memory for distortion data, if it has been
   allocated, and resets the flag to indicate that memory is no
   longer allocated.
*/

static void FreeDist (DistInfo *dist) {

	if (dist->allocated) {
	    free (dist->xcoeff);
	    free (dist->ycoeff);
	    dist->allocated = 0;
	}
}

/* This routine frees memory for photometry data. */

static void FreePhot (PhotInfo *phot) {

	if (phot->allocated) {
	    free (phot->wl);
	    free (phot->thru);
	    if (phot->pcorr != NULL) {
		free (phot->pcorr);
		phot->pcorr = NULL;
	    }
	    phot->allocated = 0;
	}
}

/* This routine frees memory for throughput of slit. */

static void FreeThroughput (ApInfo *slit) {

	if (slit->allocated) {
	    free (slit->wl);
	    free (slit->thr);
	    slit->allocated = 0;
	}
}
