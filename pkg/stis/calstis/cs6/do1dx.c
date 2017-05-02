# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <math.h>		/* fabs */

# include "xtables.h"
# include "hstio.h"

# include "stis.h"
# include "calstis6.h"
# include "hstcalerr.h"
# include "stisdef.h"
# include "stisdq.h"
# include "stispht.h"
# include "stistds.h"

static int OpenSGeo (char *, FloatHdrData *, FloatHdrData *);
static void warnDummy (char *, int, int);
static void PrSwitch6 (StisInfo6 *, char *, int);
static int MOCAdjustDisp (StisInfo6 *, DispRelation *);


/*
   Extract 1-D spectra from the input image, creating a 3-D table at
   the output. Each extver is handled separately, generating a separate
   table extension at the output. Individual spectral orders in echelle
   data are stored in contiguous rows in each output table extension. The
   output's primary header is inherited from the input primary header and
   updated with history information.


   Revision history:
   ----------------
   20 Feb 97  -  Implemented (I.Busko)
   10 Apr 97  -  Changes after code review (IB):
                 - replaced literal by LOG_LINSIZE constant.
                 - added TimeStamp6 function.
   22 Apr 97  -  Do not write empty spectrum - turned off for the moment (IB)
   25 Apr 97  -  New routine GetMAMAOff6 (IB)
   01 May 97  -  Warns if no output was written (IB)
   01 May 97  -  Warns if (command-line) A2CENTER and SPORDER do not agree (IB)
   02 May 97  -  Skip sporders that have a DUMMY pedigree in any table (IB)
   02 May 97  -  Physically writes the output's primary header (IB)
   08 May 97  -  Skip orders whose nominal A2CENTER is off-image (IB)
   09 May 97  -  Conform to new _trl standard, no need for LOG_LINSIZE (IB)
   29 May 97  -  Changed arguments in WriteRow (IB)
   08 Jul 97  -  Turn off linked list in trace structure (IB)
   17 Jul 97  -  Report heliocentric correction in HISTORY keyword (IB)
   17 Jul 97  -  Zero SHIFTA2 and MOFFSET2 if A2CENTER comes from c. line (IB)
   17 Jul 97  -  A2CENTER from command line is in physical coordinates (IB)
   28 Jul 97  -  Added logic to avoid scanning orders outside the image (IB)
   28 Jul 97  -  Re-activate trace linked list; add interpolation (IB)
   31 Jul 97  -  Skip order if row in reference table not found (IB)
   04 Aug 97  -  Fix error in A2CENTER conversion (IB)
   03 Sep 97  -  Add ATODGAIN correction (IB)
   26 Sep 97  -  Re-interpolate in trace table after crosscor (long-slit) (IB)
   08 Oct 97  -  Update swicth values in output header (IB)
   14 Nov 97  -  Add POSTARG corrections (IB)
   08 Dec 97  -  Remove POSTARG corrections (IB)
   12 Dec 97  -  Skip test for out of bounds if not echelle grating (PH)
   22 Dec 97  -  Report extraction position in _trl file (IB)
   27 Jan 98  -  PCTAB support (IB)
   03 Feb 98  -  1-indexed cl_a2center (IB)
   13 Apr 98  -  Changed arguments in WriteRow (IB)
   22 Jun 98  -  Global croscor fit mode (IB)
   23 Jun 98  -  Profile generator (IB)
   29 Jun 98  -  Rejection ranges in profile generator (IB)
   10 Jul 98  -  Fixed reported extraction box position (IB)
   20 Jul 98  -  Always get aperture description (because of A2 offset, IB)
   28 Jul 98  -  Fixed crosscor report messages (IB)
   28 Jul 98  -  Reports extraction position only if inside the data array (IB)
   17 Sep 98  -  Optimal extraction (IB)
   23 Sep 98  -  Output weights image (IB)
   28 Oct 98  -  OPR 37810: extr. algorithm taken from first order (IB)
   28 Oct 98  -  Write XTRACALG in output primary header (IB)
   19 Nov 98  -  Use SDQFLAGS for both optimal extr. and profile building (IB)
   17 Dec 98  -  Avoid geocoronal Lya (IB)
   15 Sep 99  -  Handle upper and lower case calibration swith values (IB)
   20 Oct 99  -  Close auxiliary trace lists (IB)
   17 Dec 99  -  Flux normalization factor in opt. extr. rejection (IB)
   06 Jan 00  -  Scatterd light correction algorithm (IB)
   17 Feb 00  -  IMSET selection (IB)
   14 Apr 00  -  Verbosity control (IB)
   14 Jun 00  -  Turn global crosscor off if first-order (IB)
   13 Jul 00  -  Redefined extrsize as double (IB)
   25 Jul 00  -  Fixed updating of switches in output header (IB)
   04 Aug 00  -  UpdHdrSwitch moved to its own source file (IB)
   17 Aug 00  -  Fix flux computation in hires data (IB)
   22 Sep 00  -  Turn off fluxcorr if dummy pedigree in pht entry (IB)
   01 Nov 00  -  Profile offsets (IB)
   08 Nov 00  -  Do not read MOFFTAB (B)
   01 Dec 00  -  Subsampling in profile builder (IB)
   16 Feb 01  -  Ignore small blemish in profile builder/opt. extraction (IB)
   12 Apr 01  -  OPR 43663: skip IMSET (on 1st order) if crosscor fails (IB)
   10 May 01  -  Move GetDisp6 out of the loop on sporder (IB)
   10 May 01  -  Add MOCAdjustDisp (IB, function provided by P.Hodge)
   08 Jan 02  -  Time-dependent sensitivity (IB)
   05 Feb 02  -  MSM/blaze correction (IB)
   17 Apr 02  -  BLZSHIFT keyword added to output extension header (IB)
   05 Aug 02  -  Smooth background (IB)
   07 Nov 02  -  Fix problem with Ly alpha detection in prism data (IB)
   17 Dec 02  -  Slitless extraction (IB)
   17 Jan 03  -  Change calling sequence of HelioFactor6, update v_helio (PEH)
   24 Feb 03  -  Include units for v_helio in comment field for keyword (PEH)
   16 Jun 03  -  Added function to get CTE corr info from CCD table (PB)
    1 Jul 04  -  Initialize ypos to 0.0 (PEH)
   21 Jul 04  -  If ctecorr was done, set the header keyword to COMPLETE;
                 replace UpdHdrSwitch was a call to Put_KeyS (PEH)
   27 Dec 04  -  Add explicit braces to avoid ambiguity in the
                 "Print background switch status" section (PEH)
   08 Apr 05  -  Initialize slit.gac_allocated to 0; call GetGAC6 (PEH)
   21 Apr 05  -  Rename variable stpos to xoffset, and apply it directly
                 as an offset rather than assuming the nominal location
                 was 511.5 (PEH)
   12 Sep 05  -  Don't divide xoffset by ltm[0], i.e. interpret it as an
                 offset in reference pixels rather than in image pixels (PEH)
   10 Apr 06  -  Modify the section to determine wpos:  get the wavelength
                 at pixel 512 (zero indexed, if npts=1024), and correct
                 for on-board Doppler shift (PEH)
   11 May 06  -  The value written to keyword BLZSHIFT is now the average
                 blaze shift over all orders in the output table (PEH)
    3 Nov 08  -  Call checkImsetOK to get imset_ok, and skip imset if F;
                 pass extver instead of o_ext to CreaTable, since it's really
                 the input imset, not the output.
   21 Oct 11  -  Set the output flux and error to zero if the relevant row in
		 the phottab has pedigree = DUMMY (PEH)
*/

int Do1Dx (StisInfo6 *sts, Hdr *phdr) {

/* argument:
StisInfo6 *sts    i: calibration switches and info
*/

	int status;

	TblDesc otable;		/* pointer to output table descriptor */
	RowContents row_contents; /* row contents (values written into) */

	ProfTblDesc optable;	/* pointer to profile table descriptor */

	SingleGroup in;		/* input data */
	IODescPtr im;		/* descriptor for output primary header */
	SingleGroup outw;	/* output weight image in optimal extr. */
	SingleGroup holdraw;	/* raw input image holding buffer */

	XtractInfo *extract;	/* list of extraction parameter records */
	XtractInfo *extract_o;	/* one record extracted from extract */
	XtractInfo *extract_a;	/* one record extracted from extract */
	SpTrace *trace;		/* list of spectrum traces */
	SpTrace *trace1;	/* auxiliary list of spectrum traces */
	SpTrace *trace2;	/* auxiliary list of spectrum traces */
        SpTrace *trace_y;       /* spectrum trace interpolated at y */
	PhotInfo phot;		/* for conversion to absolute flux */
	PhotInfo photb;		/* for MSM/blaze correction */
	PhotInfo photc;		/* for extraction box height correction */
	TdsInfo tds;		/* for time-dependent sensitivity correction */
	DispRelation *disp;	/* dispersion relation */
	DispRelation *disp_y;	/* dispersion relation interpolated at y */
        CoordInfo *coords;      /* list of plate scale records */
        CoordInfo *coord_o;     /* one record extracted from coords */
	ApInfo slit;		/* description of slit */
	InangInfo iac;		/* incidence-angle correction coeff */
	IntensArray inta;	/* intensity array for optimal extraction */
	ProfileArray *profa;	/* profile array list for opt. extraction */
	FloatHdrData ssgx;	/* small-scale distortion in X */
	FloatHdrData ssgy;	/* small-scale distortion in Y */
        CTICorrInfo cti;        /* for CCD CTE correction */

	int extver;		/* loop index for extension version number */
	int o_ext;		/* extension number in output file */
	int o_row;		/* current row in current output extension*/
	int minorder, maxorder;	/* min & max spectral order */
	int maxsearch;		/* maximum range in cross correlation */
	int a2phys;		/* nominal A2CENTER in image coordinates */
	int trc_status;         /* flag for trace status */
	int trc_status1;        /* flag for trace status */
	int trc_status2;        /* flag for trace status */
	int sdc_status;         /* flag for coordinate info status */
	int opt_status;         /* flag for optimal extraction */
	int outw_exists;        /* flag for weights image */
	int phu_ready;		/* is the output image PHU written already ?*/
	int phuw_ready;		/* is the weights image PHU written already?*/
	int ext_ready;		/* is current extension's header written ? */
	int success;		/* at least one order successfully extrac. ? */
	int skipping;		/* skipping this order ? */
	int imset_ok;		/* value of header keyword IMSET_OK */
	int first_extver;	/* true when processing first imset */
	double extrsize;	/* size of spectrum extraction box */
	int i, j;
	double bksize[2];       /* sizes of background boxes */
	double bkoffset[2];     /* offsets of background boxes */
	double avcrscroff;	/* average cross correlation offset */
	double delta = 0.0;	/* X offset (arcsec) from slit used to measure
				   dispersion coefficients */
	double cl_a2center;	/* read from command line (for local use) */
	int norder;		/* index into croscor arrays */
	int ccstatus;		/* status of croscor finding process */
	int ilow_end;		/* image region to be restored when applying */
	int ihigh_end;		/* the scattered light correction algorithm */
	int extver1, extver2;	/* IMSETs to be processed */
	int mref;		/* MSM/blaze correction reference order */
	double ypos = 0.;	/* MSM/blaze correction reference position */
	double wpos = 0.;	/* MSM/blaze correction reference wavelength */
	double ddisp;		/* MSM/blaze correction reference dispersion*/
	int bcorr_l;		/* MSM/blaze correction local switch */
	double blazeshift;	/* MSM/blaze correction shift value */
	double sum_blazeshift;	/* sum of blaze shift values */
	double n_blazeshift;	/* number of blaze shift values in the sum */
	double d1, d2, radvel;
	double hold;
	int dum, warn;
	char str[50];		/* temporary string area */

	int AbsFlux6 (StisInfo6 *, RowContents *, PhotInfo *,
                      PhotInfo *, ApInfo *, TdsInfo *, CTICorrInfo *,
                      double, double, int, double *, double*);
	void AddOffsets6 (StisInfo6 *, ApInfo *);
	void AdjustDisp6  (StisInfo6 *, DispRelation *, double, InangInfo *);
	int AllocOutArrays (RowContents *);
	int AllocProfile (StisInfo6 *, int, double);
	int BuildOptProf (StisInfo6 *, ProfileArray **);
	int CreaProfTable (StisInfo6 *, ProfTblDesc *);
	int CreaTable (StisInfo6 *, int, TblDesc *);
	int CrossCorr (StisInfo6 *, SpTrace *, SingleGroup *, FloatHdrData *,
                       FloatHdrData *, int, int, int);
	int DefineBackRegions (StisInfo6 *, SpTrace *, CoordInfo *,
                               SingleGroup *, FloatHdrData *, FloatHdrData *,
                               int *, int *);
	void FindLya (StisInfo6 *, SingleGroup *, ApInfo *, double *,
                      int, int *, int *, int *, int *);
	void FreeDisp6 (DispRelation **);
	void FreeInang6 (InangInfo *);
	void FreeIntensity (IntensArray *);
	void FreeOutArrays (RowContents *);
	void FreePhot6 (PhotInfo *);
	void FreeProfile (StisInfo6 *, int);
	void FreeProfileArray (ProfileArray **);
	void FreeTds (TdsInfo *);
	void FreeTrace6 (SpTrace **);
	void FreeThroughput6 (ApInfo *);
	void FreeXtract (XtractInfo **);
	void FreeCoord6 (CoordInfo **);
	int GetAbsPhot6 (StisInfo6 *, int, PhotInfo *, int, int *);
	int GetApDes6 (StisInfo6 *, ApInfo *);
	int GetApOffset6 (StisInfo6 *, float *, char *, double *);
	int GetApThr6 (StisInfo6 *, ApInfo *);
	int GetGAC6 (StisInfo6 *, ApInfo *);
        int GetCCDTab6 (StisInfo6 *, CTICorrInfo *);
	int GetDisp6 (StisInfo6 *, DispRelation **);
	int GetGrpInfo6 (StisInfo6 *, Hdr *);
	int GetInang6 (StisInfo6 *, RefTab *, int, InangInfo *);
	int GetIntens (StisInfo6 *, int, IntensArray *);
	int GetPCT6 (StisInfo6 *, PhotInfo *, double, int);
	int GetProfile (StisInfo6 *, int, ProfileArray **);
        int GetSDC6 (StisInfo6 *, CoordInfo **, int *, int *);
	int GetTds (char *, char *, TdsInfo *);
        int GetTrace6 (StisInfo6 *, int, SpTrace **);
	int GetXtract (StisInfo6 *, XtractInfo **, int *, int *);
	int GCrossCorr (StisInfo6 *, SingleGroup *, int, int, double *, int,
	                double *);
	double HelioFactor6 (StisInfo6 *, double *);
	int HistoryAddHeader (char *, char *);
	void InitProfTblDesc (ProfTblDesc *);
	void InitRowContents (RowContents *);
	void InitTblDesc (TblDesc *);
	int InterpDisp6 (DispRelation **, double, DispRelation **);
	int InterpTrace6 (SpTrace **, double, SpTrace **);
	int Lee (StisInfo6 *, RowContents *);
	void Message6 (StisInfo6 *, int);
        int ReturnCoord6 (CoordInfo **, int, CoordInfo **);
	int ReturnXtract (XtractInfo **, int, XtractInfo **);
	int SelectAlg (StisInfo6 *, XtractInfo *);
	void SetRanges (StisInfo6 *, RowContents *);
	int SmoothBack (StisInfo6 *, RowContents *);
	void TimeStamp6 (int, char *, char *);
	int doubleUpdateHeader (char *, char *, double, char *, int);
	int intUpdateHeader (char *, char *, int, char *, int);
	int strUpdateHeader (char *, char *, char *, char *, int);
	void Wave (StisInfo6 *, DispRelation *, RowContents *);
	int WriteProfRow (StisInfo6 *sts, ProfTblDesc *, RowContents *, int *);
	int WriteRow (StisInfo6 *sts, TblDesc *, RowContents *, int *);
	int X1DSpec (StisInfo6 *, SpTrace *, XtractInfo *, double,
                     SingleGroup *, SingleGroup *, FloatHdrData *,
                     FloatHdrData *, IntensArray *, RowContents *);

	/* Set flags to indicate that memory has not been allocated yet. */
	phot.allocated  = 0;
	photb.allocated = 0;
	photc.allocated = 0;
	slit.allocated  = 0;
	slit.gac_allocated  = 0;
	tds.allocated   = 0;
	iac.allocated   = 0;
	inta.allocated  = 0;
	extract     = NULL;
	extract_o   = NULL;
	extract_a   = NULL;
	trace       = NULL;
	trace1      = NULL;
	trace2      = NULL;
	trace_y     = NULL;
	disp        = NULL;
	disp_y      = NULL;
        coords      = NULL;
        coord_o     = NULL;
	phot.pcorr  = NULL;
	photb.pcorr = NULL;
	photc.pcorr = NULL;
	profa       = NULL;

     	/* Get aperture description. */
	if ((status = GetApDes6 (sts, &slit)))
	    return (status);

	/* Get aperture throughput info and grating-aperture correction. */
	if (sts->fluxcorr == PERFORM) {
	    if ((status = GetApThr6 (sts, &slit)))
		return (status);
	    if (sts->gaccorr == PERFORM) {
		if ((status = GetGAC6 (sts, &slit)))
		    return (status);
	    }
	}

	/* Get time-dependent sensitivity info. */
	if (sts->tdscorr == PERFORM) {
	    status = GetTds (sts->tdstab.name, sts->opt_elem, &tds);
            if (status == DUMMY || status == ROW_NOT_FOUND)
                sts->tdscorr = OMIT;
            else if (status)
		return (status);
	}

	/* Read the small-scale distortion data into memory. */
	if (sts->sgeocorr == PERFORM) {
	    /* Read the small-scale distortion data into memory. */
	    if ((status = OpenSGeo (sts->sdstfile.name, &ssgx, &ssgy)))
		return (status);
	}

	/* Get output coordinate parameters for all spectral orders.
           This is necessary for the application of the H&S scattered light
           correction algorithm. In this algorithm the background
           extraction boxes read from the extraction table are ignored
           and the image coordinates taken from the sdc table are used
           instead to define a image section where the background analysis
           will take place. Notice that we always have to open this file,
           even if the algorithm is not going to be run. This is because
           at this point we don't know the selected algorithm yet. Values
           assigned to minorder and maxorder will be overwritten at the
           next call to GetXtract.
         */
	if (GetSDC6 (sts, &coords, &minorder, &maxorder))
	    return (status);

	/* Get 1-D extraction parameters for all spectral orders.
           The extraction table is the main driver that controls the
           scanning of spectral orders in echelle data.
        */
	if ((status = GetXtract (sts, &extract, &minorder, &maxorder)))
	    return (status);
        if (sts->x1d_o == DUMMY) {
	    printf ("ERROR    DUMMY pedigree entry in extraction table.\n");
	    return (status);
	}

        /* Get CTE correction parameters from CDD table. */
        if (sts->detector == CCD_DETECTOR && sts->ctecorr == PERFORM &&
            (status = GetCCDTab6(sts, &cti))) {
            if (status == COLUMN_NOT_FOUND) {
                printf(
                "Warning  Column not found in CCDTAB. Skipping CTECORR\n");
                sts->ctecorr = OMIT;
            } else
                return (status);
        }

	/* Print extraction and trace reference file info. This is done
           before any actual use of the data, because if an error happens
           during reference file processing, the file name is documented
           in the _trl file before the program dies.
        */
	if (sts->verbose == 1 || sts->verbose == 2) {
	    Message6 (sts, XTRAC_INFO);
	    printf ("\n");
	}
	sts->echelle = (maxorder > 1);

	/* Turn global crosscor mode off if first-order. */

	if (!(sts->echelle))
	    sts->cc_global = 0;

	/* Override spectral order info with command-line parameter. */
	if (sts->sporder != NO_ORDER) {
	    minorder = sts->sporder;
	    maxorder = sts->sporder;
	}

	/* Alloc arrays to store croscor offsets. */
	sts->cc_size  = maxorder - minorder + 1;
	sts->cc_off   = (double *) calloc (sts->cc_size, sizeof (double));
	sts->cc_spord = (int *)    calloc (sts->cc_size, sizeof (int));
	sts->cc_rej   = (int *)    calloc (sts->cc_size, sizeof (int));

	/* Loop over IMSETs in input image, or process just the
           selected one.
        */
	phu_ready = 0;
	phuw_ready = 0;
	o_ext = 0;
	if (sts->imset == 0) {
	    extver1 = 1;
	    extver2 = sts->nimages;
	} else {
	    extver1 = sts->imset;
	    extver2 = sts->imset;
	}
	first_extver = 1;		/* true */
	for (extver = extver1;  extver <= extver2;  extver++) {

	    photb.mref  = 0;

	    if (sts->verbose == 1 || sts->verbose == 2)
	        PrGrpBegin ("imset", extver);

	    status = checkImsetOK (sts->input, extver, &imset_ok);
	    if (status < 0)	/* imset 'extver' not present in input file */
		continue;
	    else if (status > 0)
                return status;
	    if (!imset_ok) {
		printf ("Warning  imset %d skipped (IMSET_OK = F)\n",
			extver);
		continue;
	    }

	    /* Bump output extension number. */
	    o_ext++;

	    /* Initalize I/O areas. */
	    initSingleGroup (&in);
	    initSingleGroup (&outw);
	    initSingleGroup (&holdraw);
	    outw_exists = 0;
	    InitTblDesc (&otable);
	    InitRowContents (&row_contents);
	    if (sts->do_profile)
	        InitProfTblDesc (&optable);

	    /* Read current IMSET. It may contain more than one order. */
	    getSingleGroup (sts->input, extver, &in);
	    if (hstio_err())
		return (OPEN_FAILED);
	    if (sts->verbose == 1 || sts->verbose == 2)
	        printf ("         Input read into memory.\n");

	    /* Get keyword values from image extension header. */
	    if ((status = GetGrpInfo6 (sts, &in.sci.hdr)))
		return (status);


	    /* Add A1 offset entered from command line (this should be
               used for slitless data only).
	    */
	    if (sts->xoffset != 0.) {
		sts->msm_offset[0] += sts->xoffset;
	        if (sts->verbose == 1 || sts->verbose == 2) {
		    printf (
	"         Offset of %g low-res pixels added in dispersion direction\n",
			sts->xoffset);
		}
	    }

	    /* Alloc memory for output data. */
	    otable.array_size = in.sci.data.nx;
	    row_contents.npts = in.sci.data.nx;
	    if ((status = AllocOutArrays (&row_contents))) {
		printf ("ERROR    Cannot allocate output arrays.\n");
	        return (status);
	    }
	    if (sts->do_profile) {
	        /* This gets extrsize from first entry in extract table. */
	        if (sts->extrsize == NO_SIZE)
	            extrsize = extract->extrsize;
	        else
	            extrsize = sts->extrsize;

                /* Translate into image pixels and add 1 to account
                   for extra profile pixel.
                */
	        /* optable.array_size += 1; avoid that for now */
	        /* optable.array_size  = extrsize * sts->ltm[1]; */

	        /* notice that the definition in here MUST be the same
                   used in AllocProfile to define the subpixeled profile
                   array size. This code duplication cannot be avoided
                   due to the convoluted way the optimal extraction code
                   was developed.
                */
	        optable.array_size  = (int)sts->subscale *
                                      (extrsize * sts->ltm[1] + PROF_MARGIN);

	        /* This is to accomodate the new profile format where the
                   entire uncompressed profile array is written.
	        */
	        /* optable.array_size  = (int)(extrsize * sts->ltm[1] *
                                         in.sci.data.nx);  */

	        optable.array_size_off  = in.sci.data.nx;
	    }

	    /* Perform actions triggered by the fact A2CENTER was entered
               from the command line. This must be done for each individual
               IMSET because values such as ltm and ltv may change from
               IMSET to IMSET.
               - don't add SHIFTA2 to A2CENTER;
               - don't add MOFFSET2 to A2CENTER;
               - disable POSTARG2 (spatial) correction;
               - convert A2CENTER from physical to reference coordinates,
                 since the original code was developed using A2CENTER
                 always in reference coordinates.
            */
	    if (sts->cl_a2center != NO_POSITION) {
	        sts->msm_offset[1]    = 0.0;
	        sts->dither_offset[1] = 0.0;
	        sts->pos_targ[1]      = 0.0;
	        cl_a2center = (sts->cl_a2center - sts->ltv[1]) /
                               sts->ltm[1];
	    } else
	        cl_a2center = sts->cl_a2center;

	    /* These effectively disable any POSTARG correction (12/8/97) */
	    sts->pos_targ[0] = 0.0;
	    sts->pos_targ[1] = 0.0;

	    /* Get heliocentric correction factor. */
	    if (sts->heliocorr == PERFORM)
		sts->hfactor = HelioFactor6 (sts, &radvel);

	    /* Set pixel rejection flag to be used with global crosscor.
               We must first run SelectAlg to properly select the
               extraction algorithm. Notice that this only works because
               of the new requirement imposed by OPR 37810.
            */
	    if ((status = ReturnXtract (&extract, minorder, &extract_a)))
		return (status);
	    if ((status = SelectAlg (sts, extract_a))) {
	        FreeXtract (&extract_a);
		return (status);
	    }
	    FreeXtract (&extract_a);
            /* Ignore small blemish (dust mote) - 16/2/01 IB */
	    if (sts->optimal || sts->do_profile)
	        sts->sdqflags = sts->sdqflags_orig & ~(SMALLBLEM);
	    else
	        sts->sdqflags = 0;

	    /*  This is being moved from the interior of the sporder loop
                to here. The reason is that the dispersion solution
                coefficients shouldn't be selected based on sporder.
            */

	    /*  Read dispersion solution. For echelle data, apply
                correction in the A2 direction.
            */
	    if (sts->dispcorr == PERFORM) {
	        if ((status = GetDisp6 (sts, &disp)))
	            return (status);
	        if (sts->echelle) {

	            /* Could this be a problem for the MSM/blaze correction ? */

	            if ((status = MOCAdjustDisp (sts, disp)))
	                return (status);
	        }
	    }

	    /* In order to apply the MSM/blaze correction, we must get
	       some information up front. The first item to get is the
	       reference spectral order number. It is gotten by scanning
	       the _pht table. The reference order number must be passed
	       to the global crosscor function so it can find precisely
	       the A2 position of the reference order in the input image.
	    */
	    mref = 0;
	    bcorr_l = OMIT;
	    warn = 1;
	    if (sts->fluxcorr == PERFORM && sts->echelle) {
	        for (i = minorder; i <= maxorder; i++) {
	            status = GetAbsPhot6 (sts, i, &photb, 0, &warn);
	            if (status == ROW_NOT_FOUND) {
	                FreePhot6 (&photb);
	                continue;
	            } else if (status) {
	                FreePhot6 (&photb);
	                return (status);
	            } else {
	                mref = photb.mref;
	                bcorr_l = photb.blazecorr;
	                FreePhot6 (&photb);
	                break;
	            }
	        }
	        if (status)
	            return (status);
	    }

	    blazeshift = sts->blazeshift;	/* could be NO_VALUE */
	    sum_blazeshift = 0.;	/* initial values for taking average */
	    n_blazeshift = 0.;

	    /* Now this is a kludge that has to do with the Ly alpha
               finding algorithm. The original design of this algorithm
               breaks when faced with prism data, since these do not carry
               wavelength info in their CD WCS header keywords. Therefore
               an actual wavelength scale must be provided to the Ly alpha
               finding routine. Since we are going to compute one anyways
               for the blaze shift algorithm (which applies only to echelle
               data), we fake here that there *is* a reference spectral
               order even in the prism case, so to force the wavelength
               computation step to run in that case too.
             */
	    if (strcmp (sts->opt_elem, "PRISM") == 0) {
	        mref = 1;
		ypos = cl_a2center;
	    }

	    /* Get global crosscor offset. This only makes sense if
               more than one order is being extracted.
            */
	    if (minorder != maxorder) {
	        if ((status = GCrossCorr (sts, &in, minorder, maxorder,
                                          &(sts->gcrscroff), mref, &ypos))) {
	            printf ("ERROR    Cannot compute global offset.\n");
	            return (status);
	        }
	    } else
	        sts->gcrscroff = 0.0;

	    /* This image acts as a holding buffer that stores the raw
               pixels of the input image when the scatterd light algorithm
               is selected. The algorithm smooths the background in the input
               image and this image is restored to its former pristine
               condition with data taken from this buffer, since there might
               be some overlap from neighboring orders. The error array must
               be saved as well.
            */
	    if (sts->scatter) {
	        if (allocSingleGroup (&holdraw, in.sci.data.nx,
                                      in.sci.data.ny, True) == ALLOCATION_PROBLEM)
	            return (OUT_OF_MEMORY);
	        for (i = 0; i < holdraw.sci.data.nx; i++) {
	            for (j = 0; j < holdraw.sci.data.ny; j++) {
	                Pix (holdraw.sci.data,i,j) = Pix (in.sci.data,i,j);
	                Pix (holdraw.err.data,i,j) = Pix (in.sci.data,i,j);
	            }
	        }
	    }

	    /* Here we have to open the trace table (_1dt), as well
	       as poke into the dispersion solution structure (disp)
	       in order to retrieve the remaining information necessary
	       to run the MSM blaze correction algorithm. Note that this
	       code was put in place *way* after the initial requirement
	       that the code should not provide for cross-referencing
	       between spectral orders, so we end up with some more code
	       duplication (^@%$&^%#% !)

               We are also using the wavelength array generated at this
               step to find the Ly alpha avoidance region in prism data.
               This is done by forcing mref = 1 for prism data.
	    */

	    wpos = 0.0;
	    if (sts->fluxcorr == PERFORM && mref > 0) {

	        /* Get the trace for the reference order */
	        if ((status = GetTrace6 (sts, mref, &trace)))
	            return (status);

	        /* Add all pertinent offsets. */
		AddOffsets6 (sts, &slit);
	        ypos += trace->a2center + sts->total_offset[1];

	        /* Get the dispersion for the reference order. Note that
	           if no dispersion solution is activated, no phot
	           correction will take place, so the wpos and
	           ddisp values will never get used.
	        */

	        if (sts->dispcorr == PERFORM) {
		    double dshift;	/* on-board Doppler shift (low-res) */
		    double ix;		/* pixel number, could be a fraction */
		    double p, q;	/* weights for linear interpolation */
	            if ((status = GetInang6 (sts, &sts->inangtab, mref, &iac)))
	                return (status);
	            if ((status = GetApOffset6 (sts, slit.ap_offset,
                                                disp->ref_aper, &delta)))
	                return (status);
	            if ((status = InterpDisp6 (&disp, ypos, &disp_y)))
	                return (status);
                    AdjustDisp6 (sts, disp_y, delta, &iac);
	            row_contents.sporder = mref;
	            row_contents.npts = in.sci.data.nx;
	            Wave (sts, disp_y, &row_contents);
		    dshift = orbitalDopp (sts->expstart, sts->expend,
			sts->doppzero, sts->orbitper, sts->doppmag);
	            i = row_contents.npts / 2;
	            if (row_contents.npts > 1024) {	/* high-res? */
			/* add 0.5 to get the point corresponding to the
			   center of low-res pixel 512
			*/
			ix = (float)i + 0.5 + 2. * dshift;
		    } else {
			ix = (float)i + dshift;
		    }
		    i = (int)ix;			/* truncate */
		    q = ix - (float)i;
		    p = 1. - q;
		    wpos = p * row_contents.wave[i] +
	                   q * row_contents.wave[i+1];
		    ddisp = row_contents.wave[i] - row_contents.wave[i-1];
	            if (row_contents.npts > 1024)
	                ddisp *= 2.;	/* dispersion per low-res pixel */
	            if (sts->heliocorr == PERFORM) {
		        wpos  /= sts->hfactor;
		        ddisp /= sts->hfactor;
	            }
	            FreeDisp6 (&disp_y);
	            FreeInang6 (&iac);
	        }
	        FreeTrace6 (&trace);
	    }

	    /* Reset output table row counter, average crosscor and
               extension header flag.
            */
	    o_row = 0;
	    avcrscroff = 0.0;
	    ext_ready = 0;

	    /* Reset flags that control loop termination. */
	    success  = 0;
	    skipping = 0;

	    /* Do for each spectral order in input image. */

	    for (row_contents.sporder =  minorder;
                 row_contents.sporder <= maxorder;
                 row_contents.sporder++) {

	        /* If last order was skipped but at least one previous
                   order was successfully extracted, assume orders are
                   falling off from the detector and finish current IMSET.
                */
	        if (success == 1 && skipping == 1)
	            break;

	        /* This flag controls skipping by dummy reference entries. */
	        sts->x1d_o = PERFORM;

		/* Get the XtractInfo record for this order. */
		if ((status = ReturnXtract (&extract, row_contents.sporder,
                                            &extract_o)))
		    return (status);

	        /* Select extraction algorithm for this order. If no
                   reference files exist for implementing the chosen
                   algorithm, skip current spectral order.

                   This was superseded by OPR 37810. The algorithm must
                   be taken not from the current spectral order entry, but
                   always from the first matching entry. But we leave the
                   code logic in place just in case (OPRs do come and go...)
                   The use of minorder to extract the extraction record
                   presumes a sorted (by sporder) extraction table.
                */
		if ((status = ReturnXtract (&extract, minorder, &extract_a)))
		    return (status);
	        if ((status = SelectAlg (sts, extract_a))) {
	            if ((status == CAL_FILE_MISSING)) {
	                status = 0;
	                printf ("Warning  Skipping order.\n");
	                FreeXtract (&extract_a);
	                FreeXtract (&extract_o);
	                continue;
	            } else
		        return (status);
	        }
	        FreeXtract (&extract_a);

                /* Extract the CoordInfo record for this order. This is
                   used for two things: (i) to get the appropriate plate scale
                   in the A2 direction (CDELT2). The plate scale is used
                   to scale any eventual POSTARG2 value to reference pixel
                   units. Plate scale may be a function of spectral order,
                   thus we have to read it here and not as a global parameter.
                   And (ii) to get the rectified image size in the A2 direction
                   built by calstis7 for the current spectral order. This
                   is used by the scaterred light correction algorithm.
                   If no matching row exists, skip silently.

                   This code is being called only if scattered light
                   correction is active, since POSTARG correction was
                   turned off. Do not forget to call FreeCoord6 (&coord_o)
                   only if scattered light correction is on.
                 */
	        if (sts->scatter) {
	            sdc_status = ReturnCoord6 (&coords, row_contents.sporder,
                                               &coord_o);
	            if (sdc_status) {
	                status = 0;
	                FreeXtract (&extract_o);
	                FreeCoord6 (&coord_o);
	                continue;
	            }
	        }

	        /* Get spectrum trace info for this order. If no matching
                   row in trace table, silently skip to next order. Notice
                   that this implies opening and closing the trace table
                   every time. This shouldn't be a burden in elapsed time
                   if the extraction and trace table are reasonably matched.
                   If any matching trace table entry is DUMMY, skip and warn.
                */
	        trc_status = GetTrace6 (sts, row_contents.sporder, &trace);

	        if (trc_status == NO_TRACE) {
	            status = 0;
	            FreeTrace6 (&trace);
	            FreeXtract (&extract_o);
	            if (sts->scatter)
	                FreeCoord6 (&coord_o);
	            continue;
	        } else if (trc_status == ERROR_TRACE)
	            return (trc_status);
	        else if (sts->x1d_o == DUMMY) {
	            PrGrpBegin ("order", row_contents.sporder);
	            warnDummy ("SPTRCTAB", row_contents.sporder, 1);
	            printf ("\n");
	            FreeTrace6 (&trace);
	            FreeXtract (&extract_o);
	            if (sts->scatter)
	                FreeCoord6 (&coord_o);
	            continue;
	        }

	        if (sts->verbose == 1 || sts->verbose == 2)
	            PrGrpBegin ("order", row_contents.sporder);
	        else
	            printf ("Extracting order %d\n", row_contents.sporder);

                /* This is to warn in case the command-line-entered A2CENTER
                   value is too off from the respective nominal position
                   of the also command-line-entered SPORDER. If such happens,
                   the wavelength scale will not correspond to the desired
                   A2CENTER position.
                */
                if (minorder == maxorder && sts->echelle &&
                    sts->cl_a2center != NO_POSITION &&
                    sts->dispcorr == PERFORM) {
	            trc_status1 = GetTrace6 (sts, minorder-1, &trace1);
	            trc_status2 = GetTrace6 (sts, minorder+1, &trace2);
	            if (trc_status1 != ERROR_TRACE &&
	                trc_status2 != ERROR_TRACE) {
	                d1 = fabs (cl_a2center - trace1->a2center);
	                d2 = fabs (cl_a2center - trace2->a2center);
	                if (trc_status1 == NO_TRACE) d1 = 1.E4;
	                if (trc_status2 == NO_TRACE) d2 = 1.E4;
                        d1 = (d1 < d2) ? d1 : d2;
	                d2 = fabs (cl_a2center - trace->a2center);
                        if (d2 > d1)
	                    printf (
            "Warning  Extraction position corresponds to another SPORDER.\n");
	                FreeTrace6 (&trace1);
	                FreeTrace6 (&trace2);
	            }
                }

		/* Compute total offset = "dither" + wavecal + aperture. */
		AddOffsets6 (sts, &slit);

	        /* Compute nominal spectrum trace position by adding offsets
                   to either table or user-supplied value.
                 */
	        if (sts->cl_a2center == NO_POSITION)
	            sts->nm_a2center = trace->a2center;
	        else
	            sts->nm_a2center = cl_a2center;

	        sts->nm_a2center += sts->total_offset[1];

	        /* POS TARG correction in spatial direction. The plate
                   scale read from table is in arcsec per reference pixel,
                   thus the POS TARG correction can be directly added to
                   the nominal center position at this point.
                */
	        /* sts->nm_a2center += sts->pos_targ[1] / coord_o->cdelt2;   */

	        /* Interpolate in the trace table. */
	        if ((status = InterpTrace6 (&trace, sts->nm_a2center, &trace_y)))
	            return (status);
	        sts->nm_a2center = trace_y->a2center;

	        /* Check if nominal center is within the data array
                   boundaries. Skip order if off-limits. This should take
                   care of both sub-arrays and wandering off-detector
                   spectral orders.
                */
		if (sts->echelle) {	/* added by PEH 1997 Dec 12 */
		    a2phys = (int)(sts->nm_a2center * sts->ltm[1] +
                             sts->ltv[1]);
		    if (a2phys < 0 || a2phys >= in.sci.data.ny) {
			 printf (
      "Warning  Extraction position outside data array. Skipping order.\n");
			skipping = 1;
			FreeTrace6 (&trace);
			FreeTrace6 (&trace_y);
	                FreeXtract (&extract_o);
	                if (sts->scatter)
	                    FreeCoord6 (&coord_o);
			continue;
		    } else
			skipping = 0;
		} else {
		    skipping = 0;
		}

	        if (sts->verbose == 1 || sts->verbose == 2)
	            PrSwitch6  (sts, "x1dcorr", PERFORM);

	        if (!sts->cc_global) {

	            /* Refine trace position by cross correlation. If
                       supplied, use search range from command line.
                       Only perform crosscor if the range is meaningful,
                       otherwise just use nominal value for spectrum
                       trace position.
                    */
	            if (sts->maxsearch == NO_RANGE)
	                maxsearch = extract_o->maxsearch;
	            else
	                maxsearch = sts->maxsearch;
	            ccstatus = 0;

	            /* If 1st order, look for geocoronal Lya. */
	            if (minorder == maxorder) {
	                FindLya (sts, &in, &slit,
	                        row_contents.wave, row_contents.npts,
	                        &(sts->avoid1a), &(sts->avoid2a),
	                        &(sts->avoid1b), &(sts->avoid2b));
	            } else {
	                sts->avoid1a = 0;
	                sts->avoid2a = 0;
	                sts->avoid1b = 0;
	                sts->avoid2b = 0;
	            }

 	            if (maxsearch != NO_RANGE) {

	                /* Do crosscorr. */
	                ccstatus = CrossCorr (sts, trace, &in, &ssgx, &ssgy,
                                              maxsearch, sts->avoid1a,
                                              sts->avoid2a);

	                /* Locate where in offset array is the current
                           sporder. Offset array is populated by global
                           croscor offset analysis.
                        */
	                norder = 0;
	                for (i = 0; i < sts->cc_size; i++) {
	                    if (row_contents.sporder == sts->cc_spord[i]) {
	                        norder = i;
	                        break;
	                    }
	                }

	                /* Croscorr may hard-fail or soft-fail. A soft failure
                           means that signal was located but its offset is too
                           far from the average offset derived by the croscor
                           global analysis. This information is stored in the
                           cc_rej vector for each individual order.
                        */
	                if (ccstatus || sts->cc_rej[norder] == 1) {

	                    /* This is the new behavior when failing:
                               use result of global analysis.
                            */
	                    if (sts->verbose == 1 || sts->verbose == 2)

	                        printf (
	    "Warning  Cross correlation to locate spectrum failed.\n");
	                     /* Prints this only if soft-fail. */
	                    if (!ccstatus && (sts->verbose == 1))
	                        printf (
	    "Warning  The offset found was: %g pixels\n", sts->crscroff);
	                     /* Global offset exists only in echelle mode. */
	                    if (sts->echelle) {
	                        if (sts->verbose == 1 || sts->verbose == 2)
	                            printf (
	                       "Warning  Using global offset instead.\n");
	                        sts->crscroff = sts->gcrscroff;
	                    } else {
	                        sts->crscroff = 0.0;
				/* This was put in place on 12 Apr 01 to
                                   comply with OPR 43663 (IB).
                                */
	                        printf ("ERROR    Cannot extract.\n");
				skipping = 1;
				FreeTrace6 (&trace);
				FreeTrace6 (&trace_y);
				FreeXtract (&extract_o);
				if (sts->scatter)
				    FreeCoord6 (&coord_o);
				continue;
	                    }

	                    /* This was the old behavior when failing. */
	                    /* printf (
	    "Warning  Spectrum extracted at default Y position: %g\n",
                        sts->nm_a2center * sts->ltm[1] + sts->ltv[1] + 1.0);
	                    */
	                } else if (sts->verbose == 2)
		            printf (
            "         Cross correlation offset: %g pixels\n", sts->crscroff);
	            } else
	                sts->crscroff = 0.0;
	        } else

	            /* Always use the global offset. */
	            sts->crscroff = sts->gcrscroff;

	        /* Compute actual average offset for this IMSET. */
	        avcrscroff += sts->crscroff;

	        /* Update output table with croscorr offset. */
	        row_contents.cc_offset = sts->crscroff;

	        /* Update spectrum trace position with crosscor result. */
	        sts->cc_a2center = sts->nm_a2center + sts->crscroff;

	        /* Re-interpolate in the trace table with updated A2
                   position. This is not really necessary for
                   closely-spaced echelle extractions, but in long-slit
                   mode with a large search box the A2 position before
                   crosscor may be too off from the true position.
                */
	        if (!sts->echelle) {
	            FreeTrace6 (&trace_y);
	            if ((status = InterpTrace6 (&trace, sts->cc_a2center,
                                                &trace_y)))
	                return (status);
	        }

	        /* Get intensity array for optimal extraction. If no entry
                   for current spectral order can be found in OSPECTAB,
                   reset extraction algorithm to default.
                */
	        if (sts->optimal) {
	            opt_status = GetIntens (sts, row_contents.sporder, &inta);
	            if (opt_status == ROW_NOT_FOUND) {
	                FreeIntensity (&inta);
	                status = 0;
	                sts->optimal = 0;
	                printf (
          "Warning  Extraction algorithm set to UNWEIGHTED for this order.\n");
	            } else if (opt_status != 0) {
	                FreeIntensity (&inta);
	                return (opt_status);
	            }

	            /* If intensity array size does not match input
                       image X size, reset extraction algorithm to default.
                    */
	            if (inta.allocated) {
	                if (inta.nelem != in.sci.data.nx) {
	                    sts->optimal = 0;
	                    printf (
          "Warning  Intensity array size different from image X size.\n");
	                    printf (
          "Warning  Extraction algorithm set to UNWEIGHTED for this order.\n");
	                 }
	            }
	        }

	        /* Get profiles for optimal extraction. If no entry for
                   current spectral order can be found in OPROFTAB,
                   reset extraction algorithm to default.
                */
	        if (sts->optimal) {
	            opt_status = GetProfile (sts, row_contents.sporder,
                                             &profa);
	            if (opt_status == ROW_NOT_FOUND) {
	                FreeProfileArray (&profa);
	                status = 0;
	                sts->optimal = 0;
	                printf (
          "Warning  Extraction algorithm set to UNWEIGHTED for this order.\n");
	            } else if (opt_status != 0) {
	                FreeProfileArray (&profa);
	                return (opt_status);
	            }
	        }

	        /* Define actual extraction box size to use. The basic
                   value comes from XTRACTAB, but it can be overriden
                   by command line setting or optimal profile file.
                   The profile Y size takes precedence over everything
                   else. Notice that the profile size is already in
                   image (physical) pixels, while extrsize is defined
                   in reference pixel units. Thus the profile file and
                   the image being processed must have consistent
                   binning parameters.
                */
	        if (sts->extrsize == NO_SIZE)
	            extrsize = extract_o->extrsize;
	        else
	            extrsize = sts->extrsize;

	        /* If optimal extraction, the extraction size must be
                   computed consistently with the profile array size
                   definition used in Allocprofile. This code duplication
                   cannot be avoided due to the convoluted way the optimal
                   extraction code was developed.
                */
	        if (sts->optimal)
                    extrsize = ((profa->npts / sts->subscale) - PROF_MARGIN) /
                                sts->ltm[1];

	        /* Alloc memory for profile storage. This is used both
                   for optimal extractions and for profile building.
                */
	        if ((status = AllocProfile (sts, in.sci.data.nx, extrsize)))
	            return (status);

	        /* Print background switch status. The switch is
                   set to PERFORM in case of optimal extraction.
                */
	        if (sts->backcorr == PERFORM || sts->optimal) {
	            if (sts->verbose == 1 || sts->verbose == 2)
	                PrSwitch6 (sts, "backcorr", PERFORM);
	        } else {
	            if (sts->verbose == 1 || sts->verbose == 2)
	                PrSwitch6 (sts, "backcorr", OMIT);
		}

	        /* Perform actions specific to the optimal extraction
                   algorithm. This can't be done before since only here
                   we can be sure that optimal extraction will take place
                   at least for one spectral order.
                */
	        if (sts->optimal) {

	            /* Initialize output weights image. This is performed
                       only once for each IMSET.
                    */
	            if (!outw_exists) {
	                if (allocSingleGroup (&outw, in.sci.data.nx,
                                                     in.sci.data.ny, True) == ALLOCATION_PROBLEM)
	                    return (OUT_OF_MEMORY);
	                outw_exists = 1;
	            }

	            /* Generate profile array. */
	            if ((status = BuildOptProf (sts, &profa)))
	                return (status);

	            /* Compute constants. */
	            if (sts->atodgain > 0.0)
	                sts->gain = 1.0 / sts->atodgain;
	            else
	                sts->gain = 1.0 / sts->ccdgain;
	            sts->rn2  = sts->readnse * sts->crsplit;
	            sts->rn2 *= sts->rn2;

	            /* Set pixel rejection flags to whatever is selected
                       in the image header (SDQFLAGS keyword).
                    */
	            /* Ignore small blemish (dust mote) - 16/2/01 IB */
	            sts->sdqflags = sts->sdqflags_orig & ~(SMALLBLEM);

	        } else {

	            /* Unweighted algorithm does not work with dq flags,
                       but profile builder does.
                    */
	            /* Ignore small blemish (dust mote) - 16/2/01 IB */
	            if (sts->do_profile)
	                sts->sdqflags = sts->sdqflags_orig & ~(SMALLBLEM);
	            else
	                sts->sdqflags = 0;
	        }

	        /* If scattered light correction is active, re-define
                   background extraction regions.
                */
	        if (sts->scatter) {
	            if ((status = DefineBackRegions (sts, trace_y, coord_o, &in,
                            &ssgx, &ssgy, &ilow_end, &ihigh_end))) {

	                if (status == ERROR_RETURN) {
	                    status = 0;
	                    printf (
"Warning  Skipping order. No usable background regions could be found.\n");
	                    FreeXtract (&extract_o);
	                    FreeCoord6 (&coord_o);
	                    continue;
	                } else

	                    return (status);
	            }
	            FreeCoord6 (&coord_o);
	        }

	        /* Extract spectrum, or build profile. */
	        if ((status = X1DSpec(sts, trace_y, extract_o, extrsize,
                                      &in, &outw, &ssgx, &ssgy, &inta,
                                      &row_contents)))
	            return (status);

	        /* Background smoothing must be done right
                   after it is extracted.
                 */
	        if ((status = SmoothBack (sts, &row_contents)))
	            return (status);

	        /* Filter background if correcting for scatterd light. */
	        if (sts->scatter && sts->lfilter > 1) {
	            if ((status = Lee (sts, &row_contents)))
	                return (status);
	        }

	        /* Update trailer file and header. */
	        if (sts->backcorr == PERFORM || sts->optimal) {
	            if (sts->verbose == 1 || sts->verbose == 2)
	                PrSwitch6 (sts, "backcorr", COMPLETE);
		    Put_KeyS (phdr, "backcorr", "COMPLETE", "");
	        }

	        if (sts->verbose == 1 || sts->verbose == 2)
	            PrSwitch6 (sts, "x1dcorr", COMPLETE);
		Put_KeyS (phdr, "x1dcorr", "COMPLETE", "");
	        if (sts->verbose == 2)
	            TimeStamp6 (sts->printtime, "X1DCORR complete",
                                sts->rootname);

	        /* Get rid of interpolated trace. */
	        FreeTrace6 (&trace_y);

	        /* Free memory associated with optimal extraction. */
	        if (sts->optimal) {
	            FreeIntensity (&inta);
	            FreeProfileArray (&profa);
	        }

	        /* Report extraction position on trailer file, but only if
                   it is inside the data array.
                */
	        if (!sts->do_profile) {
	            hold = sts->cc_a2center * sts->ltm[1] + sts->ltv[1] + 1.0;
	            if (hold > 0.0 && hold < in.sci.data.ny) {
	                if (sts->verbose == 1 || sts->verbose == 2)

	                    printf (
                "         Spectrum extracted at y position = %g\n",hold);
	            }
	        }
	        if (sts->verbose == 2 && !sts->do_profile) {
	            /* See if there are command-line overrides. */
	            for (i = 0; i < 2; i++) {
	                if (sts->bksize[i] == NO_SIZE)
	                    bksize[i] = extract_o->bksize[i];
	                else
	                    bksize[i] = sts->bksize[i];
	                if (sts->bkoffset[i] == NO_SIZE)
	                    bkoffset[i] = extract_o->bkoffset[i];
	                else
	                    bkoffset[i] = sts->bkoffset[i];
	            }
	            printf ("         Extraction box height = %g\n",
                            extrsize * sts->ltm[1]);
printf ("         Background box 1 height = %g offset %g from A2CENTER\n",
                            bksize[0]   * sts->ltm[1],
                            bkoffset[0] * sts->ltm[1]);
printf ("         Background box 2 height = %g offset %g from A2CENTER\n",
                            bksize[1]   * sts->ltm[1],
                            bkoffset[1] * sts->ltm[1]);
	        }
	        if (!sts->do_profile) {
	            if (sts->verbose == 1 || sts->verbose == 2)
	                printf ("\n");
	        }
	        /* Generate wavelength and flux arrays. */
	        if (sts->dispcorr == PERFORM) {

	            if (sts->verbose == 1 || sts->verbose == 2)
	                PrSwitch6 (sts, "dispcorr", PERFORM);

	            /* Print reference file info. */
	            if (sts->verbose == 1 || sts->verbose == 2)
	                Message6 (sts, DISP_INFO);

	            /* Read reference data for current order. */

	            /*  This is being removed from here to outside the loop
                        in sporder. The reason is that the dispersion
                        solution coefficients shouldn't be selected based on
                        sporder. However, we keep this comment-out code in
                        here for the time beign, just in case.

	            if ((status = GetDisp6 (sts, row_contents.sporder, &disp))) {
	                if (status == ROW_NOT_FOUND) {
	                    FreeDisp6 (&disp);
	                    FreeTrace6 (&trace);
	                    FreeProfile (sts, in.sci.data.nx);
	                    FreeXtract (&extract_o);
	                    continue;
	                } else
	                    return (status);
	            }
	            if (sts->x1d_o == DUMMY) {
	                warnDummy ("DISPTAB", row_contents.sporder, 1);
	                FreeDisp6 (&disp);
	                FreeTrace6 (&trace);
	                FreeProfile (sts, in.sci.data.nx);
	                FreeXtract (&extract_o);
	                continue;
	            }

	            */

	            if ((status = GetInang6 (sts, &sts->inangtab,
                                             row_contents.sporder, &iac))) {
	                /* Just skip order if row not found. */
	                if (status == ROW_NOT_FOUND) {
	                    FreeInang6 (&iac);
	                    FreeTrace6 (&trace);
	                    FreeProfile (sts, in.sci.data.nx);
	                    FreeXtract (&extract_o);
	                    continue;
	                } else
	                    return (status);
	            }
	            if (sts->x1d_o == DUMMY) {
	                warnDummy ("INANGTAB", row_contents.sporder, 1);
	                FreeInang6 (&iac);
	                FreeTrace6 (&trace);
	                FreeProfile (sts, in.sci.data.nx);
	                FreeXtract (&extract_o);
	                continue;
	            }

	            /* Get offset of slit from that used to measure the
                       dispersion relation.
                    */
	            if ((status = GetApOffset6 (sts,
                            slit.ap_offset, disp->ref_aper, &delta))) {
	                /* Just skip order if row not found. */
	                if (status == ROW_NOT_FOUND) {
	                    FreeInang6 (&iac);
	                    FreeTrace6 (&trace);
	                    FreeProfile (sts, in.sci.data.nx);
	                    FreeXtract (&extract_o);
	                    continue;
	                } else
	                    return (status);
	            }
	            if (sts->x1d_o == DUMMY) {
	                warnDummy ("APDESTAB", row_contents.sporder, 1);
	                FreeInang6 (&iac);
	                FreeTrace6 (&trace);
	                FreeProfile (sts, in.sci.data.nx);
	                FreeXtract (&extract_o);
	                continue;
	            }
	            if (sts->verbose == 2 && !sts->do_profile)
	                printf ("         Delta = %.6g arcsec.\n", delta);

	            /* Interpolate the dispersion relation. */
	            if ((status = InterpDisp6 (&disp, sts->cc_a2center, &disp_y)))
	                return (status);

	            /* Modify the dispersion coefficients using IAC and MSM. */
	            AdjustDisp6 (sts, disp_y, delta, &iac);

	            /* Compute wavelengths. */
	            Wave (sts, disp_y, &row_contents);

	            if (sts->verbose == 1 || sts->verbose == 2)
	                PrSwitch6 (sts, "dispcorr", COMPLETE);
		    Put_KeyS (phdr, "dispcorr", "COMPLETE", "");

	            if (sts->verbose == 2)
	                TimeStamp6 (sts->printtime, "DISPCORR complete",
                                    sts->rootname);
	            if (!sts->do_profile) {
	                if (sts->verbose == 1 || sts->verbose == 2)
	                    printf ("\n");
	            }
	            if (sts->heliocorr == PERFORM) {
	                if (sts->verbose == 1 || sts->verbose == 2) {
	                    PrSwitch6 (sts, "helcorr", PERFORM);
	                    PrSwitch6 (sts, "helcorr", COMPLETE);
	                }
			Put_KeyS (phdr, "helcorr", "COMPLETE", "");
	                if (sts->verbose == 2)
	                    TimeStamp6 (sts->printtime, "HELCORR complete",
                                        sts->rootname);
	                if (!sts->do_profile) {
	                    if (sts->verbose == 1 || sts->verbose == 2)
	                        printf ("\n");
	                }
	            } else {
	                if (sts->verbose == 1 || sts->verbose == 2)
	                    PrSwitch6 (sts, "helcorr", OMIT);
	                if (!sts->do_profile) {
	                    if (sts->verbose == 1 || sts->verbose == 2)
	                        printf ("\n");
	                }
	            }

	            /* Convert to absolute flux. */
	            if (sts->fluxcorr == PERFORM && !sts->do_profile) {

	                /* Print reference file info. */
	                if (sts->verbose == 1 || sts->verbose == 2)
	                    Message6 (sts, FLUX_INFO);

	                if ((status = GetAbsPhot6 (sts,
                                row_contents.sporder, &phot, 1, &warn))) {
	                    /* Just skip order if row not found. */
	                    if (status == ROW_NOT_FOUND) {
	                        FreePhot6 (&phot);
	                        FreeInang6 (&iac);
	                        FreeDisp6 (&disp_y);
	                        FreeTrace6 (&trace);
	                        FreeProfile (sts, in.sci.data.nx);
	                        FreeXtract (&extract_o);
	                        printf ("ERROR    Skipping order %d\n",
                                         row_contents.sporder);
	                        continue;
	                    } else
	                        return (status);
	                }

	                /* If DUMMY pedigree was found, skip flux
                           calibration but continue with extraction.
                           (behavior added 22 Sep 00).
                        */
	                if (sts->x1d_o == DUMMY) {
	                    warnDummy ("PHOTTAB", row_contents.sporder, 0);
			    for (i = 0; i < row_contents.npts; i++) {
				row_contents.flux[i]  = 0.0F;
				row_contents.error[i] = 0.0F;
			    }
	                    if (sts->verbose == 1 || sts->verbose == 2) {
	                        PrSwitch6 (sts, "fluxcorr", OMIT);
	                        printf ("\n");
	                    }
	                } else {
	                    if (sts->verbose == 1 || sts->verbose == 2) {
	                        PrSwitch6 (sts, "fluxcorr", PERFORM);
                                if (sts->detector == CCD_DETECTOR) {
                                    if (sts->ctecorr == PERFORM) {
                                        PrSwitch6 (sts, "ctecorr", PERFORM);
                                        Message6 (sts, CCD_INFO);
                                    } else
                                        PrSwitch6 (sts, "ctecorr", OMIT);
                                }
                            }

	                    /* This is an inefficient way to initialize the
                               auxiliary phot structure. This is a temporary
                               solution to the problem of PCT interpolation.
                            */
	                    dum = GetAbsPhot6 (sts, row_contents.sporder,
                                               &photc, 1, &warn);

		            /* Get PCT info. A zeroed height means to get
                               the photometry correction for the maximum
                               box height available in the table.
                            */
		            if ((status = GetPCT6 (sts, &phot, 0.0, WARN))) {
	                        FreePhot6 (&phot);
			        return (status);
	                    }
		            if ((status = GetPCT6 (sts, &photc, extrsize,
                                                   NO_WARN))) {
	                        FreePhot6 (&phot);
	                        FreePhot6 (&photc);
			        return (status);
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

	                    /* Do fluxcorr */

	                    status = AbsFlux6 (sts, &row_contents,
	                                       &phot, &photc, &slit, &tds,
                                               &cti, sts->hfactor,
                                               sts->atodgain, sts->optimal,
                                               sts->profile_rejf, &blazeshift);
                            if (status) {
	                        FreePhot6 (&phot);
	                        FreePhot6 (&photc);
			        return (status);
                            }
			    if (phot.blazecorr == PERFORM) {
				sum_blazeshift += blazeshift;
				n_blazeshift += 1.;
			    }
			    if (sts->gaccorr == PERFORM) {
				printf (
			"         GACTAB correction applied for %s %s %d\n",
				sts->opt_elem, sts->aperture, sts->cenwave);
			    }
	                    if (sts->verbose == 1 || sts->verbose == 2) {
				if (sts->detector == CCD_DETECTOR &&
				    sts->ctecorr == PERFORM) {
				    PrSwitch6 (sts, "ctecorr", COMPLETE);
				    Put_KeyS (phdr, "ctecorr", "COMPLETE", "");
				}
	                        PrSwitch6 (sts, "fluxcorr", COMPLETE);
                            }
			    Put_KeyS (phdr, "fluxcorr", "COMPLETE", "");
	                    if (sts->verbose == 2) {
	                        TimeStamp6 (sts->printtime, "FLUXCORR complete",
                                           sts->rootname);
	                        printf ("\n");
	                    }
	                    FreePhot6 (&phot);
	                    FreePhot6 (&photc);
	                }
		    } else {
	                if (!sts->do_profile) {
	                    if (sts->verbose == 1 || sts->verbose == 2) {
	                        PrSwitch6 (sts, "fluxcorr", OMIT);
	                        printf ("\n");
	                    }
	                }
	            }

	            /* Free memory associated with wavelength computation */
	            FreeInang6 (&iac);
	            FreeDisp6 (&disp_y);
	        } else {

	            /* With no wavelength calibration, neither flux nor
                       heliocentric correction can be done.
                    */
	            if (sts->verbose == 1 || sts->verbose == 2)
	                PrSwitch6 (sts, "dispcorr", OMIT);
	            if (!sts->do_profile) {
	                if (sts->verbose == 1 || sts->verbose == 2) {
	                    printf ("\n");
	                    PrSwitch6 (sts, "fluxcorr", OMIT);
	                    printf ("\n");
	                }
	            }
	            if (sts->verbose == 1 || sts->verbose == 2)
	                PrSwitch6 (sts, "helcorr", OMIT);
	            if (!sts->do_profile) {
	                if (sts->verbose == 1 || sts->verbose == 2)
	                    printf ("\n");
	            }
	        }

	        /* No small geom. correction implemented yet. */
	        if (sts->verbose == 1 || sts->verbose == 2)
	             PrSwitch6 (sts, "sgeocorr", OMIT);
	        if (!sts->do_profile) {
	            if (sts->verbose == 1 || sts->verbose == 2)
	                printf ("\n");
	        }

	        /* If this is the first row to be output in this extension,
                   first create the output table extension. But if this is
                   also the first extension to be created, first create
                   the output primary header. This complicated logic avoids
                   the creation of empty headers and/or extensions.
                */
	        if (o_row == 0) {
                    if (first_extver) {
	                im = openOutputImage (sts->output, "", 0, phdr,
                                              0, 0, FITSBYTE);
	                if (hstio_err())
	                    return (OPEN_FAILED);
	                closeImage (im);
	                phu_ready = 1;
			first_extver = 0;		/* false */
	            }
	            if (!sts->do_profile) {
	                if ((status = CreaTable (sts, extver, &otable)))
	                    return (status);
	            } else {
	                if ((status = CreaProfTable (sts, &optable)))
	                    return (status);
	            }
	            ext_ready = 1;
	        }

	        /* Output current row. */
	        if (ext_ready) {
	            if (!sts->do_profile) {

	                if (WriteRow (sts, &otable, &row_contents, &o_row)) {
	                    FreeTrace6 (&trace);
	                    FreeProfile (sts, in.sci.data.nx);
	                    FreeXtract (&extract_o);
	                    continue;
	                }

	            } else {

	                /* Set flags for profile rejection. */
	                if (strlen (sts->rejranges) > 0)
	                    SetRanges (sts, &row_contents);

	                if (WriteProfRow (sts, &optable, &row_contents,
                                          &o_row)) {
	                    FreeTrace6 (&trace);
	                    FreeProfile (sts, in.sci.data.nx);
	                    FreeXtract (&extract_o);
	                    continue;
	                }
	            }
	        }

	        success = 1;
	        if (sts->verbose == 1 || sts->verbose == 2)
	            PrGrpEnd ("order", row_contents.sporder);
	        if (!sts->do_profile) {
	            if (sts->verbose == 1 || sts->verbose == 2)
	                printf ("\n");
	        }
	        /* Free memory for this spectral order. */
	        FreeTrace6 (&trace);
	        FreeProfile (sts, in.sci.data.nx);
	        FreeXtract (&extract_o);

	        /* Restore input image contents that were smoothed out
                   by scattered light correction algorithm.
                */
	        if (sts->scatter) {
	            for (i = 0; i < holdraw.sci.data.nx; i++) {
	                for (j = ilow_end; j <= ihigh_end; j++) {
	                    Pix (in.sci.data,i,j) = Pix (holdraw.sci.data,i,j);
	                    Pix (in.err.data,i,j) = Pix (holdraw.err.data,i,j);
	                }
	            }
	        }
	    }

	    /* End loop over spectral orders. */

	    if (sts->dispcorr == PERFORM)
	        FreeDisp6 (&disp);

	    if (o_row == 0)
	        printf ("Warning  No rows were written; no table created.\n");

	    if (sts->verbose && sts->trace_rotation != 0.)
	        printf ("         trace was rotated by = %.6g degree.\n",
                       sts->trace_rotation);

	    /* Update current extension header with average
               crosscor, blazeshift, and close current table
               extension.
            */
	    if (o_row > 0 && !sts->do_profile) {
	        c_tbhadd (otable.tp, "CRSCROFF", avcrscroff/o_row);
	        c_tbhpcm (otable.tp, "CRSCROFF",
                "offset from 1-D extraction cross-corr.");
                if (bcorr_l == PERFORM || blazeshift != NO_VALUE) {
		    if (n_blazeshift > 0.)
			blazeshift = sum_blazeshift / n_blazeshift;
	            c_tbhadd (otable.tp, "BLZSHIFT", blazeshift);
	            c_tbhpcm (otable.tp, "BLZSHIFT",
				"average blaze shift (pixels)");
	        }
		if (sts->heliocorr == PERFORM) {
		    c_tbhadd (otable.tp, "V_HELIO", radvel);
		    c_tbhpcm (otable.tp, "V_HELIO",
				"heliocentric radial velocity (km/s)");
		}
	        c_tbtclo (otable.tp);
	    }

	    /* Close profile table. */
	    if (sts->do_profile)
	        c_tbtclo (optable.tp);

	    /* Output weights image. */
	    if (outw_exists) {

	        if (FileExists (sts->outw) == 0) {

	            /* If no primary header was yet created for the output
                       weights image, create it as a plain copy of the
                       input image header. The logic here avoids creation of
                       a output file with primary header but no extensions.
                    */
	            if (!phuw_ready) {
	                im = openOutputImage (sts->outw, "", 0, phdr,
                                              0, 0, FITSBYTE);
	                if (hstio_err())
	                    return (OPEN_FAILED);
	                closeImage (im);
	                phuw_ready = 1;
	            }
	            if (putSingleGroup (sts->outw, extver, &outw, 0))
	                return (OPEN_FAILED);
	            printf ("         Weights image %s written to disk.\n",
                            sts->outw);
	        } else {

                   /* If a file with same name exists, warn and do not
                      write anything. This can be misleading to users
                      running the task several times and not paying
                      attention to the _trl output.
	           */
	           printf (
    "Warning  Weights image %s already exists. Nothing written this time.\n",
                 sts->outw);
	        }
	    }

	    /* Free memory for this IMSET. */
	    FreeOutArrays (&row_contents);
	    freeSingleGroup (&in);
	    freeSingleGroup (&outw);
	    if (sts->scatter)
	        freeSingleGroup (&holdraw);
	    if (sts->verbose == 1 || sts->verbose == 2)
	        PrGrpEnd ("imset",  extver);
	    if (!sts->do_profile) {
	        if (sts->verbose == 1 || sts->verbose == 2)
	            printf ("\n");
	        }
	}

	/* Free memory for all IMSETs. */
	FreeThroughput6 (&slit);
	FreeTds (&tds);
	FreeXtract (&extract);
	FreeCoord6 (&coords);
	free (sts->cc_off);
	free (sts->cc_spord);
	free (sts->cc_rej);

	if (sts->sgeocorr == PERFORM) {
	    freeFloatHdrData (&ssgx);
	    freeFloatHdrData (&ssgy);
	}

	/* Update headers. */
	if (phu_ready) {
	    if ((status = intUpdateHeader (sts->output, "NEXTEND", o_ext,
                                           "number of extensions", 0)))
	        return (status);
	    if (!(sts->do_profile)) {
	        if (sts->optimal) {
	            if ((status = strUpdateHeader (sts->output, "XTRACALG",
                            OPTIMAL, "extraction algorithm", 0)))
	                return (status);
	        } else {
	            if ((status = strUpdateHeader (sts->output, "XTRACALG",
                            UNWEIGHTED, "extraction algorithm", 0)))
	                return (status);
	        }
	    }
	    if (sts->dispcorr == PERFORM && sts->heliocorr == PERFORM) {
                sprintf (str, "Heliocentric correction = %g km/s", radvel);
	        if ((status = HistoryAddHeader (sts->output, str)))
	            return (status);
	    }
	}

	return (0);
}



/* This routine modifies the dispersion coefficients in-place to account
   for the displacement of the image on the detector from the location
   that was used for measuring the coefficients.
   This only needs to be applied for echelle data.
*/

static int MOCAdjustDisp (StisInfo6 *sts, DispRelation *disp) {

	int mref;		/* spectral order number of reference order */
	double yref;		/* Y location of reference order */
	double a4corr;		/* correction factor */
	double ydiff;		/* Y position offset */
	double a2center;	/* from the trace for order mref */
	double r_shifta2;	/* shifta2 converted to reference pixel size */
	int status;

	SpTrace *trace;		/* spectral trace for spectral order mref */
	int GetTrace6 (StisInfo6 *, int, SpTrace **);
	void FreeTrace6 (SpTrace **);

	trace = NULL;

	/* Get parameters from dispersion table. */
	if ((status = GetMOC (&sts->disptab, sts->opt_elem, sts->cenwave,
                              &mref, &yref, &a4corr)))
	    return (status);
	if (a4corr == 0.)
	    return (0);		/* nothing further to do */

	/* Get the trace for spectral order mref, and get its a2center. */
	if ((status = GetTrace6 (sts, mref, &trace)))
	    return (status);
	a2center = trace->a2center;

	/* Here's what we're doing:
		ydiff = ypos - yref
		shifta2 = ypos - a2center, and shifta2 is called msm_offset[1]
	   so:
		ydiff = shifta2 + a2center - yref
	*/
	r_shifta2 = sts->msm_offset[1] / sts->ltm[1];	/* reference pixels */
	ydiff = r_shifta2 + a2center - yref;

	/* Adjust dispersion coefficients using a4corr. */
	disp->coeff[0] -= ydiff * sts->cenwave * a4corr;
	disp->coeff[4] += ydiff * a4corr;

	FreeTrace6 (&trace);

	return (0);
}



/*
   Open the two images in a small-scale geometric correction reference file.

   The small-scale geometric correction is NOT yet implemented in calstis6.
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


static void PrSwitch6 (StisInfo6 *sts, char *text, int cal_switch) {

	if (!sts->do_profile)
	    PrSwitch (text, cal_switch);
}



static void warnDummy (char *text, int sporder, int skipmsg) {

	printf ("ERROR    DUMMY pedigree in order %d in %s.", sporder, text);
	if (skipmsg)
	    printf ("Skipping order.");
	printf ("\n\n");
}
