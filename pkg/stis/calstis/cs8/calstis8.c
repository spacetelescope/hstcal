/* calstis8 -- Sum repeatobs data

   This file contains:
	CalStis8
	StisInit8
	GetKeyInfo8
	GetSpOrder
	NotMultiOrder
	SumGrps
	SumOrder
	FindExtver
	GetGrpInfo8
	PutHdrInfo8
	SquareErr
	SqrtErr
	RptSum
*/

# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <math.h>

# include "c_iraf.h"
# include "hstio.h"

# include "stis.h"
# include "calstis8.h"
# include "err.h"
# include "stisdef.h"

static int FindExtver (int *, int, int, int);
static int GetKeyInfo8 (StisInfo8 *, Hdr *);
static int GetGrpInfo8 (SingleGroup *, double *, double *, int *);
static int GetSpOrder (StisInfo8 *, int *, int *, int *);
static void NotMultiOrder (StisInfo8 *, int *, int *, int *);
static int PutHdrInfo8 (SingleGroup *, int, double, double, int);
static int RptSum (SingleGroup *, SingleGroup *);
static void SqrtErr (SingleGroup *);
static void SquareErr (SingleGroup *);
static void StisInit8 (StisInfo8 *);
static int SumGrps (StisInfo8 *, int *, int, int);
static int SumOrder (StisInfo8 *, int *, int, int);

/* This routine adds together individual imsets.

   Phil Hodge, 1997 Oct 20:
	Pass sts->ncombine to PutHdrInfo8, and update ncombine in output
	primary header.

   Phil Hodge, 1998 Jan 22:
	Add GetSpOrder, and other changes to allow summing echelle data.
	Pass ncombine instead of sts->ncombine to PutHdrInfo8, and update
	ncombine in output SCI extension header instead of primary header.

   Phil Hodge, 1998 Oct 5:
	Change status value 1001 to GENERIC_ERROR_CODE.

   Phil Hodge, 1998 Nov 12:
	A misleading error message was printed (OPR 37635).  In GetSpOrder,
	modify the section that checks the value of sporder gotten from the
	header to allow the possibility of sporder = 1 for echelle data.
	This would happen if the input data were flat fielded but not 2-D
	rectified, for example.

   Phil Hodge, 1999 Sept 30:
	In SumOrder, divide by the sum of the exposure times, rather than
	by the exposure time of the last imset.  This is the case where
	the data have been flux calibrated, so they're being averaged
	instead of summed, and the exposure time is used as the weight.

   Phil Hodge, 2008 Dec 12:
	Get keyword IMSET_OK, and skip image sets for which this is false.
*/

int CalStis8 (char *input, char *output, int printtime, int verbose) {

	int status;

	StisInfo8 sts;	/* calibration switches, reference files, etc. */

	IODescPtr im;		/* descriptor for input image */
	Hdr phdr;		/* primary header for input image */
	int *orders;		/* array of sporder for each input imset */
	int minorder, maxorder;	/* range of sporder numbers */

	PrBegin (8);

	if (printtime)
	    TimeStamp ("CALSTIS-8 started", "");

	/* Initialize structure containing calstis information. */
	StisInit8 (&sts);

	/* Copy command-line arguments into sts. */
	strcpy (sts.input, input);
	strcpy (sts.output, output);
	sts.printtime = printtime;
	sts.verbose = verbose;
	sts.statcorr = PERFORM;

	PrFileName ("input", sts.input);
	PrFileName ("output", sts.output);

	initHdr (&phdr);

	/* Check whether the output file already exists. */
	if ((status = FileExists (sts.output)))
	    return (status);

	/* Open input image in order to read its primary header. */
	im = openInputImage (sts.input, "", 0);
	if (hstio_err())
	    return (OPEN_FAILED);
	getHeader (im, &phdr);		/* get primary header */
	if (hstio_err())
	    return (OPEN_FAILED);
	closeImage (im);

	/* Get keyword values from primary header. */
	if ((status = GetKeyInfo8 (&sts, &phdr)))
	    return (status);

	/* Read the input file to get SPORDER for each imset. */
	if ((orders = malloc (sts.nimages * sizeof(int))) == NULL)
	    return (OUT_OF_MEMORY);
	if ((status = GetSpOrder (&sts, orders, &minorder, &maxorder)))
	    return (status);

	freeHdr (&phdr);

	/* Print information about this image. */
	PrHdrInfo (sts.obsmode, sts.aperture, sts.opt_elem, sts.det);

	if (sts.printtime)
	    TimeStamp ("Begin processing", sts.rootname);

	/* Sum imsets. */
	if ((status = SumGrps (&sts, orders, minorder, maxorder)))
	    return (status);

	free (orders);
	printf ("\n");
	PrEnd (8);

	if (sts.printtime)
	    TimeStamp ("CALSTIS-8 completed", sts.rootname);

	return (0);
}

/* Initialize the calstis8 structure.  This includes information about the
   input and output images, calibration files, and flags to specify which
   calibration steps are to be done.
*/

static void StisInit8 (StisInfo8 *sts) {

	/* Assign default values. */

	sts->input[0] = '\0';
	sts->output[0] = '\0';
	sts->rootname[0] = '\0';
	sts->obsmode[0] = '\0';
	sts->aperture[0] = '\0';
	sts->opt_elem[0] = '\0';
	sts->det[0] = '\0';
	sts->echelle = 0;
	sts->nimages = 1;
	sts->nrptexp = 1;
	sts->sdqflags = 32767;
}

/* This routine gets keyword values and the calibration switch STATFLAG
   from the primary header.
*/

static int GetKeyInfo8 (StisInfo8 *sts, Hdr *phdr) {

/* arguments:
StisInfo8 *sts  io: calibration switches and info
Hdr *phdr       i: primary header
*/

	int status;
	int sdqflags;			/* "serious" data quality flags */
	int nextend;			/* number of FITS extensions */
	int nrptexp;			/* number of exposures */
	int no_def = 0;			/* missing keyword is fatal error */
	int use_default = 1;		/* use default if keyword is missing */

	/* Have the data been converted from counts to absolute flux? */
	if ((status = GetSwitch (phdr, "FLUXCORR", &sts->fluxcorr)))
	    return (status);

	if ((status = Get_KeyS (phdr, "ROOTNAME",
                                no_def, "", sts->rootname, STIS_CBUF)))
	    return (status);

	if ((status = Get_KeyS (phdr, "OBSMODE",
                                no_def, "", sts->obsmode, STIS_CBUF)))
	    return (status);

	if ((status = Get_KeyS (phdr, "APERTURE",
                                no_def, "", sts->aperture, STIS_CBUF)))
	    return (status);

	if ((status = Get_KeyS (phdr, "OPT_ELEM",
                                no_def, "", sts->opt_elem, STIS_CBUF)))
	    return (status);

	sts->echelle = sts->opt_elem[0] == 'E' || sts->opt_elem[0] == 'e';

	if ((status = Get_KeyS (phdr, "DETECTOR",
                                no_def, "", sts->det, STIS_CBUF)))
	    return (status);

	if ((status = Get_KeyI (phdr, "SDQFLAGS", use_default, 32767, &sdqflags)))
	    return (status);
	sts->sdqflags = sdqflags;

	/* Find out how many extensions there are in this file. */
	if ((status = Get_KeyI (phdr, "NEXTEND",
                                use_default, EXT_PER_GROUP, &nextend)))
	    return (status);

	/* Convert number of extensions to number of SingleGroups. */
	sts->nimages = nextend / EXT_PER_GROUP;

	/* Get NRPTEXP and compare with nimages. */
	if (!sts->echelle) {
	    if ((status = Get_KeyI (phdr, "NRPTEXP", no_def, 0, &nrptexp)))
		return (status);
	    if (nrptexp != sts->nimages) {
		printf (
		"Warning  NEXTEND implies %d exposures, but NRPTEXP is %d.\n",
			sts->nimages, nrptexp);
	    }
	}

	if (sts->nimages <= 1) {
	    printf ("Warning  NEXTEND indicates there is only one imset.\n");
	    return (NOTHING_TO_DO);
	}

	return (0);
}

/* This routine reads the input file to get SPORDER for each imset,
   assigning the values to the elements of the orders array.  minorder
   and maxorder are the minimum and maximum values of SPORDER.
*/

static int GetSpOrder (StisInfo8 *sts,
		int *orders, int *minorder, int *maxorder) {

/* arguments:
StisInfo8 *sts  io: calibration switches and info; nrptexp is updated
int orders[]     o: array of spectral order numbers, sts->nimages elements
int *minorder    o: minimum value of sporder
int *maxorder    o: maximum value of sporder
*/

	int status;

	IODescPtr im;		/* descriptor for input image */
	Hdr hdr;		/* header for a SCI extension */
	int sporder;		/* value of keyword in an extension */
	int delta_sporder;	/* +1 or -1 indicates direction of change */
	int diff;		/* difference between adjacent sporders */
	int i;
	int use_default = 1;	/* use default if missing keyword */

	if (!sts->echelle) {
	    NotMultiOrder (sts, orders, minorder, maxorder);
	    return (0);
	}

	/* Loop over all imsets in input image. */
	for (i = 0;  i < sts->nimages;  i++) {

	    /* read extension header */
	    initHdr (&hdr);
	    im = openInputImage (sts->input, "SCI", i+1);
	    if (hstio_err())
		return (OPEN_FAILED);
	    getHeader (im, &hdr);
	    if (hstio_err())
		return (OPEN_FAILED);
	    closeImage (im);

	    /* Get sporder for current imset, and save in array. */
	    if ((status = Get_KeyI (&hdr, "SPORDER", use_default, -1, &sporder)))
		return (status);
	    if (sporder <= 1) {
		if (i == 0) {			/* first imset? */
		    /* could be prior to 2-D rectification */
		    NotMultiOrder (sts, orders, minorder, maxorder);
		    return (0);
		} else if (sporder == -1) {
		    printf ("ERROR    Keyword SPORDER missing, imset %d.\n",
				i+1);
		    return (KEYWORD_MISSING);
		} else {
		    printf (
	"ERROR    SPORDER = %d is invalid for echelle data; imset %d.\n",
			sporder, i+1);
		    return (HEADER_PROBLEM);
		}
	    }
	    orders[i] = sporder;

	    if (i == 0) {
		*minorder = sporder;
		*maxorder = sporder;
	    }
	    if (sporder < *minorder)
		*minorder = sporder;
	    if (sporder > *maxorder)
		*maxorder = sporder;

	    freeHdr (&hdr);
	}

	if (sts->nimages <= 1)
	    return (0);

	/* Are the spectral orders increasing or decreasing? */
	delta_sporder = (orders[1] > orders[0]) ? 1 : -1;

	/* Find out how many imsets to add together for each sporder by
	   counting the number of times sporder changes in the opposite
	   direction from the usual.
	*/
	sts->nrptexp = 1;
	for (i = 1;  i < sts->nimages;  i++) {
	    diff = orders[i] - orders[i-1];
	    if (diff * delta_sporder < 0)
		sts->nrptexp++;
	}

	return (0);
}

/* This assigns one to each element of the orders array and to minorder
   and maxorder.  Also, nrptexp in sts is set to nimages.
*/

static void NotMultiOrder (StisInfo8 *sts,
		int *orders, int *minorder, int *maxorder) {

	int i;

	for (i = 0;  i < sts->nimages;  i++)
	    orders[i] = 1;
	*minorder = 1;
	*maxorder = 1;
	sts->nrptexp = sts->nimages;
}

/* This routine sums all the input imsets.  The output will normally be
   just one imset, but for echelle data the output can be multiple imsets,
   one for each spectral order.
*/

static int SumGrps (StisInfo8 *sts, int *orders, int minorder, int maxorder) {

	int status;

	int sporder;		/* current spectral order number */
	int oextver;		/* output imset number */

	oextver = 0;
	for (sporder = minorder;  sporder <= maxorder;  sporder++) {

	    oextver++;
	    PrGrpBegin ("imset", oextver);

	    if ((status = SumOrder (sts, orders, sporder, oextver)))
		return (status);

	    PrGrpEnd ("imset", oextver);
	}

	return (0);
}

/* This routine sums the input imsets for one spectral order. */

static int SumOrder (StisInfo8 *sts, int *orders, int sporder, int oextver) {

	int status;

	SingleGroup x, y;	/* first imset and Nth imset */
	double exptime;		/* exposure time of current imset */
	double sumexptime = 0.;	/* accumulated exposure time */
	double expend;		/* expend of last imset read */
	char *message;		/* for printtime info */
	int iextver;		/* imset number in input image */
	int nrptexp;		/* actual number of imsets combined */
	int imset_ok;		/* value of header keyword IMSET_OK */
	int done;		/* loop termination flag */

	int doStat (SingleGroup *, short);

	initSingleGroup (&x);
	initSingleGroup (&y);

	if (sts->printtime) {
	    if ((message = calloc (STIS_LINE+1, sizeof (char))) == NULL)
		return (OUT_OF_MEMORY);
	}

	/* Find and read the first image that has IMSET_OK = T. */
	iextver = 0;
	done = 0;
	while (!done) {

	    /* Find the first imset for this sporder. */
	    iextver++;
	    iextver = FindExtver (orders, sts->nimages, iextver, sporder);
	    if (iextver < 0) {
		printf ("Warning  No data for spectral order %d.\n", sporder);
		return (0);
	    }

	    getSingleGroup (sts->input, iextver, &x);
	    if (hstio_err())
		return (OPEN_FAILED);

	    /* get from x */
	    if ((status = GetGrpInfo8 (&x, &exptime, &expend, &imset_ok)))
		return (status);

	    if (imset_ok) {
		if (sts->printtime) {
		    sprintf (message, "imset %d read", iextver);
		    TimeStamp (message, sts->rootname);
		} else if (sts->verbose) {		/* don't do both */
		    printf ("         imset %d read\n", iextver);
		}
		nrptexp = 1;				/* one image so far */
		sumexptime = exptime;
		/* Weight x by the exposure time, if we have flux. */
		if (sts->fluxcorr == COMPLETE) {
		    if ((status = multk2d (&x, (float)exptime)))
			return (status);
		}
		/* Square the errors to convert to variance. */
		SquareErr (&x);				/* operate on x */
		done = 1;
	    }
	}

	done = 0;
	while (!done) {

	    /* Find the next imset for this sporder. */
	    iextver++;
	    iextver = FindExtver (orders, sts->nimages, iextver, sporder);

	    if (iextver > 0) {

		getSingleGroup (sts->input, iextver, &y);
		if (hstio_err())
		    return (OPEN_FAILED);

		/* Update exposure time info from y. */
		if ((status = GetGrpInfo8 (&y, &exptime, &expend, &imset_ok)))
		    return (status);
		if (!imset_ok) {
		    freeSingleGroup (&y);
		    continue;
		}

		sumexptime += exptime;

		/* Weight y by the exposure time. */
		if (sts->fluxcorr == COMPLETE) {
		    if ((status = multk2d (&y, (float)exptime)))
			return (status);
		}
		SquareErr (&y);				/* operate on y */

		/* Add current imset to sum (i.e. add y to x).  This differs
		    from add2d in that RptSum adds variances, rather than
		    adding errors in quadrature.
		*/
		if ((status = RptSum (&x, &y)))
		    return (status);

		if (sts->printtime) {
		    sprintf (message, "imset %d added", iextver);
		    TimeStamp (message, sts->rootname);
		} else if (sts->verbose) {
		    printf ("         imset %d added\n", iextver);
		}
		freeSingleGroup (&y);
		nrptexp++;

	    } else {
		done = 1;
	    }
	}

	/* Take the square root of variance to convert back to errors. */
	SqrtErr (&x);

	/* Divide by the sum of the exposure times. */
	if (sts->fluxcorr == COMPLETE) {
	    if (sumexptime <= 0.) {
		printf ("ERROR    Sum of EXPTIME = %.6g\n", sumexptime);
		return (GENERIC_ERROR_CODE);
	    }
	    if ((status = multk2d (&x, (float)(1./sumexptime))))
		return (status);
	}

	/* Compute statistics and update keywords in output headers. */
	if (sts->statcorr == PERFORM) {
	    printf ("\n");
	    PrSwitch ("statflag", PERFORM);
	    if ((status = doStat (&x, sts->sdqflags)))
		return (status);
	    PrSwitch ("statflag", COMPLETE);
	    if (sts->printtime)
		TimeStamp ("STATFLAG complete", sts->rootname);
	}

	if (nrptexp != sts->nrptexp) {
	    printf ("Warning  Order %d has %d imset", sporder, nrptexp);
	    if (nrptexp > 1)
		printf ("s");
	    printf (" instead of %d.\n", sts->nrptexp);
	}

	/* Update header info in output.  Note that nrptexp is passed
	   rather than sts->nrptexp, so the value can be different for
	   each output imset.
	*/
	if ((status = PutHdrInfo8 (&x, sts->statcorr,
                                   sumexptime, expend, nrptexp)))
	    return (status);

	/* Update CAL_VER and FILENAME, then write output file. */
	UCalVer (x.globalhdr);
	UFilename (sts->output, x.globalhdr);
	putSingleGroup (sts->output, oextver, &x, 0);
	if (hstio_err())
	    return (GENERIC_ERROR_CODE);
	freeSingleGroup (&x);

	if (sts->printtime)
	    TimeStamp ("Output written to disk", sts->rootname);

	if (sts->printtime)
	    free (message);
	return (0);
}

/* This routine searches the orders array for the current spectral order
   number.  The function value will be the EXTVER number if the spectral
   order was found, and it will be -1 if it was not found.

   Note that the EXTVER number is one indexed, but the orders array is
   zero indexed.  Therefore, if the spectral order is found at orders[n],
   then n+1 will be returned.

   Searching starts at array index startextver - 1.
*/

static int FindExtver (int *orders, int nimages, int startextver, int sporder) {

/* arguments:
int orders[]      i: array of spectral order numbers
int nimages       i: size of orders array
int startextver   i: EXTVER number at which seach is to begin
int sporder       i: spectral order number to be found in orders
*/

	int i;

	if (startextver < 1)
	    startextver = 1;

	for (i = startextver - 1;  i < nimages;  i++) {
	    if (orders[i] == sporder)
		return (i + 1);
	}

	return (-1);
}

/* This routine gets the exposure time and exposure end time from
   the current SCI extension.  The accumulated exposure time is then
   incremented by the current exposure time.
*/

static int GetGrpInfo8 (SingleGroup *in, double *exptime, double *expend,
		int *imset_ok) {

/* arguments:
SingleGroup *in   i: current imset
double *exptime   o: exposure time for current imset
double *expend    o: last exposure end time read from input file
int *imset_ok     o: 1 is OK, 0 means the current imset either has zero exposure
                     time or has constant pixel values
*/

	int status;
	int no_default = 0;	/* missing keyword is fatal error */
	int use_def = 1;	/* use default if keyword is missing */
	Bool value;		/* value of IMSET_OK keyword */

	if ((status = Get_KeyD (&in->sci.hdr, "EXPTIME",
                                no_default, 0., exptime)))
	    return (status);
	if ((status = Get_KeyD (&in->sci.hdr, "EXPEND",
                                no_default, 0., expend)))
	    return (status);

	/* Get IMSET_OK, which will be false if the current imset was
	   flagged as having zero exptime or constant pixel values.
	*/
	if ((status = Get_KeyB (&in->sci.hdr, "IMSET_OK", use_def, True,
			&value)) != 0)
	    return status;
	if (value)
	    *imset_ok = 1;
	else
	    *imset_ok = 0;

	return (0);
}

/* This routine adds history info and updates RPTCORR and NEXTEND in the
   the primary header.  It also updates EXPTIME, EXPEND, and NCOMBINE in
   the SCI extension header.
*/

static int PutHdrInfo8 (SingleGroup *out,
		int statcorr, double sumexptime, double expend, int nrptexp) {

/* arguments:
SingleGroup *out  i: current imset
int statcorr      i: statistics flag
double *exptime   i: accumulated exposure time
double *expend    i: last exposure end time read from input file
int nrptexp       i: the number of imsets that were combined
*/

	int status;

	/* Set the switch to indicate that rptcorr has been done. */
	if ((status = Put_KeyS (out->globalhdr, "RPTCORR", "COMPLETE",
                                "add individual repeat observations")))
	    return (status);

	/* Write history records.  Note that these will only be written
	   to the header once, even though we call this for every output
	   imset, because the primary header is not rewritten after it
	   has been written once.
	*/

	addHistoryKw (out->globalhdr, "RPTCORR complete");
	if (hstio_err())
	    return (HEADER_PROBLEM);

	if (statcorr == PERFORM) {
	    addHistoryKw (out->globalhdr, "Statistics computed after rptcorr.");
	    if (hstio_err())
		return (HEADER_PROBLEM);
	}

	/* Update NEXTEND in primary header, to indicate only one imset. */
	if ((status = Put_KeyI (out->globalhdr, "NEXTEND", EXT_PER_GROUP,
                                "number of extensions")))
	    return (status);

	/* Update NCOMBINE in the SCI extension header, to tell how many
	   imsets were combined into this one output imset.
	*/
	if ((status = Put_KeyI (&out->sci.hdr, "NCOMBINE", nrptexp,
                                "number of imsets combined")))
	    return (status);

	/* Update exposure time info in SCI extension header. */
	if ((status = Put_KeyD (&out->sci.hdr, "EXPTIME", sumexptime,
                                "exposure time")))
	    return (status);
	if ((status = Put_KeyD (&out->sci.hdr, "EXPEND", expend,
                                "exposure end time")))
	    return (status);

	return (0);
}

/* convert error to variance */

static void SquareErr (SingleGroup *x) {

	int i, j;

	for (j = 0;  j < x->err.data.ny;  j++) {
	    for (i = 0;  i < x->err.data.nx;  i++) {
		Pix (x->err.data, i, j) =
			Pix (x->err.data, i, j) * Pix (x->err.data, i, j);
	    }
	}
}

/* convert variance to error */

static void SqrtErr (SingleGroup *x) {

	int i, j;

	for (j = 0;  j < x->err.data.ny;  j++) {
	    for (i = 0;  i < x->err.data.nx;  i++) {
		Pix (x->err.data, i, j) = sqrt (Pix (x->err.data, i, j));
	    }
	}
}


/* Add two SingleGroup triplets, leaving the result in the first.

   (*a) += (*b)

   The science data arrays are added together; the error arrays are
   added; the data quality arrays are ORed.

   This differs from add2d in that RptSum assumes the error arrays
   contain variance rather than standard deviations, so those values
   will simply be added.
*/

static int RptSum (SingleGroup *a, SingleGroup *b) {

/* arguments:
SingleGroup *a   io: input data; output sum
SingleGroup *b   i: second input data
*/

	int i, j;
	short dqa, dqb, dqab;	/* data quality for a, b, combined */

	if (a->sci.data.nx != b->sci.data.nx ||
	    a->sci.data.ny != b->sci.data.ny)
	    return (SIZE_MISMATCH);

	/* science data */
	for (j = 0;  j < a->sci.data.ny;  j++) {
	    for (i = 0;  i < a->sci.data.nx;  i++) {
		Pix (a->sci.data, i, j) =
			Pix (a->sci.data, i, j) + Pix (b->sci.data, i, j);
	    }
	}

	/* error array (actually contains variance) */
	for (j = 0;  j < a->err.data.ny;  j++) {
	    for (i = 0;  i < a->err.data.nx;  i++) {
		Pix (a->err.data, i, j) =
			Pix (a->err.data, i, j) + Pix (b->err.data, i, j);
	    }
	}

	/* data quality */
	for (j = 0;  j < a->dq.data.ny;  j++) {
	    for (i = 0;  i < a->dq.data.nx;  i++) {
		dqa = DQPix (a->dq.data, i, j);
		dqb = DQPix (b->dq.data, i, j);
		dqab = dqa | dqb;
		DQSetPix (a->dq.data, i, j, dqab);
	    }
	}

	return (0);
}
