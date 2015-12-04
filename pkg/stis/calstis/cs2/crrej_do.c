# include 	<stdio.h>
# include 	<stdlib.h>
# include 	<string.h>
# include 	"c_iraf.h"
# include 	"hstio.h"

# include	"stis.h"
# include	"cs2.h"
# include	"calstis2.h"
# include	"stisdef.h"

static int countImsets (IRAFPointer);

/*  crrej_do -- Perform the cosmic ray rejection for STIS images

  Description:
  ------------
  This is mostly a file bookkeeping routine for the cosmic ray rejection task.
  It takes care of input/output files open/close, check for dimensions, read/
  write data from/to files, allocate memory spaces etc.

  Date		Author		Description
  ----		------		-----------
  06-May-1996  J.-C. Hsu	Adapt from the SPP code crrej_do.x
  03-Sep-1999  Phil Hodge	Remove niter from calling sequence.
  10-Feb-2000  Phil Hodge	Add exp_range to the calling sequence of
				cr_scaling; and update EXPSTART, etc.;
				replace put_key routines with putKey[];
				check for out of memory after calloc;
				remove iperr (since it was never used);
				return int instead of void,
				and replace exit with return.
  12-Dec-2008  Phil Hodge	Add a call to checkImsetOK.
  22-May-2012  Phil Hodge	Add function countImsets, and allocate memory
				instead of using MAX_FILES for ipsci, etc.
*/
int crrej_do (IRAFPointer tpin, char *outfile, clpar *par, int newpar[],
		float sigma[])
{

	int		max_files;	/* size for arrays ipsci, etc. */
	IODescPtr	*ipsci;		/* science image descriptor */
	IODescPtr	*ipdq;		/* data quality image descriptor */
	float		*skyval;	/* background DN values */
	float		*efac;		/* exposure factors */
	float		*tfac;		/* temperature factors */
	float		*noise;		/* readout noise */
	float		*gain;		/* A-to-D gain factors */
	int		*grp;
	char		**imgname;

	float		exptot[MAX_NEXP];

	float           temptot;
	float		texpt;
	/* exp_range[0] is the minimum value of EXPSTART, and
	   exp_range[1] is the maximum value of EXPEND (both are MJD). */
	double		exp_range[2];
	int		nimgs;
	int		niter;
	SingleGroup	sg;
	
	int		dim_x, dim_y;	/* image dimensions */
	int		i, j, n;	/* loop indices */
	float		*efacsum, *work;
	int		nrej;		/* total number of rejected pixels */
	float		skysum;		/* total sky level */
	int		status;
	int		imset_ok;

	int 	crrej_check (IRAFPointer, clpar *, int [],  char **,
				int [], IODescPtr [],
				IODescPtr [], float [], float [], int *,
				int *, int *);
	int 	cr_scaling (char *, IODescPtr [], char **,
				int, float [], float [], double []);
        int	o_cal2_in (clpar *, int [], int, float,   int *, float []);
	int 	crrej_sky (char *, IODescPtr [], IODescPtr [], int,   float []);
	void 	cr_history (SingleGroup *, clpar *);
	int 	crrej_init (IODescPtr [], clpar *, int, int, int,
				float [], float [], float [], float [],
				FloatTwoDArray *, FloatTwoDArray *,   float *);
	int 	crrej_loop (IODescPtr [], IODescPtr [], char **,
				int [], int, clpar *, int, int, int, float [],
				float [], float [], float [], float [],
				FloatTwoDArray *, FloatTwoDArray *, float *,
				ShortTwoDArray *, int *);

/* -------------------------------- begin ---------------------------------- */

	max_files = countImsets (tpin);
	if (max_files < 0)
	    return 2;
	printf ("Total number of input image sets = %d\n", max_files);

	ipsci = calloc (max_files, sizeof(IODescPtr));
	ipdq = calloc (max_files, sizeof(IODescPtr));
	skyval = calloc (max_files, sizeof(float));
	efac = calloc (max_files, sizeof(float));
	tfac = calloc (max_files, sizeof(float));
	noise = calloc (max_files, sizeof(float));
	gain = calloc (max_files, sizeof(float));
	grp = calloc (max_files, sizeof(int));
	imgname = calloc (max_files, sizeof(char *));
	if (ipsci == NULL || ipdq == NULL || skyval == NULL ||
	    efac == NULL || tfac == NULL || noise == NULL ||
	    gain == NULL || grp == NULL || imgname == NULL) {
	    printf ("ERROR    out of memory in crrej_do\n");
	    return 2;
	}
	for (n = 0;  n < max_files;  n++) {
	    imgname[n] = calloc (STIS_FNAME+1, sizeof(char));
	    if (imgname[n] == NULL) {
		printf ("ERROR    out of memory in crrej_do\n");
		return 2;
	    }
	}

	/* open input files and temporary files, check the parameters */
	if (crrej_check (tpin, par, newpar, imgname, grp, ipsci, ipdq,
			  noise, gain, &dim_x, &dim_y, &nimgs))
	    return (2);
	
	/* calculate the scaling factors due to different exposure time */
	strcpy (par->expname, "EXPTIME");
	if (cr_scaling (par->expname, ipsci, imgname, nimgs,
			efac, tfac, exp_range))
	    return (2);

	/* calculate the total exposure time */
	exptot[0] = 0.;
	temptot = 0.;
	for (n = 0; n < nimgs; ++n) {
	    exptot[0] += efac[n];
	    temptot += tfac[n]*efac[n];
	}

	/* read in the parameters */
        if (o_cal2_in (par, newpar, nimgs, exptot[0],   &niter, sigma))
	    return (2);

	/* for the case of all images have zero exposure time, use equal
	   exposure time of 1. */
	if (exptot[0] == 0.) {
            for (n = 0; n < nimgs; ++n) {
                efac[n] = 1.;
	    }
	    texpt = (float) nimgs;
	} else {
	    texpt = exptot[0];
	}


	/* allocate array space */
	efacsum = calloc (dim_x*dim_y, sizeof(float));
	work    = calloc (nimgs*dim_x, sizeof(float));
	if (efacsum == NULL || work == NULL) {
	    printf ("ERROR    out of memory in crrej_do\n");
	    return (2);
	}

	/* calculate the sky levels */
	if (crrej_sky (par->sky, ipsci, ipdq, nimgs,   skyval))
	    return (2);
	if (par->verbose) {
	    for (n = 0; n < nimgs; ++n)
		printf ("sky of '%s[sci,%d]' is %0.3f DN\n", imgname[n],
			grp[n], skyval[n]);
	}

	/* use the first input image to set up the data structure */
	initSingleGroup (&sg);

	/* Close these, because we're about to reopen them. */
	if (ipsci[0] != NULL)
	    closeImage (ipsci[0]);
	if (ipdq[0] != NULL)
	    closeImage (ipdq[0]);

	getSingleGroup (imgname[0], 1, &sg);
	if (hstio_err()) {
	    printf ("ERROR    %s\n", hstio_errmsg());
	    return (2);
	}
	n = 0;
	while ((status = checkImsetOK (imgname[0], n, &imset_ok)) == 0) {
	    if (imset_ok) {
		ipsci[0] = openInputImage (imgname[0], "SCI", n+1);
		ipdq[0] = openInputImage (imgname[0], "DQ", n+1);
		break;
	    }
	    n++;
	}

	/* compute the initial pixel values to be used to compared against all
	   images. */
	if (crrej_init (ipsci, par, nimgs, dim_x, dim_y,
			noise, gain, efac, skyval,
			&sg.sci.data, &sg.err.data, work))
	    return (2);

	/* do the iterative cosmic ray rejection calculations */
	if (crrej_loop (ipsci, ipdq, imgname, grp, nimgs, par, niter,
			dim_x, dim_y,
			sigma, noise, gain, efac, skyval,
			&sg.sci.data, &sg.err.data, efacsum, &sg.dq.data,
			&nrej))
	    return (2);

	/* calculate the total sky and force it to be non-negative */
	skysum = 0.;
	for (n = 0; n < nimgs; ++n) {
	    skysum += skyval[n];
	}

	/* write to the output image */
	for (j = 0; j < dim_y; ++j) {
	    for (i = 0; i < dim_x; ++i) {
		PPix(&sg.sci.data,i,j) = PPix(&sg.sci.data,i,j)*texpt + skysum;
		PPix(&sg.err.data,i,j) *= texpt;
	    }
	}

	/* update the exposure time of the output image */
	putKeyF (sg.globalhdr, "TEXPTIME", exptot[0], "");
	putKeyF (sg.globalhdr, "SKYSUM", skysum, "Total sky level (DN)");
	if (exptot[0] > 0.) putKeyF (sg.globalhdr, "REJ_RATE",
		(float)nrej/exptot[0], "Cosmic ray impact rate (pixels/sec)");
	putKeyF (&sg.sci.hdr, "EXPTIME", exptot[0], "");
	if ((exptot[0] > 0.) && (temptot > 0.)) {
	    putKeyF(&sg.sci.hdr, "OCCDHTAV", temptot/exptot[0], "");
	}

	/* update exposure start and stop times */
	putKeyD (sg.globalhdr, "TEXPSTRT", exp_range[0], "");
	putKeyD (sg.globalhdr, "TEXPEND", exp_range[1], "");
	putKeyD (&sg.sci.hdr, "EXPSTART", exp_range[0], "");
	putKeyD (&sg.sci.hdr, "EXPEND", exp_range[1], "");

	/* record parameters to the output file */
	cr_history (&sg, par);
	putKeyI (&sg.sci.hdr, "NCOMBINE", nimgs, "");
	UFilename (outfile, sg.globalhdr);

	putSingleGroup (outfile, 1, &sg, 0);
	if (hstio_err()) {
	    printf ("ERROR    %s\n", hstio_errmsg());
	    return (2);
	}
	freeSingleGroup (&sg);
	
	/* deallocate memory */
	free (efacsum);
	free (work);

	free (ipsci);
	free (ipdq);
	free (skyval);
	free (efac);
	free (tfac);
	free (noise);
	free (gain);
	free (grp);
	for (n = 0;  n < max_files;  n++)
	    free (imgname[n]);
	free (imgname);

	return (0);
}

/* This function returns a count of the total number of image sets in the
   set of input files.  The function value will be -1 if any of the input
   images can't be opened.
*/

static int countImsets (IRAFPointer tpin) {

	IODescPtr ip;		/* file handle for an input file */
	Hdr prihdr;		/* primary header */
	char fname[STIS_FNAME+1];
	int n;			/* loop index */
	int nextend;		/* number of extensions */
	int nimsets = 0;	/* total number of image sets */
	int status;

	initHdr (&prihdr);

	/* rewind the image template */
	c_imtrew (tpin);
	
	/* loop over all input files */
	for (n = 0;  n < c_imtlen(tpin);  n++) {

	    /* read the next input image name in the template list */
	    c_imtgetim (tpin, fname, STIS_FNAME);

	    /* open the primary header */
	    ip = openInputImage (fname, "", 0);
	    if (hstio_err()) {
		printf ("ERROR    HSTIO error %s\n", hstio_errmsg());
		return -1;
	    }
	
	    getHeader (ip, &prihdr);
	    closeImage (ip);

	    if ((status = getKeyI (&prihdr, "NEXTEND", &nextend)) != 0)
		nextend = 0;
	    nimsets += (nextend / EXT_PER_GROUP);
	}
	freeHdr (&prihdr);

	return nimsets;
}
