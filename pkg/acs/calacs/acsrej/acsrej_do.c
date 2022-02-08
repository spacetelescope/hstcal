# include   <stdio.h>
# include   <stdlib.h>
# include   <string.h>

#include "hstcal.h"
# include   "hstio.h"

# include   "acs.h"
# include   "acsrej.h"
# include   "hstcalerr.h"
# include   "rej.h"

static void closeSciDq (int, IODescPtr [], IODescPtr [], IODescPtr [], clpar *);


/*  acsrej_do -- Perform the cosmic ray rejection for ACS images

  Description:
  ------------
  This is mostly a file bookkeeping routine for the cosmic ray rejection task.
  It takes care of input/output files open/close, check for dimensions, read/
  write data from/to files, allocate memory spaces etc.

  Date          Author          Description
  ----          ------          -----------
  06-May-1996   J.-C. Hsu       Adapt from the SPP code crrej_do.x
  05-Aug-1998   W. Hack         Modified to handle ACS data
  11-Feb-1999   W. Hack         EXPTIME now in Pri. Hdr.
  18-Mar-1999   W.J. Hack       Revised to read EXPTIMEs from Primary headers
                                for cr_scaling using tpin directly
  07-Jul-1999   W.J. Hack       EXPNAME values modified for all extensions
                                in a SingleGroup.
  14-Sep-1999   W.J. Hack       Cleaned up SHADCORR usage. Added check for max
                                number of files here.
  24-Feb-2003   W.J. Hack       Updated REJ_RATE to use 'texpt' as a safe
                                value when EXPTIME=0 for all members.
  02-Sep-2015   P.L. Lim        Initialize nrej and texpt like CALWF3.
                                nrej is used to calculate REJ_RATE.
  15-Sep-2015   P.L. Lim        Corrected unit verbiage in header comments
                                (sky is electrons, not DN).
                                Added doc. Cleaned up codes.
  27-Feb-2019   M.D. DeLaPena   Read the IMAGETYP keyword from Primary of first image
                                as this info is needed by acsrej_loop.
  30-Apr-2021   M.D. DeLaPena   Read the FLASHDUR keyword from Primary of all the input
                                images in order to compute the cumulative value and update
                                FLASHDUR in the blv_tmp/blc_tmp file(s).
*/
int acsrej_do (IRAFPointer tpin, char *outfile, char *mtype, clpar *par,
               int newpar[]) {
    extern int  status;

    IODescPtr   ipsci[MAX_FILES];   /* science image descriptor */
    IODescPtr   iperr[MAX_FILES];   /* error image descriptor */
    IODescPtr   ipdq[MAX_FILES];    /* data quality image descriptor */
    float       skyval[MAX_FILES];  /* background values */
    float       efac[MAX_FILES];    /* exposure factors */
    multiamp    noise;              /* readout noise */
    multiamp    gain;               /* A-to-D gain factors */
    float       exptot;
    float       texpt;
    int         nimgs;
    SingleGroup sg;
    int         niter = 0;
    float       sigma[MAX_ITER];

    Hdr         phdr;               /* primary header */
    int         extver;             /* Current extension being processed */
    int         numext;             /* Number of extensions in each image */
    int         nextend;            /* Number of output extensions */
    char        imgname[MAX_FILES][CHAR_FNAME_LENGTH];
    char        fimage[CHAR_FNAME_LENGTH];  /* Name of first image in list */
    char        root[CHAR_FNAME_LENGTH];    /* ROOTNAME for output CRJ file */
    char        uroot[CHAR_FNAME_LENGTH];   /* Upper case version of rootname */
    char        *shadrefname;
    char        imagetyp[ACS_CBUF];

    int         ext[MAX_FILES];
    int         dim_x, dim_y;       /* image dimensions */
    int         i, j, n;            /* loop indices */
    float       *efacsum, *work;
    int         nrej;               /* total number of rejected pixels */
    float       skysum;             /* total sky level */
    int         logit;
    RefImage    shadref;
    int         shadswitch;
    double      expend, expstart;
    int         non_zero;           /* number of input images with EXPTIME>0 */
    int         found;
    char        imgdefault[CHAR_FNAME_LENGTH];  /* name of first input image with EXPTIME > 0. */
    float       cumFlashDur;
    float       cumDarktime;

    int     GetSwitch (Hdr *, char *, int *);
    int     UpdateSwitch (char *, int, Hdr *, int *);

    void    InitRefImg (RefImage *);
    int     ImgHistory (const RefImage *, Hdr *);
    int     ImgPedigree (RefImage *);

    int     acsrej_check (IRAFPointer, int, int, clpar *, int [],
                          char [][CHAR_FNAME_LENGTH], int [], IODescPtr [],
                          IODescPtr [], IODescPtr [], multiamp *, multiamp *,
                          int *, int *, int, float [], float *, float *);
    int     cr_scaling (char *, IRAFPointer, float [], int *, double *,
                        double *);
    int     rejpar_in(clpar *, int [], int, float, int *, float []);
    void    acsrej_sky (char *, IODescPtr [], IODescPtr [], int, short,
                        float []);
    void    cr_history (SingleGroup *, clpar *, int);
    int     acsrej_init (IODescPtr [], IODescPtr [], clpar *, int, int, int,
                         float [], float [], SingleGroup *, float *);
    int     acsrej_loop (IODescPtr [], IODescPtr [], IODescPtr [],
                         char [][CHAR_FNAME_LENGTH], int [], int, clpar *, int, int,
                         int, float [], multiamp, multiamp, float [], float [],
                         FloatTwoDArray *, FloatTwoDArray *, float *,
                         ShortTwoDArray *, int *, char *, char *);
    int     PutKeyFlt (Hdr *, char *, float, char *);
    int     PutKeyDbl (Hdr *, char *, double, char *);
    int     PutKeyStr (Hdr *, char *, char *, char *);
    int     GetKeyStr (Hdr *, char *, int, char *, char *, int);
    int     PutKeyInt (Hdr *, char *, int, char *);
    int     GetKeyInt (Hdr *, char *, int, int, int *);
    void    UFilename (char *, Hdr *);
    void    UMemType (char *, Hdr *);
    void    UExpname (char *, Hdr *);
    int     LoadHdr (char *, Hdr *);
    void    UpperAll (char *, char *, int);
    void    TimeStamp (char *, char *);
    void    WhichError (int);
    void    PrSwitch (char *, int);
    void    FindAsnRoot (char *, char *);
    void    initmulti (multiamp *);

    nrej = 0;

    /* -------------------------------- begin ------------------------------- */
    /* Initialize necessary structures */
    InitRefImg (&shadref);
    root[0] = '\0';
    uroot[0] = '\0';
    initmulti (&noise);
    initmulti (&gain);

    numext = 0;
    nextend = 0;

    /* Since CR-SPLIT images are in separate files, we need to
       combine the same chip's exposure from each file.  Therefore
       we will loop over each extension in the first image, determine
       what chip that corresponds to, and get the same chip from the
       rest of the images (which could be in any arbitrary
       extension in each of the images). */

    /* First, let's determine how many extensions/chips in each file */
    c_imtgetim (tpin, fimage, CHAR_FNAME_LENGTH);

    if (LoadHdr (fimage, &phdr))
        return (status = ERROR_RETURN);

    if (GetKeyInt (&phdr, "NEXTEND", NO_DEFAULT, 0, &nextend) == 0)
        numext = nextend / EXT_PER_GROUP;
    else
        numext = 1;

    /* While the primary header is available, get the image type which is now 
       needed in acsrej_loop.c.  The interest here are the specific cases of 
       combining BIAS or DARK images to generate calibration images, so only 
       need to read the IMAGETYP keyword from the first image. */
    if (GetKeyStr (&phdr, "IMAGETYP", NO_DEFAULT, "", imagetyp, ACS_CBUF)) {
        trlkwerr ("IMAGETYP", fimage);
        return(status = KEYWORD_MISSING);
    }

    shadswitch = 0;
    /* Check to see if SHADCORR was set to PERFORM in image header */
    if (GetSwitch (&phdr, "SHADCORR", &shadswitch) )
        return(status);

    /* If shadcorr was set either by the user on the command line
       or in the image header, initialize shadcorr processing. */
    if (par->shadcorr == PERFORM || shadswitch == PERFORM) {

        /* Use par->shadcorr as switch for performing shading correction */
        par->shadcorr = PERFORM;
        shadrefname = calloc(CHAR_FNAME_LENGTH, sizeof(char));

        if (GetKeyStr (&phdr, "SHADFILE", NO_DEFAULT, "", shadrefname, CHAR_FNAME_LENGTH) )
            return(status);
        strcpy (shadref.name, shadrefname);

        /* Read in PEDIGREE and DESCRIPTION for SHADFILE */
        if (ImgPedigree (&shadref))
            return (status);

        /* If a DUMMY shadfile was specified, turn off shadcorr */
        if (shadref.goodPedigree == DUMMY)
            par->shadcorr = OMIT;

        free (shadrefname);
    }

    freeHdr (&phdr);

    /* Initialize efac */
    for (n = 0; n < MAX_FILES; n++)
        efac[n] = 1.0;

    /* calculate the scaling factors due to different exposure time */
    strcpy (par->expname, "EXPTIME");
    if (cr_scaling (par->expname, tpin, efac, &nimgs, &expend, &expstart)) {
        WhichError (status);
        return (status);
    }

    /* make sure there is more than one image to process */
    if (nimgs < 2) {
        trlmessage ("Needs more than one input image.");
        return (status = NOTHING_TO_DO);
    }

    /* calculate the total exposure time */
    exptot = 0.;
    texpt = 0.;
    non_zero = 0;
    for (n = 0; n < nimgs; ++n) {
        exptot += efac[n];
        /* Count how many inputs have non-zero(valid) EXPTIME */
        if (efac[n] > 0.)
            non_zero++;
    }

    /* for the case of all images have zero exposure time (BIAS),
       use equal exposure time of 1. */
    if (exptot == 0.) {
		for (n = 0; n < nimgs; ++n) {
			efac[n] = 1.;
		}
		texpt = (float) nimgs;
		non_zero = nimgs;
	} else {
		texpt = exptot;
	}

	/* Now, start the extver loop. */
	for (extver = 1; extver <= numext; extver++) {
		if (par->printtime) {
			TimeStamp ("Start cosmic ray rejection", "");
		}

		/* open input files and temporary files, check the parameters */
		if (acsrej_check (tpin, extver, numext, par, newpar, imgname, ext,
						  ipsci, iperr, ipdq, &noise, &gain, &dim_x, &dim_y,
						  nimgs, efac, &cumFlashDur, &cumDarktime)) {
			WhichError (status);
			return(status);
		}

		/* Now that we have read in SHADCORR, report if it will be performed */
		PrSwitch ("shadcorr", par->shadcorr);

		/* read in the parameters.
		   NOTE: newpar is useless after this */
		if ( rejpar_in (par, newpar, nimgs, exptot, &niter, sigma) )
			return(status);

		/* allocate array space */
		efacsum = calloc (dim_x*dim_y, sizeof(float));
		work    = calloc (nimgs*dim_x, sizeof(float));

		/* calculate the sky levels */
		acsrej_sky (par->sky, ipsci, ipdq, nimgs, par->badinpdq, skyval);
		if (status != ACS_OK) {
			WhichError (status);
			return (status);
		}
		if (par->verbose) {
			for (n = 0; n < nimgs; n++) {
				sprintf (MsgText, "sky of '%s[sci,%d]' is %0.3f electrons",
						 imgname[n], ext[n], skyval[n]);
				trlmessage (MsgText);
			}
		}

		/* use the first input image to set up the data structure */
		initSingleGroup (&sg);

		/* Find the first image in the input list which has an
		   EXPTIME > 0. To use for initializing the output SingleGroup.
		*/
		found = 0;
		n = 0;
		/* By default, simply use the first one, so initialize accordingly. */
		strcpy (imgdefault, imgname[0]);
		do {
			if (efac[n] > 0.) {
				strcpy(imgdefault, imgname[n]);
				found = 1;
			}
			n++;
		} while (found == 0);

		getSingleGroup (imgdefault, extver, &sg);

		if (non_zero > 1) {
			/* compute the initial pixel values to be used to compare against
			   all images. */
			if (non_zero < nimgs) {
				trlwarn ("Some input exposures had EXPTIME = 0.");
			}
			if (acsrej_init (ipsci, ipdq, par, nimgs, dim_x, dim_y, efac,
							 skyval, &sg, work)) {
				WhichError(status);
				closeSciDq(nimgs, ipsci, iperr, ipdq, par);
				return (status);
			}

			if (par->printtime)
				TimeStamp ("Calculated initial guess for extension", "");

			/* do the iterative cosmic ray rejection calculations */
			if (acsrej_loop (ipsci, iperr, ipdq, imgname, ext, nimgs, par,
							 niter, dim_x, dim_y, sigma, noise, gain, efac,
							 skyval, &sg.sci.data, &sg.err.data, efacsum,
							 &sg.dq.data, &nrej, shadref.name, imagetyp)) {
				WhichError(status);
				closeSciDq(nimgs, ipsci, iperr, ipdq, par);
				return (status);
			}
		} else {
			trlwarn ("Cosmic-ray rejection NOT performed!");
			if (non_zero > 0) {
				trlwarn ("Some input exposures had EXPTIME = 0.");
				trlwarn ("Output product will not be cosmic-ray cleaned!");
			}
        } /* End if(non_zero) block */

        /* must close all images, now that we are done reading them */
        closeSciDq(nimgs, ipsci, iperr, ipdq, par);

        /* calculate the total sky (electrons)... */
        skysum = 0.;
        for (n = 0; n < nimgs; ++n) {
            skysum += skyval[n];
        }
        /* ... and force it to be non-negative */
        if (skysum < 0.) skysum = 0.;

        if (par->printtime){
            if (non_zero > 1) {
                TimeStamp ("Finished detecting cosmic rays on extension", "");
            } else {
                TimeStamp ("Done checking this extension", "");
            }
        }

        /* write to the output image */
        if (non_zero > 0) {
            for (j = 0; j < dim_y; ++j) {
                for (i = 0; i < dim_x; ++i) {
                    /* SCI: Convert e/s to electrons, add sky back in */
                    PPix(&sg.sci.data, i, j) = (PPix(&sg.sci.data, i, j) *
                                                texpt) + skysum;
                }
            }
            /* Set at least one pixel to a different value to insure
               that an image array actually gets produced.  This was done
               so the result of a test where the input ERR are all zero will
               produce an actual output ERR array. */
            PPix(&sg.err.data, 0, 0) = -1.;
        } else {
            /* DQ values here do not make sense, but I guess it is not
               a real issue because the output is dummy anyway. */
            for (j = 0; j < dim_y; ++j) {
                for (i = 0; i < dim_x; ++i) {
                    PPix(&sg.sci.data, i, j) = par->fillval;  /* always 0 */
                    PPix(&sg.err.data, i, j) = 0.;
                    /* Set DQ value to one which will always be
                       considered BAD.
                       1 = Reed-Solomon decoding error */
                    PPix(&sg.dq.data, i, j) = 1;
                }
            }
            /* Set at least one pixel to a different value to insure
               that an image array actually gets produced.
               8 = Masked by aperture feature */
            PPix(&sg.err.data, 0, 0) = -1.;
            PPix(&sg.dq.data, 0, 0) = 8;
        }

        /* record parameters to the output file */

        /* update the exposure time of the output images.
           exptot is just direct sum of all EXPTIMES.
           texpt is the actual sum used in calculations.
        */
        PutKeyFlt (sg.globalhdr, "TEXPTIME", exptot, "");
        PutKeyFlt (sg.globalhdr, "EXPTIME", exptot, "");
        PutKeyDbl (sg.globalhdr, "EXPSTART", expstart,
                   "computed exposure start time (Modified Julian Date)");
        PutKeyDbl (sg.globalhdr, "EXPEND", expend,
                   "computed exposure end time (Modified Julian Date)");
        PutKeyFlt (sg.globalhdr, "SKYSUM", skysum,
                   "Total sky level (electrons)");
        PutKeyFlt (sg.globalhdr, "REJ_RATE", (float)nrej / texpt,
                   "Cosmic ray impact rate (pixels/sec)");

        /* Cannot update the comment portion of the FITS card for an existing keyword
           with the high-level interface so using lower-level functions for the update.
        */
        PutKeyFlt (sg.globalhdr, "FLASHDUR", cumFlashDur, "");
        FitsKw kw = findKw(sg.globalhdr, "FLASHDUR");
        putKwComm(kw, "Cumulative exposure time in seconds");

        PutKeyFlt (sg.globalhdr, "DARKTIME", cumDarktime, "");
        kw = findKw(sg.globalhdr, "DARKTIME");
        putKwComm(kw, "Cumulative dark time in seconds");

        if (par->shadcorr) {
            logit = 0;
            if (UpdateSwitch ("SHADCORR", par->shadcorr, sg.globalhdr, &logit))
                return (status);
            PrSwitch ("shadcorr", COMPLETE);
            /* Records SHADFILE information in header comments... */
            if (logit) {
                if (ImgHistory (&shadref, sg.globalhdr))
                    return (status);
            }
        }

        cr_history (&sg, par, nextend);
        PutKeyInt (&sg.sci.hdr, "NCOMBINE", nimgs, "");
        UFilename (outfile, sg.globalhdr);
        UMemType (mtype, sg.globalhdr);
        FindAsnRoot (outfile, root);
        UpperAll (root, uroot, strlen(root) + 1);
        UExpname (root, &sg.sci.hdr);
        UExpname (root, &sg.err.hdr);
        UExpname (root, &sg.dq.hdr);
        PutKeyStr (sg.globalhdr, "ROOTNAME", uroot,
                   "Rootname of the observation set");

        /* Output CHIP to the same EXTVER as the CHIP ID */
        putSingleGroup (outfile, extver, &sg, 0);
        freeSingleGroup (&sg);

        if (par->printtime)
            TimeStamp ("Finished writing out extension", "");

        /* deallocate memories */
        free (efacsum);
        free (work);
    } /* End of extver loop */

    /* Set status to a value which will be understood by
       CALACS to turn off subsequent processing. */
    if (non_zero == 0)
        status = NO_GOOD_DATA;

    return (status);
}


/* Helper function to clean up image pointers... */
static void closeSciDq(int nimgs, IODescPtr ipsci[], IODescPtr iperr[],
                       IODescPtr ipdq[], clpar *par) {
    int n;

    /* must close all images, now that we are done reading them */
    for (n = 0; n < nimgs; ++n) {
        closeImage (ipsci[n]);
        closeImage (iperr[n]);
        closeImage (ipdq[n]);
    }
}
