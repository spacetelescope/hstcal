# include   <stdio.h>
# include   <stdlib.h>
# include   <string.h>

# include "hstio.h"

# include   "wf3.h"
# include   "wf3rej.h"
# include   "hstcalerr.h"
# include   "wf3info.h"
# include   "rej.h"

static void closeSciDq (int, IODescPtr [], IODescPtr [], clpar *);

/*  rej_do -- Perform the cosmic ray rejection for WFC3 images

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
    14-Sep-1999   W.J. Hack       Cleaned up SHADCORR usage. Added check for max
                                  number of files here.
    29-Aug-2000   H.A. Bushouse   Revised for WFC3 use.
    18-Jun-2002   H. Bushouse	  Modified to compute correct total dark time
				  of images being combined and update the
				  EXPSTART keyword to reflect that time
				  (following CALACS changes).
    27-Oct-2003   H. Bushouse	  Modified to avoid divide-by-zero when total
				  exptime for all members is zero. Also put in
				  various upgrades to handle one or more input
				  images with EXPTIME = 0.
				  (CALACS changes).
    16-May-2008   H. Bushouse     Pass detector to cr_history, for identifying
                                  WFC3 IR images.
    26-Aug-2008   H. Bushouse     Set non_zero=nimgs when resetting efrac for 
                                  the case where all inputs have 0 exptime.
                                  This will allow bias images to process.
    14-Dec-2011   H. Bushouse     Upgraded to pass new bunit array to and from
				  all functions that need it, in order to
				  handle input data that are in count rates.
				  (PR 69969; Trac #814).
  30-Aug-12       M. Sosey        Checks the value of EXPFLAG in all the input image and if
                            any one image contains something other than NORMAL, it reports
                            the value as MIXED in the output crj header. PR #72001

  11-Aug-15      M. Sosey  REJ_RATE was being reported wrong because nrej wasn't initialized
*/

int rej_do (IRAFPointer tpin, char *outfile, char *mtype, clpar *par,
	    int newpar[], int detector) {

    extern int  status;

    IODescPtr   ipsci[MAX_FILES];   /* science image descriptor */
    IODescPtr   ipdq[MAX_FILES];    /* data quality image descriptor */
    float       skyval[MAX_FILES];  /* background DN values */
    float       efac[MAX_FILES];    /* exposure factors */
    DataUnits	bunit[MAX_FILES];   /* image data units */
    multiamp    noise;              /* readout noise */
    multiamp    gain;               /* A-to-D gain factors */
    float       exptot;
    float       texpt;
    int         nimgs;
    SingleGroup sg;
    int         niter = 0;
    float       sigma[MAX_ITER];
    
    Hdr         phdr;               /* primary header */
    int         extver;             /* Current extension being processed*/
    int         numext;             /* Number of extensions in each image */
    int         nextend;            /* Number of output extensions */
    char        imgname[MAX_FILES][SZ_FNAME+1];
    char        fimage[SZ_FNAME+1];  /* Name of first image in list */
    char        root[SZ_FNAME+1];    /* ROOTNAME for output CRJ file */
    char        uroot[SZ_FNAME+1];   /* Upper case version of rootname */
    char        *shadrefname;  

    int         ext[MAX_FILES];
    int         dim_x, dim_y;       /* image dimensions */
    int         i, j, n;            /* loop indices */
    float       *efacsum, *work;
    int         nrej;               /* total number of rejected pixels */
    float       skysum;             /* total sky level */
    int         logit;
    RefImage    shadref;
    int         shadswitch;
    double      expend, expstart;   /* exposure end & start times */
    int		non_zero;	    /* number of input images with EXPTIME>0 */
    int		found;
    char	imgdefault[SZ_FNAME+1]; /* name of first input image with
					   EXPTIME>0 */

    int     GetSwitch (Hdr *, char *, int *);
    int     UpdateSwitch (char *, int, Hdr *, int *);

    void    InitRefImg (RefImage *);
    int     ImgHistory (RefImage *, Hdr *);
    int     ImgPedigree (RefImage *);

    int     rej_check (IRAFPointer, int, int, clpar *, int [],
		       char [][SZ_FNAME+1], int [], IODescPtr [], IODescPtr [],
		       multiamp *, multiamp *, int *, int *, int, char []);
    int     cr_scaling (char *, IRAFPointer, float [], int *, double *,
			double *, DataUnits []);
    int     rejpar_in (clpar *, int [], int, float,   int *, float []);
    void    rej_sky (char *, IODescPtr [], IODescPtr [], int, short, float [],
		     DataUnits [], float []);
    void    cr_history (SingleGroup *, clpar *, int, int);
    int     rej_init (IODescPtr [], IODescPtr [], clpar *, int, int, int,
                multiamp, multiamp, float [], float [], DataUnits [],
		SingleGroup *, float *);
    int     rej_loop (IODescPtr [], IODescPtr [], char [][SZ_FNAME+1],
                int [], int, clpar *, int, int, int, float [], multiamp,
		multiamp, float [], float [], DataUnits [], FloatTwoDArray *,
		FloatTwoDArray *, float *, ShortTwoDArray *, int *, char *);
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
    char    expflagFinal[]="NORMAL"; /*24 char keyword max in header */
    
    nrej=0;
    
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
        extension in each of the images).  
    */

    /* First, let's determine how many extensions/chips in each file */
    c_imtgetim (tpin, fimage, SZ_FNAME);

    if (LoadHdr (fimage, &phdr) )
        return (status = ERROR_RETURN);
        
    if (GetKeyInt (&phdr, "NEXTEND", NO_DEFAULT, 0, &nextend) == 0) 
        numext = nextend / EXT_PER_GROUP;
    else 
        numext = 1;

    /* Check to see if SHADCORR was set to PERFORM in image header */
    shadswitch = 0;
    if (GetSwitch (&phdr, "SHADCORR", &shadswitch) ) 
        return(status);

    /* If shadcorr was set either by the user on the command line
        or in the image header, initialize shadcorr processing.  */
    if (par->shadcorr == PERFORM || shadswitch == PERFORM) {
    
        /* Use par->shadcorr as switch for performing shading correction */
        par->shadcorr = PERFORM;
        shadrefname = calloc(SZ_FNAME+1, sizeof(char));
    
        if (GetKeyStr (&phdr, "SHADFILE", NO_DEFAULT, "", shadrefname,SZ_FNAME))
            return(status);
        strcpy (shadref.name, shadrefname);

        /* Read in PEDIGREE and DESCRIPTION for SHADFILE */
        if (ImgPedigree (&shadref) )
            return (status);

        /* If a DUMMY shadfile was specified, turn off shadcorr */
        if (shadref.goodPedigree == DUMMY) 
            par->shadcorr = OMIT;
        free (shadrefname);
    }

    freeHdr (&phdr);

    /* Initialize efac */
    for (n = 0; n < MAX_FILES; n++) efac[n] = 1.0;
    
    /* Calculate the scaling factors due to different exposure time */
    strcpy (par->expname, "EXPTIME");
    if (cr_scaling (par->expname, tpin, efac, &nimgs, &expend, &expstart,
		    bunit)) {
        WhichError (status);
        return (status);
    }

    /* Make sure there is more than one image to process */
    if (nimgs < 2) {
        trlmessage ("Needs more than one input image.");
        return (status = NOTHING_TO_DO);
    }

    /* Calculate the total exposure time */
    exptot = 0.;
    texpt = 0.;
    
    non_zero = 0;
    for (n = 0; n < nimgs; ++n) {
        exptot += efac[n];
	/* Count how many inputs have non-zero (valid) EXPTIME */
	if (efac[n] > 0.) non_zero++;
    }

    /* For the case of all images have zero exposure time, use equal 
        exposure time of 1. */
    if (exptot == 0.) {
        for (n = 0; n < nimgs; ++n) {
            efac[n] = 1.;
        }
        texpt = (float) nimgs;
	non_zero = nimgs;
    } else {
        texpt = exptot;
    }

    /* Now, start the loop over all extensions in first input image. */
    for (extver = 1; extver <= numext; extver++) { 

        if (par->printtime) {
            TimeStamp ("Start cosmic ray rejection","");
        }

        /* Open input files and temporary files, check the parameters */
        if (rej_check (tpin, extver, numext, par, newpar, imgname, ext,
            ipsci, ipdq, &noise, &gain, &dim_x, &dim_y, nimgs, expflagFinal)) {
            WhichError (status);
            return(status);
        }

        /* Now that we have read in SHADCORR, report if it will be performed */
        PrSwitch ("shadcorr", par->shadcorr);

        /* Read in the parameters */
        if (rejpar_in (par, newpar, nimgs, exptot, &niter, sigma) )
            return(status);

        /* Allocate array space */
        efacsum = calloc (dim_x*dim_y, sizeof(float));
        work    = calloc (nimgs*dim_x, sizeof(float));

        /* Calculate the sky levels */
        rej_sky (par->sky, ipsci, ipdq, nimgs, par->badinpdq, efac,
		 bunit, skyval);
        if (status != WF3_OK) {
            WhichError (status);
            return (status);
        }
        if (par->verbose) {
            for (n = 0; n < nimgs; n++) {
                 sprintf (MsgText, "sky of '%s[sci,%d]' is %0.3f DN",
			  imgname[n], ext[n], skyval[n]);
                 trlmessage (MsgText);
            }
        }

        /* Use the first input image to set up the data structure */
        initSingleGroup (&sg);

	/* Find the first image in the input list that has an EXPTIME>0 to
	** use for initializing the output SingleGroup */
	found = 0;
	n = 0;
	/* By default, simply use the first one, so initialize accordingly */
	strcpy (imgdefault, imgname[0]);
	do {
	   if (efac[n] > 0.) {
	       strcpy (imgdefault, imgname[n]);
	       found = 1;
	   }
	   n++;
	} while (found == 0);

        getSingleGroup (imgdefault, extver, &sg);

	if (non_zero > 1) {

	    if (non_zero < nimgs)
		trlwarn ("Some input exposures had EXPTIME = 0.");

            /* Compute the initial pixel values to be used to compare against
	    ** all images. */
            if (rej_init (ipsci, ipdq, par, nimgs, dim_x, dim_y, 
			  noise, gain, efac, skyval, bunit, &sg, work)) {
                WhichError(status);
                closeSciDq(nimgs, ipsci, ipdq, par);
                return (status);
            }

            if (par->printtime)
                TimeStamp ("Calculated initial guess for extension", "");

            /* Do the iterative cosmic ray rejection calculations */
            if (rej_loop (ipsci, ipdq, imgname, ext, nimgs, par, niter, dim_x,
		      dim_y, sigma, noise, gain, efac, skyval, bunit,
		      &sg.sci.data, &sg.err.data, efacsum, &sg.dq.data, &nrej,
		      shadref.name)){
                WhichError(status);
                closeSciDq(nimgs, ipsci, ipdq, par);
                return (status);
            }

	} else {

	    trlwarn ("Cosmic-ray rejection NOT performed!");
	    if (non_zero > 0) {
		trlwarn ("Some input exposures had EXPTIME = 0.");
		trlwarn ("Output product will not be cosmic-ray cleaned!");
	    } /*else {
		trlwarn ("ALL input exposures had EXPTIME = 0.");
		trlwarn ("Output product will be BLANK!");
	    } */
	} /* End if(non_zero) block */
        
        /* Must close all images, now that we are done reading them */
        closeSciDq (nimgs, ipsci, ipdq, par);

        /* Calculate the total sky ... */ 
        skysum = 0.;
        for (n = 0; n < nimgs; ++n) {
             skysum += skyval[n];
        }
        /* ... and force it to be non-negative */
        if (skysum < 0.) skysum = 0.;
        
        if (par->printtime) {
	    if (non_zero > 1) {
                TimeStamp ("Finished detecting cosmic rays on extension", "");
	    } else {
		TimeStamp ("Done checking this extension","");
	    }
	}

        /* Write to the output image */
	if (non_zero > 0) {
            for (j = 0; j < dim_y; ++j) {
                for (i = 0; i < dim_x; ++i) {
                    PPix(&sg.sci.data,i,j) = PPix(&sg.sci.data,i,j)*texpt +
					     skysum;
                    PPix(&sg.err.data,i,j) *= texpt;
                }
            }
	} else {
            for (j = 0; j < dim_y; ++j) {
                for (i = 0; i < dim_x; ++i) {
                    PPix(&sg.sci.data,i,j) = par->fillval;
                    PPix(&sg.err.data,i,j) = 0.;
		    /* Set DQ value to one which will always be considered BAD*/
		    PPix(&sg.dq.data,i,j) = 1;
                }
            }
	    /* Set at least one pixel to a different value to insure that an
	    ** image array actually gets written by HSTIO */
	    PPix(&sg.err.data,0,0) = -1.;
	    PPix(&sg.dq.data,0,0) = 8;
	}

        /* Update the exposure time of the output images */
        PutKeyFlt (sg.globalhdr, "TEXPTIME", exptot, "");
        PutKeyFlt (sg.globalhdr, "SKYSUM", skysum, "Total sky level (DN)");
        PutKeyDbl (sg.globalhdr, "EXPSTART", expstart, 
		   "computed exposure start time (Modified Julian Date)");
        PutKeyDbl (sg.globalhdr, "EXPEND", expend,
		   "exposure end time (Modified Julian Date)");
           
	/* Upgraded to use 'texpt' as a safe value when EXPTIME=0 for
	** all members. */
        PutKeyFlt (sg.globalhdr, "REJ_RATE", (float)nrej/texpt, 
		   "Cosmic ray impact rate (pixels/sec)");
        PutKeyFlt (sg.globalhdr, "EXPTIME", exptot, "");
        PutKeyStr (sg.globalhdr, "EXPFLAG", expflagFinal,
                "Exposure flag for combined dataset");
        if (par->shadcorr) {
            logit = 0;
            if (UpdateSwitch ("SHADCORR", par->shadcorr, sg.globalhdr, &logit))
                return (status);
            PrSwitch ("shadcorr", COMPLETE);

            if (logit) {
                /*Records SHADFILE information in header comments... */
                if (ImgHistory (&shadref, sg.globalhdr))
                    return (status);
            }
        }

        /* Record parameters to the output file */
        cr_history (&sg, par, nextend, detector);
        PutKeyInt (&sg.sci.hdr, "NCOMBINE", nimgs, "");
        UFilename (outfile, sg.globalhdr);
        UMemType (mtype, sg.globalhdr);
        FindAsnRoot (outfile, root);
        UpperAll (root, uroot, strlen(root)+1 );
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

        /* deallocate memory */
        free (efacsum);
        free (work);
    }

    /* Set status to a value that will be understood by CALWF3 to
    ** turn off subsequent processing */
    if (non_zero == 0) status = NO_GOOD_DATA;

    return (status);

}

/* Helper function to clean up image pointers... */

static void closeSciDq (int nimgs, IODescPtr ipsci[], IODescPtr ipdq[],
			clpar *par) {

    int n;
    
    /* Must close all images, now that we are done reading them */
    for (n = 0; n < nimgs; ++n) {
        closeImage (ipsci[n]);
        closeImage (ipdq[n]);

    }
}

