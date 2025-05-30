# include	<stdio.h>
# include	<string.h>
# include	"hstio.h"
# include	"c_iraf.h"
# include	"xtables.h"

# include	"stis.h"
# include	"cs2.h"
# include	"calstis2.h"

/*  o_cal2_in -- Read CL parameters for the task crrej.

  Description:
  ------------
  Reads CL parameters and does necessary checkings
  
  Input parameters from crrej reference table:
  -------------------------------------------

  "skysub"		Sky levels subtraction scheme
  "crsigmas"            Rejection thresholds
  "crradius"            Radius (in pixels) to propagate the cosmic ray
  "crthresh"            Propagation factor
  "initgues"		Scheme of computing initial-guess image 
  "scalense"		multiplicative noise in percents
  "badinpdq"		Data quality pset
  "crmask"		flag CR-rejected pixels in input files?

  Input parameters from input image primary header:
  ------------------------------------------------


  Date		Author		Description
  ----		------		-----------
  06-May-1996  J.-C. Hsu	adapt from the SPP code crrej_in.x 
  10-Feb-2000  Phil Hodge	Return int instead of void,
				and replace exit with return.
   5-Feb-2004  Phil Hodge	Define exp_in, meanexp, mindiff and diff
				to be double rather than float, and
				rearrange the section where mindiff and
				diff are compared when searching for a
				matching row in the crrejtab.  This is to
				work around the "x86 excess precision"
				problem for Linux machines.
*/
int o_cal2_in (clpar *par, int newpar[], int nimgs, float exptot,   int *niter,
float sigma[])
{
	IRAFPointer	tp;
	IRAFPointer	colptr, colptr1;
	int		i, nrows, nmatch, row;
	int		crsplit_in, crsplit, maxcrsplit;
	double		exp_in, meanexp, mindiff, diff;

	int 		strtor (char *, float []);
        void 		PrRefInfo (char *, char *, char *, char *, char *);

/* -------------------------------- begin ---------------------------------- */

	crsplit_in = nimgs;

	exp_in = exptot / (double) crsplit_in;
	par->meanexp = exp_in;

	/* if all parameters are specified by the user, no need to open the 
		reference CRREJ table */
	if (newpar[0] < MAX_PAR) {
	    tp = c_tbtopn (par->tbname, IRAF_READ_ONLY, 0);
            if (c_iraferr() != 0) {
                trlerror("CRREJTAB table '%s' does not exist", par->tbname);
                return (2);
            }
	    nrows = c_tbpsta (tp, TBL_NROWS);

	    /* read the columns CRSPLIT and MEANEXP */
	    c_tbcfnd1 (tp, "crsplit", &colptr);
	    if (colptr == 0) {
	        trlerror("column CRSPLIT does not exist in CRREJTAB");
	        return (2);
	    }
	    c_tbcfnd1 (tp, "meanexp", &colptr1);
	    if (colptr1 == 0) {
	        trlerror("column MEANEXP does not exist in CRREJTAB");
	        return (2);
	    }

	    /* find the largest value in the CRSPLIT column */
	    for (i = 1; i <= nrows; i++) {
	        c_tbegti (tp, colptr, i, &crsplit);
	        if (i == 1) maxcrsplit = crsplit;
	        if (crsplit > maxcrsplit)
		    maxcrsplit = crsplit;
	    }
	    if (crsplit_in >= maxcrsplit) crsplit_in = maxcrsplit;

	    /* find the matching row in CRREJTAB to use */
	    nmatch = 0;
	    for (i = 1; i <= nrows; i++) {
	        c_tbegti (tp, colptr, i, &crsplit);
	        c_tbegtd (tp, colptr1, i, &meanexp);
	        diff = meanexp - exp_in;
	        if (crsplit_in == crsplit && diff >= 0.) {
		    nmatch++;
		    if (nmatch == 1 || diff <= mindiff) {
			row = i;
			mindiff = diff;
		    }
	        }
	    }
	    if (nmatch == 0) {
	        trlerror(" No matching CRSPLIT and MEANEXP in CRREJTAB");
	        return (2);
	    }

	    /* read the sigmas parameter */ 
	    if (newpar[CRSIGMAS] == 0) {
	        c_tbcfnd1 (tp, "crsigmas", &colptr);
	        if (colptr == 0) {
	            trlerror("column CRSIGMAS does not exist in CRREJTAB");
	            return (2);
	        }
	        c_tbegtt (tp, colptr, row, par->sigmas, STIS_LINE);
		if (c_iraferr()) {
		    trlerror("can't read CRSIGMAS from row %d", row);
		    return (2);
		}
	    }

	    /* read other parameters */
	    if (newpar[SKYSUB] == 0) {
	        c_tbcfnd1 (tp, "skysub", &colptr);
	        if (colptr == 0) {
	            trlerror("column SKYSUB does not exist in CRREJTAB");
	            return (2);
	        }
	        c_tbegtt (tp, colptr, row, par->sky, STIS_LINE);

		/* do not subtract sky if spectroscopic mode */
		if (strcmp (par->obstype, "SPECTROSCOPIC") == 0) 
		    strcpy (par->sky, "none");
	    }

	    /* CR propagation parameters */
	    if (newpar[CRRADIUS] == 0) {
	        c_tbcfnd1 (tp, "crradius", &colptr);
	        if (colptr == 0) {
	            trlerror("column CRRADIUS does not exist in CRREJTAB");
	            return (2);
	        }
	        c_tbegtr (tp, colptr, row, &par->rej);
	    }
	    if (newpar[CRTHRESH] == 0) {
	        c_tbcfnd1 (tp, "crthresh", &colptr);
	        if (colptr == 0) {
	            trlerror("column CRTHRESH does not exist in CRREJTAB");
	            return (2);
	        }
	        c_tbegtr (tp, colptr, row, &par->psigma);
	    }

            /* figure out how to do initial comparison image */
	    if (newpar[INITGUES] == 0) {
	        c_tbcfnd1 (tp, "initgues", &colptr);
	        if (colptr == 0) {
	            trlerror("column INITGUES does not exist in CRREJTAB");
	            return (2);
	        }
	        c_tbegtt (tp, colptr, row, par->initial, STIS_LINE);
	    }

	    /* read the noise model */
	    if (newpar[SCALENSE] == 0) {
	        c_tbcfnd1 (tp, "scalense", &colptr);
	        if (colptr == 0) {
	            trlerror("column SCALENSE does not exist in CRREJTAB");
	            return (2);
	        }
	        c_tbegtr (tp, colptr, row, &par->scalenoise);
	    }

	    if (newpar[BADINPDQ] == 0) {
	        c_tbcfnd1 (tp, "badinpdq", &colptr);
	        if (colptr == 0) {
	            trlerror("column BADINPDQ does not exist in CRREJTAB");
	            return (2);
	        }
	        c_tbegts (tp, colptr, row, &par->badbits);
	    }

	    if (newpar[CRMASK] == 0) {
	        c_tbcfnd1 (tp, "crmask", &colptr);
	        if (colptr == 0) {
	            trlerror("column CRMASK does not exist in CRREJTAB");
	            return (2);
	        }
	        c_tbegti (tp, colptr, row, &par->mask);
	    }

	    c_tbtclo (tp);
	}

        PrRefInfo ("crrejtab", par->tbname, "", "", "");

	/* parse the sigmas string into numbers */
        *niter = strtor (par->sigmas, sigma);
        if (*niter > MAX_ITER) {
            trlmessage("max number of iterations is exceeded");
	    return (2);
	}
        if (*niter <= 0) {
            trlmessage("number of iterations is zero");
	    return (2);
	}

	/* other fixed (for now) parameters */
	par->crval = (short) 8192;
	par->fillval = 0.;

	/* print out which parameter are used */
	if (par->verbose) {
	    trlmessage("\n number of imsets = %d", nimgs);
	    trlmessage(" ref table used: %s", par->tbname);
	    trlmessage(" initial guess method: %s", par->initial);
	    trlmessage(" total exposure time = %0.1f", exptot);
	    trlmessage(" sigmas used: %s", par->sigmas);
	    trlmessage(" sky subtraction used: %s", par->sky);
	    trlmessage(" rejection radius = %0.1f", par->rej);
	    trlmessage(" propagation threshold = %0.1f", par->psigma);
	    trlmessage(" scale noise = %0.1f%%", par->scalenoise);
	    trlmessage(" input bad bits value = %d", par->badbits);
	    trlmessage(" reset crmask = %d (1 = yes)", par->mask);
	}

	return (0);
}
