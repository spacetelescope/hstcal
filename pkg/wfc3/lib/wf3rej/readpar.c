# include   <stdio.h>
# include   <string.h>
# include   <stdlib.h>
# include   <ctype.h>
#include "hstcal.h"
# include   "hstio.h"
# include   "xtables.h"

# include   "wf3.h"
# include   "wf3rej.h"
# include   "hstcalerr.h"
# include   "wf3dq.h"
# include   "rej.h"

static int strtor (char *, float []);

/*  rejpar_in -- Read parameters either from user input or table.

  Description:
  ------------
  Reads CL parameters and does necessary checkings
  
  Input parameters from crrej reference table:
  -------------------------------------------
  Col. Name     Parameter
  "skysub"      sky         Sky levels subtraction scheme
  "crsigmas"    sigmas      Rejection thresholds
  "crradius"    radius      Radius (in pixels) to propagate the cosmic ray
  "crthresh"    thresh      Propagation factor
  "initgues"    initgues    Scheme of computing initial-guess image
  "scalense"    scalense    multiplicative noise in percents
  "badinpdq"    badinpdq    Data quality pset
  "crmask"      mask        flag CR-rejected pixels in input files?

  Input parameters from input image primary header:
  ------------------------------------------------


  Date          Author          Description
  ----          ------          -----------
  24-Sep-1998   W.J. Hack       Initial ACS Version	
  29-Aug-2000	H.A. Bushouse	Initial WFC3 Version
*/

int rejpar_in (clpar *par, int newpar[], int nimgs, float exptot, int *niter,
	       float sigma[]) {

    extern int status;

    IRAFPointer     tp;
    IRAFPointer     colptr, colptr1;
    int             i, nrows, nmatch, row;
    int             crsplit_in, crsplit, maxcrsplit;
    float           exp_in, meanexp, mindiff, diff;
    char            maskstr[SZ_CBUF+1];

    void    PrRefInfo (char *, char *, char *, char *, char *);
    void    WhichError (int);

/* -------------------------------- begin ---------------------------------- */

    crsplit_in = nimgs;
    maxcrsplit=0;
    row=0;
    mindiff=0.0f;
    diff=0.0f;
    
    exp_in = exptot / (float) crsplit_in;
    par->meanexp = exp_in;

    /* if all parameters are specified by the user, no need to open the 
	    reference CRREJ table */
    if (newpar[0] < MAX_PAR) {

        tp = c_tbtopn (par->tbname, IRAF_READ_ONLY, 0);
        if (c_iraferr() != 0) {
            trlerror("CRREJTAB table '%s' does not exist", par->tbname);
            return (status = TABLE_ERROR);
        }
        nrows = c_tbpsta (tp, TBL_NROWS);

        /* read the columns CRSPLIT and MEANEXP */
        c_tbcfnd1 (tp, "crsplit", &colptr);
        if (colptr == 0) {
            trlerror("column CRSPLIT does not exist in CRREJTAB");
            return (status = COLUMN_NOT_FOUND);
        }
        c_tbcfnd1 (tp, "meanexp", &colptr1);
        if (colptr1 == 0) {
            trlerror("column MEANEXP does not exist in CRREJTAB\n");
            return (status = COLUMN_NOT_FOUND);
        }
        nmatch = 0;

        /* find the largest value in the CRSPLIT column */
        for (i = 1; i <= nrows; i++) {
            c_tbegti (tp, colptr, i, &crsplit);
            if (i == 1) maxcrsplit = crsplit;
            if (crsplit > maxcrsplit)
            maxcrsplit = crsplit;
        }
        if (crsplit_in >= maxcrsplit) crsplit_in = maxcrsplit;

        /* find the matching row in CRREJTAB to use */
        for (i = 1; i <= nrows; i++) {
            c_tbegti (tp, colptr, i, &crsplit);
            c_tbegtr (tp, colptr1, i, &meanexp);
            diff = meanexp - exp_in;
            if (crsplit_in == crsplit && diff >= 0.) {
                nmatch++;
                if (nmatch == 1) mindiff = diff;
                if (diff <= mindiff) {
                    row = i;
                    mindiff = diff;
                }
            }
        }
        if (nmatch == 0) {
            trlerror(" No matching CRSPLIT and MEANEXP in CRREJTAB");
            return (status = ROW_NOT_FOUND);
        }

        /* read the sigmas parameter */ 
        if (newpar[CRSIGMAS] == 0) {
            c_tbcfnd1 (tp, "crsigmas", &colptr);
            if (colptr == 0) {
                trlerror("column CRSIGMAS does not exist in CRREJTAB");
                return (status = COLUMN_NOT_FOUND);
            }
            c_tbegtt (tp, colptr, row, par->sigmas, CHAR_LINE_LENGTH);
        }

        /* read other parameters */
        if (newpar[SKYSUB] == 0) {
            c_tbcfnd1 (tp, "skysub", &colptr);
            if (colptr == 0) {
                trlerror("column SKYSUB does not exist in CRREJTAB");
                return (status = COLUMN_NOT_FOUND);
            }
            c_tbegtt (tp, colptr, row, par->sky, SZ_FITS_REC);
        }

        /* CR propagation parameters */
        if (newpar[CRRADIUS] == 0) {
            c_tbcfnd1 (tp, "crradius", &colptr);
            if (colptr == 0) {
                trlerror("column CRRADIUS does not exist in CRREJTAB");
                return (status = COLUMN_NOT_FOUND);
            }
            c_tbegtr (tp, colptr, row, &par->radius);
        }
        if (newpar[CRTHRESH] == 0) {
            c_tbcfnd1 (tp, "crthresh", &colptr);
            if (colptr == 0) {
                trlerror("column CRTHRESH does not exist in CRREJTAB");
                return (status = COLUMN_NOT_FOUND);
            }
            c_tbegtr (tp, colptr, row, &par->thresh);
        }

            /* figure out how to do initial comparison image */
        if (newpar[INITGUES] == 0) {
            c_tbcfnd1 (tp, "initgues", &colptr);
            if (colptr == 0) {
                trlerror("column INITGUES does not exist in CRREJTAB");
                return (status = COLUMN_NOT_FOUND);
            }
            c_tbegtt (tp, colptr, row, par->initgues, SZ_FITS_REC);
        }

        /* read the noise model */
        if (newpar[SCALENSE] == 0) {
            c_tbcfnd1 (tp, "scalense", &colptr);
            if (colptr == 0) {
                trlerror("column SCALENSE does not exist in CRREJTAB");
                return (status = COLUMN_NOT_FOUND);
            }
            c_tbegtr (tp, colptr, row, &par->scalense);
        }

        if (newpar[BADINPDQ] == 0) {
            c_tbcfnd1 (tp, "badinpdq", &colptr);
            if (colptr == 0) {
                trlerror("column BADINPDQ does not exist in CRREJTAB");
                return (status = COLUMN_NOT_FOUND);
            }
            c_tbegts (tp, colptr, row, &par->badinpdq);
        }

        if (newpar[CRMASK] == 0) {
            c_tbcfnd1 (tp, "crmask", &colptr);
            if (colptr == 0) {
                trlerror("column CRMASK does not exist in CRREJTAB");
                return (status = COLUMN_NOT_FOUND);
            }
            c_tbegti (tp, colptr, row, &par->mask);
        }

        c_tbtclo (tp);
    }

    PrRefInfo ("crrejtab", par->tbname, "", "", "");

    /* parse the sigmas string into numbers */
    *niter = strtor (par->sigmas, sigma);
    if (status != WF3_OK) {
        WhichError (status);
        return (status);
    }
    if (*niter > MAX_ITER) {
        trlerror("No more than %d iterations permitted.", MAX_ITER);
        return (status = ERROR_RETURN);
    }
    if (*niter <= 0) {
        trlerror("Number of iterations is ZERO.");
        return (status = ERROR_RETURN);
    }

    /* other fixed (for now) parameters */
    par->crval = (short) DATAREJECT;
    par->fillval = 0.;

    /* print out which parameter are used */
    if (par->verbose) {
        trlmessage("\n number of images = %d", nimgs);
        trlmessage(" CRREJ ref table used: %s", par->tbname);
        trlmessage(" initial guess method: %s", par->initgues);
        trlmessage(" total exposure time = %0.1f", exptot);
        trlmessage(" sigmas used: %s", par->sigmas);
        trlmessage(" sky subtraction used: %s", par->sky);
        trlmessage(" rejection radius = %0.1f", par->radius);
        trlmessage(" propagation threshold = %0.1f", par->thresh);
        // Double-escape for trlmessage. Direct trlmessage calls don't require double-escaping.
        // This will be a direct trlmessage call in the future.
        trlmessage(" scale noise = %0.1f%%%%", par->scalense);
        trlmessage(" input bad bits value = %d", par->badinpdq);

        if (par->mask == 1) {
            strcpy (maskstr,"YES");
        } else {
            strcpy (maskstr, "NO");
        }
        trlmessage(" reset crmask = %s\n", maskstr);
    }
    return (status);
}

/*  strtor -- convert a string of real numbers into a real array 

    Description:
    ------------
    If the input string is blank(s), this routine will return 0.  If there are
    characters other than digits, decimal point, comma, semi-colon, slash or 
    blank in the input string, this routine will issue an error message.
    A null (or blank) string is valid and the function value will be zero
    (i.e. no value specified). However, a null substring (e.g. repeated
    separators) is an error.
    
    NOTE: This function sets 'status' upon error, but returns a 
        different variable.

    Date            Author          Description
    ----            ------          -----------
    09-May-1996     J.-C. Hsu       Adapt from the SPP code strtor.x
    26-Apr-2010     Howard Bushouse Add support for embedded blanks via
                                    updates from STIS library strtor.c.
*/

static int strtor (char *str, float arr[]) {

    extern int status;

    /* indexes in the string; the substring to be copied to tmp (from
       which a numerical value will be read) runs from ip0 to ip
    */
    int	    ip0, ip;
    int	    n, i;
    int	    done;
    char    tmp[100];

    n = 0;
    ip0 = 0;
    ip = 0;

    /* Initialize value to allow for proper error-checking
        after this function returns to the calling routine */
    status = WF3_OK;

    while (str[ip0] == ' ') {
           ip0++;
           ip++;
    }
    if (str[ip0] == '\0')
        return (0);

    done = 0;
    while (!done) {
        if (str[ip] == ',' || str[ip] == ';' || str[ip] == '/' || 
            str[ip] == ' ' || str[ip] == '\0') {
            for (i = 0; i < (ip-ip0); ++i)
                tmp[i] = str[ip0+i];
            
            tmp[ip-ip0] = '\0';
            if (!(isdigit (tmp[0]) || tmp[0] == '-' || tmp[0] == '.')) {
                trlerror("illegal input string '%s'", str);
                status = INVALID_VALUE;
                return (0);
            }
            arr[n] = (float) atof(tmp);
            ++n;
            if (str[ip] == '\0')
                break;
            ip++;		/* increment past the separator */
            ip0 = ip;

            while (str[ip0] == ' ') {
		   ip0++;
		   ip++;
	    }
        } else {
	    ++ip;
	}
        if (str[ip0] == '\0') 
	    break;
    }

    return (n);
}
