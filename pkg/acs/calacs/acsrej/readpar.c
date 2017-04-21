# include   <stdio.h>
# include   <string.h>
# include   <stdlib.h>
# include <ctype.h>
# include   "hstio.h"
# include   "xtables.h"

# include   "acs.h"
# include   "acsrej.h"
# include   "err.h"
# include   "acsdq.h"
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
*/
int rejpar_in (clpar *par, int newpar[], int nimgs, float exptot, int *niter, float sigma[])
{
    extern int status;

    IRAFPointer     tp;
    IRAFPointer     colptr, colptr1;
    int             i, nrows, nmatch, row;
    int             crsplit_in, crsplit, maxcrsplit;
    float           exp_in, meanexp, mindiff, diff;
    char            maskstr[ACS_CBUF+1];

    void    PrRefInfo (char *, char *, char *, char *, char *);
    void    WhichError (int);

/* -------------------------------- begin ---------------------------------- */

    crsplit_in = nimgs;

    exp_in = exptot / (float) crsplit_in;
    par->meanexp = exp_in;

    /* if all parameters are specified by the user, no need to open the 
	    reference CRREJ table */
    if (newpar[0] < MAX_PAR) {

        tp = c_tbtopn (par->tbname, IRAF_READ_ONLY, 0);
            if (c_iraferr() != 0) {
                sprintf (MsgText,"CRREJTAB table '%s' does not exist", par->tbname);
                trlerror (MsgText);
                return (status = TABLE_ERROR);
            }
        nrows = c_tbpsta (tp, TBL_NROWS);

        /* read the columns CRSPLIT and MEANEXP */
        c_tbcfnd1 (tp, "crsplit", &colptr);
        if (colptr == 0) {
            trlerror ("column CRSPLIT does not exist in CRREJTAB");
            return (status = COLUMN_NOT_FOUND);
        }
        c_tbcfnd1 (tp, "meanexp", &colptr1);
        if (colptr1 == 0) {
            trlerror ("column MEANEXP does not exist in CRREJTAB\n");
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
            trlerror (" No matching CRSPLIT and MEANEXP in CRREJTAB");
            return (status = ROW_NOT_FOUND);
        }

        /* read the sigmas parameter */ 
        if (newpar[CRSIGMAS] == 0) {
            c_tbcfnd1 (tp, "crsigmas", &colptr);
            if (colptr == 0) {
                trlerror ("column CRSIGMAS does not exist in CRREJTAB");
                return (status = COLUMN_NOT_FOUND);
            }
            c_tbegtt (tp, colptr, row, par->sigmas, ACS_LINE);
        }

        /* read other parameters */
        if (newpar[SKYSUB] == 0) {
            c_tbcfnd1 (tp, "skysub", &colptr);
            if (colptr == 0) {
                trlerror ("column SKYSUB does not exist in CRREJTAB");
                return (status = COLUMN_NOT_FOUND);
            }
            c_tbegtt (tp, colptr, row, par->sky, ACS_FITS_REC);
        }

        /* CR propagation parameters */
        if (newpar[CRRADIUS] == 0) {
            c_tbcfnd1 (tp, "crradius", &colptr);
            if (colptr == 0) {
                trlerror ("column CRRADIUS does not exist in CRREJTAB");
                return (status = COLUMN_NOT_FOUND);
            }
            c_tbegtr (tp, colptr, row, &par->radius);
        }
        if (newpar[CRTHRESH] == 0) {
            c_tbcfnd1 (tp, "crthresh", &colptr);
            if (colptr == 0) {
                trlerror ("column CRTHRESH does not exist in CRREJTAB");
                return (status = COLUMN_NOT_FOUND);
            }
            c_tbegtr (tp, colptr, row, &par->thresh);
        }

            /* figure out how to do initial comparison image */
        if (newpar[INITGUES] == 0) {
            c_tbcfnd1 (tp, "initgues", &colptr);
            if (colptr == 0) {
                trlerror ("column INITGUES does not exist in CRREJTAB");
                return (status = COLUMN_NOT_FOUND);
            }
            c_tbegtt (tp, colptr, row, par->initgues, ACS_FITS_REC);
        }

        /* read the noise model */
        if (newpar[SCALENSE] == 0) {
            c_tbcfnd1 (tp, "scalense", &colptr);
            if (colptr == 0) {
                trlerror ("column SCALENSE does not exist in CRREJTAB");
                return (status = COLUMN_NOT_FOUND);
            }
            c_tbegtr (tp, colptr, row, &par->scalense);
        }

        if (newpar[BADINPDQ] == 0) {
            c_tbcfnd1 (tp, "badinpdq", &colptr);
            if (colptr == 0) {
                trlerror ("column BADINPDQ does not exist in CRREJTAB");
                return (status = COLUMN_NOT_FOUND);
            }
            c_tbegts (tp, colptr, row, &par->badinpdq);
        }

        if (newpar[CRMASK] == 0) {
            c_tbcfnd1 (tp, "crmask", &colptr);
            if (colptr == 0) {
                trlerror ("column CRMASK does not exist in CRREJTAB");
                return (status = COLUMN_NOT_FOUND);
            }
            c_tbegti (tp, colptr, row, &par->mask);
        }

        c_tbtclo (tp);
    }

    PrRefInfo ("crrejtab", par->tbname, "", "", "");
    
    /* parse the sigmas string into numbers */
    *niter = strtor (par->sigmas, sigma);
    if (status != ACS_OK) {
        WhichError (status);
        return (status);
    }
    
    if (*niter > MAX_ITER) {
        sprintf (MsgText,"No more than %d iterations permitted.", MAX_ITER);
        trlerror (MsgText);
        return (status = ERROR_RETURN);
    }
    if (*niter <= 0) {
        trlerror ("Number of iterations is ZERO.");
        return (status = ERROR_RETURN);
    }

    /* other fixed (for now) parameters */
    par->crval = (short) DATAREJECT;
    par->fillval = 0.;
    
    /* print out which parameter are used */
    if (par->verbose) {
        sprintf (MsgText,"\n number of images = %d", nimgs);
        trlmessage (MsgText);
        sprintf (MsgText," CRREJ ref table used: %s", par->tbname);
        trlmessage (MsgText);
        sprintf (MsgText," initial guess method: %s", par->initgues);
        trlmessage (MsgText);
        sprintf (MsgText," total exposure time = %0.1f", exptot);
        trlmessage (MsgText);
        sprintf (MsgText," sigmas used: %s", par->sigmas);
        trlmessage (MsgText);
        sprintf (MsgText," sky subtraction used: %s", par->sky);
        trlmessage (MsgText);
        sprintf (MsgText," rejection radius = %0.1f", par->radius);
        trlmessage (MsgText);
        sprintf (MsgText," propagation threshold = %0.1f", par->thresh);
        trlmessage (MsgText);
        sprintf (MsgText," scale noise = %0.1f%%", par->scalense);
        trlmessage (MsgText);
        sprintf (MsgText," input bad bits value = %d", par->badinpdq);
        trlmessage (MsgText);
        
        if (par->mask == 1) {
            strcpy (maskstr,"YES");
        } else {
            strcpy (maskstr, "NO");
        }
        sprintf (MsgText," reset crmask = %s\n", maskstr);
        trlmessage (MsgText);
    }
    return (status);
}

/*  strtor -- convert a string of real numbers into a real array 

  Description:
  ------------
  If the input string is blank(s), this routine will return 0.  If there are
  characters other than digits, decimal point, comma, semi-colon, slash or
  blank in the input string, this routine will print an error message and
  return -1.
  A null (or blank) string is valid, and the function value will be zero
  (i.e. no value specified).  However, a null substring (e.g. repeated
  separators) is an error.

  Date		Author		Description
  ----		------		-----------
  09-May-1996  J.-C. Hsu	Adapt from the SPP code strtor.x
  10-Feb-2000  Phil Hodge	Replace exit (2) with return (0).
  20-Mar-2002  Ivo Busko	Move to library, return -1 on error;
				support multiple separators.
  24-Jan-2005  Phil Hodge	Fix indexing error.
*/
int strtor (char *str, float arr[])
{
	/* indexes in the string; the substring to be copied to tmp (from
	   which a numerical value will be read) runs from ip0 to ip
	*/
	int	ip0, ip;
	int	n, i;
	int	done;
	char	tmp[100];

	n = 0;
	ip0 = 0;
	ip = 0;

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
		    printf ("illegal input string '%s'\n", str);
		    return (-1);
		}
	    	arr[n] = (float) atof(tmp);
		++n;
	        if (str[ip] == '\0')
		    break;
		ip++;			/* increment past the separator */
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
