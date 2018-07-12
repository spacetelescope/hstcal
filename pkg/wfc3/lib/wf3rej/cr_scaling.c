# include   <stdio.h>
# include   <string.h>
#include "hstcal.h"
# include "hstio.h"

# include   "wf3.h"
# include   "hstcalerr.h"
# include   "wf3info.h"
# include   "rej.h"

/*  cr_scaling -- Determine the scaling factors according to exposure times or
        other user specified scheme. */

/*
Description:
------------
If using the exposure time, the scaling factors are normalized to ratios 
relative to the max exposure. 

    Date            Author      Description
    ----            ------      -----------
    24-Sep-1998     W.J. Hack   Initial Version
    18-Mar-1999     W.J. Hack   Revised to read EXPTIMEs from Primary headers
                                using image-template list directly
    20-Oct-1999     W.J. Hack   Revised to compute number of good input images
                                and insure they are less than MAX_FILES.
    14-Apr-2000     W.J. Hack   Revised to also return final EXPEND appropriate
                                for output CR-combined product
    18-Jun-2002     H. Bushouse Modified to fine initial EXPSTART and return it
				(following CALACS changes).
    27-Oct-2003     H. Bushouse Upgrades for handling inputs with EXPTIME=0
				(following CALACS changes).
    14-Dec-2011     H. Bushouse Upgraded to read BUNIT keyword value from each
				input image. (PR 69969; Trac #814)
    07-May-2012 M. Sosey, updated the BUNIT read to be case insensitive
                (trac #887)
*/

int cr_scaling (char *expname, IRAFPointer tpin, float efac[], int *nimgs,
		double *expend, double *expstart, DataUnits bunit[]) {

    extern int status;

    Hdr         hdr;
    int         nzero, k;
    char        fdata[CHAR_FNAME_LENGTH + 1], units[12];
    IODescPtr   ip;
    int         numimgs;        /* How many good input images are there? */

    double      end, keyend, keystart, start;

    int         GetKeyFlt (Hdr *, char *, int, float, float *);
    int         GetKeyDbl (Hdr *, char *, int, double, double *);
    int         streq_ic (char *, char *);/* strings equal? (case insensitive)*/

    /* -------------------------------- begin ------------------------------- */

    /* Rewind the image template pointer */
    c_imtrew(tpin);

    *nimgs = c_imtlen(tpin);
    end = 0.0;
    keyend = 0.0;
    start = 1e+10;
    keystart = 0.0;
 
    /* Check to make sure there are not too many images to work with... */
    if (*nimgs > MAX_FILES) {
        trlerror("There are too many input images to combine. "); 
        return(status = NOTHING_TO_DO);
    }

    /* if the parameter scaling is null, all images have equal weight. 
        If no keyword name is given for the exposure time, assume equal
        weights of 1 for all images.
    */
    if (expname[0] == '\0') {
        return (status);
    }

    /* Use exposure time as scaling factor */
    nzero = 0;	 
 
    /* loop all input files counting how many usable inputs there are */
    numimgs = 0;
    for (k = 0; k < *nimgs; ++k) {

        /* read the next input image name in the template list */
        c_imtgetim (tpin, fdata, CHAR_FNAME_LENGTH);

        /* open the primary header */
        ip = openInputImage (fdata, "", 0);
        if (hstio_err()) {
            sprintf (MsgText, "Cannot open data file '%s'", fdata);
            trlerror (MsgText);
            return (status = OPEN_FAILED);
        }

        initHdr (&hdr);

        /* read in primary header from image */
        getHeader (ip, &hdr);

        if (GetKeyFlt (&hdr, expname, USE_DEFAULT, 0., &efac[k]) != 0) {
            sprintf (MsgText,
		     "cannot read '%s' from the primary header of '%s'",
		      expname, fdata);
            trlerror (MsgText);
            freeHdr (&hdr);
            return(status = KEYWORD_MISSING);
        }
        
        if (efac[k] < 0.) {
            sprintf (MsgText, "exposure time of file '%s' is negative", fdata);
            trlerror (MsgText);
            freeHdr (&hdr);
            return(status = INVALID_VALUE);
        }
        if (efac[k] == 0.) {
            nzero++;
        }
        
        numimgs++;
        if (GetKeyDbl (&hdr, "EXPEND", USE_DEFAULT, 0., &keyend) != 0) {
            sprintf (MsgText,
		     "cannot read 'EXPEND' from the primary header of '%s'",
		     fdata);
            trlerror (MsgText);
            freeHdr (&hdr);
            return(status = KEYWORD_MISSING);
        }
        if (GetKeyDbl (&hdr, "EXPSTART", USE_DEFAULT, 0., &keystart) != 0) {
            sprintf (MsgText,
		     "cannot read 'EXPSTART' from the primary header of '%s'",
		     fdata);
            trlerror (MsgText);
            freeHdr (&hdr);
            return(status = KEYWORD_MISSING);
        }
        
        end = (keyend > end) ? keyend: end;
	start = (keystart < start) ? keystart : start;
        closeImage (ip);
        freeHdr (&hdr);

	/* Open first SCI extension header to read BUNIT keyword */
        ip = openInputImage (fdata, "SCI", 1);
        if (hstio_err()) {
            sprintf (MsgText, "Cannot open data file '%s'", fdata);
            trlerror (MsgText);
            return (status = OPEN_FAILED);
        }

        initHdr (&hdr);

        /* read in the extension header from image */
        getHeader (ip, &hdr);

        /* Read the BUNIT keyword */
        units[0] = '\0';
        if (getKeyS (&hdr, "BUNIT", units)) {
            trlkwerr ("BUNIT", fdata);
            return (status = 1);
        }

        /* Check the BUNIT keyword value for validity */
        if (streq_ic (units,"COUNTS/S") || streq_ic (units, "ELECTRONS/S"))
            bunit[k] = COUNTRATE;
        else if  (streq_ic (units,"COUNTS") || streq_ic (units, "ELECTRONS"))
            bunit[k] = COUNTS;
        else {
            sprintf (MsgText, "%s value for BUNIT does not match", units);
            trlerror (MsgText);
            return (status = INVALID_VALUE);
        }


        closeImage (ip);
        freeHdr (&hdr);
    }
    
    if (nzero > 0 && nzero < *nimgs) {
        trlwarn ("Some (but not all) input imsets have zero exposure time.");
	trlwarn ("Final product will be compromised!");

	/* This type of error will need to be handled differently in order
	** to allow pipeline processing of this type of dataset.
        return (status = INVALID_VALUE);
	*/
    }
    
    /* Only return the number of valid input images,
    ** initial EXPSTART and final EXPEND value */
    *nimgs = numimgs;
    *expend = end;
    *expstart = start;
    
    return (status);
}
