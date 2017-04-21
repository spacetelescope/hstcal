# include   <stdio.h>
# include   "hstio.h"

# include   "acs.h"
# include   "err.h"
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
    14-Mar-2002     W.J. Hack   Added computation of cumulative DARKTIME
    4-Apr-2002      W.J. Hack   added initialization of 'totd'
   24-Apr-2002      W.J. Hack   removed darktime altogether, find initial EXPSTART
*/
int cr_scaling (char *expname, IRAFPointer tpin, float efac[], int *nimgs, double *expend, double *expstart)
{
    extern int status;

    Hdr         prihdr;
    int         nzero, k;
    char        fdata[ACS_FNAME + 1];
    IODescPtr   ip;
    int         numimgs;        /* How many good input images are there? */

    double     end, keyend, keystart, start;

    int         GetKeyFlt (Hdr *, char *, int, float, float *);
    int         GetKeyDbl (Hdr *, char *, int, double, double *);
    /* -------------------------------- begin ---------------------------------- */

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
        c_imtgetim (tpin, fdata, ACS_FNAME);

        /* open the primary header */
        ip = openInputImage (fdata, "", 0);
        if (hstio_err()) {
            sprintf (MsgText, "Cannot open data file '%s'", fdata);
            trlerror (MsgText);
            return (status = OPEN_FAILED);
        }

        initHdr (&prihdr);

        /* read in primary header from image */
        getHeader (ip, &prihdr);

        if (GetKeyFlt (&prihdr, expname, USE_DEFAULT, 0., &efac[k]) != 0) {
            sprintf (MsgText, "cannot read '%s' from the primary header of '%s'", expname, fdata);
            trlerror (MsgText);
            freeHdr (&prihdr);
            return(status = KEYWORD_MISSING);
        }
        
        if (efac[k] < 0.) {
            sprintf (MsgText, "exposure time of file '%s' is negative", fdata);
            trlerror (MsgText);
            freeHdr (&prihdr);
            return(status = INVALID_VALUE);
        }
        if (efac[k] == 0.) {
            nzero++;
        }
        
        numimgs++;
        if (GetKeyDbl (&prihdr, "EXPEND", USE_DEFAULT, 0., &keyend) != 0) {
            sprintf (MsgText, "cannot read 'EXPEND' from the primary header of '%s'", fdata);
            trlerror (MsgText);
            freeHdr (&prihdr);
            return(status = KEYWORD_MISSING);
        }
        if (GetKeyDbl (&prihdr, "EXPSTART", USE_DEFAULT, 0., &keystart) != 0) {
            sprintf (MsgText, "cannot read 'EXPSTART' from the primary header of '%s'", fdata);
            trlerror (MsgText);
            freeHdr (&prihdr);
            return(status = KEYWORD_MISSING);
        }
        
        end = (keyend > end) ? keyend: end;
        start = (keystart < start) ? keystart : start;
        closeImage (ip);
        freeHdr (&prihdr);
    }
    
    if (nzero > 0 && nzero < *nimgs) {
        trlwarn ("Some (but not all) input imsets have zero exposure time.");
        trlwarn ("Final product will be compromised!");
        
        /* This type of error will need to be handled differently in order
            to allow pipeline processing of this type of dataset. 
        return (status = INVALID_VALUE);
        */
    }
    
    /* Only return the number of valid input images,
        initial EXPSTART and final EXPEND value
    */
    *nimgs = numimgs;
    *expend = end;
    *expstart = start;
    
    return (status);
}
