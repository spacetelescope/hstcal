# include   <stdio.h>
# include   <string.h>
#include "hstcal.h"
# include "hstio.h"

# include   "wf3.h"
# include   "wf3rej.h"
# include   "rej.h"
# include   "hstcalerr.h"
# include   "wf3info.h"


static int getnsegn (Hdr *, char *, multiamp *, multiamp *);
static int checkgn (multiamp, char *);
static int getampxy (Hdr *, int, int, char *, int, int, int *, int *);

/*  rej_check -- check input files of crrej

  Description:
  ------------
  Open and check input images/masks to have consistent dimensions and number
  of groups.
  
  Date      Author          Description
  ----      ------          -----------
  22-Sep-98 W. Hack         Initial version using multiamp noise,gain
  13-Sep-99 W.J. Hack       Simplified logic to only read in values from
                            first image where possible.
  20-Oct-99 W.J. Hack       getACSnsegn function revised to use 'for' loop
                            to read keywords.  
  29-Aug-00 H.A. Bushouse   Revised for WFC3 use.
  24-Jun-09 H. Bushouse     Fixed ampx/ampy values in getampxy for IR subarrays.
  30-Aug-12 M. Sosey        Checks the value of EXPFLAG in all the input image and if
                            any one image contains something other than NORMAL, it reports
                            the value as MIXED in the output crj header. PR #72001
*/

int rej_check (IRAFPointer tpin, int extver, int ngrps, clpar *par,
	       int newpar[], char imgname[][CHAR_FNAME_LENGTH+1], int grp[],
	       IODescPtr ipsci[], IODescPtr ipdq[], multiamp *noise,
	       multiamp *gain, int *dim_x, int *dim_y, int nimgs, char expflagFinal[]) {

    extern int status;

    IODescPtr   ip;
    Hdr         prihdr, scihdr;         /* header structures */
    char        fdata[CHAR_FNAME_LENGTH+1];
    char        det[SZ_CBUF+1];
    int         detector;
    multiamp    gn, ron;
    char        ccdamp[NAMPS+1], ccdamp0[NAMPS+1];
    char        mixed[]="MIXED";
    char        expflag[24];
    int         k, n;
    int         i;
    int         chip;
    int         ampx, ampy;

    int         GetKeyInt (Hdr *, char *, int, int, int *);
    int         GetKeyStr (Hdr *, char *, int, char *, char *, int);
    int         DetCCDChip (char *, int, int, int *);
    int         streq_ic (char *, char *);/* strings equal? (case insensitive)*/
    void        initmulti (multiamp *);

    /* -------------------------------- begin ------------------------------- */

    /* rewind the image template */
    c_imtrew (tpin);

    /* initialize local variables */
    ccdamp[0] = '\0';
    ccdamp0[0] = '\0';
    expflag[0]='\0';
    chip = 0;
    detector = UNKNOWN_DETECTOR;
    
    /* loop over all input files */
    for (k = 0; k < nimgs; ++k) {

        /* read the next input image name in the template list */
        c_imtgetim (tpin, fdata, CHAR_FNAME_LENGTH);
        
        /* open the primary header */
        ip = openInputImage (fdata, "", 0);
        if (hstio_err()) {
            sprintf (MsgText, "Cannot open data file '%s'", fdata);
            trlerror (MsgText);
            return (status = OPEN_FAILED);
        }
        strcpy (imgname[k], fdata);

        initHdr (&prihdr);

        /* Get primary header for keyword information */
        getHeader (ip, &prihdr);

    /*Read the EXPFLAG keyword from all images */
    if (GetKeyStr (&prihdr, "EXPFLAG", NO_DEFAULT, "",expflag, 24)){
        trlkwerr ("EXPFLAG",fdata);
        return (status=  KEYWORD_MISSING);
    }
    
    /* compare the expflag values */
    if (strcmp(expflag,expflagFinal) != 0){ /*they are not the same value*/
        if(strcmp(expflagFinal,mixed) !=0) {/*the final keyword is not already mixed*/
            strcpy(expflagFinal,mixed);
        }
    }
    

	/* Read the CCDAMP keyword from all images */
        if (GetKeyStr (&prihdr, "CCDAMP", NO_DEFAULT, "", ccdamp, NAMPS)) {
            trlkwerr ("CCDAMP", fdata);
            return(status = KEYWORD_MISSING);
        }

        /* The following keywords are read from first image only */        
        if (k == 0) {

	    n = extver;

            /* Save the CCDAMP value for the first image */
            strcpy (ccdamp0, ccdamp);

            /* read the CRREJ reference table name, if necessary */
            if (newpar[0] < MAX_PAR && par->tbname[0] == '\0') {
                if (GetKeyStr (&prihdr, "CRREJTAB", NO_DEFAULT, "", par->tbname,
			       CHAR_FNAME_LENGTH)) {
                    trlkwerr ("CRREJTAB", fdata);
                    return(status = KEYWORD_MISSING);
                }
            }

	    /* Detector keyword */
            if (GetKeyStr (&prihdr, "DETECTOR", NO_DEFAULT, "", det, SZ_CBUF)) {
                trlkwerr ("DETECTOR", fdata);
                return(status = KEYWORD_MISSING);
            }

	    /* Check detector keyword value */
            if (strcmp (det, "IR") == 0) {
                detector = IR_DETECTOR;
            } else if (strcmp (det, "UVIS") == 0) {
                detector = CCD_DETECTOR;
            } else {
                sprintf (MsgText, "DETECTOR = %s is unknown.", det);
                trlwarn (MsgText);
                detector = UNKNOWN_DETECTOR;
            }

            /* Read in READNOISE and ATODGN keywords */
            initmulti (&ron);
            initmulti (&gn);
            if (getnsegn (&prihdr, fdata, &ron, &gn) )
                return (status);
            if (checkgn (gn, fdata)) {
                return (status);
            }
        
        } else {
        
	    /* These checks are performed for all images after the first */

            /* Make sure each image has the same value of CCDAMP */
            if (!streq_ic (ccdamp,ccdamp0)) {
                sprintf (MsgText, "%s uses different CCDAMP", fdata);
                trlerror (MsgText);
                return (status = INVALID_VALUE);
            }

            /* Determine which extension corresponds to desired chip 
            ** for the remainder of the images */
            if (DetCCDChip(fdata, chip, ngrps, &n) ) {
                return (status);
            }
                
        }
        
        grp[k] = n;
        ipsci[k] = openInputImage (fdata, "SCI", n);
        ipdq[k]  = openInputImage (fdata, "DQ",  n);

        /* Use the first image's attributes to compare with the rest of 
        ** the files */
        if (k == 0) {
            *dim_x = getNaxis1(ipsci[k]);
            *dim_y = getNaxis2(ipsci[k]);

            /* Open sci hdr to get CCDHIP ID */
            initHdr (&scihdr);
            getHeader (ipsci[k], &scihdr);

            /* Which CHIP corresponds to the input extension?
            ** This chip ID will be used for the remainder of
            ** the images in the list 
            */			
            if (GetKeyInt (&scihdr, "CCDCHIP", USE_DEFAULT, 1, &chip))
                chip = 1;
            
            /* Now we need to read in CCDTAB file to get AMP regions */
            if (getampxy (&prihdr, detector, chip, ccdamp, *dim_x, *dim_y,
			  &ampx, &ampy) )
                return (status);

            freeHdr (&scihdr);
        }

        /* Verify the image size to be the same as the first image */
        if (getNaxis1(ipsci[k]) != *dim_x || getNaxis2(ipsci[k]) != *dim_y){
            sprintf (MsgText,
	   "file '%s[sci,%d]' does not have the same size as the first image\n",
	    imgname[k], grp[k]);
            trlerror (MsgText);
            return (status = SIZE_MISMATCH);
        }

        /* Keep memory clean... */
        closeImage (ip);
        freeHdr (&prihdr);

    } /* End loop over k (images in input list) */	

    /* Save noise and gain values for use in rest of WF3REJ */
    for (i = 0; i < NAMPS; i++) {
        if (gn.val[i] == 0. ) {
	    noise->val[i] = 0.;
        } else {
	    noise->val[i] = ron.val[i]/gn.val[i];
        }
        gain->val[i]  = gn.val[i];
    }
    noise->colx = ampx;
    noise->coly = ampy;
    gain->colx  = ampx;
    gain->coly  = ampy;

    /* Store the chip ID that is being processed */
    noise->chip = chip;
    gain->chip  = chip;
    noise->detector = detector;
    gain->detector  = detector;						

    return (status);
}


static int getnsegn (Hdr *hdr, char *fdata, multiamp *ron, multiamp *gn) {

    extern int status;
    
    char   *nsekeyw[NAMPS] = { "READNSEA", "READNSEB", "READNSEC", "READNSED"};
    char   *gnkeyw[NAMPS] = { "ATODGNA", "ATODGNB", "ATODGNC", "ATODGND"};
    int    i;
    
    int    GetKeyFlt (Hdr *, char *, int, float, float *);

    /* Get the readout noise values from the SCI header */
    for (i = 0; i < NAMPS; i++){    
        if (GetKeyFlt (hdr, nsekeyw[i], USE_DEFAULT, 0., &ron->val[i]) != 0) {
            trlkwerr (nsekeyw[i], fdata);
            return (status = KEYWORD_MISSING);
        }
    }

    /* Get the A-to-D gain values from the SCI header*/
    for (i = 0; i < NAMPS; i++){    
        if (GetKeyFlt (hdr, gnkeyw[i], USE_DEFAULT, 0., &gn->val[i]) != 0) {
            trlkwerr (gnkeyw[i], fdata);
            return (status = KEYWORD_MISSING);
        }
    }

    return (status);
}

static int checkgn (multiamp gn, char *fdata) {

    extern int status;
    
    int i;
    int ampset = 0;
    
    for (i = 0; i < NAMPS; i++) {
        if (gn.val[i] > 0.) ampset++;
    }   
    
    if (ampset == 0) {
        sprintf (MsgText, "All ATODGN values in file %s are 0.\n", fdata);
        trlerror (MsgText);
        return (status = INVALID_VALUE);
    }

    return (status);

}

static int getampxy (Hdr *hdr, int det, int chip, char *ccdamp, int dimx,
		     int dimy, int *ampx, int *ampy) {

    extern int status;
    WF3Info wf3rej;
    char tabname[CHAR_LINE_LENGTH+1];

    void WF3Init (WF3Info *);
    int GetKeys (WF3Info *, Hdr *);
    int GetKeyStr (Hdr *, char *, int, char *, char *, int);
    int GetCCDTab (WF3Info *, int, int);

    WF3Init (&wf3rej);
    
    tabname[0] = '\0';

    /* Copy input values into wf3 structure */
    wf3rej.detector = det;
    wf3rej.chip = chip;
    strcpy(wf3rej.ccdamp, ccdamp);

    /* We need to read in CCDGAIN, CCDOFST[A-D], and BINAXIS[1,2] 
    ** from hdr as selection criteria for row from CCDTAB. 
    ** Get keyword values from primary header using same function 
    ** used in WF3CCD.
    */
    if (GetKeys (&wf3rej, hdr)) {
        freeHdr (hdr);
        return (status);
    }

    /* Read CCDTAB to get AMPX, AMPY values */
    if (GetKeyStr (hdr, "CCDTAB", USE_DEFAULT, "", tabname, SZ_FITS_REC))
        return (status);
    strcpy (wf3rej.ccdpar.name, tabname);

    if (GetCCDTab (&wf3rej, dimx, dimy) ) {
        return (status);
    }

    /* Assign AMPX and AMPY values read from table */
    *ampx = wf3rej.ampx;
    *ampy = wf3rej.ampy;

    /* For IR images the AMPX/AMPY values in CCDTAB are only correct for
    ** full-frame exposures, therefore reset them in case we're processing
    ** IR subarray images. They're always equal to half the image size. */
    if (wf3rej.detector == IR_DETECTOR) {
	*ampx = dimx / 2;
	*ampy = dimy / 2;
    }

    return (status);
}
