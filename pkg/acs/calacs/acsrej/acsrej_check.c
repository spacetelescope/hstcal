# include   <stdio.h>
# include   <string.h>
#include "hstcal.h"
# include   "hstio.h"

# include   "acs.h"
# include   "acsrej.h"
# include   "rej.h"
# include   "hstcalerr.h"
# include   "acsinfo.h"
# include   "getacskeys.h"


static int getACSnsegn (Hdr *, char *, multiamp *, multiamp *);
static int checkgn (multiamp, char *);
static int getACSampxy (Hdr *, int, int, char *, int, int, int *, int *);

/*  acsrej_check -- check input files of acsrej

  Description:
  ------------
  Open and check input images/masks to have consistent dimensions and number
  of groups.

  Date        Author     Description
  ----        ------     -----------
  22-Sep-1998 W.J. Hack  Initial version using multiamp noise, gain.
  13-Sep-1999 W.J. Hack  Simplified logic to only read in values from
                         first image where possible.
  20-Oct-1999 W.J. Hack  getACSnsegn function revised to use 'for' loop
                         to read keywords.
  19-Oct-2015 P.L. Lim   Calculations all done in electrons now.
  30-Apr-2021 M.D. DeLaPena Compute the cumulative flash duration over all images.
  24-May-2021 M.D. DeLaPena Compute the cumulative value of the DARKTIME keywords and
                         reset the CRJ/CRC DARKTIME to this cumulativate valuve.
*/

int acsrej_check (IRAFPointer tpin, int extver, int ngrps, clpar *par,
                  int newpar[],
                  char imgname[][CHAR_FNAME_LENGTH], int grp[],
                  IODescPtr ipsci[], IODescPtr iperr[], IODescPtr ipdq[],
                  multiamp *noise, multiamp *gain, int *dim_x, int *dim_y,
                  int nimgs, float efac[MAX_FILES], float *cumFlashDur, float *cumDarktime) {
    /*
      Parameters:

      tpin    i: Pointer of image list.
      extver  i: EXTVER to process. For WFC, 1 is CD and 2 is AB.
      ngrps   i: Total number of EXTVER.
      par, newpar  i: User specified parameters.
      imagename  o: Array of image names.
      grp     o: Array of EXTVER for each input image.
      ipsci   o: Array of pointers to SCI extension of the given EXTVER,
                 each pointer is an input image. Unit now in electrons.
      iperr   o: Array of pointers to ERR extension of the given EXTVER,
                 each pointer is an input image. Unit now in electrons.
      ipdq    o: Array of pointers to DQ extension of the given EXTVER,
                 each pointer is an input image.
      noise   o: Calibrated readnoise. Noise value in electrons
                 for each amp is read from primary header keyword READNSE[AMP],
                 where [AMP] can be A, B, C, or D.
      gain    o: Calibrated gain in e/DN from primary header keyword
                 ATODGN[AMP], where [AMP] can be A, B, C, or D.
      dim_x, dim_y  o: Image dimension taken from the first input image.
                       All images must have the same dimension.
      nimgs   i: Number of input images.
      efac    i: Exposure times (seconds) for input images.
      cumFlashDur o: Cumulative flash duration over all images.
      cumDarktime o: Cumulative darktime over all images.
    */

    extern int status;

    IODescPtr   ip;
    Hdr         prihdr, scihdr;         /* header structures */
    char        fdata[CHAR_FNAME_LENGTH];
    char        det[ACS_CBUF], flashcur[ACS_CBUF], flashcur0[ACS_CBUF];
    int         detector, flshcorr, flshcorr0;
    float       flashdur, flashlevel, flashlevel0;
    float       darktime;
    multiamp    gn, ron;
    char        ccdamp[NAMPS+1], ccdamp0[NAMPS+1];
    int         k, n;
    int         i;
    int         chip;
    int         ampx, ampy;
    float       totalFlashDuration, totalDarktime;

    int         GetKeyInt (Hdr *, char *, int, int, int *);
    int         GetKeyFlt (Hdr *, char *, int, float, float *);
    int         GetKeyStr (Hdr *, char *, int, char *, char *, int);
    int         GetSwitch (Hdr *, char *, int *);
    int         DetCCDChip (char *, int, int, int *);
    int         streq_ic (char *, char *);  /* str equal? (case insensitive) */
    void        initmulti (multiamp *);
    /* -------------------------------- begin ------------------------------- */

    /* rewind the image template */
    c_imtrew (tpin);

    /* initialize local variables */
    ccdamp[0] = '\0';
    ccdamp0[0] = '\0';
    det[0] = '\0';
    flashcur[0] = '\0';
    flashcur0[0] = '\0';
    chip = 0;
    flashdur = 0.0;
    flashlevel = 0.0;
    flashlevel0 = 0.0;
    totalFlashDuration = 0.0;
    darktime = 0.0;
    totalDarktime = 0.0;

    /* loop through all input files */
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

        /* Read in keywords here that are required from all images */

        /* make sure the same CCDAMP is used for all observations */
        if (GetKeyStr (&prihdr, "CCDAMP", NO_DEFAULT, "", ccdamp, NAMPS)) {
            trlkwerr ("CCDAMP", fdata);
            return(status = KEYWORD_MISSING);
        }

        if (GetKeyFlt (&prihdr, "DARKTIME", USE_DEFAULT, 0.0, &darktime)) {
            trlkwerr ("DARKTIME", fdata);
            return(status = KEYWORD_MISSING);
        }
        totalDarktime += darktime;

        /* Post-flash keywords */
        if (GetSwitch (&prihdr, "FLSHCORR", &flshcorr)) {
            trlkwerr ("FLSHCORR", fdata);
            return (status = KEYWORD_MISSING);
        }
        if (flshcorr == PERFORM) {
            if (GetKeyStr (&prihdr, "FLASHCUR", NO_DEFAULT, "",
                           flashcur, ACS_CBUF)) {
                trlkwerr ("FLASHCUR", fdata);
                return(status = KEYWORD_MISSING);
            }
            if (GetKeyFlt (&prihdr, "FLASHDUR", USE_DEFAULT, 0.0, &flashdur)) {
                trlkwerr ("FLASHDUR", fdata);
                return(status = KEYWORD_MISSING);
            }
            flashlevel = flashdur / efac[k];
            totalFlashDuration += flashdur;
        } else {
            flshcorr = OMIT; /* OMIT/COMPLETE = flsh is not an issue anymore */
        }

        /* Read in keywords from first image */
        if (k == 0) {
            /* read the CRREJ reference table name from the first file,
               if necessary */
            if (newpar[0] < MAX_PAR && par->tbname[0] == '\0') {
                if(GetKeyStr (&prihdr, "CRREJTAB", NO_DEFAULT, "", par->tbname, CHAR_FNAME_LENGTH) ) {
                    trlkwerr ("CRREJTAB", fdata);
                    return(status = KEYWORD_MISSING);
                }
            }

            /* Detector */
            if (GetKeyStr (&prihdr, "DETECTOR", NO_DEFAULT, "",
                           det, ACS_CBUF)) {
                trlkwerr ("DETECTOR", fdata);
                return(status = KEYWORD_MISSING);
            }

            if (strcmp (det, "HRC") == 0) {
                detector = HRC_CCD_DETECTOR;
            } else if (strcmp (det, "WFC") == 0) {
                detector = WFC_CCD_DETECTOR;
            } else {
                sprintf (MsgText, "DETECTOR = %s is unknown.", det);
                trlwarn (MsgText);
                detector = UNKNOWN_DETECTOR;
            }

            /* Remember the values for the first image */
            strcpy(ccdamp0, ccdamp);
            flshcorr0 = flshcorr;
            if (flshcorr == PERFORM) {
                strcpy(flashcur0, flashcur);
                flashlevel0 = flashlevel;
            }

            /* Start by opening the inputted extension of the first image */
            n = extver;

            /* Read in READNOISE and ATODGN values from first image */
            initmulti (&ron);
            initmulti (&gn);
            if (getACSnsegn (&prihdr, fdata, &ron, &gn))
                return (status);
            if (checkgn (gn, fdata)) {
                return (status);
            }

        /* END CHECK ON FIRST (k=0) IMAGE HEADER INFO */
        } else {
            /* Make sure each image has the same value of CCDAMP */
            if (!streq_ic(ccdamp, ccdamp0)) {
                sprintf (MsgText,
                         "%s uses different CCDAMP, expected %s but got %s",
                         fdata, ccdamp0, ccdamp);
                trlerror (MsgText);
                return (status = INVALID_VALUE);
            }

            /* Make sure each image has the same post-flash properties.

               Notes from Norman Grogin:

               The global mode-subtraction will not solve the varying FLASH
               issue, as the flash illumination is far from flat across the
               field of view.

               Throw an error if it detects that: (the image stack does not all
               have the same value of FLASH --- zero or otherwise) &&
               (FLASHCORR != COMPLETE).

               Since the comparison is done in e/s, the comparison of FLASHDUR
               should instead be FLASHDUR/EXPTIME for each exposure.
               For example, two frames would be OK to CRREJ together if
               FLASHCORR != COMPLETE and (one had twice the exposure *and*
               twice the FLASHDUR of the other). But would not be OK if
               FLASHDUR were the same but the EXPTIME were different.

               If FLASHCORR == COMPLETE when the ACSREJ is performed, then
               no worry about the FLASHDUR. The needed information will live on
               in the ERR arrays. */
            if (flshcorr != flshcorr0) {
                sprintf (MsgText,
                         "%s uses different FLSHCORR, expected %d but got %d",
                         fdata, flshcorr0, flshcorr);
                trlerror (MsgText);
                return (status = INVALID_VALUE);
            }
            if (flshcorr == PERFORM) {
                if (!streq_ic(flashcur, flashcur0)) {
                    sprintf (MsgText,
                         "%s uses different FLASHCUR, expected %s but got %s",
                             fdata, flashcur0, flashcur);
                    trlerror (MsgText);
                    return (status = INVALID_VALUE);
                }
                if (flashlevel != flashlevel0) {
                    sprintf (MsgText,
   "%s uses different flash levels (FLASHDUR/EXPTIME), expected %f but got %f",
                             fdata, flashlevel0, flashlevel);
                    trlerror (MsgText);
                    return (status = INVALID_VALUE);
                }
            }

            /* Determine which extension corresponds to desired chip
               for the remainder of the images */
            if (DetCCDChip(fdata, chip, ngrps, &n)) {
                return (status);
            }
        }

        grp[k] = n;
        ipsci[k] = openInputImage (fdata, "SCI", n);
        iperr[k] = openInputImage (fdata, "ERR", n);
        ipdq[k]  = openInputImage (fdata, "DQ",  n);

        /* use the first image's attributes to compare with the rest of
         the files */
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
            if (getACSampxy (&prihdr, detector, chip, ccdamp, *dim_x, *dim_y,
                             &ampx, &ampy) )
                return (status);

            freeHdr (&scihdr);
        }

        /* verify the image size to be the same as the first image */
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
    *cumFlashDur = totalFlashDuration;
    *cumDarktime = totalDarktime;

    /* Record noise and gain values for use in rest of ACSREJ */
    for (i = 0; i < NAMPS; i++) {
        if (gn.val[i] == 0. ) {
            noise->val[i] = 0.;
        } else {
            noise->val[i] = ron.val[i];  /* electrons */
        }
        gain->val[i] = gn.val[i];  /* e/DN */
    }
    noise->colx = ampx;
    noise->coly = ampy;
    gain->colx = ampx;
    gain->coly = ampy;

    /* store the chip ID which is being processed */
    noise->chip = chip;
    gain->chip = chip;
    noise->detector = detector;
    gain->detector = detector;


    return (status);
}


static int getACSnsegn (Hdr *hdr, char *fdata, multiamp *ron, multiamp *gn) {

    extern int status;

    char    *nsekeyw[NAMPS] = { "READNSEA", "READNSEB", "READNSEC", "READNSED" };
    char    *gnkeyw[NAMPS] = { "ATODGNA", "ATODGNB", "ATODGNC", "ATODGND" };
    int     i;

    int     GetKeyFlt (Hdr *, char *, int, float, float *);

    /* get the readout noise values from the SCI header */
    for (i = 0; i < NAMPS; i++){
        if (GetKeyFlt (hdr, nsekeyw[i], USE_DEFAULT, 0., &ron->val[i]) != 0) {
            trlkwerr (nsekeyw[i], fdata);
            return (status = KEYWORD_MISSING);
        }
    }
    /* get the A-to-D gain values from the SCI header*/
    for (i = 0; i < NAMPS; i++){
        if (GetKeyFlt (hdr, gnkeyw[i], USE_DEFAULT, 0., &gn->val[i]) != 0) {
            trlkwerr (gnkeyw[i], fdata);
            return (status = KEYWORD_MISSING);
        }
    }
    return (status);
}

static int checkgn (multiamp gn, char *fdata) {
    extern int  status;

    int         i;
    int         ampset = 0;

    for (i = 0; i < NAMPS; i++){
        if (gn.val[i] > 0.) ampset++;
    }

    if (ampset == 0) {
        sprintf (MsgText, "All ATODGN values in file %s are 0.\n", fdata);
        trlerror (MsgText);
        return (status = INVALID_VALUE);
    }

    return (status);

}
static int getACSampxy (Hdr *hdr, int det, int chip, char *ccdamp, int dimx, int dimy, int *ampx, int *ampy) {

    extern int status;
    ACSInfo acsrej;
    char tabname[CHAR_LINE_LENGTH];

    void ACSInit (ACSInfo *);
    int GetKeyStr (Hdr *, char *, int, char *, char *, int);
    int GetCCDTab (ACSInfo *, int, int);

    ACSInit (&acsrej);

    tabname[0] = '\0';

    acsrej.detector = det;
    acsrej.chip = chip;
    strcpy(acsrej.ccdamp, ccdamp);

    /* Now, we need to read in CCDGAIN, CCDOFST[A-D], and BINAXIS[1,2]
        from hdr as selection criteria for row from CCDTAB.
        Get keyword values from primary header using same function
        used in ACSCCD.
    */
    if (getAndCheckACSKeys (&acsrej, hdr)) {
        freeHdr (hdr);
        return (status);
    }

    /* Now, read CCDTAB to get AMPX, AMPY values */
    if (GetKeyStr (hdr, "CCDTAB", USE_DEFAULT, "", tabname, ACS_FITS_REC))
        return (status);
    strcpy (acsrej.ccdpar.name, tabname);

    if (GetCCDTab (&acsrej, dimx, dimy) ) {
        return (status);
    }

    /* Now, assign AMPX and AMPY values read from table */
    *ampx = acsrej.ampx;
    *ampy = acsrej.ampy;

    return (status);
}
