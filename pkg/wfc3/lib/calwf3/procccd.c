# include <stdio.h>
# include <stdlib.h>	/* calloc */
# include <string.h>

#include "hstcal.h"
# include "hstio.h"

# include "wf3.h"
# include "calwf3.h"
# include "hstcalerr.h"
# include "wf3corr.h"
# include "wf3asn.h"	/* Contains association table structures */

/* ProcessCCD: This routine controls the overall flow of processing
   for CCD images.

   H.Bushouse, 2003 Apr 25:
   Modified to treat CRCORR and RPTCORR the same for WFC3 UVIS images.

   H.Bushouse, 2003 Oct 27:
   Modified to handle condition of no good data in WF3Rej (Same as CALACS).

   H.Bushouse, 2008 Oct 1:
   Fixed handling of EXPSCORR=PERFORM so that WF32D gets called for
   all images, and fixed save_tmp setting so that blv_tmp files get
   deleted after EXPSCORR processing.

   H.Bushouse, 2009 Jun 24:
   Added logic to always use CRCORR=PERFORM internally for both CRJ and
   RPT associations, instead of CRCORR for one and RPTCORR for the other.

   H.Bushouse, 2010 Apr 30: calwf3 v2.1
   Modified flow of processing so that CRCORR (wf3rej) executes BEFORE
   the individual exposures are processed with WF32D, so that DQ updates
   to blv_tmp files made by wf3rej will appear in final flt files.
   (PR 64963; Trac #545)

   H.Bushouse, 2010 May 10: calwf3 v2.1
   Fixed logic in new processing flow so that blv_tmp's don't get
   deleted before they're needed in wf32d and so blv_tmp's get copied to
   final flt if not performing wf32d.

   H.Bushouse, 2011 Jun 16: calwf3 v2.4
   Modified logic involved in handling error returns from WF3Rej so that
   WF32d processing still takes place for individual exposures if
   EXPSCORR=PERFORM. (PR 68593; Trac #722)

   M. Sosey, 2014 September: added CTE correction  and
   indented the code to make it easier to parse. onecpu was added
   to the function signature so that the user can have more control over
   the number of threads/cpus that are used during parallel processing for CTE correction.
   (see #1193 UVIS2.0 update)

   M. Sosey, 2015: fixed some logic associated with crcorr processing and
   removing the tmp files

   M. Sosey 2016: updated the trailer file concatination; the trailer files for
   association files are not being concatinated correctly.

 */

int ProcessCCD (AsnInfo *asn, WF3Info *wf3hdr, int *save_tmp, int printtime, int onecpu) {

    extern int status;

    RefFileInfo sciref;		    /* ref file keywords and names */
    CCD_Switch sci_sw;		    /* all cal switches for science file */
    CCD_Switch wf3ccd_sci_sw;  	/* WF3CCD switches for science file */
    CCD_Switch wf32d_sci_sw;   	/* WF32D  switches for science file */
    CCD_Switch wf3cte_sci_sw;   /* CTE switches for science file */

    int prod;			        /* which CR-split/Rptobs product
                                           are we working with?	*/
    int posid;			        /* counter for input products	*/
    int expid;			        /* counter for input exposures 	*/
    int newpreface = NO;		/* switch for keeping previous comments
                                           for all remaining trailer files */

    char *wf3rej_input;	        /* list of input names for WF3REJ */
    char *wf3rej_cte_input;     /* listof input CTE corrected names for WF3REJ*/
    char *wf3rej_msgtext;	    /* string for list of input filenames */
    int  nchars;		        /* Number of chars for input string to WF3REJ */

    int  WF3Rej_0   (char *, char *, char *, int, int, int);
    int  CopyFFile  (char *, char *);
    void SetCCDSw   (CCD_Switch *, CCD_Switch *, CCD_Switch *,CCD_Switch *);
    void ResetCCDSw (CCD_Switch *, CCD_Switch *);
    void ResetSwitch (CCD_Switch *, CCD_Switch *);
    void WF3Defaults (WF3Info *);
    void InitRefFile (RefFileInfo *);
    void FreeRefFile (RefFileInfo *);
    int  CCDRefInit (WF3Info *, CCD_Switch *, RefFileInfo *);
    int  WF3cte (char *, char *, CCD_Switch *, RefFileInfo *, int, int, int);
    int  WF3ccd (char *, char *, CCD_Switch *, RefFileInfo *, int, int);
    int  WF32d (char *, char *,CCD_Switch *, RefFileInfo *, int, int);
    int  GetAsnMember (AsnInfo *, int, int, int, WF3Info *);
    int  GetSingle (AsnInfo *, WF3Info *);
    int  CheckCorr (AsnInfo *, WF3Info *);
    int  InsertWF3Suffix (WF3Info *);
    int  updateAsnTable (AsnInfo *, int, int);

    /* initial value; */
    posid=0;
    wf3rej_input=NULL;
    wf3rej_cte_input=NULL;
    wf3rej_msgtext=NULL;

    /* Reset RPTCORR setting from ASN table to use CRCORR for UVIS */
    if (asn->rptcorr == PERFORM) {
        asn->crcorr = PERFORM;
        asn->rptcorr = OMIT;
    }


    /* LOOP OVER THE PRODUCTS/POSITIONS FOR EACH CR-SPLIT/REPEAT-OBS SET. */
    for (prod = 0; prod < asn->numprod; prod++) {
        if (asn->verbose){
            sprintf (MsgText, "CALWF3: processing UVIS product %d, spmems is %i, total products is %i", prod,asn->spmems[posid], asn->numprod);
            trlmessage (MsgText);
        }

        /* PROCESS THIS PARTIAL/SINGLE/FULL PRODUCT... */
        for (posid = 1; posid <= asn->numsp; posid++) {
            if (asn->verbose) {
                sprintf (MsgText,"CALWF3: processing posid = %d, num sub-products=%i", posid, asn->numsp);
                trlmessage (MsgText);
            }

            /*  Allocate space for WF3REJ NON-CTE input image list */
            if (asn->crcorr == PERFORM ) {

                nchars = asn->spmems[posid] * (CHAR_FNAME_LENGTH+1);
                wf3rej_input = (char *) calloc( nchars + 1, sizeof(char));
                wf3rej_input[0] = '\0';

                nchars = asn->spmems[posid] * (CHAR_FNAME_LENGTH+1);
                wf3rej_cte_input = (char *) calloc( nchars + 1, sizeof(char));
                wf3rej_cte_input[0] = '\0';

            }


            /*IS THIS A BUG? WITH 2 MEMBERS AND 1 PRODUCT IN THE ASN TABLE THE CODE IN WF3TABLE WHERE SPMEMS
              IS POPULATED ALWAYS SETS IT TO 1. EXPID ALWAYS ENDS UP AS 1 AND POSID INCREASES , BUT SPMEMS[POSID] IS ALWAYS 1
              SO THE THE FLT IMAGES GET CALIBRATED WITH THEIR OWN REFERENCE FILES INSTEAD OF THE INTENDED FIRST IMAGE REFERENCE
              FILES. IF THERE WERE 2 PRODUCTS IN THE TABLE WOULD THIS BE DIFFERENT? IT ACTUALLY IS THE PREFERRED PROCESS NOW. */
            for (expid = 1; expid <= asn->spmems[posid]; expid++) {

                /** FROM THIS POINT ON, WE ARE WORKING WITH ONE IMAGE AT A TIME...
                  READ IN MEMBER KEYWORD INFO FROM --FIRST-wf3rej- EXP ** IMAGE HEADER */

                if (expid == 1) {
                    /* (RE-)INITIALIZE THE LISTS OF REFERENCE FILE
                     ** KEYWORDS AND NAMES. */
                    InitRefFile (&sciref);

                    /* ASSIGN DEFAULT VALUES IN WF3. */
                    WF3Defaults (wf3hdr);

                    if (asn->process == SINGLE ) {
                        if (GetSingle (asn, wf3hdr) )
                            return (status);
                    } else {
                        if (GetAsnMember (asn, prod, posid, expid,
                                    wf3hdr))
                            return (status);
                    }

                    /* INITIALIZE VARIABLES, AND GET INFO (CALIBRATION
                     ** SWITCHES AND REFERENCE FILE NAMES) FROM PRIMARY
                     ** HEADERS. */
                    if (CCDRefInit (wf3hdr, &sci_sw, &sciref))
                        return (status);

                    /* MAKE SURE 'SAVE_TMP' NO LONGER
                     ** SET TO 'DUMMY'.
                     **
                     ** IF SAVE_TMP WAS SET FROM THE COMMAND-LINE, THEN
                     ** IGNORE WHAT IS SET IN THE HEADER KEYWORD
                     ** EXPSCORR. */

                    if (*save_tmp == DUMMY) {
                        *save_tmp = NO;
                    }

                    if (asn->verbose) {
                        trlmessage
                            ("CALWF3: Got reference file information");
                    }
                    /* Store the trailer file comments into preface */
                    newpreface = YES;

                } else { /*EXPID ==1 END IF*/

                    if (asn->process == SINGLE ) {

                        if (GetSingle (asn, wf3hdr))
                            return (status);
                    } else {
                        if (GetAsnMember (asn, prod, posid, expid,wf3hdr))
                            return (status);
                    }

                    /* CONSTRUCT OUTPUT AND TEMPORARY FILE NAMES. */
                    if (InsertWF3Suffix (wf3hdr))
                        return (status);
                } /*END ELSE FOR EXPID != 1*/


                /* CHECK ASN TABLE SETTINGS FOR CRCORR/RPTCORR
                 ** AGAINST IMAGE HEADER VALUES.  IF THEY ARE NOT
                 ** CONSISTENT, PRINT OUT A ERROR MESSAGE AND QUIT.
                 ** ONLY NEED TO PERFORM THIS CHECK ON THE FIRST IMAGE
                 ** IN THE ASSOCIATION. WJH 21 MAY 1999 */
                if (prod == 0 && posid == 1 && expid == 1) {
                    if (CheckCorr(asn, wf3hdr))
                        return (status);
                }

                /* COPY SWITCH VALUES FOR WF3CCD AND WF32D AND WF3CTE ARGUMENT
                 ** LISTS. */
                SetCCDSw (&sci_sw, &wf3ccd_sci_sw, &wf3cte_sci_sw, &wf32d_sci_sw);

                if (asn->verbose) {
                    trlmessage ("CALWF3: Processing switches set... ");
                }

                /* READY TO PROCESS INPUT IMAGE NOW... */
                /* DO WE HAVE ANYTHING TO DO WITH THIS IMAGE? 	*/
                if (wf3hdr->sci_basic_2d  != PERFORM &&
                        wf3hdr->sci_basic_ccd != PERFORM &&
                        wf3hdr->sci_crcorr    != PERFORM &&
                        wf3hdr->sci_rptcorr   != PERFORM &&
                        wf3hdr->sci_basic_cte  != PERFORM) {

                    trlwarn("No calibration switch was set to PERFORM.");
                    FreeRefFile (&sciref);
                    return (status = NOTHING_TO_DO);

                }

                /* COPY ALL MESSAGES UP TO HERE INTO PREFACE BUFFER
                   TO BE PREPENDED TO ALL INPUT FILE TRAILER FILES */
                if (newpreface == YES) {
                    InitTrlPreface();
                }
                SetTrlPrefaceMode (YES);

                if (wf3hdr->sci_basic_cte == PERFORM ) {

                    /*correct for CTE issues and THEN complete the calibration
                      also run through the pipeline completely without CTE corr after,
                      so in this case, two sets of fully calibrated output images will be created

                      This step processes both chips together
                      and has some of its own calibration files so it's
                      outside the single chip loop below

                      This also occurs before DoCCD which populates
                      all the useful data structures :/

                      So wfc3 version of cte correction occurs
                      as the first thing, instead of like ACS which
                      happens after acsccd and the blevcorr step.

                     */


                    if ( WF3cte(wf3hdr->rawfile, wf3hdr->rac_tmp, &sci_sw, &sciref, printtime, asn->verbose, onecpu) )
                        return (status);

                    if (wf3hdr->sci_basic_ccd == PERFORM) {

                        /* THIS IS CALWF3CCD; NOTE WE USE WF3CCD_SCI_SW. */

                        /*DO WF3CCD STUFF ON THE RAC IMAGE */
                        if (WF3ccd (wf3hdr->rac_tmp, wf3hdr->blc_tmp,
                                    &wf3ccd_sci_sw, &sciref, printtime,
                                    asn->verbose))
                                return (status);

                            /* RESET SWITCHES FOR NON-CTE PROCESSING RUN */
                            ResetSwitch(&sci_sw, &wf3ccd_sci_sw);

                    } else {

                        /* COPY THE INPUT FILE TO BLC_TMP. WE NEED TO DO THIS
                           BECAUSE WF3REJ MAY (DEPENDING ON THE CRMASK COLUMN
                           IN THE CRREJTAB) MODIFY THE DQ FLAGS IN ITS INPUT
                           FILE. */
                        if (CopyFFile (wf3hdr->rac_tmp, wf3hdr->blc_tmp))
                            return (status);
                    }



                    /* WE ARE NOW WORKING WITH BLC_TMP HERE... */
                    /* COPY BLC_TMP NAME TO WF3REJ_INPUT STRING */
                    if (wf3hdr->sci_crcorr == PERFORM ||  wf3hdr->sci_rptcorr == PERFORM) {
                        strcat (wf3rej_cte_input, wf3hdr->blc_tmp);
                        if (expid < asn->spmems[posid]) {
                            /* DON'T ADD A COMMA TO THE END OF THE LAST FILENAME */
                            strcat(wf3rej_cte_input, ",");
                        }

                        /* ALSO ADD BLC_TMP and RAC_TMP TO LIST OF FILES TO BE DELETED */
                        strcpy (asn->product[prod].subprod[posid].exp[expid].blc_tmp,
                                wf3hdr->blc_tmp);
                        strcpy (asn->product[prod].subprod[posid].exp[expid].rac_tmp,
                                wf3hdr->rac_tmp);

                    }

                }  /*END CTE PROCESSING TO BLC_TMP LEVEL*/


                /*ALWAYS REPEAT THE PROCESS WITHOUT THE CTE CORRECTION SO BOTH ARE PRODUCED*/
                if (wf3hdr->sci_basic_ccd == PERFORM) {

                    if (WF3ccd (wf3hdr->rawfile, wf3hdr->blv_tmp,
                                &wf3ccd_sci_sw, &sciref, printtime,
                                asn->verbose))
                        return (status);

                        /* RESET SWITCHES */
                        ResetSwitch(&sci_sw, &wf3ccd_sci_sw);

                } else {

                    /* Copy the input file to blv_tmp. We need to do this
                     ** because wf3rej may (depending on the crmask column
                     ** in the crrejtab) modify the DQ flags in its input
                     ** file. */
                    if (CopyFFile (wf3hdr->rawfile, wf3hdr->blv_tmp))
                        return (status);
                }


                /* WE ARE NOW WORKING WITH BLV_TMP HERE... */
                /* COPY BLV_TMP NAME TO WF3REJ_INPUT STRING */
                if (wf3hdr->sci_crcorr == PERFORM ||  wf3hdr->sci_rptcorr == PERFORM) {
                    strcat (wf3rej_input, wf3hdr->blv_tmp);
                    if (expid < asn->spmems[posid]) {
                        /* DON'T ADD A COMMA TO THE END OF THE LAST FILENAME */
                        strcat(wf3rej_input, ",");
                    }

                    /* ALSO ADD BLV_TMP TO LIST OF FILES TO BE DELETED */
                    strcpy (asn->product[prod].subprod[posid].exp[expid].blv_tmp,
                            wf3hdr->blv_tmp);
                }



                /* IF WE ARE NOT PERFORMING CRCORR OR RPTCORR, THEN
                   FINISH PROCESSING INDIVIDUAL FILES WITH WF32D. */

                if (wf3hdr->sci_crcorr != PERFORM && wf3hdr->sci_rptcorr != PERFORM) {
                    if (asn->debug) {
                        trlmessage ("UVIS crcorr/rptcorr not set to Perform...");
                    }

                    /*DO THIS WITHOUT CTE PERFORM*/
                    if (wf3hdr->sci_basic_2d == PERFORM) {
                        SetTrlPrefaceMode (NO);

                        /* BASIC 2-D PROCESSING (FLAT FIELD, ETC). */
                        if (WF32d (wf3hdr->blv_tmp, wf3hdr->fltfile,
                                    &wf32d_sci_sw, &sciref, printtime,
                                    asn->verbose))
                            return(status);


                    } else {
                        /* REMEMBER BLV_TMP AS FINAL OUTPUT...*/
                        if (CopyFFile (wf3hdr->blv_tmp, wf3hdr->fltfile))
                            return (status);
                    }


                    /* PROCESS CTE DATA */
                    if (wf3hdr->sci_basic_2d == PERFORM && wf3hdr->sci_basic_cte == PERFORM) {
                        SetTrlPrefaceMode (NO);

                        /* Basic 2-D processing (flat field, etc). */
                        if (WF32d (wf3hdr->blc_tmp, wf3hdr->flcfile,
                                    &wf32d_sci_sw, &sciref, printtime,
                                    asn->verbose))
                            return (status);

                    } else {/* REMEMBER BLC_TMP AS FINAL OUTPUT...*/
                        if ( wf3hdr->sci_basic_cte == PERFORM){
                            if (CopyFFile (wf3hdr->blc_tmp, wf3hdr->flcfile))
                                return (status);
                        }
                    }/*END PROCESS CTE WITH WF32D*/

                    if (*save_tmp == NO){
                        if (wf3hdr->sci_basic_cte == PERFORM) {
                            remove (wf3hdr->blc_tmp);
                            remove (wf3hdr->rac_tmp);
                        }
                        remove (wf3hdr->blv_tmp);
                    }
                }

            }/* END EXPID  LOOP OVER INDIVIDUAL EXPOSURES,GO ON TO CREATE PRODUCTS */

            /* Reset the trailer file preface to NULL since
            ** it has already been copied into trailer files */
            ResetTrlPreface();

            /*** DO CRCORR/RPTCORR PROCESSING ***/
            if (wf3hdr->sci_crcorr == PERFORM ||  wf3hdr->sci_rptcorr == PERFORM) {

                /* COSMIC RAY REJECTION, FOLLOWED BY BASIC 2-D PROCESSING. */
                if (asn->verbose) {
                    sprintf (MsgText,
                            "CALWF3: Now processing position %d from product %d",
                            posid, prod);
                    trlmessage (MsgText);
                }

                /* COPY SWITCH VALUES FOR WF3CCD AND WF32D ARGUMENT LISTS. */
                SetCCDSw (&sci_sw, &wf3ccd_sci_sw, &wf3cte_sci_sw, &wf32d_sci_sw);

                /* RESET DQICORR, SINCE IT WILL ALREADY HAVE BEEN DONE */
                wf32d_sci_sw.dqicorr = COMPLETE;

                if (asn->debug || asn->verbose) {
                    sprintf (MsgText,"Non-cte input to WF3REJ is:");
                    trlmessage (MsgText);


                    /* NEED TO ALLOCATE MEMORY FOR A SEPARATE STRING TO
                     ** BE LONG ENOUGH TO HOLD ALL INPUT NAMES WHEN
                     ** PRINTING IT OUT. CAUSES PIPELINE PROBLEMS
                     ** OTHERWISE. HAB 20-JUN-2004 */
                    wf3rej_msgtext = (char *) calloc(strlen(wf3rej_input)+25,sizeof(char));
                    sprintf (wf3rej_msgtext, "%s", wf3rej_input);
                    trlmessage (wf3rej_msgtext);
                    free (wf3rej_msgtext);
                    sprintf (MsgText,"Output from WF3REJ is: %s",
                            wf3hdr->crj_tmp);
                    trlmessage (MsgText);
                }

                /* SET UP THE CORRECT VALUE OF ASN_MTYP FOR WF3REJ
                   OUTPUT IMAGE. EXP-CRJ*/
                if (strncmp (wf3hdr->asn_table, wf3hdr->crj_tmp, 9) ) {
                    /* CRJ_TMP IS A SUB-PRODUCT */
                    strcpy (wf3hdr->mtype,
                            asn->product[prod].subprod[posid].mtype);
                } else {

                    /* CRJ_TMP LEADS TO A PRODUCT WITH SAME ROOTNAME AS
                     ** ASN_TABLE */
                    strcpy (wf3hdr->mtype, asn->product[prod].mtype);
                }

                /* REJECT COSMIC RAYS FOR NON-CTE DATA AND UPDATE THE OUTPUT SPT */
                if (WF3Rej_0 (wf3rej_input, wf3hdr->crj_tmp,
                            wf3hdr->mtype, printtime, asn->verbose,1)) {
                    if (status == NO_GOOD_DATA || status == NOTHING_TO_DO) {
                        /* Set CRCORR to skipped so that we don't try to
                           apply WF32d to crj_tmp product */
                        wf3hdr->sci_crcorr = SKIPPED;
                        /* Reset status to good now that we have dealt
                         ** with this condition */
                        status = WF3_OK;
                    } else {
                        return (status);
                    }
                }

                /* FREE UP MEMORY USED BY WF3REJ_INPUT FOR THIS
                 ** SUBPRODUCT */
                free(wf3rej_input);

                /*THE ASN TABLE IS ONLY UPDATED FOR NON-CTE DATA*/
                if (updateAsnTable (asn, prod, posid))
                    return (status);


                /*NOW DO CR-REJECTION FOR THE CTE DATA*/
                if (wf3hdr->sci_basic_cte == PERFORM){

                    if (asn->debug || asn->verbose) {
                        /* NEED TO ALLOCATE MEMORY FOR A SEPARATE STRING TO
                         ** BE LONG ENOUGH TO HOLD ALL INPUT NAMES WHEN
                         ** PRINTING IT OUT. CAUSES PIPELINE PROBLEMS
                         ** OTHERWISE. HAB 20-JUN-2004 */

                        wf3rej_msgtext = (char *) calloc(strlen(wf3rej_cte_input)+25,sizeof(char));
                        sprintf (wf3rej_msgtext, "%s", wf3rej_cte_input);
                        trlmessage (wf3rej_msgtext);
                        free(wf3rej_msgtext);
                    }

                        /*the mtype was already set in the non-cte loop to EXP-CRJ*/

                    /* REJECT COSMIC RAYS. */
                    if (WF3Rej_0 (wf3rej_cte_input, wf3hdr->crc_tmp,
                                wf3hdr->mtype, printtime, asn->verbose, 0)) {
                        if (status == NO_GOOD_DATA || status == NOTHING_TO_DO) {
                            /* SET CRCORR TO SKIPPED SO THAT WE DON'T TRY TO
                               APPLY WF32D TO CRC_TMP PRODUCT */
                            wf3hdr->sci_crcorr = SKIPPED;
                            /* RESET STATUS TO GOOD NOW THAT WE HAVE DEALT
                             ** WITH THIS CONDITION */
                            status = WF3_OK;
                        } else {
                            return (status);
                        }
                    }
                    
                    /* FREE UP MEMORY USED BY WF3REJ_INPUT FOR THIS
                     ** SUBPRODUCT */
                    free(wf3rej_cte_input);

                }/*END REJ FOR CTE DATA*/

                trlmessage("Starting WF32D");

                /* MOVE ON TO WF32D */

                if (wf3hdr->sci_basic_2d == PERFORM &&  wf3hdr->sci_crcorr == PERFORM) {

                   /* FLATFIELD THE SUMMED, COSMIC RAY REJECTED IMAGE. */

                    if (WF32d (wf3hdr->crj_tmp, wf3hdr->crjfile,
                                &wf32d_sci_sw, &sciref, printtime, asn->verbose))
                        return (status);

                    printf("\n**** wf3hdr->sci_basic_cte is %i *****\n\n",wf3hdr->sci_basic_cte);

                    if (wf3hdr->sci_basic_cte == PERFORM){
                        /* Flatfield the summed, cosmic ray rejected CTE image. */
                        if (WF32d (wf3hdr->crc_tmp, wf3hdr->crcfile,
                                    &wf32d_sci_sw, &sciref, printtime, asn->verbose))
                            return (status);
                    }

                } else if ( wf3hdr->sci_crcorr == PERFORM){

                    if (wf3hdr->sci_basic_cte == PERFORM) {
                        /* REMEMBER CR-COMBINED IMAGE AS FINAL OUTPUT NAME,
                         ** SINCE WF32D WAS NOT RUN. */
                        if (CopyFFile (wf3hdr->crc_tmp, wf3hdr->crcfile))
                            return (status);
                    }

                    if (CopyFFile (wf3hdr->crj_tmp, wf3hdr->crjfile))
                        return (status);
                }

                /* REMEMBER WHAT CRC_TMP FILES TO DELETE */
                if (wf3hdr->sci_basic_cte == PERFORM) {
                    strcpy (asn->product[prod].subprod[posid].crc_tmp,
                        wf3hdr->crc_tmp);
                }

                /* REMEMBER WHAT CRJ_TMP FILES TO DELETE */
                strcpy (asn->product[prod].subprod[posid].crj_tmp,   wf3hdr->crj_tmp);

                trlmessage("Starting EXPSCORR loops...");

                /* IF EXPSCCORR = PERFORM, FINISH PROCESSING INDIVIDUAL EXPOSURES*/
                if (sci_sw.expscorr == PERFORM) {

                    for (expid=1; expid <= asn->spmems[posid]; expid++) {

                        /* GET INPUT & OUTPUT FILE NAMES FOR THIS EXP */
                        if (asn->process == SINGLE) {
                            if (GetSingle (asn, wf3hdr))
                                return (status);
                        } else {
                            if (GetAsnMember (asn, prod, posid, expid,
                                        wf3hdr))
                                return (status);
                        }
                        if (InsertWF3Suffix (wf3hdr))
                            return (status);

                        if (wf3hdr->sci_basic_2d == PERFORM){
                            SetTrlPrefaceMode(YES);

                            /* RESET SWITCHES*/
                            ResetCCDSw (&wf3ccd_sci_sw, &wf32d_sci_sw);

                            /* Basic 2-D processing (flat field, etc). */
                            if (WF32d (wf3hdr->blv_tmp, wf3hdr->fltfile,
                                        &wf32d_sci_sw, &sciref, printtime,
                                        asn->verbose))
                                return (status);

                            if (wf3hdr->sci_basic_cte == PERFORM){
                                if (WF32d (wf3hdr->blc_tmp, wf3hdr->flcfile,
                                        &wf32d_sci_sw, &sciref, printtime,
                                        asn->verbose))
                                    return (status);
                             }

                        } else {

                            /* SAVE BLV_TMP AS FINAL OUTPUT */
                            if (CopyFFile (wf3hdr->blv_tmp,
                                        wf3hdr->fltfile))
                                return (status);

                            if (wf3hdr->sci_basic_cte == PERFORM){
                                if (CopyFFile (wf3hdr->blc_tmp,
                                        wf3hdr->flcfile))
                                    return (status);
                            }
                        }

                    } /* END LOOP OVER INDIVDUAL EXPOSURES */
                } /* END OF EXPSCORR PROCESSING */


                /* DELETE INTERMEDIATE FILES */
                if (*save_tmp == NO) {
                    for (expid=1; expid <= asn->spmems[posid]; expid++) {
                        remove (asn->product[prod].subprod[posid].exp[expid].blv_tmp);
                    }
                }

                if (wf3hdr->sci_basic_cte == PERFORM) {
                    /* Delete intermediate files */
                    if (*save_tmp == NO) {
                        for (expid=1; expid <= asn->spmems[posid]; expid++) {
                            remove (asn->product[prod].subprod[posid].exp[expid].blc_tmp);
                            remove (asn->product[prod].subprod[posid].exp[expid].rac_tmp);
                        }
                    }
                 }

             } /* End of CRCORR/RPTCORR processing*/

        }  /* End of processing for SINGLE/PARTIAL product, END POSID */
        /* DONE WITH THE TMP FILES..DELETE AS ORDERED */

        if (asn->process != SINGLE) {
            /* We only have something to delete if process != SINGLE */
            if ((*save_tmp == NO ) && (wf3hdr->sci_crcorr == PERFORM ||
                        wf3hdr->sci_rptcorr == PERFORM)) {
                for (posid = 1; posid <= asn->numsp; posid++) {
                    remove (asn->product[prod].subprod[posid].crj_tmp);
                }
                if (wf3hdr->sci_basic_cte == PERFORM){
                    for (posid = 1; posid <= asn->numsp; posid++){
                        remove (asn->product[prod].subprod[posid].crc_tmp);
                    }
                }
            }
        } /*END REMOVE TMP FILES*/

    } /* END LOOP OVER CCD PRODUCTS HERE...	END PROD	*/

    /* DONE WITH LISTS OF REFERENCE FILE KEYWORDS AND NAMES. */
    FreeRefFile (&sciref);
    return (status);
}
