# include <stdio.h>
# include <stdlib.h>	/* calloc */
# include <string.h>

# include "hstio.h"

# include "wf3.h"
# include "calwf3.h"
# include "wf3err.h"
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
 */

int ProcessCCD (AsnInfo *asn, WF3Info *wf3hdr, int *save_tmp, int printtime, int onecpu) {

    extern int status; 

    RefFileInfo sciref;		    /* ref file keywords and names */
    CCD_Switch sci_sw;		    /* all cal switches for science file */
    CCD_Switch wf3ccd_sci_sw;  	/* WF3CCD switches for science file */
    CCD_Switch wf32d_sci_sw;   	/* WF32D  switches for science file */
    CCD_Switch wf3cte_sci_sw;   /* CTE switches for science file */
    int save_crj;			    /* don't delete crj_tmp file? */
    int save_rac;               /* don't delete the RAC file? */

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

    int  WF3Rej_0   (char *, char *, char *, int, int);
    int  CopyFFile  (char *, char *);
    void SetCCDSw   (CCD_Switch *, CCD_Switch *, CCD_Switch *,CCD_Switch *);
    void ResetCCDSw (CCD_Switch *, CCD_Switch *);
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

    save_crj = *save_tmp;		/* initial value; */	
    save_rac = *save_tmp;

    posid=0;	/* initial value; */	

    /* Reset RPTCORR setting from ASN table to use CRCORR for UVIS */
    if (asn->rptcorr == PERFORM) {
        asn->crcorr = PERFORM;
        asn->rptcorr = OMIT;
    }

    /* Loop over the products/positions for each CR-split/Repeat-Obs set. */

    for (prod = 0; prod < asn->numprod; prod++) {
        if (asn->verbose){
            sprintf (MsgText, "CALWF3: processing UVIS product %d, spmems is %i, total products is %i", prod,asn->spmems[posid], asn->numprod);	
            trlmessage (MsgText);
        }	

        /*see bug notes in next for */
        /* Process this PARTIAL/SINGLE/FULL product... */
        for (posid = 1; posid <= asn->numsp; posid++) {
            if (asn->verbose) {
                sprintf (MsgText,"CALWF3: processing posid = %d, num sub-products=%i", posid, asn->numsp);
                trlmessage (MsgText);
            }

            /*  Allocate space for WF3REJ input image list */
            if (asn->crcorr == PERFORM || asn->rptcorr == PERFORM) {
                nchars = asn->spmems[posid] * (SZ_FNAME+1);
                wf3rej_input = (char *) calloc( nchars + 1, sizeof(char));

                /* Initialize this string to NULL */
                wf3rej_input[0] = '\0';
            }

            /*  Allocate space for WF3REJ CTE input image list */
            if (asn->crcorr == PERFORM || asn->rptcorr == PERFORM) {
                if (wf3hdr->sci_basic_cte == PERFORM){
                    nchars = asn->spmems[posid] * (SZ_FNAME+1);
                    wf3rej_cte_input = (char *) calloc( nchars + 1, sizeof(char));

                    /* Initialize this string to NULL */
                    wf3rej_cte_input[0] = '\0';
                }
            }

            /*is this a bug? with 2 members and 1 product in the asn table the code in wf3table where spmems 
              is populated always sets it to 1. expid always ends up as 1 and posid increases , but spmems[posid] is always 1
              so the the flt images get calibrated with their own reference files instead of the intended first image reference
              files. If there were 2 products in the table would this be different? It actually is the preferred process now. */			
            for (expid = 1; expid <= asn->spmems[posid]; expid++) {

                /** From this point on, we are working with ONE image at a time... **/

                /* Read in Member keyword info from --FIRST-- EXP ** image header */	

                /*if (prod == 0 && posid == 1 && expid == 1) {*/
                if (expid == 1) {
                    /* (Re-)Initialize the lists of reference file
                     ** keywords and names. */
                    InitRefFile (&sciref);

                    /* Assign default values in wf3. */
                    WF3Defaults (wf3hdr);

                    if (asn->process == SINGLE ) {
                        if (GetSingle (asn, wf3hdr) )
                            return (status);
                    } else {	
                        if (GetAsnMember (asn, prod, posid, expid,
                                    wf3hdr)) 

                            return (status);
                    } 

                    /* Initialize variables, and get info (calibration
                     ** switches and reference file names) from primary
                     ** headers. */
                    if (CCDRefInit (wf3hdr, &sci_sw, &sciref))
                        return (status);


                    /* Make sure 'save_tmp' and 'save_crj' are no longer
                     ** set to 'DUMMY'. 
                     **
                     ** If save_tmp was set from the command-line, then 
                     ** ignore what is set in the header keyword
                     ** EXPSCORR. */

                    if (*save_tmp == DUMMY) {
                        *save_tmp = NO;
                        save_crj = *save_tmp;
                        save_rac = *save_tmp;
                    }

                    if (asn->verbose) {
                        trlmessage 
                            ("CALWF3: Got reference file information");
                    }
                    /* Store the trailer file comments into preface */
                    newpreface = YES;

                } else { /*expid ==1 end if*/

                    if (asn->process == SINGLE ) {

                        if (GetSingle (asn, wf3hdr))
                            return (status);
                    } else {	               
                        if (GetAsnMember (asn, prod, posid, expid,wf3hdr)) 
                            return (status);
                    } 

                    /* Construct output and temporary file names. */
                    if (InsertWF3Suffix (wf3hdr))
                        return (status);	    		
                } /*end else for expid != 1*/

                /* Check ASN table settings for CRCORR/RPTCORR 
                 ** against image header values.  If they are not
                 ** consistent, print out a ERROR message and quit.
                 ** Only need to perform this check on the first image
                 ** in the association. WJH 21 May 1999 */
                if (prod == 0 && posid == 1 && expid == 1) {
                    if (CheckCorr(asn, wf3hdr))
                        return (status);
                }                

                /* Copy switch values for WF3CCD and WF32D and WF3CTE argument
                 ** lists. */
                SetCCDSw (&sci_sw, &wf3ccd_sci_sw, &wf3cte_sci_sw, &wf32d_sci_sw);

                if (asn->verbose) {
                    trlmessage ("CALWF3: Processing switches set... ");
                }

                /* Ready to process input image now... */
                /* Do we have anything to do with this image? 	*/
                if (wf3hdr->sci_basic_2d  != PERFORM &&
                        wf3hdr->sci_basic_ccd != PERFORM &&
                        wf3hdr->sci_crcorr    != PERFORM &&
                        wf3hdr->sci_rptcorr   != PERFORM &&
                        wf3hdr->sci_basic_cte  != PERFORM) {

                    trlwarn 
                        ("No calibration switch was set to PERFORM.");
                    FreeRefFile (&sciref);
                    return (status = NOTHING_TO_DO);

                } else {

                    if (wf3hdr->sci_basic_cte == PERFORM) {

                        /*correct for CTE issues and THEN complete the calibration
                          also run through the pipeline completely without CTE corr after,
                          so in this case, two sets of fully calibrated output images will be created*/                 

                        if ( WF3cte(wf3hdr->rawfile, wf3hdr->rac_tmp, &sci_sw, &sciref, printtime, asn->verbose, onecpu) )
                            return (status);                

                        /* First, check for CTE correction, then intialilize the Input CCD image with WF3CCD.
                         ** Then, either do cosmic ray rejection followed by
                         ** WF32D, or just run WF32D. */

                        /*check to see if PCTECORR is on
                          This step processes both chips together
                          and has some of its own calibration files so it's
                          outside the single chip loop below 

                          This also occurs before DoCCD which populates
                          all the useful data structures :/

                          So wfc3 version of cte correction occurs
                          as the first thing, instead of like ACS which
                          happens after acsccd and the blevcorr step. This
                          means I have to read in all the information just
                          for the pctestep.

                          Is it better to move the population of wf3 structure
                          higher in the code or just live with it or keep it the
                          same as the other instrument pipelines for pseudo clarity
                          to others?

                         */

                        if (wf3hdr->sci_basic_ccd == PERFORM) {

                            /* This is calwf3ccd; note we use wf3ccd_sci_sw. */

                            /* Copy all messages up to here into preface buffer
                             ** to be prepended to all input file trailer files */
                            if (newpreface == YES) {
                                InitTrlPreface();
                            }
                            SetTrlPrefaceMode (YES);

                            /*do wf3ccd stuff on the RAC image inside wf3cte instead */
                            if (WF3ccd (wf3hdr->rac_tmp, wf3hdr->blc_tmp,
                                        &wf3ccd_sci_sw, &sciref, printtime,
                                        asn->verbose)) 		
                                return (status);

                            /* Reset switches */
                            trlmessage("resetting ccd switches");
                            ResetCCDSw (&wf3ccd_sci_sw, &wf32d_sci_sw);
                        }
                        trlmessage("Finished resetting ccd switch");
                        
                        /* We now are working with blc_tmp here... */
                        /* Copy (cat) blc_tmp name to wf3rej_input string */
                        if (wf3hdr->sci_crcorr == PERFORM || 
                                wf3hdr->sci_rptcorr == PERFORM) {
                                strcat (wf3rej_cte_input, wf3hdr->blc_tmp);
                                trlmessage("strcat blc_tmp");
                                if (expid < asn->spmems[posid]) {
                                    trlmessage("strcat comma");
                                    /* Don't add a comma to the end of the last
                                     ** filename */
                                    strcat(wf3rej_cte_input, ",");
                                }
                            trlmessage("Updating delete file list for blc_tmp");
                            /* Also add blc_tmp to list of files to be deleted */
                            strcpy (asn->product[prod].subprod[posid].exp[expid].blc_tmp,
                                wf3hdr->blc_tmp);
                        }
                    }
                    
         
                    /*ALWAYS REPEAT THE PROCESS WITHOUT THE CTE CORRECTION SO BOTH ARE PRODUCED*/
                    trlmessage("Should be starting non-cte ccd reduction now");
                    if (wf3hdr->sci_basic_ccd == PERFORM) {
                        trlmessage("Starting  non-cte reduction");

                        /* This is calwf3ccd; note we use wf3ccd_sci_sw. */

                        /* Copy all messages up to here into preface buffer
                         ** to be prepended to all input file trailer files */

                        trlmessage("Calling wf3ccd for original raw file and no CTE");
                        if (WF3ccd (wf3hdr->rawfile, wf3hdr->blv_tmp,
                                    &wf3ccd_sci_sw, &sciref, printtime,
                                    asn->verbose)) 		
                            return (status);

                        /* Reset switches */
                        ResetCCDSw (&wf3ccd_sci_sw, &wf32d_sci_sw);

                        /* We now are working with blv_tmp here... */
                        /* Copy (cat) blv_tmp name to wf3rej_input string */
                        if (wf3hdr->sci_crcorr == PERFORM || 
                                wf3hdr->sci_rptcorr == PERFORM) {
                            strcat (wf3rej_input, wf3hdr->blv_tmp);
                            if (expid < asn->spmems[posid]) {
                                /* Don't add a comma to the end of the last
                                 ** filename */
                                strcat(wf3rej_input, ",");
                            }

                            /* Also add blv_tmp to list of files to be deleted */
                            strcpy (asn->product[prod].subprod[posid].exp[expid].blv_tmp,
                                    wf3hdr->blv_tmp);
                        }


                    } else {

                        /* Copy the input file to blv_tmp. We need to do this
                         ** because wf3rej may (depending on the crmask column
                         ** in the crrejtab) modify the DQ flags in its input
                         ** file. */
                        if (CopyFFile (wf3hdr->rawfile, wf3hdr->blv_tmp))
                            return (status);
                    }



                    /* If we are not performing CRCORR or RPTCORR, then
                     ** finish processing  CTE corrected individual files with WF32D. */ 
                    if (wf3hdr->sci_crcorr != PERFORM && wf3hdr->sci_rptcorr != PERFORM) { 
                        if (asn->debug) {
                            trlmessage ("UVIS crcorr/rptcorr not set to Perform...");
                        }

                        if (wf3hdr->sci_basic_2d == PERFORM && wf3hdr->sci_basic_cte == PERFORM) {

                            SetTrlPrefaceMode (NO);

                            /* Basic 2-D processing (flat field, etc). */
                            trlmessage("Starting WF32d with blc_tmp file");
                            if (WF32d (wf3hdr->blc_tmp, wf3hdr->flcfile, 
                                        &wf32d_sci_sw, &sciref, printtime,
                                        asn->verbose)) 
                                return (status);				
                            /* Remember blc_tmp as final output...*/
                            if (CopyFFile (wf3hdr->blc_tmp, wf3hdr->flcfile))
                                return (status);

                        } 
                            
                        /*DO this with or without CTE perform*/
                        if (wf3hdr->sci_basic_2d == PERFORM ) {

                            SetTrlPrefaceMode (NO);

                            /* Basic 2-D processing (flat field, etc). */
                            if (WF32d (wf3hdr->blv_tmp, wf3hdr->fltfile, 
                                        &wf32d_sci_sw, &sciref, printtime,
                                        asn->verbose)) 
                                 return(status);
                                 
                             /* Remember blv_tmp as final output...*/
                            if (CopyFFile (wf3hdr->blv_tmp, wf3hdr->fltfile))
                                return (status);
                            

                        } 

                        
                    } 

                    /* Only delete these if not doing CRCORR/RPTCORR */
                    if ((*save_tmp != YES) &&
                            (wf3hdr->sci_crcorr != PERFORM) &&
                            (wf3hdr->sci_rptcorr != PERFORM)){
                        if (wf3hdr->sci_basic_cte == PERFORM) {
                            remove (wf3hdr->blc_tmp);
                            remove (wf3hdr->blv_tmp);
                        } else {
                            remove (wf3hdr->blv_tmp);
                        }

                    }
                }  
                
            }   /* End loop over individual EXPosures,go on to create products */

            /* Reset the trailer file preface to NULL since
             ** it has already been copied into trailer files */
            ResetTrlPreface();


            /* Do CRCORR/RPTCORR processing  */
            if (wf3hdr->sci_crcorr == PERFORM || 
                    wf3hdr->sci_rptcorr == PERFORM) {

                /* Cosmic ray rejection, followed by basic 2-D
                 ** processing. */
                if (asn->verbose) {
                    sprintf (MsgText,
                            "CALWF3: Now process position %d from product %d",
                            posid, prod);
                    trlmessage (MsgText);
                }

                /* Copy switch values for WF3CCD and WF32D argument
                 ** lists. */
                SetCCDSw (&sci_sw, &wf3ccd_sci_sw, &wf3cte_sci_sw, &wf32d_sci_sw);	

                /* Reset DQICORR, since it will already have been done */
                wf32d_sci_sw.dqicorr = COMPLETE;

                if (asn->debug || asn->verbose) {
                    sprintf (MsgText,"Input to WF3REJ is:");
                    trlmessage (MsgText);

                    /* Need to allocate memory for a separate string to
                     ** be long enough to hold all input names when
                     ** printing it out. Causes pipeline problems
                     ** otherwise. HAB 20-Jun-2004 */

                    wf3rej_msgtext = calloc(strlen(wf3rej_input)+25,
                            sizeof(char));
                    sprintf (wf3rej_msgtext, "%s", wf3rej_input);
                    trlmessage (wf3rej_msgtext);
                    free (wf3rej_msgtext);
                    sprintf (MsgText,"Output to WF3REJ is: %s",
                            wf3hdr->crj_tmp);
                    trlmessage (MsgText);
                }

                /* Set up the correct value of ASN_MTYP for WF3REJ
                 ** output image. */                
                if (strncmp (wf3hdr->asn_table, wf3hdr->crj_tmp, 9) ) {

                    /* crj_tmp is a sub-product */
                    strcpy (wf3hdr->mtype,
                            asn->product[prod].subprod[posid].mtype);
                } else {

                    /* crj_tmp leads to a product with same rootname as
                     ** asn_table */
                    strcpy (wf3hdr->mtype, asn->product[prod].mtype);
                }   

                /* Reject cosmic rays. */
                if (WF3Rej_0 (wf3rej_input, wf3hdr->crj_tmp,
                            wf3hdr->mtype, printtime, asn->verbose)) {
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

                if (updateAsnTable (asn, prod, posid))
                    return (status);

                /* Free up memory used by wf3rej_input for this
                 ** subproduct */
                free(wf3rej_input); 


                /*NOW DO IT FOR THE CTE DATA*/
                if (asn->debug || asn->verbose) {
                    sprintf (MsgText,"Input to WF3REJ is:");
                    trlmessage (MsgText);

                    /* Need to allocate memory for a separate string to
                     ** be long enough to hold all input names when
                     ** printing it out. Causes pipeline problems
                     ** otherwise. HAB 20-Jun-2004 */

                    wf3rej_msgtext = calloc(strlen(wf3rej_cte_input)+25,
                            sizeof(char));
                    sprintf (wf3rej_msgtext, "%s", wf3rej_cte_input);
                    trlmessage (wf3rej_msgtext);
                    free (wf3rej_msgtext);
                    sprintf (MsgText,"Output to WF3REJ is: %s",
                            wf3hdr->crc_tmp);
                    trlmessage (MsgText);
                }

                /* Set up the correct value of ASN_MTYP for WF3REJ
                 ** output image. */                
                if (strncmp (wf3hdr->asn_table, wf3hdr->crc_tmp, 9) ) {

                    /* crc_tmp is a sub-product */
                    strcpy (wf3hdr->mtype,
                            asn->product[prod].subprod[posid].mtype);
                } else {

                    /* crc_tmp leads to a product with same rootname as
                     ** asn_table */
                    strcpy (wf3hdr->mtype, asn->product[prod].mtype);
                }   

                /* Reject cosmic rays. */
                if (WF3Rej_0 (wf3rej_cte_input, wf3hdr->crc_tmp,
                            wf3hdr->mtype, printtime, asn->verbose)) {
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

                if (updateAsnTable (asn, prod, posid))
                    return (status);

                /* Free up memory used by wf3rej_input for this
                 ** subproduct */
                free(wf3rej_cte_input); 


                /* MOVE ON TO WF32D */

                if (wf3hdr->sci_basic_2d == PERFORM && 
                        wf3hdr->sci_crcorr == PERFORM) {

                    if (wf3hdr->sci_basic_cte == OMIT){
                        /* Flatfield the summed, cosmic ray rejected image. */
                        if (WF32d (wf3hdr->crj_tmp, wf3hdr->crjfile,
                                    &wf32d_sci_sw, &sciref, printtime, asn->verbose))
                            return (status);
                    } else {
                        /* Flatfield the summed, cosmic ray rejected CTE image. */
                        if (WF32d (wf3hdr->crc_tmp, wf3hdr->crcfile,
                                    &wf32d_sci_sw, &sciref, printtime, asn->verbose))
                            return (status);

                    }

                } else if (wf3hdr->sci_crcorr == PERFORM) {
                    if (wf3hdr->sci_basic_cte == OMIT) {

                        /* Remember CR-combined image as final output name,
                         ** since WF32D was not run. */
                        if (CopyFFile (wf3hdr->crj_tmp, wf3hdr->crjfile))
                            return (status);

                        /* Remember what crj_tmp files to delete */
                        strcpy (asn->product[prod].subprod[posid].crj_tmp,
                                wf3hdr->crj_tmp);

                        /* If EXPSCORR=PERFORM, finish processing individual
                           exposures with Wf32d. */
                        if (sci_sw.expscorr == PERFORM) {

                            for (expid=1; expid <= asn->spmems[posid]; expid++) {

                                /* Get input & output file names for this exp */
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

                                if (wf3hdr->sci_basic_2d == PERFORM) {
                                    SetTrlPrefaceMode (NO);

                                    /* Reset switches */
                                    ResetCCDSw (&wf3ccd_sci_sw, &wf32d_sci_sw);

                                    /* Basic 2-D processing (flat field, etc). */
                                    if (WF32d (wf3hdr->blv_tmp, wf3hdr->fltfile,
                                                &wf32d_sci_sw, &sciref, printtime,
                                                asn->verbose))
                                        return (status);

                                } else {

                                    /* Save blv_tmp as final output */
                                    if (CopyFFile (wf3hdr->blv_tmp,
                                                wf3hdr->fltfile))
                                        return (status);
                                }

                            } /* End loop over indivdual exposures */

                        } /* End of EXPSCORR processing */

                        /* Delete intermediate files */
                        if (*save_tmp != YES) {
                            for (expid=1; expid <= asn->spmems[posid]; expid++) {
                                if (asn->debug) {
                                    sprintf (MsgText, "Deleting %s...",
                                            asn->product[prod].subprod[posid].exp[expid].blv_tmp);
                                    trlmessage (MsgText);
                                }
                                remove (asn->product[prod].subprod[posid].exp[expid].blv_tmp);
                            }
                        }
                    } else { /*process CTE images*/

                        /* Remember CR-combined image as final output name,
                         ** since WF32D was not run. */
                        if (CopyFFile (wf3hdr->crc_tmp, wf3hdr->crcfile))
                            return (status);

                        /* Remember what crc_tmp files to delete */
                        strcpy (asn->product[prod].subprod[posid].crc_tmp,
                                wf3hdr->crc_tmp);

                        /* If EXPSCORR=PERFORM, finish processing individual
                           exposures with Wf32d. */
                        if (sci_sw.expscorr == PERFORM) {

                            for (expid=1; expid <= asn->spmems[posid]; expid++) {

                                /* Get input & output file names for this exp */
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

                                if (wf3hdr->sci_basic_2d == PERFORM) {
                                    SetTrlPrefaceMode (NO);

                                    /* Reset switches */
                                    ResetCCDSw (&wf3ccd_sci_sw, &wf32d_sci_sw);
                                    trlmessage("Reset CCDsw");
                                    /* Basic 2-D processing (flat field, etc). */
                                    if (WF32d (wf3hdr->blc_tmp, wf3hdr->fltfile,
                                                &wf32d_sci_sw, &sciref, printtime,
                                                asn->verbose))
                                        return (status);

                                } else {

                                    /* Save blc_tmp as final output */
                                    if (CopyFFile (wf3hdr->blc_tmp,
                                                wf3hdr->fltfile))
                                        return (status);
                                }
                                trlmessage("Finished WF32D");

                            } /* End loop over indivdual exposures */

                        } /* End of EXPSCORR processing */

                        /* Delete intermediate files */
                        if (*save_tmp != YES) {
                            for (expid=1; expid <= asn->spmems[posid]; expid++) {
                                if (asn->debug) {
                                    sprintf (MsgText, "Deleting %s...",
                                            asn->product[prod].subprod[posid].exp[expid].blc_tmp);
                                    trlmessage (MsgText);
                                }
                                remove (asn->product[prod].subprod[posid].exp[expid].blc_tmp);
                            }
                        }                

                    } /*end of CTE CRCORR processing*/

                } /* End of CRCORR/RPTCORR processing */	

            }  /* End of processing for SINGLE/PARTIAL product */

            /* Done with the tmp files.... */
            if (asn->process != SINGLE) {

                /* We only have something to delete if process != SINGLE */
                if ((save_crj != YES ) && (wf3hdr->sci_crcorr == PERFORM ||
                            wf3hdr->sci_rptcorr == PERFORM)) {
                    for (posid = 1; posid <= asn->numsp; posid++) {
                        remove (asn->product[prod].subprod[posid].crc_tmp);
                    }
                }
            } /*end remove tmp files*/

        } /* End loop over CCD products here...		*/

        /* Done with lists of reference file keywords and names. */
        FreeRefFile (&sciref);
    }
        return (status);
    }

