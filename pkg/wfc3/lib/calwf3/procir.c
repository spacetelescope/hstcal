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
# include "trlbuf.h"

/* ProcessIR: This routine controls the overall flow of processing
   for IR images.

   H.Bushouse, 2008 May 16:
   Added logic and supporting functionality for calling WF3Rej_0 to
   combine Repeat-Obs exposures into a crj product.

   H.Bushouse, 2011 Nov 23:
   Modified the logic used to call wf3rej for rptcorr processing to use
   wf3ir_sci_sw.rptcorr instead of wf3hdr->sci_rptcorr, because the
   former is based on the setting of the RPTCORR header keyword.
   (PR 69952; Trac 807)

   M.Sosey, 2012 Jan 24:
   Updated to check for the number of images in wf3hdr as a check
   against observations that have RPTCORR=PERFORM, but are not in an
   association. This was causing a crash when attempting to run wf3rej
   on the non-existent association.

   H.Bushouse, 2012 Mar 21:
   Updated to set sub-product status to PRESENT after running wf3rej and
   report RPTCORR processing status via trlmessage when wf3rej is not run.
   (PR 70922; Trac #869)
   
   M. Sosey 2015 August 17:
   Reformatted because indentation was driving me nuts.

 */

int ProcessIR (AsnInfo *asn, WF3Info *wf3hdr, int printtime) {

    extern int status;

    RefFileInfo sciref;		/* ref file keywords and names */
    IR_Switch sci_sw;		/* all cal switches for science file */
    IR_Switch wf3ir_sci_sw;		/* WF3IR switches for science file */

    int prod;			/* which CR-split/Rptobs product 
                                 ** are we working with?	*/
    int posid;			/* counter for input products	*/
    int expid;			/* counter for input exposures */
    int newpreface = NO;		/* switch for keeping previous comments
                                         ** for all remaining trailer files */

    char *wf3rej_input;     /* list of input names for WF3REJ */
    char *wf3rej_msgtext;   /* string for list of input filenames */
    int  nchars;            /* Number of chars for input string to WF3REJ */

    /* CALWF3 function definitions */
    int  WF3Rej_0    (char *, char *, char *, int, int,int);
    int  CopyFFile   (char *, char *);
    void SetIRSw     (IR_Switch *, IR_Switch *);
    void WF3Defaults (WF3Info *);
    void InitRefFile (RefFileInfo *);
    void FreeRefFile (RefFileInfo *);
    int  IRRefInit   (WF3Info *, IR_Switch *, RefFileInfo *);
    int  GetAsnMember(AsnInfo *, int, int, int, WF3Info *);
    int  GetSingle   (AsnInfo *, WF3Info *);
    int  InsertWF3Suffix (WF3Info *);
    int  WF3ir (char *, char *, IR_Switch *, RefFileInfo *, int, int,
            float, float, int);
    int  updateAsnTable (AsnInfo *, int, int);
    void PrSwitch (char *, int);

    wf3rej_input=NULL;
    wf3rej_msgtext=NULL;
    
    /* Loop over the products/positions for each
     ** Repeat-Obs or Repeat-Obs/DITHER set.
     ** PRODID/prod starts at 0, while numprod starts at 1... */
    for (prod = 0; prod < asn->numprod; prod++) {
        if (asn->verbose) {
            sprintf (MsgText, "CALWF3: processing IR product %d", prod);
            trlmessage (MsgText);
        }

        /* Loop over members of a position/product;
         ** individual Repeat-obs exposures */
        /* Process SINGLE/PARTIAL/FULL product */

        for (posid=1; posid <= asn->numsp; posid++) {
            if (asn->verbose) {
                sprintf (MsgText,
                        "CALWF3: processing posid = %d", posid);
                trlmessage (MsgText);
            }

            /*  Allocate space for WF3REJ input image list */
            if (asn->rptcorr == PERFORM) {
                nchars = asn->spmems[posid] * (CHAR_FNAME_LENGTH+1);
                wf3rej_input = (char *) calloc( nchars + 1, sizeof(char));

                /* Initialize this string to NULL */
                wf3rej_input[0] = '\0';
            }
            for (expid=1; expid <= asn->spmems[posid]; expid++) { 

                /* From this point on, we are working with ONE image
                 ** at a time... */

                /* Read in Member keyword info from --FIRST-- EXP
                 ** image header */	
                if (prod == 0 && posid == 1 && expid == 1) {

                    /* (Re-)Initialize the lists of reference file
                     ** keywords and names */
                    InitRefFile (&sciref);

                    /* Assign default values in wf3. */
                    WF3Defaults (wf3hdr);

                    if (asn->process == SINGLE) {
                        if (GetSingle (asn, wf3hdr))
                            return (status);
                    } else {	
                        if (GetAsnMember (asn, prod, posid, expid,
                                    wf3hdr)) 
                            return (status);
                    } 				

                    /* Initialize variables, and get info (calibration
                     ** switches and reference file names) from primary
                     ** headers. */
                    if (IRRefInit (wf3hdr, &sci_sw, &sciref))
                        return (status);

                    if (asn->verbose) {
                        trlmessage 
                            ("CALWF3: Got reference file information");
                    }

                    if (sci_sw.rptcorr == PERFORM && wf3hdr->nimages < 2) {
                        trlmessage
                            ("RPTCORR will be omitted because there's only one image.");
                        sci_sw.rptcorr = OMIT;
                    }

                    /* Store the trailer file comments into preface */
                    newpreface = YES;

                } else {

                    /* Read in Member keyword info from image header */
                    if (asn->process == SINGLE) {
                        if (GetSingle (asn, wf3hdr))
                            return (status);
                    } else {	
                        if (GetAsnMember (asn, prod, posid, expid,
                                    wf3hdr)) 
                            return (status);
                    } 

                    /* Construct output and temporary file names. */
                    if (InsertWF3Suffix (wf3hdr))
                        return (status);	    		
                }

                /* Copy switch values for WF3IR argument list */
                SetIRSw (&sci_sw, &wf3ir_sci_sw);

                if (asn->verbose) {
                    trlmessage ("CALWF3: Processing switches set... ");
                }

                /* Ready to process input image now...	*/
                /* Do we have anything to do with this image?  */
                if (wf3hdr->sci_basic_ir != PERFORM) {
                    trlerror("No calibration switch was set to PERFORM.");
                    FreeRefFile (&sciref);
                    return (status = NOTHING_TO_DO);
                }

                if (wf3hdr->sci_basic_ir == PERFORM) {

                    /* Copy all messages up to here into preface buffer
                     ** to be prepended to all input file trailer files */
                    if (newpreface == YES) {
                        InitTrlPreface();
                    }

                    SetTrlPrefaceMode (YES);

                    /* Apply the calibration steps */
                    if (WF3ir (wf3hdr->rawfile, wf3hdr->fltfile,
                                &wf3ir_sci_sw, &sciref, printtime,
                                asn->verbose, 0.0, 0.0, 0))
                        return (status);

                } else {

                    /* Copy the input file to fltfile. We need to do
                     ** this because WF3SUM needs it for input. */
                    trlwarn ("No processing performed on image.");
                    if (wf3hdr->sci_rptcorr != PERFORM)
                        status = NOTHING_TO_DO;
                    sprintf (MsgText, "Copying input to %s ",
                            wf3hdr->fltfile);
                    trlwarn (MsgText);

                    if (CopyFFile (wf3hdr->rawfile, wf3hdr->fltfile))
                        return (status);
                }

                /* We now are working with flt here... */
                /* Copy (cat) flt name to wf3rej_input string */
                if (wf3hdr->sci_rptcorr == PERFORM) {
                    strcat (wf3rej_input, wf3hdr->fltfile);
                    if (expid < asn->spmems[posid]) {
                        /* Don't add a comma to the end of the last
                         ** filename */
                        strcat(wf3rej_input, ",");
                    }
                }
            }    /* Finished processing EXP images */

            /* Do RPTCORR processing */
            /* if (wf3hdr->sci_rptcorr == PERFORM) { */
            if (wf3ir_sci_sw.rptcorr == PERFORM) {

                if (asn->verbose) {
                    sprintf (MsgText,
                            "CALWF3: Now process position %d from product %d",
                            posid, prod);
                    trlmessage (MsgText);
                }

                if (asn->verbose) {
                    sprintf (MsgText,"Input to WF3REJ is:");
                    trlmessage (MsgText);

                    /* Need to allocate memory for a separate string to
                     ** be long enough to hold all input names when
                     ** printing it out. Causes pipeline problems
                     ** otherwise. */
                    wf3rej_msgtext = calloc(strlen(wf3rej_input)+25,
                            sizeof(char));
                    sprintf (wf3rej_msgtext, "%s", wf3rej_input);
                    trlmessage (wf3rej_msgtext);
                    free (wf3rej_msgtext);
                    sprintf (MsgText,"Output to WF3REJ is: %s",
                            wf3hdr->crjfile);
                    trlmessage (MsgText);
                }

                /* Set up the correct value of ASN_MTYP for WF3REJ
                 ** output image. */
                if (strncmp (wf3hdr->asn_table, wf3hdr->crjfile, 9) ) {
                    /* crjfile is a sub-product */
                    strcpy (wf3hdr->mtype,
                            asn->product[prod].subprod[posid].mtype);
                } else {
                    /* crjfile leads to a product with same rootname as
                     ** asn_table */
                    strcpy (wf3hdr->mtype, asn->product[prod].mtype);
                }

                /* Reject cosmic rays. */
                if (WF3Rej_0 (wf3rej_input, wf3hdr->crjfile,
                            wf3hdr->mtype, printtime, asn->verbose,1)) {
                    if (status == NO_GOOD_DATA) {
                        /* Reset status to good now that we have dealt
                         ** with this condition */
                        status = WF3_OK;
                    } else {
                        return (status);
                    }
                }

                asn->product[prod].subprod[posid].prsnt = True;
                if (updateAsnTable (asn, prod, posid))
                    return (status);

                /* Free up memory used by wf3rej_input for this
                 ** subproduct */
                free(wf3rej_input);

            } else {
                if (asn->process == FULL) {
                    asn->product[prod].subprod[posid].prsnt = False;
                    trlmessage ("");
                    PrSwitch ("rptcorr", wf3ir_sci_sw.rptcorr);
                }

            } /* End RPTCORR processing */

        }    /* Finished processing this position's sub-product */
        }    /* Finished looping over all the products */

        /* Reset the trailer file preface to NULL since
         ** it has already been copied into trailer files */
        ResetTrlPreface();

        /* Done with lists of reference file keywords and names. */
        FreeRefFile (&sciref);

        return (status);
    }

