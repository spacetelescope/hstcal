# include <stdio.h>
# include <stdlib.h>	/* calloc */
# include <string.h>
# include "xtables.h"   /* for IRAFPointer definition */

# include "hstio.h"

# include "acs.h"
# include "calacs.h"
# include "acserr.h"
# include "acscorr.h"
# include "acsasn.h"	/* Contains association table structures */

# define NOPOSID 0


/* calacs -- integrated calacs processing

   Warren J. Hack, 1998 May 28:
	Initial version.
	
   Warren J. Hack, 1999 Jan 7:
   	Added support for EXPSCORR to control 'save_tmp', if set to PERFORM.
   Warren J. Hack, 2002 Feb 1:
    Changed so that when no calibration switches are set to PERFORM, it 
        only WARNs the user, not provide ERROR messages.
    WJH, 2002 Apr 17:   
        Removed all references to 'statcorr'. Perform 'doStat' all the time.
    WJH, 2002 Apr 18:
        Added check against 'RPTCORR' to be consistent with 'CRCORR' check
            for determining whether any calibration switches are set to 
            'PERFORM'.  If only RPTCORR is on, then it will not stop but 
            copy raw files and sum them with ACSSUM.	
    WJH, 2002 Apr 23:
        Update DQICORR for CRJ product so that doDQI is not run again in ACS2D.
    WJH, 2002 July 26:
        Corrected memory usage error in BuildSumInput. 
    WJH, 2002 Aug 19:
        Corrected a bug where RPT-OBS sub-products MTYPE was getting overwritten
        by the PROD-DTH MTYPE in dithered RPT-OBS associations.
        Fix involved using subprod.mtype instead of prod.mtype.     
    WJH, 2005 Feb 14:
        Removed any reference to updateAsnStat, since only OPUS should
        update the ASN_STAT keyword in the ASN table header. 
*/


static int ACSRej_0 (char *, char *, char *, int, int, int);


static int CopyFFile (char *, char *);
static void SetACSSw (CalSwitch *, CalSwitch *, CalSwitch *);
static void ResetACSSw (CalSwitch *, CalSwitch *);

int CalAcsRun (char *input, int printtime, int save_tmp, int verbose, int debug) {

/* arguments:
char *input	 i: name of the FITS file/table to be processed
int printtime    i: true --> print time stamps at intermediate steps
int *save_tmp     i: true --> save temporary files
int verbose      i: true --> print info during processing
int debug 	     i: true --> print debugging info during processing
*/


	extern int status;

	ACSInfo acshdr;		/* calibration switches, etc */
	AsnInfo	asn;		/* association table data    */
	char *acssum_input; /* Input list for ACSSUM */
    char acssum_output[ACS_FNAME+1];
    char acssum_mtype[SZ_STRKWVAL+1];      /* Role of exposure in association */
	char *acsdth_input; /* Input list for ACSDTH */
    char *acssum_msgtext;
	
    IRAFPointer tpin; 
    int nimgs, i;
    char sumfile[ACS_LINE];

	int prod, posid;
	
	void PrBegin (char *);
	void PrEnd (char *);
	void PrFileName (char *, char *);
	void TimeStamp (char *, char *);
/* Association table routines			*/
	void initAsnInfo (AsnInfo *);
	void freeAsnInfo (AsnInfo *);
	int LoadAsn (AsnInfo *);
	int ProcessCCD (AsnInfo *, ACSInfo *, int *, int);
	int ProcessMAMA (AsnInfo *, ACSInfo *, int);	
	int AcsSum (char *, char *, char *, int, int);
	int AcsDth (char *, char *, int, int, int);
	char *BuildSumInput (AsnInfo *, int, int);
	char *BuildDthInput (AsnInfo *, int);
	void updateAsnTable (AsnInfo *, int, int);
	void InitDthTrl (char *, char *); 
	
    /* Post error handler */
    push_hstioerr (errchk);
	    
	PrBegin ("CALACS");	/* *** CALACS -- Version ... *** */

	if (printtime)
	    TimeStamp ("CALACS started", "");

	/* Determine if input is a single file or an association
	**	table, then populate ASN structure with appropriate
	**	information to control processing.		
	*/
	initAsnInfo(&asn);
			
	if (debug) {
		trlmessage ("Initialized Association data ... ");
	}
	
	/* Copy Input filename to ASN structure	*/
	strcpy (asn.input, input);
	
	/* Print image name. */
	trlmessage ("\n");
	
	PrFileName ("input", asn.input);

	/* Set verbose flag... */
	asn.verbose = verbose;
	asn.debug = debug;
		
	/* LoadAsn will determine whether input is a single file or an 
	**	Association table.  If a single image, it will look in that images
	**	association table to see what products are associated with it, and
	**	process them accordingly.  If it is just a single image as its
	**	own output, it will proceed as a 1 element table.	
	**	Based on routines from n_getAsnTable() and n_setup() in CALNICB		
	*/
	if(LoadAsn(&asn)) {
		freeAsnInfo (&asn);
		return (status);
	}

    /* Check to see that detector is known, as it could come through with a
        value of 0 or UNKNOWN_DETECTOR.  
        WJH  2Mar99
    */
	if (asn.detector == UNKNOWN_DETECTOR || asn.detector == 0) {
		trlwarn ("Unknown detector type for observations.");
		freeAsnInfo(&asn);	
    	return (status = NOTHING_TO_DO);
	}

	if (asn.verbose) {
		sprintf (MsgText,"CALACS: Detector %s, type %d ",asn.instr, asn.detector);
		trlmessage (MsgText);
	}

/* 	Determine what detector we are working with...	*/
	if (asn.detector != MAMA_DETECTOR ) { /* Process CCD data ... */
		if (asn.verbose) {
			trlmessage ("CALACS: processing a CCD product");
		}
		if (ProcessCCD(&asn, &acshdr, &save_tmp, printtime)) { 
			if (status == NOTHING_TO_DO) {
                trlwarn ("No processing desired for CCD data.");            
            } else {
                trlerror ("Couldn't process CCD data");
            }

			freeAsnInfo(&asn);
			return (status);			
		}	

	} else { /* Process MAMA observations here */

        trlmessage("Starting to process MAMA data now...");
        
		if (ProcessMAMA(&asn, &acshdr, printtime)) { 
			if (status == NOTHING_TO_DO){
                trlwarn ("No processing desired for MAMA data.");            
            } else{
                trlerror ("Couldn't process MAMA data");
            }
			freeAsnInfo(&asn);
			return (status);			
		}
	}
    trlmessage("Finished MAMA processing...");
    
	/* 
	** Perform RPTCORR here, as needed.
	*/
	if (asn.rptcorr == PERFORM) {
/*		We need to pass a list of input fltfile images for each subproduct
		to be produced, with the subproducts as inputs to DTHCORR...
*/	
		acssum_input = NULL;
        acssum_mtype[0] = '\0';
		
        trlmessage("Starting RPTCORR now...");
        
		/* For each DTH product... */
		for (prod = 0; prod < asn.numprod; prod++) {
			/* Loop through all the input subproducts... */
			for (posid = 1; posid <= asn.numsp; posid++) {
                if (asn.spmems[posid] == 0) {
                    continue;
                } 
				acssum_input = BuildSumInput (&asn, prod, posid);
                
                strcpy(acssum_output, asn.product[prod].subprod[posid].spname);
				if (asn.verbose) {
                    /* 
                        Need to allocate memory for a separate string to be long
                        enough to hold all the input names when printing it
                        out.  Caused pipeline problems otherwise. WJH 19-Mar-04
                    */
                    acssum_msgtext = calloc(strlen(acssum_input)+25, sizeof(char));
                    sprintf (acssum_msgtext, "Input to ASCSUM is: %s",acssum_input);
                    trlmessage(acssum_msgtext);
                    free(acssum_msgtext);
                    
					sprintf (MsgText, "Output for ASCSUM is: %s",acssum_output);
					trlmessage(MsgText);
				}
                            
                strcpy(acssum_mtype, asn.product[prod].subprod[posid].mtype);
                if (acssum_mtype[0] == '\0'){
		            strcpy(acssum_mtype, asn.product[prod].mtype);
                }
                trlmessage("ACSSUM starting now...");

	    		if (AcsSum (acssum_input, acssum_output, acssum_mtype, printtime, verbose)) {
                    return (status);
                }

				/* Pass posid of 0 to indicate a PRODUCT is to be updated */
				updateAsnTable(&asn, prod, posid);
			}
		}
		free (acssum_input);
	} /* End RPTCORR Processing */
	
	/* Add DTH processing here... */
		/* For each DTH product... */
        if (asn.process == FULL){
	    if (asn.verbose) {
		    trlmessage ("CALACS: Building DTH products");
	    }
		acsdth_input = NULL;
		for (prod = 0; prod < asn.numprod; prod++) {
			/* Create empty DTH product, header only */
			/* This uses only one sub-product for the header template,
			** 	but later versions should use a function similar to
			** 	BuildSumInput to create list of subproducts as inputs...
			*/            
			acsdth_input = BuildDthInput (&asn, prod);

            /* 
            We always want to create a final concatenated trailer file for
            the entire association whether there is a product or not. So, we
            set up the trailer file based on the association file name itself. 
            */
            
            /* If desired, we could optionally use the full _drz.tra filename
                as the trailer file name, based on the output dither product name...
            if (strcmp(asn.product[prod].prodname,"") != 0) {
    	        InitDthTrl (acsdth_input, asn.product[prod].prodname);        
            } else { */
            
            InitDthTrl(acsdth_input, asn.rootname);
            
            /* End brace for optional dither product name assignment...
             } */
             
			/* Check if we have a PROD-DTH specified...*/
			if (strcmp(asn.product[prod].prodname,"") != 0) {

            	if ((asn.dthcorr == PERFORM || asn.dthcorr == DUMMY)) {
    				if (AcsDth (acsdth_input, asn.product[prod].prodname, asn.dthcorr, printtime, asn.verbose) )
					    return (status);

				    /* Pass posid of 0 to indicate a PRODUCT is to be updated */
				    updateAsnTable(&asn, prod, NOPOSID);	
            	}

			} else {
				trlwarn ("No DTH product name specified. No product created.");
			    /*status = ACS_OK; */
                
            }
		}
		free (acsdth_input);
	}
	if (asn.verbose) {
		trlmessage ("CALACS: Finished processing product ");
	}
	
	freeAsnInfo(&asn);
	
	trlmessage ("\n");
	PrEnd ("CALACS");		/* *** CALACS complete *** */

	if (printtime)
	    TimeStamp ("CALACS completed", acshdr.rootname);
	
    /* Return the final value of status, which should be ACS_OK if
        all went well or some error condition otherwise. */
	return (status);
}



char *BuildSumInput (AsnInfo *asn, int prod, int posid) {

	int nchars;
    int acssum_len;
	int i;
	char *acssum_input;
	char tmpexp[ACS_LINE];
	char tmpflt[ACS_LINE];
	int MkName (char *, char *, char *, char *, char *, int);

	/* Determine how long this string needs to be... */
	/*nchars = asn->spmems[posid] * ACS_FNAME; */
    nchars = 1;
    
    /* Keep track of individual filename lengths and total length*/
    acssum_len = 0;
	/* Now, lets search the association table for all inputs... */
	for (i=1; i <= asn->spmems[posid]; i++) {
		strcpy(tmpexp, asn->product[prod].subprod[posid].exp[i].expname);
                
		if (MkName (tmpexp, "_raw", "_flt", "", tmpflt, ACS_LINE)) {
			strcpy (tmpflt,asn->product[prod].subprod[posid].exp[i].name);  
	    	strcat (tmpflt, "_flt.fits");
		}
		/* Sum together lengths of all filenames 
            to be used as input to ACSSUM.
        */
        acssum_len += strlen(tmpflt)+1;
	}

    /* Now that we know how long it will be, 
        allocate space for input string. 
    */
	acssum_input = (char *) calloc(acssum_len + 1, sizeof(char));

	/* Go back through the association table to string all inputs together... */
	for (i=1; i <= asn->spmems[posid]; i++) {
		strcpy(tmpexp, asn->product[prod].subprod[posid].exp[i].expname);

		if (MkName (tmpexp, "_raw", "_flt", "", tmpflt, ACS_LINE)) {
			strcpy (tmpflt,asn->product[prod].subprod[posid].exp[i].name);  
	    	strcat (tmpflt, "_flt.fits");
		}

		strcat(acssum_input, tmpflt);
		if (i < (asn->spmems[posid])) {
			/* Don't add a comma to the end of the last filename*/
			strcat(acssum_input, ",");						
		}
	}
    strcat(acssum_input,"\0");
    
	return(acssum_input);
}
char *BuildDthInput (AsnInfo *asn, int prod) {

	int nchars;
	int i;
	char *acsdth_input;
	char tmpexp[ACS_LINE];
	char tmpflt[ACS_LINE];
	int MkName (char *, char *, char *, char *, char *, int);

	/* Determine how long this string needs to be... */
	nchars = asn->numsp * ACS_FNAME;
	acsdth_input = (char *) calloc( nchars + 1, sizeof(char));
	/* Initialize this string to NULL */
	acsdth_input[0] = '\0';

	/* Now, lets search the association table for all inputs... */
	for (i=1; i <= asn->numsp; i++) {
		strcpy(tmpexp, asn->product[prod].subprod[i].spname);

		/*if (MkName (tmpexp, "_crj_tmp", "_crj", "", tmpflt, ACS_LINE)) {
			strcpy (tmpflt,asn->product[prod].subprod[posid].exp[i].name);  
	    	strcat (tmpflt, "_crj.fits");
		}
		*/
		strcat(acsdth_input, tmpexp);
		if (i < (asn->numsp)) {
			/* Don't add a comma to the end of the last filename*/
			strcat(acsdth_input, ",");						
		}
	}
	return(acsdth_input);
}

int ProcessCCD (AsnInfo *asn, ACSInfo *acshdr, int *save_tmp, int printtime) {

	extern int status; 
	
	RefFileInfo sciref;			/* ref file keywords and names */
	CalSwitch sci_sw;			/* all cal switches for science file */
	CalSwitch acsccd_sci_sw;  	/* ACSCCD switches for science file */
	CalSwitch acs2d_sci_sw;   	/* ACS2D  switches for science file */
	int save_crj;				/* don't delete crj_tmp file? */

	int prod;					/* which CR-split/Rptobs product 
									are we working with?	*/
	int posid;					/* counter for input products	*/
	int expid;					/* counter for input exposures 	*/
	int newpreface = NO;		/* switch for keeping previous comments for 
									all remaining trailer files */
	
	char *acsrej_input; /* list of input names for ACSREJ */
    char *acsrej_msgtext;       /* string for list of input filenames */
	int	nchars;			/* Number of chars for input string to ACSREJ */
	
	void ACSDefaults (ACSInfo *);	
	void InitRefFile (RefFileInfo *);
	void FreeRefFile (RefFileInfo *);
	int ACSRefInit (ACSInfo *, CalSwitch *, RefFileInfo *);
	int ACSccd (char *, char *, CalSwitch *, RefFileInfo *, int, int);
	int ACS2d (char *, char *,CalSwitch *, RefFileInfo *, int, int);
	int GetAsnMember (AsnInfo *, int, int, int, ACSInfo *);
	int GetSingle (AsnInfo *, ACSInfo *);
	int CheckCorr (AsnInfo *, ACSInfo *);
	int InsertACSSuffix (ACSInfo *);
	void updateAsnTable (AsnInfo *, int, int);
    int UpdateSwitch (char *, int, Hdr *, int *);
    	
	save_crj = *save_tmp;		/* initial value; */	
	
	/* 	Loop over the products/positions for each CR-split/Repeat-Obs
		set.							*/
	for (prod = 0; prod < asn->numprod; prod++) {
		if (asn->verbose){
			sprintf (MsgText, "CALACS: processing CCD product %d", prod);	
			trlmessage (MsgText);
		}	
		/* 
		**	Process this PARTIAL/SINGLE/FULL product...
		*/
		for (posid = 1; posid <= asn->numsp; posid++) {

			if (asn->verbose) {
				sprintf (MsgText,"CALACS: processing posid = %d", posid);
				trlmessage (MsgText);
			}
			/*  Allocate space for ACSREJ input image list */
			if (asn->crcorr == PERFORM) {
				nchars = asn->spmems[posid] * ACS_FNAME;
				acsrej_input = (char *) calloc( nchars + 1, sizeof(char));
				/* Initialize this string to NULL */
				acsrej_input[0] = '\0';
			}
			
			for (expid = 1; expid <= asn->spmems[posid]; expid++) {
			/* From this point on, we are working with ONE image at a time... */
                
				/* Read in Member keyword info from --FIRST-- EXP image header */	
				/*if(prod == 0 && posid == 1 && expid == 1) { */
				if(expid == 1) {
					/* (Re-)Initialize the lists of reference file keywords and names. */
					InitRefFile (&sciref);

					/* Assign default values in acs. */
					ACSDefaults (acshdr);

					if (asn->process == SINGLE ) {
						if (GetSingle (asn, acshdr) )
							return (status);
					} else {	
						if( GetAsnMember(asn, prod, posid, expid, acshdr)) 
							return (status);
					} 
					/* Initialize variables, and get info (calibration switches and
					** reference file names) from primary headers.
					*/
					if (ACSRefInit (acshdr, &sci_sw, &sciref))
	    				return (status);

					/* Make sure 'save_tmp' and 'save_crj' are no 
						longer set to 'DUMMY'. 
						
						If save_tmp was set from the command-line, then 
						ignore what is set in the header keyword 
						EXPSCORR.
					*/
					if (*save_tmp == DUMMY) {
						*save_tmp = NO;
						save_crj = *save_tmp;
					}

					if (asn->verbose) {
						trlmessage ("CALACS: Got reference file information");
					}
					/* Store the trailer file comments into preface */
					newpreface = YES;

				} else {

					if (asn->process == SINGLE ) {
						if (GetSingle (asn, acshdr) )
							return (status);
					} else {
						if( GetAsnMember(asn, prod, posid, expid, acshdr)) 
							return (status);                        
					} 
					
					/* Construct output and temporary file names. */
						if (InsertACSSuffix (acshdr))
	    					return (status);	    									
				}

                /* Check ASN table settings for CRCORR/RPTCORR 
                    against image header values.  If they are not
                    consistent, print out a ERROR message and quit.
                    Only need to perform this check on the first image
                    in the association. WJH 21 May 1999
                */
				if(prod == 0 && posid == 1 && expid == 1) {
                    if (CheckCorr(asn, acshdr))
                        return (status);
                }                
		/* Copy switch values for ACSCCD and ACS2D argument lists. */
	    		SetACSSw (&sci_sw, &acsccd_sci_sw, &acs2d_sci_sw);
                                
				if (asn->verbose) {
						trlmessage ("CALACS: Processing switches set... ");
			    }
				/* Ready to process input image now...				*/
				/* Do we have anything to do with this image? 	*/
	    		if (acshdr->sci_basic_2d != PERFORM &&
	    			acshdr->sci_basic_ccd != PERFORM &&
	    			acshdr->sci_crcorr != PERFORM &&
                    acshdr->sci_rptcorr != PERFORM ){

	    			trlwarn ("No calibration switch was set to PERFORM.");
					FreeRefFile (&sciref);
					return (status = NOTHING_TO_DO);
					
			 	}

				/* First, intialilize the Input CCD image with ACSCCD.  Then, 
				** either do cosmic ray rejection followed by ACS2D, or just
				** run ACS2D.  
				*/
	    		
				/* First do atodcorr, dqicorr, and blevcorr
				** (or copy the name). 
				*/
	    		if (acshdr->sci_basic_ccd == PERFORM) {
				/* This is calacsccd; note that we use acsccd_sci_sw.
				*/
					/* Copy all messages up to here into preface buffer
						to be prepended to all input file trailer files */
					if (newpreface == YES) {
						InitTrlPreface();
					}

					SetTrlPrefaceMode (YES);
                    
		    		if (ACSccd (acshdr->rawfile, acshdr->blv_tmp,
					&acsccd_sci_sw, &sciref, printtime, asn->verbose)) 		
						return (status);
											
					/* Reset switches */
		    		ResetACSSw (&acsccd_sci_sw, &acs2d_sci_sw);

	    		} else {
				/* Copy the input file to blv_tmp.  We need to do this
				** because acsrej may (depending on the crmask column
				** in the crrejtab) modify the DQ flags in its input file.
				*/
			    	if (CopyFFile (acshdr->rawfile, acshdr->blv_tmp))
			    		return (status);
	        	}
				/* We now are working with blv_tmp here... */
				/* Copy (cat) blv_tmp name to acsrej_input string */
				if (acshdr->sci_crcorr == PERFORM) {
					strcat(acsrej_input, acshdr->blv_tmp);
					if (expid < asn->spmems[posid]) {
						/* Don't add a comma to the end of the last filename*/
						strcat(acsrej_input, ",");					
					}
					/* Also, add blv_tmp to list of files to be deleted */
					strcpy(asn->product[prod].subprod[posid].exp[expid].blv_tmp,
					acshdr->blv_tmp);
				}
				
				if (asn->debug){
					trlmessage("Run ACS2D for CCD next...");
				}
				
				/* If we are not performing CRCORR or
                **    EXPSCORR is set to PERFORM, then finish
				**	processing individual EXP files with ACS2D...
				*/ 
				if (acshdr->sci_crcorr != PERFORM || sci_sw.expscorr == PERFORM) { 
					if (asn->debug) {
						trlmessage("CRCORR not set to Perform or EXPSCORR turned on...");
					}
					
                    if (acshdr->sci_basic_2d == PERFORM) {

					    SetTrlPrefaceMode (NO);
                        if (asn->copy_input == PERFORM) {
                            if (CopyFFile (acshdr->fltfile, acshdr->blv_tmp))
	    		    		    return (status);
                            remove (acshdr->fltfile);
                        }
					
	    			    /* Basic 2-D processing (flat field, etc). */
	        			if (ACS2d (acshdr->blv_tmp, acshdr->fltfile, 
							&acs2d_sci_sw, &sciref, printtime, asn->verbose)) 
			  					return (status);				
					} else {
                        /* If we are not performing CRCORR, then copy
                            blv_tmp file into a final _flt file.
                        */
                        if (acshdr->sci_crcorr != PERFORM) {
    			    	    if (CopyFFile (acshdr->blv_tmp, acshdr->fltfile))
	    		    		    return (status);
                        }
	      			}
				} /* Finished processing blv_tmp file through ACS2D */
                
				/* Now, delete _blv_tmp files ONLY if not doing CRCORR */
				if ( (*save_tmp != YES) && (acshdr->sci_crcorr != PERFORM) )
		    		remove (acshdr->blv_tmp);
					
			}    	/* End loop over individual EXPosures, 
			 	 	**	go on to create products			*/
			
			/* Reset the trailer file preface to NULL since
				it has already been copied into trailer files */
			ResetTrlPreface();

			if (acshdr->sci_crcorr == PERFORM) {
			/* Cosmic ray rejection, followed by basic 2-D processing.		*/
				if (asn->verbose) {
					sprintf (MsgText,"CALACS: Now process position %d from product %d",
						posid, prod);
					trlmessage (MsgText);
				}

				/* Copy switch values for ACSCCD and ACS2D argument lists. */
	    		SetACSSw (&sci_sw, &acsccd_sci_sw, &acs2d_sci_sw);	
                
                /* Reset DQICORR, since it will already have been done...
                    WJH  23-Apr-2002
                */
                acs2d_sci_sw.dqicorr = COMPLETE;
                
				if (asn->debug){
					sprintf(MsgText,"Input to ACSREJ is:");
					trlmessage(MsgText);
                    /* 
                        Need to allocate memory for a separate string to be long
                        enough to hold all the input names when printing it
                        out.  Caused pipeline problems otherwise. WJH 19-Mar-04
                    */
                    acsrej_msgtext = calloc(strlen(acsrej_input)+25, sizeof(char));
                    sprintf (acsrej_msgtext, "%s",acsrej_input);
                    trlmessage(acsrej_msgtext);
                    free(acsrej_msgtext);

					sprintf(MsgText,"Output to ACSREJ is: %s", acshdr->crj_tmp);
					trlmessage(MsgText);
				}

                /* Set up the correct value of ASN_MTYP for ACSREJ output image. */                
                if ( strncmp(acshdr->asn_table, acshdr->crj_tmp, 9) ) {
                    /* crj_tmp is a sub-product */
                    strcpy (acshdr->mtype, asn->product[prod].subprod[posid].mtype);
                } else {
                    /* crj_tmp leads to a product with same rootname as asn_table */
                    strcpy (acshdr->mtype, asn->product[prod].mtype);
                }   
                
	    		/* Reject cosmic rays. */
	    	    if (ACSRej_0 (acsrej_input, acshdr->crj_tmp, acshdr->mtype, acshdr->newbias, printtime, asn->verbose)){ 
                    if (status == NO_GOOD_DATA){
                        /* Turn off further processing... */
                        acshdr->sci_basic_2d = SKIPPED;
                        /* Reset STATUS to good now that we have dealt with
                            this condition. */
                        status = ACS_OK;
                    } else {
                	    return (status); 
                    }
                }

				updateAsnTable (asn, prod, posid);
				
				/* Free up memory used by acsrej_input for this subproduct */
				free(acsrej_input); 
				
				if (*save_tmp != YES) {
					for (expid = 1; expid <= asn->spmems[posid]; expid++) {
						if (asn->debug) {
							sprintf(MsgText,"Deleting %s...", asn->product[prod].subprod[posid].exp[expid].blv_tmp);
							trlmessage(MsgText);					
						}
						remove (asn->product[prod].subprod[posid].exp[expid].blv_tmp);
					}
				}
		
	    		if (acshdr->sci_basic_2d == PERFORM) {
					/* Flat field the summed, cosmic ray rejected image. */
		    		if (ACS2d (acshdr->crj_tmp, acshdr->crjfile,
					&acs2d_sci_sw, &sciref, printtime, asn->verbose))
						return (status);							

	    		} else {

					/* Remember CR-combined image as final output name, since 
					**	ACS2D was not run.
					*/
			    	if (CopyFFile (acshdr->crj_tmp, acshdr->crjfile))
			    		return (status);
	    			/*	strcpy (acshdr->crj_tmp, acshdr->crjfile);
		    		save_crj = YES;                         */
	        	}
			
				/* Remember what crj_tmp files to delete */
				strcpy (asn->product[prod].subprod[posid].crj_tmp, acshdr->crj_tmp);
					
			} /* End of CRCORR processing */	
				  			  
		}  /* End of processing for SINGLE/PARTIAL product */
		
		/* Done with the tmp files.... */
		if (asn->process != SINGLE) {
			/* We only have something to delete if process != SINGLE */
			if ( (save_crj != YES ) && acshdr->sci_crcorr == PERFORM) {
				for (posid = 1; posid <= asn->numsp; posid++) {
					remove (asn->product[prod].subprod[posid].crj_tmp);
				}
			}
		}
		
	} 	/* End loop over CCD products here...		*/


	/* Done with lists of reference file keywords and names. */
	FreeRefFile (&sciref);
	
	return (status);
}


int ProcessMAMA (AsnInfo *asn, ACSInfo *acshdr, int printtime) {

	extern int status; 
	
	RefFileInfo sciref;			/* ref file keywords and names */
	CalSwitch sci_sw;			/* all cal switches for science file */
	CalSwitch acsccd_sci_sw;  	/* ACSCCD switches for science file */
	CalSwitch acs2d_sci_sw;   	/* ACS2D  switches for science file */

	int prod;					/* which CR-split/Rptobs product 
									are we working with?	*/
	int posid;					/* counter for input products	*/
	int expid;					/* counter for input exposures */
	int newpreface = NO;	    /* switch for keeping previous comments for 
									all remaining trailer files */
	
	void ACSDefaults (ACSInfo *);	
	void InitRefFile (RefFileInfo *);
	void FreeRefFile (RefFileInfo *);
	int ACSRefInit (ACSInfo *, CalSwitch *, RefFileInfo *);
	int ACS2d (char *, char *, CalSwitch *, RefFileInfo *, int, int);
	int GetAsnMember (AsnInfo *, int, int, int, ACSInfo *);
	int GetSingle (AsnInfo *, ACSInfo *);
	int InsertACSSuffix (ACSInfo *);

	if (asn->debug) {
		trlmessage ("CALACS: processing MAMA observations");
	}

	/* 	Loop over the products/positions for each
	**	Repeat-Obs or Repeat-Obs/DITHER set.
	**	PRODID/prod starts at 0, while numprod starts at 1...			
	*/
	for (prod = 0; prod < asn->numprod; prod++) {

	/* loop over members of a position/product;
		individual CR-SPLIT/Repeat-obs exposures	
	*/
	  /* Process SINGLE/PARTIAL/FULL product		*/
	  
		for (posid = 1; posid <= asn->numsp; posid++) {
			if (asn->verbose){
				sprintf (MsgText,"CALACS: processing MAMA posid = %d", posid);				  
				trlmessage (MsgText);
			}

			for (expid = 1; expid <= asn->spmems[posid]; expid++) { 

			/* Read in Member keyword info from --FIRST-- EXP image header */	
			if(prod == 0 && posid == 1 && expid == 1) {
				/* (Re-)Initialize the lists of reference file keywords and names. */
				InitRefFile (&sciref);

				/* Assign default values in acs. */
				ACSDefaults (acshdr);

				if (asn->process == SINGLE ) {
					if (GetSingle (asn, acshdr) )
						return (status);
				} else {	
					if( GetAsnMember(asn, prod, posid, expid, acshdr)) 
						return (status);
				} 				

				/* Initialize variables, and get info (calibration switches and
				** reference file names) from primary headers.
				*/
				if (ACSRefInit (acshdr, &sci_sw, &sciref))
	    			return (status);

				if (asn->verbose) {
					trlmessage ("CALACS: Got reference file information");
				}
				/* Store the trailer file comments into preface */
				newpreface = YES;
								
			} else {

				/* Read in Member keyword info from image header	*/
				if (asn->process == SINGLE ) {
					if (GetSingle (asn, acshdr) )
						return (status);
				} else {	
					if( GetAsnMember(asn, prod, posid, expid, acshdr)) 
						return (status);
				} 

				/* Construct output and temporary file names. */
					if (InsertACSSuffix (acshdr))
	    				return (status);	    		
			}
			/* From this point on, we are working with ONE image at a time... */


	/* Copy switch values for ACSCCD and ACS2D argument lists. */
	    		SetACSSw (&sci_sw, &acsccd_sci_sw, &acs2d_sci_sw);
	
	/* Ready to process input image now...				*/
				if (acshdr->sci_basic_2d == PERFORM) {
								
					/* Copy all messages up to here into preface buffer
						to be prepended to all input file trailer files */
					if (newpreface == YES) {
						InitTrlPreface();
					}
					
					/* Flat field the input MAMA image. */
		    		if (ACS2d (acshdr->rawfile, acshdr->fltfile,
					&acs2d_sci_sw, &sciref, printtime, asn->verbose)) 
						return (status);
						
	    		} else {

					/* Copy the input file to fltfile.  We need to do this
				   because ACSSUM needs it for input.
					*/
					trlwarn ("No processing performed on image.");
                    if (acshdr->sci_rptcorr != PERFORM) 
                        status = NOTHING_TO_DO;
                    
					sprintf (MsgText, "Copying input to %s ",acshdr->fltfile);
					trlwarn(MsgText);

		    		if (CopyFFile (acshdr->rawfile, acshdr->fltfile))
		    			return (status);
	        	}
			} 				/* Finished processing EXP images		*/
		}		/* Finished processing this position's sub-product */
  	} /* Finished looping over all the products */

	/* 
		Reset the trailer file preface to NULL since
		it has already been copied into trailer files 
	*/
	ResetTrlPreface();

	/* Done with lists of reference file keywords and names. */
	FreeRefFile (&sciref);
	
	return (status);
}



/* This routine copies switch values from sci_sw to
   acsccd_sci_sw and acs2d_sci_sw.
*/

static void SetACSSw (CalSwitch *sci_sw, CalSwitch *acsccd_sci_sw, CalSwitch *acs2d_sci_sw) {

/* arguments:
CalSwitch *sci_sw             i: all calibration switches for science file
CalSwitch *acsccd_sci_sw     o: acsccd switches for science file
CalSwitch *acs2d_sci_sw      o: acs2d switches for science file
*/

	/* These are the switches for calacs prior to running acsrej
	   (cosmic-ray rejection).
	*/	
	acsccd_sci_sw->crcorr = sci_sw->crcorr;
	acsccd_sci_sw->dqicorr = sci_sw->dqicorr;
	acsccd_sci_sw->atodcorr = sci_sw->atodcorr;
	acsccd_sci_sw->blevcorr = sci_sw->blevcorr;
	acsccd_sci_sw->biascorr = sci_sw->biascorr;
	acsccd_sci_sw->flashcorr = sci_sw->flashcorr;
	acsccd_sci_sw->glincorr = OMIT;
	acsccd_sci_sw->lflgcorr = OMIT;
	acsccd_sci_sw->darkcorr = OMIT;
	acsccd_sci_sw->flatcorr = OMIT;
	acsccd_sci_sw->shadcorr = OMIT;
	acsccd_sci_sw->photcorr = OMIT;
	acsccd_sci_sw->rptcorr = OMIT;
    acsccd_sci_sw->expscorr = sci_sw->expscorr;	

	/* These are the switches for acs2d for the science file. */
	acs2d_sci_sw->dqicorr = sci_sw->dqicorr;
	acs2d_sci_sw->atodcorr = sci_sw->atodcorr;
	acs2d_sci_sw->blevcorr = sci_sw->blevcorr;
	acs2d_sci_sw->glincorr = sci_sw->glincorr;
	acs2d_sci_sw->lflgcorr = sci_sw->lflgcorr;
	acs2d_sci_sw->biascorr = sci_sw->biascorr;
	acs2d_sci_sw->flashcorr = sci_sw->flashcorr;
	acs2d_sci_sw->darkcorr = sci_sw->darkcorr;
	acs2d_sci_sw->flatcorr = sci_sw->flatcorr;
	acs2d_sci_sw->shadcorr = sci_sw->shadcorr;
	acs2d_sci_sw->photcorr = sci_sw->photcorr;
	acs2d_sci_sw->rptcorr = sci_sw->rptcorr;
	acs2d_sci_sw->crcorr = sci_sw->crcorr;
    acs2d_sci_sw->expscorr = sci_sw->expscorr;	
	
}

/* This routine resets dqicorr, atodcorr, and blevcorr for acs2d.
   These steps are done prior to cosmic ray rejection, then the
   steps for acs2d are done.  It isn't essential to reset these
   switches, since acs2d will recognize that blevcorr has already
   been done and will turn off that switch locally, but dqicorr will be
   done a second time, so it's a good idea to explicitly turn that off.
*/

static void ResetACSSw (CalSwitch *acsccd_sci_sw, CalSwitch *acs2d_sci_sw) {

/* arguments:
CalSwitch *acsccd_sci_sw     i: acsccd switches for science file
CalSwitch *acs2d_sci_sw      o: acs2d switches for science file
*/

	/* These are the switches performed by acsccd prior to running
	   acsrej (cosmic-ray rejection).
	*/
	if (acsccd_sci_sw->dqicorr == PERFORM)
	    acs2d_sci_sw->dqicorr = OMIT;
	if (acsccd_sci_sw->atodcorr == PERFORM)
	    acs2d_sci_sw->atodcorr = OMIT;
	if (acsccd_sci_sw->blevcorr == PERFORM)
	    acs2d_sci_sw->blevcorr = OMIT;
	if (acsccd_sci_sw->biascorr == PERFORM)
	    acs2d_sci_sw->biascorr = OMIT;
	if (acsccd_sci_sw->flashcorr == PERFORM)
	    acs2d_sci_sw->flashcorr = OMIT;
}

# define FITS_BUFSIZE  2880	/* size of a FITS block */

/* This routine copies a FITS file. */

static int CopyFFile (char *infile, char *outfile) {

/* arguments:
char *infile    i: name of input file
char *outfile   i: name of output file
*/

	extern int status;

	FILE *ifp, *ofp;	/* for input and output files */
	void *buf;		/* buffer for copying blocks */
	int nin, nout;		/* number read and written */
	int done;

	if ((buf = calloc (FITS_BUFSIZE, sizeof(char))) == NULL)
	    return (status = OUT_OF_MEMORY);

	if ((ofp = fopen (outfile, "wb")) == NULL) {
	    sprintf (MsgText,"Can't create temporary file %s.", outfile);
	    trlerror(MsgText);
		free (buf);
	    return (status = INVALID_TEMP_FILE);
	}

	if ((ifp = fopen (infile, "rb")) == NULL) {
	    sprintf ("Can't open %s.", infile);
	    trlerror (MsgText);
		fclose (ofp);
	    remove (outfile);
	    free (buf);
	    return (status = OPEN_FAILED);
	}

	done = 0;
	while (!done) {
	    nin = fread (buf, sizeof(char), FITS_BUFSIZE, ifp);
	    if (ferror (ifp)) {
		sprintf (MsgText, "Can't read from %s (copying to %s).",
				infile, outfile);
		trlerror (MsgText);
		fclose (ofp);
		fclose (ifp);
		free (buf);
		return (status = FILE_NOT_READABLE);
	    }
	    if (feof (ifp))
		done = 1;

	    nout = fwrite (buf, sizeof(char), nin, ofp);
	    if (nout < nin) {
		sprintf (MsgText, "Can't copy %s to %s.", infile, outfile);
		trlerror (MsgText);
		fclose (ofp);
		fclose (ifp);
		free (buf);
		return (status = COPY_NOT_POSSIBLE);
	    }
	}

	fclose (ofp);
	fclose (ifp);
	free (buf);

	return (status);
}

/* This function is an interface for calling ACSRej.  It sets up the
   default parameters, copies printtime and verbose into the par structure,
   and calls ACSRej.
*/

# include "acsrej.h"

static int ACSRej_0 (char *input, char *output, char *mtype, int newbias, int printtime, int verbose) {

	extern int status;
	clpar 	par;			/* parameters used */
	int 	newpar[MAX_PAR+1];	/* user specifiable parameters */

	void rej_reset (clpar *, int []);
	int AcsRej (char *, char *, char *, clpar *, int []);

	rej_reset (&par, newpar);
	par.printtime = printtime;
	par.verbose = verbose;
    par.newbias = newbias;

	
	status = AcsRej (input, output, mtype, &par, newpar);
	
	return (status);
}


/* This function compares the values of CRCORR and RPTCORR from
     the image header with those deduced from the ASN table. */
int CheckCorr (AsnInfo *asn, ACSInfo *acshdr) {
    
    extern int status;
    
    if ( (asn->crcorr != acshdr->sci_crcorr) || (asn->rptcorr != acshdr->sci_rptcorr)) {
        trlerror ("CRCORR and/or RPTCORR values not consistent with ASN table!");
        return (status = ERROR_RETURN);
    }
    return (status);
}
