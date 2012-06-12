# include <stdio.h>
# include <stdlib.h>	/* calloc */
# include <string.h>

# include "hstio.h"

# include "wf3.h"
# include "calwf3.h"
# include "wf3err.h"
# include "wf3corr.h"
# include "wf3asn.h"	/* Contains association table structures */

# define NOPOSID 0


/* calwf3 -- integrated calwf3 processing

   Howard A. Bushouse, 2000 August 22:
	Initial version (adapted from calacs.c by W. Hack)
   H.Bushouse, 2001 May 8:
	Modified to keep in sync with changes to calacs:
	Added post-flash processing step;
	Corrected bugs dealing with dither associations that have only
	single images at each pointing.
   H.Bushouse, 2002 June 20:
	Modified to keep in sync with changes to calacs:
	Changed so that when no calibration switches are set to PERFORM, it
	only WARNs the user, not provide ERROR messages. Remove use of
	"statcorr" (always compute statistics by default). Added updateAsnStat
	routine. Added check against "RPTCORR" to be consistent with "CRCORR"
	check for determining whether any calibration switches are set to
	PERFORM. If only RPTCORR is on, then it will not stop but copy raw
	files and sum them with WF3SUM. Update DQICORR for CRJ product so that
	doDQI is not run again in WF32D.
   H.Bushouse, 2002 Nov 26:
	Modified to keep in sync with changes to calacs:
	Corrected memory usage error in BuildSumInput.
	Corrected a bug where RPT-OBS sub-products MTYPE was getting overwritten
	by the PROD-DTH MTYPE in dithered RPT-OBS associations.
	Fix involved using subprod.mtype instead of prod.mtype.
   H.Bushouse, 2003 Apr 25:
	Removed all handling of RPTCORR processing from CalWf3Run, because
	WFC3 repeat-obs images will be cr-combined (using WF3REJ) from
	within the ProcessCCD and ProcessIR routines, just like CR-SPLITs.
	Also modified SetCCDSw routine to treat RPTCORR and CRCORR the same.
   H.Bushouse, 2006 June 20:
	Modified to keep in sync with CALACS calacs.c changes:
	Removed any reference to updateAsnStat, since only OPUS should update
	the ASN_STAT keyword in the ASN table header.
   H.Bushouse, 2010 Oct 20:
	Modified BuildDthInput to skip processing for sub-products that have
	< 2 members. Also modified the DTH processing portion of CalWf3Run
	to skip processing for products that have no members. Both of these
	are for handling associations with missing members.  (PR 66366)
   H.Bushouse, 2011 Jan 14:
	Updated CopyFFile to update the FILENAME keyword in the output file.
	(PR 67225, Trac #646)
   H.Bushouse, 2012 Mar 21:
	Upgraded BuildDthInput to handle situations where no sub-product
	(i.e. _crj file) has been produced, which means the trailer file for
	the final (asn level) product must be built from the trailers of the
	individual asn members. (PR 70922, Trac #869)
   M.Sosey, 2012 May 07:
     added the option "-r" to print the current version and exit cleanly
*/

int CalWf3Run (char *input, int printtime, int save_tmp, int verbose, int debug) {

/* arguments:
char *input	i: name of the FITS file/table to be processed
int printtime	i: true --> print time stamps at intermediate steps
int save_tmp	i: true --> save temporary files
int verbose	i: true --> print info during processing
int debug	i: true --> print debugging info during processing
*/

	extern int status;

	WF3Info wf3hdr;		/* calibration switches, etc */
	AsnInfo	asn;		/* association table data    */
	char *wf3dth_input;	/* Input list for WF3DTH */

	int prod;
	
	void PrBegin (char *);
	void PrEnd (char *);
	void PrFileName (char *, char *);
	void TimeStamp (char *, char *);

	/* Association table routines			*/
	void initAsnInfo (AsnInfo *);
	void freeAsnInfo (AsnInfo *);
	int LoadAsn (AsnInfo *);
	int ProcessCCD (AsnInfo *, WF3Info *, int *, int);
	int ProcessIR  (AsnInfo *, WF3Info *, int);	
	int Wf3Dth (char *, char *, int, int, int);
	char *BuildDthInput (AsnInfo *, int);
	int updateAsnTable (AsnInfo *, int, int);
	void InitDthTrl (char *, char *);
	
	/* Post error handler */
	push_hstioerr (errchk);
	    
	PrBegin ("CALWF3");	/* *** CALWF3 -- Version ... *** */

	if (printtime)
	    TimeStamp ("CALWF3 started", "");

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
	** Association table.  If a single image, it will look in that images
	** association table to see what products are associated with it, and
	** process them accordingly.  If it is just a single image as its
	** own output, it will proceed as a 1 element table.	
	** Based on routines from n_getAsnTable() and n_setup() in CALNICB
	*/
	if (LoadAsn(&asn)) {
	   freeAsnInfo (&asn);
	   return (status);
	}

	/* Check to see that detector is known, as it could come through with
	** a value of 0 or UNKNOWN_DETECTOR.  
        ** WJH  2Mar99
	*/
	if (asn.detector == UNKNOWN_DETECTOR || asn.detector == 0) {
	    trlwarn ("Unknown detector type for observations.");
	    freeAsnInfo(&asn);	
	    return (status = NOTHING_TO_DO);
	}

	if (asn.verbose) {
	    sprintf (MsgText,"CALWF3: Detector %s, type %d ",
		     asn.instr, asn.detector);
	    trlmessage (MsgText);
	}
	
	/* Determine what detector we are working with...	*/
	if (asn.detector == CCD_DETECTOR ) { /* Process CCD data ... */

	    if (asn.verbose) {
		trlmessage ("CALWF3: processing a UVIS product");
	    }
	    if (ProcessCCD (&asn, &wf3hdr, &save_tmp, printtime)) { 
		if (status == NOTHING_TO_DO) {
		    trlwarn ("No processing desired for CCD data.");
		} else {
		    trlerror ("Couldn't process CCD data");
		}
		freeAsnInfo(&asn);
		return (status);			
	    }		

	} else { /* Process IR observations here */

	    if (asn.verbose) {
		trlmessage ("CALWF3: processing an IR product");
	    }
	    if (ProcessIR (&asn, &wf3hdr, printtime)) { 
		if (status == NOTHING_TO_DO) {
		    trlwarn ("No processing desired for IR data.");
		} else {
		    trlerror ("Couldn't process IR data");
		}
		freeAsnInfo(&asn);
		return (status);			
	    }			
	}


	/* Add DTH processing here... */
	if (asn.process == FULL) {
	    if (asn.verbose) {
		trlmessage ("CALWF3: Building DTH products");
	    }

	    /* For each DTH product... */
	    wf3dth_input = NULL;
	    for (prod = 0; prod < asn.numprod; prod++) {

		 /* Create empty DTH product, header only */
		 /* This uses only one sub-product for the header template,
		 ** but later versions should use a function similar to
		 ** BuildSumInput to create list of subproducts as inputs...
		 */

		 wf3dth_input = BuildDthInput (&asn, prod);

		 /* Skip this product if the input list is empty */
		 if (wf3dth_input == NULL) continue;

		 /* We always want to create a final concatenated trailer file
		 ** for the entire association whether there is a product or
		 ** not. So we set up the trailer file based on the association
		 ** file name itself. */

		 /* If desired, we could optionally use the full _drz.tra
		 ** filename as the trailer file name, based on the output
		 ** dither product name...
		 if (strcmp(asn.product[prod].prodname,"") != 0) {
		     InitDthTrl (wf3dth_input, asn.product[prod].prodname);
		 } else { */

		 InitDthTrl (wf3dth_input, asn.rootname);

		 /* } */

		 /* Check if we have a PROD-DTH specified */
		 if (strcmp(asn.product[prod].prodname,"") != 0) {

		     if ((asn.dthcorr == PERFORM || asn.dthcorr == DUMMY)) {
		          if (Wf3Dth (wf3dth_input,
				 asn.product[prod].prodname, asn.dthcorr,
				 printtime, asn.verbose) )
			      return (status);

		          /* Pass posid=0 to indicate a PRODUCT is to
			  ** be updated */
		          updateAsnTable(&asn, prod, NOPOSID);

		     } else {

		        trlwarn
		  	("No DTH product name specified. No product created.");
			/* status = WF3_OK; */
		     }
		 }
	         free (wf3dth_input);
	    }
	} /* End DTHCORR Processing */
	
	if (asn.verbose) {
	    trlmessage ("CALWF3: Finished processing product ");
	}
	
	freeAsnInfo(&asn);
	
	trlmessage ("\n");
	PrEnd ("CALWF3");		/* *** CALWF3 complete *** */

	if (printtime)
	    TimeStamp ("CALWF3 completed", wf3hdr.rootname);
	
	/* Return the final value of status, which should be WF3_OK if
	** all went well or some error condition otherwise. */
	return (status);

}

char *BuildSumInput (AsnInfo *asn, int prod, int posid) {

	int nchars;
	int wf3sum_len;
	int i;
	char *wf3sum_input;
	char tmpexp[SZ_LINE+1];
	char tmpflt[SZ_LINE+1];
	int MkName (char *, char *, char *, char *, char *, int);

	/* Determine how long this string needs to be... */
	/*nchars = asn->spmems[posid] * (SZ_FNAME+1);*/
	nchars = 1;

	/* Keep track of individual filename lengths and total length */
	wf3sum_len = 0;

	/* Now, lets search the association table for all inputs... */
	for (i=1; i <= asn->spmems[posid]; i++) {
	     strcpy(tmpexp, asn->product[prod].subprod[posid].exp[i].expname);

	     if (MkName (tmpexp, "_raw", "_flt", "", tmpflt, SZ_LINE)) {
		 strcpy (tmpflt, asn->product[prod].subprod[posid].exp[i].name);
		 strcat (tmpflt, "_flt.fits");
	     }

	     /* Sum together lengths of all filenames to be used as
	     ** input to WF3SUM. */
	     wf3sum_len += strlen(tmpflt)+1;
	}

	/* Now that we know how long it will be, allocate space for
	** input string */
	wf3sum_input = (char *) calloc(wf3sum_len + 1, sizeof(char));

	/* Go back through the association table to string all inputs together*/
	for (i=1; i <= asn->spmems[posid]; i++) {
	     strcpy (tmpexp, asn->product[prod].subprod[posid].exp[i].expname);

	     if (MkName (tmpexp, "_raw", "_flt", "", tmpflt, SZ_LINE)) {
		 strcpy (tmpflt,asn->product[prod].subprod[posid].exp[i].name);
		 strcat (tmpflt, "_flt.fits");
	     }
	     strcat(wf3sum_input, tmpflt);

	     /* Don't add a comma to the end of the last filename*/
	     if (i < (asn->spmems[posid])) {
		 strcat(wf3sum_input, ",");
	     }
	}
	strcat (wf3sum_input, "\0");

	return(wf3sum_input);
}


char *BuildDthInput (AsnInfo *asn, int prod) {

	int nchars;
	int i, j;
	char *wf3dth_input;
	char tmpexp[SZ_LINE+1];
	int MkName (char *, char *, char *, char *, char *, int);

	/* Determine how long this string needs to be... */
	nchars = asn->numsp * SZ_FNAME;
	wf3dth_input = (char *) calloc( nchars + 1, sizeof(char));

	/* Initialize this string to NULL */
	wf3dth_input[0] = '\0';

	/* Now, lets search the association table for all inputs... */
	for (i=1; i <= asn->numsp; i++) {

	     /* Skip CR and RPT sub-products that only have 1 member */
	     if ((asn->crcorr==PERFORM || asn->crcorr==DUMMY ||
		  asn->rptcorr==PERFORM || asn->rptcorr==DUMMY) && 
		  asn->spmems[i] < 2)
		  continue;

	     /* Check to see if a sub-product was produced */
	     if (!asn->product[prod].subprod[i].prsnt) {

		 /* If the sub-product wasn't produced, we have to build the
		 ** list of names from the individual asn members */
                 for (j=1; j <= asn->spmems[i]; j++) {
		      strcpy(tmpexp, asn->product[prod].subprod[i].exp[j].name);
		      strcat(tmpexp,"_flt.fits");
		      strcat(wf3dth_input, tmpexp);
		      if (j < asn->spmems[i])
			  strcat(wf3dth_input, ",");
		 }

	     } else {

		 /* Otherwise just use the sub-product name */
	         strcpy(tmpexp, asn->product[prod].subprod[i].spname);
	         strcat(wf3dth_input, tmpexp);
	     }

	     if (i < (asn->numsp)) {
		 /* Add a comma to the end if this is not the last filename*/
		 strcat(wf3dth_input, ",");
	     }
	}

	return(wf3dth_input);
}


/* This routine copies switch values from sci_sw to
   wf3ccd_sci_sw and wf32d_sci_sw.


   H.Bushouse 2003-04-25  Modified to treat RPTCORR exactly the same
                          as CRCORR for WFC3 UVIS use.
*/

void SetCCDSw (CCD_Switch *sci_sw, CCD_Switch *wf3ccd_sci_sw,
	       CCD_Switch *wf32d_sci_sw) {

/* arguments:
CCD_Switch *sci_sw            i: all calibration switches for science file
CCD_Switch *wf3ccd_sci_sw     o: wf3ccd switches for science file
CCD_Switch *wf32d_sci_sw      o: wf32d switches for science file
*/

	/* These are the switches for calwf3 prior to running wf3rej
	** (cosmic-ray rejection). */	
	wf3ccd_sci_sw->dqicorr   = sci_sw->dqicorr;
	wf3ccd_sci_sw->atodcorr  = sci_sw->atodcorr;
	wf3ccd_sci_sw->blevcorr  = sci_sw->blevcorr;
	wf3ccd_sci_sw->biascorr  = sci_sw->biascorr;
	wf3ccd_sci_sw->flashcorr = sci_sw->flashcorr;
	wf3ccd_sci_sw->darkcorr  = OMIT;
	wf3ccd_sci_sw->flatcorr  = OMIT;
	wf3ccd_sci_sw->shadcorr  = OMIT;
	wf3ccd_sci_sw->photcorr  = OMIT;
	wf3ccd_sci_sw->crcorr    = sci_sw->crcorr;
	wf3ccd_sci_sw->rptcorr   = sci_sw->rptcorr;
	wf3ccd_sci_sw->expscorr  = sci_sw->expscorr;	

	/* These are the switches for wf32d for the science file. */
	wf32d_sci_sw->dqicorr   = sci_sw->dqicorr;
	wf32d_sci_sw->atodcorr  = sci_sw->atodcorr;
	wf32d_sci_sw->blevcorr  = sci_sw->blevcorr;
	wf32d_sci_sw->biascorr  = sci_sw->biascorr;
	wf32d_sci_sw->flashcorr = sci_sw->flashcorr;
	wf32d_sci_sw->darkcorr  = sci_sw->darkcorr;
	wf32d_sci_sw->flatcorr  = sci_sw->flatcorr;
	wf32d_sci_sw->shadcorr  = sci_sw->shadcorr;
	wf32d_sci_sw->photcorr  = sci_sw->photcorr;
	wf32d_sci_sw->crcorr    = sci_sw->crcorr;
	wf32d_sci_sw->rptcorr   = sci_sw->rptcorr;
	wf32d_sci_sw->expscorr  = sci_sw->expscorr;	
	
}

/* This routine resets dqicorr, atodcorr, and blevcorr for wf32d.
   These steps are done prior to cosmic ray rejection, then the
   steps for wf32d are done.  It isn't essential to reset these
   switches, since wf32d will recognize that blevcorr has already
   been done and will turn off that switch locally, but dqicorr will be
   done a second time, so it's a good idea to explicitly turn that off.
*/

void ResetCCDSw (CCD_Switch *wf3ccd_sci_sw, CCD_Switch *wf32d_sci_sw) {

/* arguments:
CCD_Switch *wf3ccd_sci_sw     i: wf3ccd switches for science file
CCD_Switch *wf32d_sci_sw      o: wf32d switches for science file
*/

	/* These are the switches performed by wf3ccd prior to running
	** wf3rej (cosmic-ray rejection). */

	if (wf3ccd_sci_sw->dqicorr == PERFORM)
	    wf32d_sci_sw->dqicorr = OMIT;
	if (wf3ccd_sci_sw->atodcorr == PERFORM)
	    wf32d_sci_sw->atodcorr = OMIT;
	if (wf3ccd_sci_sw->blevcorr == PERFORM)
	    wf32d_sci_sw->blevcorr = OMIT;
	if (wf3ccd_sci_sw->biascorr == PERFORM)
	    wf32d_sci_sw->biascorr = OMIT;
	if (wf3ccd_sci_sw->flashcorr == PERFORM)
	    wf32d_sci_sw->flashcorr = OMIT;
}


/* This routine copies switch values from sci_sw to wf3ir_sci_sw.
*/

void SetIRSw (IR_Switch *sci_sw, IR_Switch *wf3ir_sci_sw) {

/* arguments:
IR_Switch *sci_sw           i: all calibration switches for science file
IR_Switch *wf3ir_sci_sw     o: wf3ir switches for science file
*/

	wf3ir_sci_sw->zsigcorr = sci_sw->zsigcorr;
	wf3ir_sci_sw->zoffcorr = sci_sw->zoffcorr;
	wf3ir_sci_sw->dqicorr  = sci_sw->dqicorr;
	wf3ir_sci_sw->blevcorr = sci_sw->blevcorr;
	wf3ir_sci_sw->noiscorr = sci_sw->noiscorr;
	wf3ir_sci_sw->darkcorr = sci_sw->darkcorr;
	wf3ir_sci_sw->nlincorr = sci_sw->nlincorr;
	wf3ir_sci_sw->flatcorr = sci_sw->flatcorr;
	wf3ir_sci_sw->unitcorr = sci_sw->unitcorr;
	wf3ir_sci_sw->photcorr = sci_sw->photcorr;
	wf3ir_sci_sw->crcorr   = sci_sw->crcorr;
	wf3ir_sci_sw->rptcorr  = sci_sw->rptcorr;

}

# define FITS_BUFSIZE  2880	/* size of a FITS block */

/* This routine copies a FITS file. */

int CopyFFile (char *infile, char *outfile) {

/* arguments:
char *infile    i: name of input file
char *outfile   i: name of output file
*/

	extern int status;

	FILE *ifp, *ofp;	/* for input and output files */
	void *buf;		/* buffer for copying blocks */
	int nin, nout;		/* number read and written */
	int done;
	IODescPtr im;		/* descriptor for primary header unit */
	Hdr phdr;		/* primary header */

	void UFilename (char *, Hdr *);

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

	/* Update the FILENAME keyword in the primary header of the 
	** output file. */

	initHdr (&phdr);
	im = openUpdateImage (outfile, "", 0, &phdr);
	UFilename (outfile, &phdr);
	putHeader (im);

	closeImage (im);
	freeHdr (&phdr);

	return (status);
}

/* This function is an interface for calling WF3Rej.  It sets up the
   default parameters, copies printtime and verbose into the par structure,
   and calls WF3Rej.
*/

# include "wf3rej.h"

int WF3Rej_0 (char *input, char *output, char *mtype, int printtime,
	      int verbose) {

	extern int status;
	clpar 	par;			/* parameters used */
	int 	newpar[MAX_PAR+1];	/* user specifiable parameters */

	void rej_reset (clpar *, int []);
	int Wf3Rej (char *, char *, char *, clpar *, int []);

	rej_reset (&par, newpar);
	par.printtime = printtime;
	par.verbose = verbose;

	status = Wf3Rej (input, output, mtype, &par, newpar);
	
	return (status);
}


/* This function compares the values of CRCORR and RPTCORR from
     the image header with those deduced from the ASN table. */

int CheckCorr (AsnInfo *asn, WF3Info *wf3hdr) {
    
    extern int status;
    
    if ((asn->crcorr  != wf3hdr->sci_crcorr) ||
	(asn->rptcorr != wf3hdr->sci_rptcorr)) {
        trlerror 
		("CRCORR and/or RPTCORR values not consistent with ASN table!");
        return (status = ERROR_RETURN);
    }
    return (status);
}

