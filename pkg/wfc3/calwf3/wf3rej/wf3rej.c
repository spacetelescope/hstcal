# include <stdio.h>
# include <string.h>
# include "hstio.h"


# include "wf3.h"
# include "wf3rej.h"
# include "hstcalerr.h"
# include "rej.h"
# include "wf3info.h"

int InitRejTrl (char *, char *); 
int MkOutName (char *, char **, char **, int, char *,int);
void WhichError (int);

/*  WF3REJ -- cosmic ray rejection with shading correction and 
    section-by-section checking for WF3

Description:
------------

Date          Author		Description
----          -------------	-----------
08-27-1998	W.J. Hack	Adapt from CALSTIS2
08-29-2000	H.A. Bushouse	Adapt from ACSREJ
10-27-2003	H. Bushouse	Always try to create new SPT file even if
REJ product not created (CALACS changes).
02-14-2005	H. Bushouse	Fix memory reallocation problem in InitRejTrl.
05-22-2008	H. Bushouse	Make use of wf3.detector to id WFC3 IR images.
06-16-2011	H. Bushouse	Fixed error return from rej_do so that original
status gets passed up for use in caller.
(PR 68593; Trac #722)
06-10-2015  M. Sosey, updates for CTE addtion, #1193 
08-01-2015 M. Sosey, removed detector dependence on InitRejTrl
08-18-2015 M. Sosey, variable size array allocation for input list and made return from
   InitRejTrl an int for error checking

*/

int Wf3Rej (char *in_list, char *output, char *mtype, clpar *par, int newpar[], int makespt) 
{

    extern int      status;
    IRAFPointer     tpin;
    int             flag;
    int		    old_status;
    WF3Info	    wf3;
    Hdr		    phdr;
    char	    det[SZ_CBUF+1];
    char	    in_name[SZ_FNAME+1];

    void        TimeStamp (char *, char *);
    void        PrBegin (char *);
    void        PrSwitch (char *, int);
    void        PrEnd (char *);
    int         rej_do(IRAFPointer, char *, char *, clpar *, int [], int);
    int         mkNewSpt (char *, char *, char *);
    int 	LoadHdr (char *, Hdr *);
    int         GetKeyStr (Hdr *, char *, int, char *, char *, int);

    /* ----------------------------- Begin ------------------------------*/

    /* Determine which detector is in use from first input image */
    
    tpin = c_imtopen (in_list);
    c_imtgetim (tpin, in_name, SZ_FNAME);
    initHdr (&phdr);
    if (LoadHdr (in_name, &phdr)) {
        sprintf (MsgText, "Could not load header from %s", in_name);
        trlerror (MsgText);
        return (status=ERROR_RETURN);
    }
    if (GetKeyStr (&phdr, "DETECTOR", NO_DEFAULT, "", det, SZ_CBUF)) {
        trlkwerr ("DETECTOR", in_name);
        return(status = KEYWORD_MISSING);
    }
    freeHdr (&phdr);
    c_imtclose (tpin);

    /* Check detector keyword value */
    if (strcmp (det, "IR") == 0) {
        wf3.detector = IR_DETECTOR;
    } else if (strcmp (det, "UVIS") == 0) {
        wf3.detector = CCD_DETECTOR;
    } else {
        sprintf (MsgText, "DETECTOR = %s is unknown.", det);
        trlwarn (MsgText);
        wf3.detector = UNKNOWN_DETECTOR;
    }
    
    /* Determine input and output trailer files, then initialize
     ** output file by combining inputs into output file */
     
    if (InitRejTrl(in_list, output)){
        WhichError(status);
        return (status);
    }
        
    /* Quit on error condition */
    if (status)
        return (status);

    PrBegin ("WF3REJ");

    /* Announce start of the task */
    if (par->printtime)
        TimeStamp ("WF3REJ started", "");

    /* Make sure the output datafile does not exist */
    flag = ckNewFile (output);
    if (flag > 0) {
        if (flag == 1) {
            sprintf (MsgText, "Output file `%s' already exists.", output);
            trlerror (MsgText);
            return (status = ERROR_RETURN);
        } else {
            sprintf (MsgText, "Can't clobber `%s'.", output);
            trlerror (MsgText);
            return (status = ERROR_RETURN);
        }
    }

    /* Open the input file template */
    tpin = c_imtopen (in_list);

    trlmessage ("\n");
    if (wf3.detector == IR_DETECTOR)
        PrSwitch ("rptcorr", PERFORM);
    else
        PrSwitch ("crcorr", PERFORM);

    /* Perform the calculation */
    if (rej_do (tpin, output, mtype, par, newpar, wf3.detector)) {
        trlmessage ("Did NOT create WF3REJ product. Trying to build SPT...");
        /* We always want to try to create a new SPT file */
        old_status = status;
        status = WF3_OK;
        if (makespt){
            if (mkNewSpt (in_list, mtype, output)) {
                return (status);
            }
        }
        /* Clean up before returning */
        c_imtclose(tpin);
        WriteTrlFile ();
        /* Pass rej_do error status to caller */
        return (old_status);
    }

    /* create new SPT file for output product */
    if (makespt){
        if (mkNewSpt (in_list, mtype, output)) {
            return(status);
        }
    }
    
    /* Update calibration switch to reflect proper execution
     **  of this processing step.  */
    if (wf3.detector == IR_DETECTOR)
        PrSwitch ("rptcorr", COMPLETE);
    else
        PrSwitch ("crcorr", COMPLETE);

    /* Close file template */
    c_imtclose (tpin);

    trlmessage ("\n");
    PrEnd ("WF3REJ");

    if (par->printtime)
        TimeStamp ("WF3REJ complete", "");

    /* Write out temp trailer file to final file */
    WriteTrlFile ();

    return (WF3_OK);
}

/* ------------------------------------------------------------------*/
/*                          InitRejTrl                               */
/* ------------------------------------------------------------------*/

int InitRejTrl (char *input, char *output) {

    extern int  status;
    IRAFPointer tpin;
    int         n;
    int         nfiles;             /* Number of files in list */

    char        *trl_in;            /* trailer filename for input */
    char        trl_out[SZ_LINE+1]; /* output trailer filename */
    char        in_name[SZ_FNAME];
    char        out_name[SZ_FNAME];


    /* Input and output suffixes for uvis. */
    char *isuffix[] = {"_blv_tmp", "_blc_tmp","_flt","_flc","_crj","_crc","_crj_tmp","_crc_tmp"};
	char *trlsuffix[] = {"", "", "", "", "", "", "", ""};
    int nsuffix=8;
    
    int MkName (char *, char *, char *, char *, char *, int);
    int MkNewExtn(char *, char *);
    void WhichError (int);

    if ((trl_in =  (char *) calloc(strlen(input)+ 1, sizeof(char)))== NULL){
        printf("\nCannot allocate memory for input string\n");
        return (status=OUT_OF_MEMORY);
    }
    
    /* Initialize TRL filenames */
    trl_in[0]  = '\0';
    trl_out[0] = '\0';
    in_name[0] = '\0';
    out_name[0]= '\0';

    /* Open the input file template */
    tpin = c_imtopen (input);
    nfiles = c_imtlen(tpin);

    /* Loop over the input file names */
    for (n = 0; n < nfiles; ++n) {
        c_imtgetim (tpin, in_name, SZ_FNAME);

        /* Start by stripping off suffix from input/output filenames */
        if (MkOutName (in_name, isuffix, trlsuffix, nsuffix , out_name, SZ_FNAME)) {
	        WhichError (status);
	        sprintf (MsgText, "Couldn't determine trailer filename for %s",in_name);
	        trlmessage (MsgText);
        }

        if (MkNewExtn (out_name,TRL_EXTN)){   
 	        WhichError (status);
	        sprintf (MsgText, "Couldn't create trailer filename for %s",out_name);
	        trlmessage (MsgText);
	    }
                
        /* Append each filename to create list of input trailer files */
        strcat(trl_in,out_name );

        /* Put a comma after all but the last filename */
        if (n < (nfiles-1)) strcat (trl_in, ",");		
    }
    
    if (MkOutName(output,isuffix,trlsuffix,nsuffix,trl_out,SZ_FNAME)){
        WhichError(status);
        sprintf(MsgText,"Couldn't create trailer filename for %s",output);
    }
    if (MkNewExtn (trl_out,TRL_EXTN)){   
 	    WhichError (status);
	    sprintf (MsgText, "Couldn't create trailer filename for %s",out_name);
	    trlmessage (MsgText);
	}

    /* Sets up temp trailer file for output and copies input
     ** trailer file into it.  */
    InitTrlFile (trl_in, trl_out);

    /* Deallocate memory */
    free(trl_in);
    c_imtclose (tpin);
    return(status);
}

