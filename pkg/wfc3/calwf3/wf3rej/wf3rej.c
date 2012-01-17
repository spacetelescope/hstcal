# include <stdio.h>
# include "hstio.h"
# include <string.h>

# include "wf3.h"
# include "wf3rej.h"
# include "wf3err.h"
# include "rej.h"
# include "wf3info.h"

static void InitRejTrl (char *, char *, int); 

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
*/

int Wf3Rej (char *in_list, char *output, char *mtype, clpar *par, int newpar[]) 
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
    if (LoadHdr (in_name, &phdr)) {
        sprintf (MsgText, "Could not load header from %s", in_name);
        trlerror (MsgText);
        return (status);
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
	if (mkNewSpt (in_list, mtype, output)) {
	    return (status);
	}
	/* Clean up before returning */
	c_imtclose(tpin);
	WriteTrlFile ();
	/* Pass rej_do error status to caller */
	return (old_status);
    }

    /* create new SPT file for output product */
    if (mkNewSpt (in_list, mtype, output)) {
        return(status);
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

static void InitRejTrl (char *input, char *output, int detector) {

    extern int  status;
    IRAFPointer tpin;
    int         n;
    int         nfiles;                     /* Number of files in list */

    char        *trl_in;                   /* trailer filename for input */
    char        trl_out[SZ_LINE+1];       /* output trailer filename */
    char        in_name[SZ_FNAME+1];
    char        out_name[SZ_FNAME+1];

    int         trl_len;

    char        isuffix[] = "_blv_tmp";
    char        osuffix[] = "_crj_tmp";

    int MkName (char *, char *, char *, char *, char *, int);
    void WhichError (int);

/* ----------------------------- Begin ------------------------------*/

    trl_in =  realloc (NULL, (SZ_LINE +1));
    trl_len = SZ_LINE + 1;

    if (trl_in == NULL) {
	printf ("Out of memory: Couldn't allocate for CRJ_TMP trailer file.");
 	status = OUT_OF_MEMORY;
	trl_len = 0;
    }	

    /* Initialize TRL filenames */
    trl_in[0]  = '\0';
    trl_out[0] = '\0';

    /* Open the input file template */
    tpin = c_imtopen (input);
    nfiles = c_imtlen(tpin);

    /* Loop over the input file names */
    for (n = 0; n < nfiles; ++n) {
        c_imtgetim (tpin, in_name, SZ_FNAME);

        /* Check which WFC3 detector we're using */
	if (n == 0) {

	    /* Reset file name suffixes for IR images */
	    if (detector == IR_DETECTOR) {
		strcpy (isuffix, "_flt");
        	strcpy (osuffix, "_crj");
	    }
	}

        /* Start by stripping off suffix from input/output filenames */
        if (MkName (in_name, isuffix, "", TRL_EXTN, out_name, SZ_FNAME)) {
            WhichError (status);
            printf ("Couldn't determine trailer filename for %s", input);
            continue;
        }

        if ((strlen(out_name) + strlen(trl_in) + 1) >= trl_len) {
            trl_len += strlen(out_name)*(nfiles-n);
            trl_in = realloc (trl_in, trl_len);
        }

        /* Append each filename to create list of input trailer files */
        strcat(trl_in, out_name);

        /* Put a comma after all but the last filename */
        if (n < (nfiles-1)) strcat (trl_in, ",");		
    }

    if (MkName (output, osuffix, "", TRL_EXTN, trl_out, SZ_LINE)) {
        WhichError (status);
        printf ("Couldn't create trailer filename for %s", output);
    }

    /* Sets up temp trailer file for output and copies input
    ** trailer file into it.  */
    InitTrlFile (trl_in, trl_out);

    /* Deallocate memory */
    free(trl_in);
    c_imtclose (tpin);
}

