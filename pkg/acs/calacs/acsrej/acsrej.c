# include <stdio.h>
# include <hstio.h>
# include <string.h>

# include "acs.h"
# include "acsrej.h"
# include "err.h"
# include "rej.h"

static void InitRejTrl (char *, char *); 

/*  ACSREJ -- cosmic ray rejection with shading correction and 
 section-by-section checking for ACS
 
 Description:
 ------------
 
 Date          Author      Description
 ----          ------      -----------
 08-27-1998    W.J. Hack   Adapt from CALSTIS2
 */

int AcsRej (char *in_list, char *output, char *mtype, clpar *par, int newpar[]) 
{
  
  extern int      status;
  IRAFPointer     tpin;
  int             flag;
  int             old_status;
  
  void        TimeStamp (char *, char *);
  void        PrBegin (char *);
  void        PrSwitch (char *, int);
  void        PrEnd (char *);
  int         acsrej_do(IRAFPointer, char *, char *, clpar *, int []);
  int         mkNewSpt (char *, char *, char *);
  /* ----------------------------- Begin ------------------------------*/
  
  /* Determine input and output trailer files, then initialize
   output file by combining inputs into output file */
  InitRejTrl(in_list, output);
  
  /* Quit on error condition */
  if (status)
    return (status);
  
  PrBegin ("ACSREJ");
  
  /*  announce start of the task */
  if (par->printtime) 
    TimeStamp ("ACSREJ version started", "");
  
  
  /* make sure the output datafile does not exist */
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
  
  /* open the input file template */
  tpin = c_imtopen (in_list);
  
  trlmessage ("\n");
  PrSwitch ("crcorr", PERFORM);
  
  /* perform the calculation */
  if (acsrej_do (tpin, output, mtype, par, newpar)) {
    trlmessage("Did NOT create ACSREJ product. Trying to build SPT...");
    /* We always want to try to create a new SPT file...*/
    old_status = status;
    status = ACS_OK;   
    if (mkNewSpt (in_list, mtype, output)) {
      return(status);
    }
    status = old_status;
    c_imtclose(tpin);
    return (status);
  }
  
  /* create new SPT file for output product */
  if (mkNewSpt (in_list, mtype, output)) {
    return(status);
  }
  
  /* Update calibration switch to reflect proper execution
   of this processing step.
   */
  PrSwitch ("crcorr", COMPLETE);
  
  /* close file template */
  c_imtclose (tpin);
  
  trlmessage ("\n");
  PrEnd ("ACSREJ");
  
  if (par->printtime)
    TimeStamp ("ACSREJ complete", "");
  
  /* Write out temp trailer file to final file */
  WriteTrlFile ();
  
  return (ACS_OK);
}

/* ------------------------------------------------------------------*/
/*                          InitRejTrl                               */
/* ------------------------------------------------------------------*/

static void InitRejTrl (char *input, char *output) {
  
  extern int  status;
  IRAFPointer tpin;
  int         n;
  int         nfiles;                     /* Number of files in list */
  
  char        *trl_in;                   /* trailer filename for input */
  char        trl_out[ACS_LINE+1];       /* output trailer filename */
  char        in_name[ACS_FNAME+1];
  char        out_name[ACS_FNAME+1];
  
  int         trl_len;
  
  char        isuffix[] = "_blv_tmp";
  char        osuffix[] = "_crj_tmp";
  char        itrlsuffix[] = "";
  char        otrlsuffix[] = "";
  
  int MkName (char *, char *, char *, char *, char *, int);
  void WhichError (int);
  
  /* check whether we're dealing with _blv_tmp or _blc_tmp files and
   * adjust isuffix and osuffix accordingly */
  if (strstr(input, isuffix) == NULL) {
    strcpy(isuffix, "_blc_tmp");
    strcpy(osuffix, "_crc_tmp");
  /*  strcpy(itrlsuffix, "_flc");
    strcpy(otrlsuffix, "_crc");  */
  }
  
  /* ----------------------------- Begin ------------------------------*/
  
  trl_in =  realloc (NULL, (ACS_LINE +1));
  trl_len = ACS_LINE + 1;
  
  if (trl_in == NULL) {
    printf ("Out of memory: Couldn't allocate for CRJ_TMP trailer file.");
    status = OUT_OF_MEMORY;
    trl_len = 0;
  }	
  /* Initialize TRL filenames */
  trl_in[0] = '\0';
  trl_out[0] = '\0';
  
  /* open the input file template */
  tpin = c_imtopen (input);
  nfiles = c_imtlen(tpin);
  
  for (n = 0; n < nfiles; ++n) {
    c_imtgetim (tpin, in_name, ACS_FNAME);
    
    /* Start by stripping off suffix from input/output filenames */
    if (MkName (in_name, isuffix, itrlsuffix, TRL_EXTN, out_name, ACS_FNAME)) {
      WhichError (status);
      printf ("Couldn't determine trailer filename for %s", input);
      continue;
    }
    
    if ( (strlen(out_name) + strlen(trl_in) + 1) >= trl_len) {
      trl_len += strlen(out_name) * (nfiles - n);
      trl_in = realloc (trl_in, trl_len);
    }
    
    /* Append each filename to create list of input trailer files */
    strcat(trl_in, out_name);
    /* But don't put a comma after the last filename */
    if (n < (nfiles-1)) strcat (trl_in, ",");		
  }
  
  if (MkName (output, osuffix, otrlsuffix, TRL_EXTN, trl_out, ACS_LINE)) {
    WhichError (status);
    printf ("Couldn't create trailer filename for %s", output);
  }
  
  /* Sets up temp trailer file for output and copies input
   trailer file into it.
   */
  InitTrlFile (trl_in, trl_out);
  
  /* Deallocate memory */
  free(trl_in);
  c_imtclose (tpin);
}
