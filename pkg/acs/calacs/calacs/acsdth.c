# include <stdio.h>
# include <stdlib.h>	/* calloc */
# include <string.h>

/* This section for EmptyGroup functions */
# if defined(VMS)
# include <stat.h>
# else
# include <sys/types.h>
# include <sys/stat.h>
# endif

# include "hstio.h"
# include "xtables.h"

# include "acs.h"
# include "calacs.h"
# include "hstcalerr.h"

/* Custom HSTIO-based routines for creating a SingleGroup without arrays
   (not used) */
/*
static int putEmptyShortHD(char *, char *, int, ShortHdrData *, int);
static int putEmptyFloatHD(char *, char *, int, FloatHdrData *, int);
static int putEmptyGroup(char *, int, SingleGroup *, int);
static int putEmptySci(char *, int, SciHdrData *, int);
static int putEmptyErr(char *, int, ErrHdrData *, int);
static int putEmptyDQ(char *, int, DQHdrData *, int);
*/

/* acsdth -- produce DTH product (empty)

 Warren J. Hack, 1998 Oct 14:
 Initial version.
 Warren J. Hack, 1999 Nov 11:
 Adds creation of SPT file for output product.
 Warren J. Hack, 2000 Jun 23:
 Modified trailer file creation to include ALL input files.
 Also modified SPT file creation to include all extensions
 from input SPT files.
 Warren J. Hack, 2002 Feb 6:
 Moved trailer file initialization outside this so that a
 trailer file can be created/appended for all ASN tables regardless
 of whether they contain actual dither members or not.
 This initialization is now handled by CalAcsRun in 'calacs.c'.
 Warren J. Hack, 2002 April 18:
 Eliminated creation of dummy '_dth.fits' product.

 int AcsDth (char *input, char *output, int dthcorr, int printtime, int verbose) {
 */

int AcsDth (char *in_list, char *output, int dthcorr, int printtime, int verbose) {

	extern int status;

	/*int option = 0;*/  /* For creating new output image */
	/*int logit;*/  /* true if we should log file name */
	char        root[ACS_CBUF+1];    /* ROOTNAME for output file */
  char        mtype[ACS_CBUF+1];  /* role of exposure in Association */
  IRAFPointer tpin;
  char        input[ACS_FNAME];  /* Name of image in list */

	void PrBegin (char *);
	void PrEnd (char *);
	int GetKeyInt (Hdr *, char *, int, int, int *);
	int GetKeyStr (Hdr *, char *, int, char *, char *, int);
	int UpdateSwitch (char *, int, Hdr *, int *);
  void TimeStamp (char *, char *);
	int LoadHdr (char *, Hdr *);
	void UCalVer (Hdr *);
	void UFilename (char *, Hdr *);
	void UMemType (char *, Hdr *);
	void UExpname (char *, Hdr *);
  int PutKeyFlt (Hdr *, char *, float, char *);
  int PutKeyInt (Hdr *, char *, int, char *);
	/*void InitDthTrl (char *, char *); */
  int         mkNewSpt (char *, char *, char *);

  /* ----------------------- Start Code --------------------------------*/

	/* Determine the names of the trailer files based on the input
   and output file names, then initialize the trailer file buffer
   with those names.
   */

  /*	InitDthTrl (in_list, output); */
	root[0] = '\0';
	sprintf(mtype,"PROD-DTH");

	/* Start the task... */
	PrBegin ("ACSDTH");

  /* open the input file template */
  tpin = c_imtopen (in_list);
  /* First, let's determine how many extensions/chips in each file */
  c_imtgetim (tpin, input, ACS_FNAME);

	if (printtime)
    TimeStamp ("ACSDTH started", "");

  sprintf(MsgText,"The task PyDrizzle needs to be run in order to generate");
  trlmessage(MsgText);
  sprintf(MsgText,"a geometrically corrected, drizzle-combined product.");
  trlmessage(MsgText);
  sprintf(MsgText,"PyDrizzle requires PyRAF. See pyraf.stsci.edu for more details.");
  trlmessage(MsgText);

  /* create new SPT file for output product */
  if (mkNewSpt (in_list, mtype, output)) {
    return(status);
  }

  c_imtclose(tpin);
	PrEnd ("ACSDTH");

	/* Write out temp trailer file to final file */
	WriteTrlFile ();

	return (status);
}


/* Custom routines copied from HSTIO and revised, to output a SingleGroup
 without any array data.
 Not used.
 */
#if false
static int putEmptyGroup(char *fname, int ever, SingleGroup *x, int option) {
  struct stat buf;
  if (option == 0) {
    if (stat(fname,&buf) == -1)
      putSingleGroupHdr(fname,x,0);
  }
  putEmptySci(fname,ever,&(x->sci),option); if (hstio_err()) return -1;
  putEmptyErr(fname,ever,&(x->err),option); if (hstio_err()) return -1;
  putEmptyDQ(fname,ever,&(x->dq),option); if (hstio_err()) return -1;
  /*        clear_err(); 		*/
  return 0;
}

static int putEmptySci(char *fname, int ever, SciHdrData *x, int option) {
  return putEmptyFloatHD(fname,"SCI",ever,x,option); }
static int putEmptyErr(char *fname, int ever, ErrHdrData *x, int option) {
  return putEmptyFloatHD(fname,"ERR",ever,x,option); }
static int putEmptyDQ(char *fname, int ever, DQHdrData *x, int option) {
  return putEmptyShortHD(fname,"DQ",ever,x,option); }

static int putEmptyFloatHD(char *fname, char *ename, int ever, FloatHdrData *x, int option)
{
  if (option == 0)
    x->iodesc = openOutputImage(fname, ename, ever, &(x->hdr),
                                x->data.tot_nx, x->data.tot_ny, FITSFLOAT);
  else if (option & Overwrite)
    x->iodesc = openUpdateImage(fname, ename, ever, &(x->hdr));
  if (hstio_err()) return -1;
  closeImage(x->iodesc);
  /*        clear_err(); 		*/
  return 0;
}

static int putEmptyShortHD(char *fname, char *ename, int ever, ShortHdrData *x, int option)
{
  if (option == 0)
    x->iodesc = openOutputImage(fname, ename, ever, &(x->hdr),
                                x->data.tot_nx, x->data.tot_ny, FITSSHORT);
  else if (option & Overwrite)
    x->iodesc = openUpdateImage(fname, ename, ever, &(x->hdr));
  if (hstio_err()) return -1;
  closeImage(x->iodesc);
  /*        clear_err(); 		*/
  return 0;
}
#endif


void InitDthTrl (char *inlist, char *output) {

	extern int status;

	IRAFPointer tpin;
	int  n, nfiles;

	char *trl_in; 	/* trailer filename for input */
	int  trl_len;
	char trl_out[ACS_LINE+1]; 	/* output trailer filename */
  char input[ACS_FNAME];  /* Name of image in list */
  char out_name[ACS_FNAME];

	char *isuffix[]={"_sfl", "_crj", "_crc", "_flt", "_flc"};
	char *osuffix[]={"_drz", "_drz", "_drc", "_drz", "_drc"};
	char *trlsuffix[]={"", "", "", "", ""};
	int nsuffix = 5;

	int MkOutName (char *, char **, char **, int, char *, int);
	int MkNewExtn (char *, char *);
	void WhichError (int);

	/* Allocate space for trailer file input list */
	trl_in = realloc(NULL, (ACS_LINE+1));
	trl_len = ACS_LINE+1;

	/* Make TRL filenames */
	trl_in[0] = '\0';
	trl_out[0] = '\0';
  out_name[0] = '\0';

	tpin = c_imtopen(inlist);
	nfiles = c_imtlen(tpin);

  for (n = 0; n < nfiles; ++n) {
    c_imtgetim (tpin, input, ACS_FNAME);
		/* Start by stripping off suffix from input/output filenames */
		if (MkOutName (input, isuffix, trlsuffix, nsuffix, out_name, ACS_LINE)) {
			WhichError (status);
			sprintf (MsgText, "Couldn't determine trailer filename for %s", input);
			trlmessage (MsgText);
		}

		/* Now, convert trailer filename extensions from '.fits' to '.trl' */
		if (MkNewExtn (out_name, TRL_EXTN) ) {
			sprintf(MsgText, "Error creating input trailer filename %s", out_name);
			trlerror (MsgText);
			WhichError (status);
		}

    if ( (strlen(out_name) + strlen(trl_in) + 1) >= trl_len) {
      trl_len += ACS_LINE;
      trl_in = realloc (trl_in, trl_len);
    }

    /* Append each filename to create list of input trailer files */
    strcat(trl_in, out_name);
    /* But don't put a comma after the last filename */
    if (n < (nfiles-1)) strcat (trl_in, ",");
    /* Reset value for the output filename for the next image */
    out_name[0] = '\0';
  }

	if (MkOutName (output, osuffix, trlsuffix, nsuffix, trl_out, ACS_LINE)) {
		WhichError (status);
		sprintf (MsgText, "Couldn't create trailer filename for %s", output);
		trlmessage (MsgText);
	}

	if (MkNewExtn (trl_out, TRL_EXTN) ) {
		sprintf(MsgText, "Error creating output trailer filename %s", trl_in);
		trlerror (MsgText);
		WhichError (status);
	}

	/* Sets up temp trailer file for output and copies input
   trailer file into it.
   */
	InitTrlFile (trl_in, trl_out);

  /* Deallocate memory */
  free(trl_in);
  c_imtclose (tpin);
}
