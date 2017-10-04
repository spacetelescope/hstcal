# include <stdio.h>
# include <string.h>
# include <ctype.h>
#include "hstcal.h"
# include "xtables.h"

# include "hstio.h"

# include "acs.h"
# include "hstcalerr.h"

/* mkNewSpt -- Create a new SPT file for the output file.

  Description:
  ------------
    New file will be based on primary header of first found input file
    given and will not have any extensions.

  Date          Author      Description
  ----          ------      -----------
  11-10-1999    W.J. Hack   Initial version
                            based on InitRejTrl and n_mkSPT from CALNICB
  11-18-2017    M.D. De La Pena  Clarified message when input *_spt.fits files not found.  
                            Generalized creation of output association SPT file.
*/

int mkNewSpt (char *in_list, char *mtype, char *output) {

/*
    arguments:
    char    *in_list            i: input filename/list to copy SPT data
    char    *mtype              i: type of exposure in association
    char    *output             o: rootname of output SPT file
*/

    extern int  status;
    IRAFPointer tpin;
    int         n;
	Hdr         header;		                /* SPT header */
	FILE        *fp;		                /* file pointer */
	IODescPtr   im;		           /* descriptor for input image */

    char        in_name[CHAR_FNAME_LENGTH+1];       /* filename of input data */
    char        in_spt[CHAR_FNAME_LENGTH+1];        /* filename of SPT source data */
    char        out_spt[CHAR_FNAME_LENGTH+1];       /* output SPT filename */
    char        rootname[CHAR_FNAME_LENGTH+1];       /* output SPT rootname */
    char        obsnum[4];                  /* Observation number string */
    ShortHdrData  stmp;
    int         nimgs, i;
    int         nx;
	int 		extnum, nextn;


    char        *isuffix[] = {"_blv_tmp","_crj_tmp","_flt","_crj","_dth","_sfl"};
    char        *osuffix[] = { "_spt","_spt","_spt","_spt","_spt","_spt"};

    int         nsuffix = 6;                /* How many suffixes to check */

    int         MkOutName (char *, char **, char **, int, char *, int);
    int         MkName (char *, char *, char *, char *, char *, int);
    void        WhichError (int);
    int         FileExists (char *);
    int         LoadHdr (char *, Hdr *);
    int         PutKeyStr(Hdr *, char *, char *, char *);
    int         PutKeyInt (Hdr *, char *, int, char *);
    int			GetKeyInt (Hdr *, char *, int, int, int *);

    Bool        doCreateOutputSPT  = False;
    Bool        isOutputSPTCreated = False;
    int         numSPT = 0;

/* ----------------------------- Begin ------------------------------*/

    out_spt[0] = '\0';
    /*
        Create output SPT filename from output data filename.
    */
    if(MkOutName (output, isuffix, osuffix, nsuffix, out_spt, CHAR_FNAME_LENGTH)){
        sprintf (MsgText, "Couldn't create output SPT filename for %s", output);
        trlerror (MsgText);
        WhichError (status);
        return (status);
    }

    // FUTURE WORK: What if there are no input files so CALACS would not make an output file anyway?
    /* In this case it is not necessary to See if an output SPT file already exists */
    if (FileExists (out_spt)) {
        sprintf (MsgText, "Output SPT filename for %s exists - cannot overwrite file.", out_spt);
        trlerror (MsgText);
        return (status);
    }

    /*
        Now, let's get the data to be copied into this new SPT file.
    */
    /* open the input file template */
    tpin = c_imtopen (in_list);
    nimgs = c_imtlen(tpin);

    /* Loop over all images in input list, and append them to output
        SPT file.
        nimgs  = Total number of SPT filenames in input list
        numSPT = Total number of actual existing input SPT files
        extnum = Total number of SPT extensions in all existing SPT files to be written to output
    */
	extnum = 0;  
    for (i = 0; i < nimgs; i++) {
        /* Get first name from input (list?) */
        in_name[0] = '\0';
        c_imtgetim (tpin, in_name, CHAR_FNAME_LENGTH);

        in_spt[0] = '\0';
        /* Create input SPT filename to look for */
        if (MkOutName (in_name, isuffix, osuffix, nsuffix, in_spt, CHAR_FNAME_LENGTH)) {
            sprintf (MsgText, "Couldn't create input SPT filename for %s", in_name);
            trlerror (MsgText);
            WhichError (status);
            return (status);
        }
	    /* Check for existence of source/input SPT file */
        /* If there are no input SPT files, then no combined/ASN SPT file will be created */
	    if ((fp = fopen (in_spt, "rb")) == NULL) {
	        sprintf (MsgText, "Cannot find file \"%s\" - processing can proceed.  Output association SPT is comprised of any found individual SPT files.\n", in_spt);
	        trlwarn  (MsgText);
            status = ACS_OK;        /* don't abort */
	        continue;				/* try the rest of the images in list */
	    } else
	        (void)fcloseWithStatus(&fp);

        /* Keep a count of the number of existing input SPT files.  Only create the output
           SPT once an input SPT is found.
        */
        numSPT++;
        if (!isOutputSPTCreated)
            doCreateOutputSPT = True;

        /* Create Primary header of new output SPT file from first found input
            image...
        */
	    /* Read the primary header of the input SPT file */
        nextn = 0;
		if (LoadHdr (in_spt, &header) )
            	    return (status);

        /* If any of the input SPT files actually contain multiple extensions, then
           the total number of output extensions is the total number of extensions in 
           each of the existing input SPT files in the most general case.
        */
		(void) GetKeyInt (&header, "NEXTEND", USE_DEFAULT, 1, &nextn);

        if (doCreateOutputSPT) {
            doCreateOutputSPT  = False;
            isOutputSPTCreated = True;
            /* Create ROOTNAME for output SPT file */
            rootname[0] = '\0';
            if (MkName (out_spt, "_spt", "", " ", rootname, CHAR_FNAME_LENGTH)) {
                sprintf (MsgText, "Couldn't create output SPT ROOTNAME for %s", out_spt);
                trlerror (MsgText);
                WhichError (status);
                return (status);
            } else {
		        sprintf(MsgText, "Created output SPT rootname %s...\n",out_spt);
                trlmessage (MsgText);
		    }

	        /* Update the FILENAME header keyword */
	        if (PutKeyStr (&header, "FILENAME", out_spt, ""))
	            return (status = 1);

	        /* Update the ASN_MTYP header keyword */
	        if (PutKeyStr(&header, "ASN_MTYP", mtype, "Role of the Exposure in the Association"))
	            return (status = 1);

	        /* NOW, update the ROOTNAME header keyword */
	        for (n = 0; n < strlen(rootname)-1; n++)
	             rootname[n] = toupper(rootname[n]);
	        if (PutKeyStr (&header, "ROOTNAME", rootname, ""))
	            return (status = 1);

            /* Update the OBSERVTN header keyword */
            strncpy (obsnum, &rootname[6], 3); obsnum[3] = '\0';
            if (putKeyS (&header, "OBSERVTN", obsnum, ""))
                return (status = 1);

	        /* Update the NEXTEND header keyword to reflect number of input
                images.  This will be updated as necessary at the end of
				this routine.
            */
	        if (PutKeyInt (&header, "NEXTEND", nimgs, ""))
	            return (status = 1);

	        /* Write the new SPT file */
            /* Open the image; this also writes the header */
            im = openOutputImage (out_spt, "", 0, &header, 0, 0, FITSBYTE);
            if (hstio_err()) {
                trlopenerr (out_spt);
                return (status = OPEN_FAILED);
            }

            /* Close the image */
            closeImage (im);
        }

        /* Write the 1-D data */
		for (nx = 0; nx < nextn; nx++){
			extnum = extnum + 1;
        	initShortHdrData(&stmp);
        	getShortHD(in_spt, "UDL", (nx+1), &stmp);
        	putShortHD (out_spt, "UDL", extnum, &stmp, 0);
        	freeShortHdrData(&stmp);
        }

	    /* Close the header here... */
        freeHdr (&header);
    }

    /* If necessary, re-open the primary header to update the NEXTEND keyword */
    if ((extnum != numSPT) || (numSPT != nimgs)) {

        im = openUpdateImage (out_spt, "", 0, &header);
        if (hstio_err()) {
            trlopenerr (out_spt);
            return (status = OPEN_FAILED);
        }

        if (PutKeyInt (&header, "NEXTEND", extnum, ""))
            return (status = 1);

        /* Write the updated header and clean up */
        putHeader(im);
        closeImage (im);
        freeHdr(&header);

        sprintf(MsgText,"Updated output SPT file to reflect %d extensions...\n",extnum);
        trlmessage(MsgText);
    }

    c_imtclose (tpin);
	/* Successful return */
	return (status);
}
