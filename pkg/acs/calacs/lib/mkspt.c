# include <stdio.h>
# include <string.h>
# include <ctype.h>
#include "hstcal.h"
# include "xtables.h"

# include "hstio.h"

# include "acs.h"
# include "hstcalerr.h"

/* mkNewSpt -- Create a new _spt.fits file for the output file.

  Description:
  ------------
    New file will be based on primary header of first found input file
    given and will not have any extensions.

  Date          Author      Description
  ----          ------      -----------
  11-10-1999    W.J. Hack   Initial version
                            based on InitRejTrl and n_mkSPT from CALNICB
  11-18-2017    M.D. De La Pena  Clarified message when input *_spt.fits files not found.  
                            Generalized creation of output association _spt.fits file.
*/

int mkNewSpt (char *in_list, char *mtype, char *output) {

/*
    arguments:
    char    *in_list            i: input filename/list to copy SPT data
    char    *mtype              i: type of exposure in association
    char    *output             o: rootname of output _spt.fits file
*/

    extern int  status;
    IRAFPointer tpin = NULL;
    int         n;
	Hdr         header;          /* _spt.fits header */
	FILE        *fp = NULL;      /* file pointer */
	IODescPtr   im = NULL;       /* descriptor for output image */
	IODescPtr   imUpdate = NULL; /* new descriptor if output image has to be updated */

    char        in_name[CHAR_FNAME_LENGTH+1];       /* filename of input data */
    char        in_spt[CHAR_FNAME_LENGTH+1];        /* filename of _spt.fits source data */
    char        out_spt[CHAR_FNAME_LENGTH+1];       /* output _spt.fits filename */
    char        rootname[CHAR_FNAME_LENGTH+1];      /* output _spt.fits rootname */
    char        obsnum[4];                          /* Observation number string */
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
        Create output _spt.fits filename from output data filename.
    */
    if(MkOutName (output, isuffix, osuffix, nsuffix, out_spt, CHAR_FNAME_LENGTH)){
        sprintf (MsgText, "Couldn't create output _spt.fits filename for %s", output);
        trlerror (MsgText);
        WhichError (status);
        return (status);
    }

    // FUTURE WORK: What if there are no input files so CALACS would not make an output file anyway?
    /* In this case it is not necessary to see if an output _spt.fits file already exists */
    if (FileExists (out_spt)) {
        sprintf (MsgText, "Output _spt.fits filename for %s exists - cannot overwrite file.", out_spt);
        trlerror (MsgText);
        return (status);
    }

    /*
        Now, let's get the data to be copied into this new _spt.fits file.
    */
    /* open the input file template */
    tpin = c_imtopen (in_list);
    nimgs = c_imtlen(tpin);

    /* Loop over all images in input list, and append them to output
        _spt.fits file.
        nimgs  = Total number of _spt.fits filenames in input list
        numSPT = Total number of actual existing input _spt.fits files
        extnum = Total number of SPT extensions in all existing _spt.fits files to be written to output
    */
	extnum = 0;  
    for (i = 0; i < nimgs; i++) {
        /* Get first name from input (list?) */
        in_name[0] = '\0';
        c_imtgetim (tpin, in_name, CHAR_FNAME_LENGTH);

        in_spt[0] = '\0';
        /* Create input _spt.fits filename */
        if (MkOutName (in_name, isuffix, osuffix, nsuffix, in_spt, CHAR_FNAME_LENGTH)) {
            sprintf (MsgText, "Couldn't create input _spt.fits filename for %s", in_name);
            trlerror (MsgText);
            WhichError (status);
            return (status);
        }

        /* Check for existence of source/input _spt.fits file */
        /* If there are no input _spt.fits files, then no combined/ASN _spt.fits file will be created */
	    if ((fp = fopen (in_spt, "rb")) == NULL) {
	        sprintf (MsgText, "Cannot find file \"%s\" - processing can proceed.  The\noutput association _spt.fits is comprised of any found individual _spt.fits files.\n", in_spt);
	        trlwarn  (MsgText);
            status = ACS_OK;        /* don't abort */
	        continue;				/* try the rest of the images in list */
	    } else
	        (void)fcloseWithStatus(&fp);

        /* Keep a count of the number of existing input _spt.fits files.  Only create the output
           _spt.fits once an input _spt.fits is found.
        */
        numSPT++;
        if (!isOutputSPTCreated)
            doCreateOutputSPT = True;

        /* Create Primary header of new output _spt.fits file from first found input
            image...
        */
	    /* Read the primary header of the input _spt.fits file */
        nextn = 0;
        if (LoadHdr (in_spt, &header)) {
            freeHdr (&header);
            return (status);
        }

        /* If any of the input _spt.fits files actually contain multiple extensions, then
           the total number of output extensions is the total number of extensions in 
           each of the existing input _spt.fits files in the most general case.
        */
		(void) GetKeyInt (&header, "NEXTEND", USE_DEFAULT, 1, &nextn);

        if (doCreateOutputSPT) {
            doCreateOutputSPT  = False;
            isOutputSPTCreated = True;
            /* Create ROOTNAME for output _spt.fits file */
            rootname[0] = '\0';
            if (MkName (out_spt, "_spt", "", " ", rootname, CHAR_FNAME_LENGTH)) {
                sprintf (MsgText, "Couldn't create output _spt.fits ROOTNAME for %s", out_spt);
                trlerror (MsgText);
                freeHdr (&header);
                WhichError (status);
                return (status);
            } else {
                sprintf(MsgText, "Created output _spt.fits rootname %s...\n",out_spt);
                trlmessage (MsgText);
		    }

            /* Update the FILENAME header keyword */
            if ((status = PutKeyStr (&header, "FILENAME", out_spt, ""))) {
                freeHdr (&header);
                return (status);
             }

            /* Update the ASN_MTYP header keyword */
            if ((status = PutKeyStr(&header, "ASN_MTYP", mtype, "Role of the Exposure in the Association"))) {
                freeHdr (&header);
                return (status);
            }

            /* NOW, update the ROOTNAME header keyword */
            for (n = 0; n < strlen(rootname)-1; n++)
                rootname[n] = toupper(rootname[n]);
            if ((status = PutKeyStr (&header, "ROOTNAME", rootname, ""))) {
                freeHdr (&header);
                return (status);
            }

            /* Update the OBSERVTN header keyword */
            strncpy (obsnum, &rootname[6], 3); obsnum[3] = '\0';
            if ((status = putKeyS (&header, "OBSERVTN", obsnum, ""))) {
                freeHdr (&header);
                return (status);
            }

            /* Update the NEXTEND header keyword to reflect number of input
               images.  This will be updated as necessary at the end of this routine.
            */
            if ((status = PutKeyInt (&header, "NEXTEND", nimgs, ""))) {
                freeHdr (&header);
                return (status);
            }

            /* Write the new _spt.fits file */
            /* Open the image; this also writes the header */
            im = openOutputImage (out_spt, "", 0, &header, 0, 0, FITSBYTE);
            if (!im || hstio_err()) {
                trlopenerr (out_spt);
                freeHdr (&header);
                return (status = OPEN_FAILED);
            }

            /* Close the image */
            closeImage (im);
            im = NULL;
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

        imUpdate = openUpdateImage (out_spt, "", 0, &header);
        if (!imUpdate || hstio_err()) {
            trlopenerr (out_spt);
            freeHdr (&header);
            return (status = OPEN_FAILED);
        }

        if ((status = PutKeyInt (&header, "NEXTEND", extnum, ""))) {
            closeImage (imUpdate);
            imUpdate = NULL;
            freeHdr (&header);
            return (status);
        }

        /* Write the updated header and clean up */
        if ((status = putHeader(imUpdate))) {
            closeImage (imUpdate);
            imUpdate = NULL;
            freeHdr (&header);
            return (status);
        }

        closeImage (imUpdate);
        imUpdate = NULL;
        freeHdr (&header);

        sprintf(MsgText,"Updated output _spt.fits file to reflect %d extensions...\n",extnum);
        trlmessage(MsgText);
    }

    c_imtclose (tpin);
	/* Successful return */
	return (status);
}

