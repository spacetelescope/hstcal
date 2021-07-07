# include <ctype.h>
# include <stdio.h>
# include <string.h>
#include "hstcal.h"
# include "xtables.h"

# include "hstio.h"

# include "acs.h"
# include "hstcalerr.h"

/* mkNewSpt -- Create a new *_spt.fits file for the output file.

  Description:
  ------------
    New file will be based on primary header of first input file
    given and will not have any extensions.

  Date          Author      Description
  ----          ------      -----------
  11-10-1999    W.J. Hack   Initial version
                            based on InitRejTrl and n_mkSPT from CALNICB

*/

int mkNewSpt (char *in_list, char *mtype, char *output) {

/*
    arguments:
    char    *in_list            i: input filename/list to copy *_spt.fits data
    char    *mtype              i: type of exposure in association
    char    *output             o: rootname of output *_spt.fits file
*/

    extern int  status;
    IRAFPointer tpin = NULL;
    int         n;
	Hdr         header;		                /* *_spt.fits header */
	FILE        *fp = NULL;		                /* file pointer */
	IODescPtr   im  = NULL;		           /* descriptor for input image */

    char        in_name[CHAR_FNAME_LENGTH+1];       /* filename of input data */
    char        in_spt[CHAR_FNAME_LENGTH+1];        /* filename of *_spt.fits source data */
    char        out_spt[CHAR_FNAME_LENGTH+1];       /* output *_spt.fits filename */
    char        rootname[CHAR_FNAME_LENGTH+1];       /* output *_spt.fits rootname */
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
    int 	GetKeyInt (Hdr *, char *, int, int, int *);

/* ----------------------------- Begin ------------------------------*/

    out_spt[0] = '\0';
    /*
        Create output *_spt.fits filename from output data filename.
    */
    if(MkOutName (output, isuffix, osuffix, nsuffix, out_spt, CHAR_FNAME_LENGTH)){
        sprintf (MsgText, "Couldn't create output *_spt.fits filename for %s", output);
        trlerror (MsgText);
        WhichError (status);
        return (status);
    }

    /* See if an output *_spt.fits file already exists */
    /* Low-level code issues a warning if the file already exists, and it cannot be deleted. */
    if (FileExists (out_spt)) {
        sprintf (MsgText, "Solution is to delete %s file or to set environment variable 'imclobber' to 'yes'.", out_spt);
        trlerror (MsgText);
        return (status);
    }

    /*
        Now, let's get the data to be copied into this new *_spt.fits file.
    */
    /* open the input file template */
    tpin = c_imtopen (in_list);
    nimgs = c_imtlen(tpin);

    /* Loop over all images in input list, and append them to output
        *_spt.fits file.
    */
	extnum = 0;
    for (i = 0; i < nimgs; i++) {
        /* Get first name from input (list?) */
        in_name[0] = '\0';
        c_imtgetim (tpin, in_name, CHAR_FNAME_LENGTH);

        in_spt[0] = '\0';
        /* Create input *_spt.fits filename */
        if (MkOutName (in_name, isuffix, osuffix, nsuffix, in_spt, CHAR_FNAME_LENGTH)) {
            sprintf (MsgText, "Couldn't create input *_spt.fits filename for %s", in_name);
            trlerror (MsgText);
            WhichError (status);
            return (status);
        }
	    /* Check for existence of source/input *_spt.fits file */
        /* If the input *_spt.fits file corresponding to the first *_raw.fits file in the list is missing, this */
        /* routine will fail as the header from the first *_spt.fits is used for the header of the output *_spt.fits file. */
        /* This routine needs to be re-worked to use the header from any found input *_spt.fits file. */
	    if ((fp = fopen (in_spt, "rb")) == NULL) {
			if (i == 0) {
	            sprintf (MsgText, "Cannot find the first input file \"%s\" which is needed to create the output file header.\n", in_spt);
	            trlerror (MsgText);
            } else {
	            sprintf (MsgText, "Cannot find input file \"%s\", but processing can proceed.  Output\nassociation *_spt.fits is comprised of any found individual *_spt.fits files.\n", in_spt);
	            trlwarn (MsgText);
            }
            status = ACS_OK;        /* don't abort */
	        continue;				/* try the rest of the images in list */
	    } else
	        (void)fcloseWithStatus(&fp);

        /* Create Primary header of new output *_spt.fits file from first input
            image...
        */
	    /* Read the primary header of the input *_spt.fits file */
        nextn = 0;
		if (LoadHdr (in_spt, &header)) {
			freeHdr (&header);
			return (status);
		}

		GetKeyInt (&header, "NEXTEND", USE_DEFAULT, 1, &nextn);

        if (i == 0) {
            /* Create ROOTNAME for output *_spt.fits file */
            rootname[0] = '\0';
            if (MkName (out_spt, "_spt", "", " ", rootname, CHAR_FNAME_LENGTH)) {
                sprintf (MsgText, "Couldn't create output *_spt.fits ROOTNAME for %s", out_spt);
                trlerror (MsgText);
                WhichError (status);
                freeHdr (&header);
                return (status);
            } else {
		        sprintf(MsgText, "Created output *_spt.fits rootname %s...\n",out_spt);
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
                images.
            */
	        if ((status = PutKeyInt (&header, "NEXTEND", nimgs, ""))) {
				freeHdr (&header);
	            return (status);
			}

            sprintf(MsgText,"Updated output *_spt.fits file to reflect %d extensions...\n",nimgs);
            trlmessage(MsgText);

	        /* Write the new *_spt.fits file */
            /* Open the image; this also writes the header */
            im = openOutputImage (out_spt, "", 0, &header, 0, 0, FITSBYTE);
            if (hstio_err()) {
                trlopenerr (out_spt);
                freeHdr (&header);
                return (status = OPEN_FAILED);
            }
            /* Close the image */
            closeImage (im);
            im = NULL;
        }

        /* Copy input *_spt.fits files into output file */
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

    c_imtclose (tpin);
    tpin = NULL;
	/* Successful return */
	return (status);
}
