# include <ctype.h>
# include <stdio.h>
# include <string.h>
# include "xtables.h"

# include "hstio.h"

# include "acs.h"
# include "hstcalerr.h"

/* mkNewSpt -- Create a new SPT file for the output file.

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

    char        in_name[ACS_FNAME+1];       /* filename of input data */
    char        in_spt[ACS_FNAME+1];        /* filename of SPT source data */
    char        out_spt[ACS_FNAME+1];       /* output SPT filename */
    char        rootname[ACS_FNAME+1];       /* output SPT rootname */
    char        obsnum[4];                  /* Observation number string */
    ShortHdrData  stmp;
    int         nimgs, i;
    int         nx;
	int 		extnum,nextn;


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
        Create output SPT filename from output data filename.
    */
    if(MkOutName (output, isuffix, osuffix, nsuffix, out_spt, ACS_FNAME)){
        sprintf (MsgText, "Couldn't create output SPT filename for %s", output);
        trlerror (MsgText);
        WhichError (status);
        return (status);
    }

    /* See if an output SPT file already exists */
    if (FileExists (out_spt)) {
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
    */
	extnum = 0;
    for (i = 0; i < nimgs; i++) {
        /* Get first name from input (list?) */
        in_name[0] = '\0';
        c_imtgetim (tpin, in_name, ACS_FNAME);

        in_spt[0] = '\0';
        /* Create input SPT filename to look for */
        if (MkOutName (in_name, isuffix, osuffix, nsuffix, in_spt, ACS_FNAME)) {
            sprintf (MsgText, "Couldn't create input SPT filename for %s", in_name);
            trlerror (MsgText);
            WhichError (status);
            return (status);
        }
	    /* Check for existence of source/input SPT file */
	    if ((fp = fopen (in_spt, "rb")) == NULL) {
	        sprintf (MsgText, "Can't find input file \"%s\"", in_spt);
	        trlwarn  (MsgText);
            status = ACS_OK;        /* don't abort */
	        continue;				/* try the rest of the images in list */
	    } else
	        fclose (fp);

        /* Create Primary header of new output SPT file from first input
            image...
        */
	    /* Read the primary header of the input SPT file */
        nextn = 0;
		if (LoadHdr (in_spt, &header) )
            	    return (status);
		if (GetKeyInt (&header, "NEXTEND", USE_DEFAULT, 1, &nextn)){
			nextn = 1;
		}
        if (i == 0) {
            /* Create ROOTNAME for output SPT file */
            rootname[0] = '\0';
            if (MkName (out_spt, "_spt", "", " ", rootname, ACS_FNAME)) {
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
                images.
            */
	        if (PutKeyInt (&header, "NEXTEND", nimgs, ""))
	            return (status = 1);

            sprintf(MsgText,"Updated output SPT file to reflect %d extensions...\n",nimgs);
            trlmessage(MsgText);

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

        /* Uncomment this section to copy input SPT files into output
            when HSTIO is fixed to work with 1-D data... */
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
	/* Successful return */
	return (status);
}
