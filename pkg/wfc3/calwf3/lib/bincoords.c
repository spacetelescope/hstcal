# include <stdio.h>
# include "hstio.h"
# include "wf3.h"		/* For USE_DEFAULT */
# include "err.h"

# define NUM_KEYWORDS 10

static int addkeyd (Hdr *, char *, double, char *);

/* This routine updates the coordinate parameters in the three extension
   headers to account for a linear transformation, a change in scale and
   a shift of the origin.  The parameters are gotten from inhdr and
   written to scihdr, errhdr, and dqhdr.  inhdr may be the same as one
   of the output headers.

   The bin factor must be positive.  The offset may be positive, negative,
   or zero.  When an output image is a binned subset of an input image
   (like blkavg), we will have bin > 1, and offset >= 0.  When an output
   image is expanded (like blkrep), we will have bin < 1, and offset can
   be <= 0.

   Note that offset is in units of the smaller pixels.  That is, when
   binning (bin > 1), offset is in units of the input pixels.  When
   unbinning (bin < 1), offset is in units of the output pixels.

   The CD matrix is mapped to the 1-D array cd in the following order:
	CD1_1, CD1_2, CD2_1, CD2_2

   Phil Hodge, 1998 Mar 9:
	Use addkeyd instead of PutKeyDbl, in order to suppress bogus warning
	messages about keywords being added to the header.

   Howard Bushouse, 2002 Mar 19:
	Created the BinCoordsIR routine as a copy of BinCoords, to use with
	WFC3 IR images.
*/

int BinCoords (Hdr *inhdr, double *block, double *offset,
		Hdr *scihdr, Hdr *errhdr, Hdr *dqhdr) {

/* arguments:
Hdr inhdr         i: header from which to get the coordinate parameters
double block[2]   i: number of input pixels for one output pixel
double offset[2]  i: offset of binned image in units of unbinned pixels
Hdr scihdr, errhdr, dqhdr   o: headers to receive modified coord parameters
*/

	extern int status;

	double ltm[2], ltv[2];
	double cd[4];		/* cd1_1, cd1_2, cd2_1, cd2_2 */
	double crpix[2];
	int status_arr[NUM_KEYWORDS];
	int i;
    
	int GetKeyDbl (Hdr *, char *, int, double, double *);
	int BinUpdate (double *, double *, double *, double *, double *,
		       double *);
        
	/* Initialize array of status values for keywords */
	for (i = 0; i < NUM_KEYWORDS; i++ ) 
	     status_arr[i] = 0;

	/* Get the coordinate parameters from the input header. */
	status_arr[0] = GetKeyDbl (inhdr, "CRPIX1", USE_DEFAULT, 0., &crpix[0]);
	status_arr[1] = GetKeyDbl (inhdr, "CRPIX2", USE_DEFAULT, 0., &crpix[1]);
	status_arr[2] = GetKeyDbl (inhdr, "CD1_1", USE_DEFAULT, 1., &cd[0]);
	status_arr[3] = GetKeyDbl (inhdr, "CD1_2", USE_DEFAULT, 0., &cd[1]);
	status_arr[4] = GetKeyDbl (inhdr, "CD2_1", USE_DEFAULT, 0., &cd[2]);
	status_arr[5] = GetKeyDbl (inhdr, "CD2_2", USE_DEFAULT, 1., &cd[3]);
	status_arr[6] = GetKeyDbl (inhdr, "LTM1_1", USE_DEFAULT, 1., &ltm[0]);
	status_arr[7] = GetKeyDbl (inhdr, "LTM2_2", USE_DEFAULT, 1., &ltm[1]);
	status_arr[8] = GetKeyDbl (inhdr, "LTV1", USE_DEFAULT, 0., &ltv[0]);
	status_arr[9] = GetKeyDbl (inhdr, "LTV2", USE_DEFAULT, 0., &ltv[1]);

	for (i = 0; i < NUM_KEYWORDS; i++) {
	     if (status_arr[i]) {
		 status = status_arr[i];
		 return (status);
	     }
	} 
        
	/* Modify the coordinate parameters. */
	if (BinUpdate (block, offset, ltm, ltv, cd, crpix))
	    return (status);

	/* Add or update the MWCS keywords in all three extensions. */
	status_arr[0] = addkeyd (scihdr, "CRPIX1", crpix[0], "X ref pixel");
	status_arr[1] = addkeyd (scihdr, "CD1_1", cd[0], "");
	status_arr[2] = addkeyd (scihdr, "CD2_1", cd[2], "");
	status_arr[3] = addkeyd (scihdr, "LTV1", ltv[0],
			"linear transformation vector, X component");
	status_arr[4] = addkeyd (scihdr, "LTM1_1", ltm[0],
			"reciprocal of sampling rate in X");

	status_arr[5] = addkeyd (scihdr, "CRPIX2", crpix[1], "Y ref pixel");
	status_arr[6] = addkeyd (scihdr, "CD1_2", cd[1], "");
	status_arr[7] = addkeyd (scihdr, "CD2_2", cd[3], "");
	status_arr[8] = addkeyd (scihdr, "LTV2", ltv[1],
			"linear transformation vector, Y component");
	status_arr[9] = addkeyd (scihdr, "LTM2_2", ltm[1],
			"reciprocal of sampling rate in Y");

	for (i = 0; i < NUM_KEYWORDS; i++) {
	     if (status_arr[i]) {
		 status = status_arr[i];
		 return (status);
	     }
	} 

	status_arr[0] = addkeyd (errhdr, "CRPIX1", crpix[0], "X ref pixel");
	status_arr[1] = addkeyd (errhdr, "CD1_1", cd[0], "");
	status_arr[2] = addkeyd (errhdr, "CD2_1", cd[2], "");
	status_arr[3] = addkeyd (errhdr, "LTV1", ltv[0],
			"linear transformation vector, X component");
	status_arr[4] = addkeyd (errhdr, "LTM1_1", ltm[0],
			"reciprocal of sampling rate in X");

	status_arr[5] = addkeyd (errhdr, "CRPIX2", crpix[1], "Y ref pixel");
	status_arr[6] = addkeyd (errhdr, "CD1_2", cd[1], "");
	status_arr[7] = addkeyd (errhdr, "CD2_2", cd[3], "");
	status_arr[8] = addkeyd (errhdr, "LTV2", ltv[1],
			"linear transformation vector, Y component");
	status_arr[9] = addkeyd (errhdr, "LTM2_2", ltm[1],
			"reciprocal of sampling rate in Y");

	for (i = 0; i < NUM_KEYWORDS; i++) {
	     if (status_arr[i]) {
		 status = status_arr[i];
		 return (status);
	     }
	} 

	status_arr[0] = addkeyd (dqhdr, "CRPIX1", crpix[0], "X ref pixel");
	status_arr[1] = addkeyd (dqhdr, "CD1_1", cd[0], "");
	status_arr[2] = addkeyd (dqhdr, "CD2_1", cd[2], "");
	status_arr[3] = addkeyd (dqhdr, "LTV1", ltv[0],
			"linear transformation vector, X component");
	status_arr[4] = addkeyd (dqhdr, "LTM1_1", ltm[0],
			"reciprocal of sampling rate in X");

	status_arr[5] = addkeyd (dqhdr, "CRPIX2", crpix[1], "Y ref pixel");
	status_arr[6] = addkeyd (dqhdr, "CD1_2", cd[1], "");
	status_arr[7] = addkeyd (dqhdr, "CD2_2", cd[3], "");
	status_arr[8] = addkeyd (dqhdr, "LTV2", ltv[1],
			"linear transformation vector, Y component");
	status_arr[9] = addkeyd (dqhdr, "LTM2_2", ltm[1],
			"reciprocal of sampling rate in Y");

	for (i = 0; i < NUM_KEYWORDS; i++) {
	     if (status_arr[i]) {
		 status = status_arr[i];
		 return (status);
	     }
	} 

	return (status);
}

int BinCoordsIR (Hdr *inhdr, double *block, double *offset,
		 Hdr *scihdr, Hdr *errhdr, Hdr *dqhdr, Hdr *smplhdr,
		 Hdr *intghdr) {

/* arguments:
Hdr inhdr         i: header from which to get the coordinate parameters
double block[2]   i: number of input pixels for one output pixel
double offset[2]  i: offset of binned image in units of unbinned pixels
Hdr scihdr, errhdr, dqhdr, smplhdr, intghdr   o: headers to receive modified coord parameters
*/

	extern int status;

	double ltm[2], ltv[2];
	double cd[4];		/* cd1_1, cd1_2, cd2_1, cd2_2 */
	double crpix[2];
	int status_arr[NUM_KEYWORDS];
	int i;
    
	int GetKeyDbl (Hdr *, char *, int, double, double *);
	int BinUpdate (double *, double *, double *, double *, double *,
		       double *);
        
	/* Initialize array of status values for keywords */
	for (i = 0; i < NUM_KEYWORDS; i++ ) 
	     status_arr[i] = 0;

	/* Get the coordinate parameters from the input header. */
	status_arr[0] = GetKeyDbl (inhdr, "CRPIX1", USE_DEFAULT, 0., &crpix[0]);
	status_arr[1] = GetKeyDbl (inhdr, "CRPIX2", USE_DEFAULT, 0., &crpix[1]);
	status_arr[2] = GetKeyDbl (inhdr, "CD1_1", USE_DEFAULT, 1., &cd[0]);
	status_arr[3] = GetKeyDbl (inhdr, "CD1_2", USE_DEFAULT, 0., &cd[1]);
	status_arr[4] = GetKeyDbl (inhdr, "CD2_1", USE_DEFAULT, 0., &cd[2]);
	status_arr[5] = GetKeyDbl (inhdr, "CD2_2", USE_DEFAULT, 1., &cd[3]);
	status_arr[6] = GetKeyDbl (inhdr, "LTM1_1", USE_DEFAULT, 1., &ltm[0]);
	status_arr[7] = GetKeyDbl (inhdr, "LTM2_2", USE_DEFAULT, 1., &ltm[1]);
	status_arr[8] = GetKeyDbl (inhdr, "LTV1", USE_DEFAULT, 0., &ltv[0]);
	status_arr[9] = GetKeyDbl (inhdr, "LTV2", USE_DEFAULT, 0., &ltv[1]);

	for (i = 0; i < NUM_KEYWORDS; i++) {
	     if (status_arr[i]) {
		 status = status_arr[i];
		 return (status);
	     }
	} 
        
	/* Modify the coordinate parameters. */
	if (BinUpdate (block, offset, ltm, ltv, cd, crpix))
	    return (status);

	/* Add or update the MWCS keywords in all three extensions. */
	status_arr[0] = addkeyd (scihdr, "CRPIX1", crpix[0], "X ref pixel");
	status_arr[1] = addkeyd (scihdr, "CD1_1", cd[0], "");
	status_arr[2] = addkeyd (scihdr, "CD2_1", cd[2], "");
	status_arr[3] = addkeyd (scihdr, "LTV1", ltv[0],
			"linear transformation vector, X component");
	status_arr[4] = addkeyd (scihdr, "LTM1_1", ltm[0],
			"reciprocal of sampling rate in X");

	status_arr[5] = addkeyd (scihdr, "CRPIX2", crpix[1], "Y ref pixel");
	status_arr[6] = addkeyd (scihdr, "CD1_2", cd[1], "");
	status_arr[7] = addkeyd (scihdr, "CD2_2", cd[3], "");
	status_arr[8] = addkeyd (scihdr, "LTV2", ltv[1],
			"linear transformation vector, Y component");
	status_arr[9] = addkeyd (scihdr, "LTM2_2", ltm[1],
			"reciprocal of sampling rate in Y");

	for (i = 0; i < NUM_KEYWORDS; i++) {
	     if (status_arr[i]) {
		 status = status_arr[i];
		 return (status);
	     }
	} 

	status_arr[0] = addkeyd (errhdr, "CRPIX1", crpix[0], "X ref pixel");
	status_arr[1] = addkeyd (errhdr, "CD1_1", cd[0], "");
	status_arr[2] = addkeyd (errhdr, "CD2_1", cd[2], "");
	status_arr[3] = addkeyd (errhdr, "LTV1", ltv[0],
			"linear transformation vector, X component");
	status_arr[4] = addkeyd (errhdr, "LTM1_1", ltm[0],
			"reciprocal of sampling rate in X");

	status_arr[5] = addkeyd (errhdr, "CRPIX2", crpix[1], "Y ref pixel");
	status_arr[6] = addkeyd (errhdr, "CD1_2", cd[1], "");
	status_arr[7] = addkeyd (errhdr, "CD2_2", cd[3], "");
	status_arr[8] = addkeyd (errhdr, "LTV2", ltv[1],
			"linear transformation vector, Y component");
	status_arr[9] = addkeyd (errhdr, "LTM2_2", ltm[1],
			"reciprocal of sampling rate in Y");

	for (i = 0; i < NUM_KEYWORDS; i++) {
	     if (status_arr[i]) {
		 status = status_arr[i];
		 return (status);
	     }
	} 

	status_arr[0] = addkeyd (dqhdr, "CRPIX1", crpix[0], "X ref pixel");
	status_arr[1] = addkeyd (dqhdr, "CD1_1", cd[0], "");
	status_arr[2] = addkeyd (dqhdr, "CD2_1", cd[2], "");
	status_arr[3] = addkeyd (dqhdr, "LTV1", ltv[0],
			"linear transformation vector, X component");
	status_arr[4] = addkeyd (dqhdr, "LTM1_1", ltm[0],
			"reciprocal of sampling rate in X");

	status_arr[5] = addkeyd (dqhdr, "CRPIX2", crpix[1], "Y ref pixel");
	status_arr[6] = addkeyd (dqhdr, "CD1_2", cd[1], "");
	status_arr[7] = addkeyd (dqhdr, "CD2_2", cd[3], "");
	status_arr[8] = addkeyd (dqhdr, "LTV2", ltv[1],
			"linear transformation vector, Y component");
	status_arr[9] = addkeyd (dqhdr, "LTM2_2", ltm[1],
			"reciprocal of sampling rate in Y");

	for (i = 0; i < NUM_KEYWORDS; i++) {
	     if (status_arr[i]) {
		 status = status_arr[i];
		 return (status);
	     }
	} 

	status_arr[0] = addkeyd (smplhdr, "CRPIX1", crpix[0], "X ref pixel");
	status_arr[1] = addkeyd (smplhdr, "CD1_1", cd[0], "");
	status_arr[2] = addkeyd (smplhdr, "CD2_1", cd[2], "");
	status_arr[3] = addkeyd (smplhdr, "LTV1", ltv[0],
			"linear transformation vector, X component");
	status_arr[4] = addkeyd (smplhdr, "LTM1_1", ltm[0],
			"reciprocal of sampling rate in X");

	status_arr[5] = addkeyd (smplhdr, "CRPIX2", crpix[1], "Y ref pixel");
	status_arr[6] = addkeyd (smplhdr, "CD1_2", cd[1], "");
	status_arr[7] = addkeyd (smplhdr, "CD2_2", cd[3], "");
	status_arr[8] = addkeyd (smplhdr, "LTV2", ltv[1],
			"linear transformation vector, Y component");
	status_arr[9] = addkeyd (smplhdr, "LTM2_2", ltm[1],
			"reciprocal of sampling rate in Y");

	for (i = 0; i < NUM_KEYWORDS; i++) {
	     if (status_arr[i]) {
		 status = status_arr[i];
		 return (status);
	     }
	} 

	status_arr[0] = addkeyd (intghdr, "CRPIX1", crpix[0], "X ref pixel");
	status_arr[1] = addkeyd (intghdr, "CD1_1", cd[0], "");
	status_arr[2] = addkeyd (intghdr, "CD2_1", cd[2], "");
	status_arr[3] = addkeyd (intghdr, "LTV1", ltv[0],
			"linear transformation vector, X component");
	status_arr[4] = addkeyd (intghdr, "LTM1_1", ltm[0],
			"reciprocal of sampling rate in X");

	status_arr[5] = addkeyd (intghdr, "CRPIX2", crpix[1], "Y ref pixel");
	status_arr[6] = addkeyd (intghdr, "CD1_2", cd[1], "");
	status_arr[7] = addkeyd (intghdr, "CD2_2", cd[3], "");
	status_arr[8] = addkeyd (intghdr, "LTV2", ltv[1],
			"linear transformation vector, Y component");
	status_arr[9] = addkeyd (intghdr, "LTM2_2", ltm[1],
			"reciprocal of sampling rate in Y");

	for (i = 0; i < NUM_KEYWORDS; i++) {
	     if (status_arr[i]) {
		 status = status_arr[i];
		 return (status);
	     }
	} 

	return (status);
}

/* This routine will update the keyword if it exists, and it will silently
   add the keyword if it doesn't exist.  This differs from PutKeyDbl only in
   that it doesn't print a warning if the keyword doesn't already exist.
*/

static int addkeyd (Hdr *hd, char *keyword, double value, char *comment) {

/* arguments:
Hdr *hd           i: pointer to header to be updated
char *keyword     i: name of keyword
double value      i: value to be updated or added
char *comment     i: comment to add, if keyword doesn't exist
*/

	extern int status;
	FitsKw key;		/* location of keyword in header */

	key = findKw (hd, keyword);
	if (key == NotFound) {
	    addDoubleKw (hd, keyword, value, comment);
	} else {
	    putDoubleKw (key, value);
	}

	if (hstio_err())
	    return (status = ERROR_RETURN);

	return (status);
}
