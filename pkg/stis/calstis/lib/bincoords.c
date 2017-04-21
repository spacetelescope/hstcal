# include <stdio.h>
# include "c_iraf.h"
# include "hstio.h"

# include "stis.h"
# include "hstcalerr.h"
# include "stisdef.h"

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
	Use addkeyd instead of PutKeyD, in order to suppress bogus warning
	messages about keywords being added to the header.

   Phil Hodge, 1998 Oct 5:
	Change status value 1001 to HEADER_PROBLEM.
*/

int BinCoords (Hdr *inhdr, double *block, double *offset,
		Hdr *scihdr, Hdr *errhdr, Hdr *dqhdr) {

/* arguments:
Hdr inhdr         i: header from which to get the coordinate parameters
double block[2]   i: number of input pixels for one output pixel
double offset[2]  i: offset of binned image in units of unbinned pixels
Hdr scihdr, errhdr, dqhdr   o: headers to receive modified coord parameters
*/

	int status;

	double ltm[2], ltv[2];
	double cd[4];		/* cd1_1, cd1_2, cd2_1, cd2_2 */
	double crpix[2];
	int use_def = 1;	/* use default if missing keyword */

	/* Get the coordinate parameters from the input header. */

	status = Get_KeyD (inhdr, "CRPIX1", use_def, 0., &crpix[0]);
	status = Get_KeyD (inhdr, "CRPIX2", use_def, 0., &crpix[1]);
	status = Get_KeyD (inhdr, "CD1_1", use_def, 1., &cd[0]);
	status = Get_KeyD (inhdr, "CD1_2", use_def, 0., &cd[1]);
	status = Get_KeyD (inhdr, "CD2_1", use_def, 0., &cd[2]);
	status = Get_KeyD (inhdr, "CD2_2", use_def, 1., &cd[3]);

	status = Get_KeyD (inhdr, "LTM1_1", use_def, 1., &ltm[0]);
	status = Get_KeyD (inhdr, "LTM2_2", use_def, 1., &ltm[1]);
	status = Get_KeyD (inhdr, "LTV1", use_def, 0., &ltv[0]);
	status = Get_KeyD (inhdr, "LTV2", use_def, 0., &ltv[1]);

	if (status)
	    return (status);

	/* Modify the coordinate parameters. */
	if ((status = BinUpdate (block, offset, ltm, ltv, cd, crpix)))
	    return (status);

	/* Add or update the MWCS keywords in all three extensions. */

	status = addkeyd (scihdr, "CRPIX1", crpix[0], "X ref pixel");
	status = addkeyd (scihdr, "CD1_1", cd[0], "");
	status = addkeyd (scihdr, "CD2_1", cd[2], "");
	status = addkeyd (scihdr, "LTV1", ltv[0],
			"linear transformation vector, X component");
	status = addkeyd (scihdr, "LTM1_1", ltm[0],
			"reciprocal of sampling rate in X");

	status = addkeyd (scihdr, "CRPIX2", crpix[1], "Y ref pixel");
	status = addkeyd (scihdr, "CD1_2", cd[1], "");
	status = addkeyd (scihdr, "CD2_2", cd[3], "");
	status = addkeyd (scihdr, "LTV2", ltv[1],
			"linear transformation vector, Y component");
	status = addkeyd (scihdr, "LTM2_2", ltm[1],
			"reciprocal of sampling rate in Y");

	if (status)
	    return (status);

	status = addkeyd (errhdr, "CRPIX1", crpix[0], "X ref pixel");
	status = addkeyd (errhdr, "CD1_1", cd[0], "");
	status = addkeyd (errhdr, "CD2_1", cd[2], "");
	status = addkeyd (errhdr, "LTV1", ltv[0],
			"linear transformation vector, X component");
	status = addkeyd (errhdr, "LTM1_1", ltm[0],
			"reciprocal of sampling rate in X");

	status = addkeyd (errhdr, "CRPIX2", crpix[1], "Y ref pixel");
	status = addkeyd (errhdr, "CD1_2", cd[1], "");
	status = addkeyd (errhdr, "CD2_2", cd[3], "");
	status = addkeyd (errhdr, "LTV2", ltv[1],
			"linear transformation vector, Y component");
	status = addkeyd (errhdr, "LTM2_2", ltm[1],
			"reciprocal of sampling rate in Y");

	if (status)
	    return (status);

	status = addkeyd (dqhdr, "CRPIX1", crpix[0], "X ref pixel");
	status = addkeyd (dqhdr, "CD1_1", cd[0], "");
	status = addkeyd (dqhdr, "CD2_1", cd[2], "");
	status = addkeyd (dqhdr, "LTV1", ltv[0],
			"linear transformation vector, X component");
	status = addkeyd (dqhdr, "LTM1_1", ltm[0],
			"reciprocal of sampling rate in X");

	status = addkeyd (dqhdr, "CRPIX2", crpix[1], "Y ref pixel");
	status = addkeyd (dqhdr, "CD1_2", cd[1], "");
	status = addkeyd (dqhdr, "CD2_2", cd[3], "");
	status = addkeyd (dqhdr, "LTV2", ltv[1],
			"linear transformation vector, Y component");
	status = addkeyd (dqhdr, "LTM2_2", ltm[1],
			"reciprocal of sampling rate in Y");

	if (status)
	    return (status);

	return (0);
}

/* This routine will update the keyword if it exists, and it will silently
   add the keyword if it doesn't exist.  This differs from Put_KeyD only in
   that it doesn't print a warning if the keyword doesn't already exist.
*/

static int addkeyd (Hdr *hd, char *keyword, double value, char *comment) {

/* arguments:
Hdr *hd           i: pointer to header to be updated
char *keyword     i: name of keyword
double value      i: value to be updated or added
char *comment     i: comment to add, if keyword doesn't exist
*/

	FitsKw key;		/* location of keyword in header */

	key = findKw (hd, keyword);
	if (key == NotFound) {
	    addDoubleKw (hd, keyword, value, comment);
	} else {
	    putDoubleKw (key, value);
	}

	if (hstio_err())
	    return (HEADER_PROBLEM);

	return (0);
}
