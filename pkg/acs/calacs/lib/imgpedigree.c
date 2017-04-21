# include <stdio.h>
# include <string.h>
# include "hstio.h"
# include "acs.h"
# include "err.h"

/* This routine opens the reference image header and get pedigree and
   descrip.  If the input name is null (or the first character is a blank),
   or if the file cannot be opened, the file will be flagged as not
   existing, but the error return code will be zero (OK).
   It is not an error if pedigree and/or descrip is not found in the header.
   If pedigree is present, if the first five letters are DUMMY, then
   goodPedigree will be set to zero; otherwise goodPedigree will be set
   to one.

   Note that for images the pedigree and descrip keywords are read from
   the primary header.

   Phil Hodge, 1998 Jan 14:
	Use GotFileName instead of explicitly checking ref->name[0];
	initialize ref->goodPedigree to GOOD_PEDIGREE.
   Warren Hack, 1998 May 26:
	Revised for ACS.
*/

int ImgPedigree (RefImage *ref) {

	extern int status;

	FitsKw key;		/* location of keyword in header */
	IODescPtr im;		/* descriptor for primary header unit */
	Hdr phdr;		/* primary header */
	int GotFileName (char *);

	initHdr (&phdr);
	ref->goodPedigree = GOOD_PEDIGREE;	/* initial value */

	if (!GotFileName (ref->name)) {
	    ref->exists = EXISTS_NO;
	    return (status);
	}

	/* Open the primary header of the reference file. */
	im = openInputImage (ref->name, "", 0);
	if (hstio_err()) {
	    ref->exists = EXISTS_NO;
	    clear_hstioerr();
        status = 0;
	    return (status);
	}
	ref->exists = EXISTS_YES;
	getHeader (im, &phdr);
	if (hstio_err())
	    return (status = HEADER_PROBLEM);

	/* Get pedigree and descrip.  If either or both are missing,
	   that's not an error in this case.
	*/
	key = findKw (&phdr, "PEDIGREE");
	if (key == NotFound) {
	    ref->pedigree[0] = '\0';
	} else {
	    getStringKw (key, ref->pedigree, ACS_FITS_REC);
	    if (hstio_err()) {
		trlerror ("Trying to get PEDIGREE.");
		return (status = HEADER_PROBLEM);
	    }
	}

	key = findKw (&phdr, "DESCRIP");
	if (key == NotFound) {
	    ref->descrip[0] = '\0';
	} else {
	    getStringKw (key, ref->descrip, ACS_FITS_REC);
	    if (hstio_err()) {
		trlerror ("Trying to get DESCRIP.");
		return (status = HEADER_PROBLEM);
	    }
	}

	/* Is this a dummy reference file? */
	if (strncmp (ref->pedigree, "DUMMY", 5) == 0)
	    ref->goodPedigree = DUMMY_PEDIGREE;	/* dummy, so pedigree is bad */
	else
	    ref->goodPedigree = GOOD_PEDIGREE;		/* pedigree is good */

	/* Done with this image for the time being. */
	closeImage (im);
	freeHdr (&phdr);

	return (status);
}
