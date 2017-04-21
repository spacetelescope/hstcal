# include <stdio.h>
# include <string.h>
# include "xtables.h"
# include "c_iraf.h"	/* clear_cvoserr */
# include "acs.h"
# include "err.h"

/* This routine opens the reference table header and gets pedigree and
   descrip.  If the input name is null (or the first character is a blank),
   or if the file cannot be opened, the file will be flagged as not
   existing, but the error return code will be zero (OK).
   It is not an error if pedigree and/or descrip is not found in the header.
   If pedigree is present, if the first five letters are DUMMY, then
   goodPedigree will be set to zero; otherwise goodPedigree will be set
   to one.

   Note that for tables the pedigree and descrip keywords are read from
   the header of the bintable extension, not the primary header.

   Phil Hodge, 1997 Nov 13:
	Don't print error message if the table can't be opened.

   Phil Hodge, 1998 Jan 14:
	Use GotFileName instead of explicitly checking ref->name[0];
	initialize ref->goodPedigree to GOOD_PEDIGREE.

   Warren Hack, 1998 May 26:
	Revised for ACS
*/

int TabPedigree (RefTab *ref) {

	extern int status;

	IRAFPointer tp;		/* for the reference table */
	int GotFileName (char *);

	ref->goodPedigree = GOOD_PEDIGREE;	/* initial value */

	if (!GotFileName (ref->name)) {
	    ref->exists = EXISTS_NO;
	    return (status);
	}

	/* Open the reference table. */
	tp = c_tbtopn (ref->name, IRAF_READ_ONLY, 0);
	if ((status = c_iraferr())) {
	    ref->exists = EXISTS_NO;
	    clear_cvoserr();
	    return (status);
	}
	ref->exists = EXISTS_YES;

	/* Get pedigree and descrip.  If either or both are missing,
	   that's not an error in this case.
	*/
	c_tbhgtt (tp, "PEDIGREE", ref->pedigree, ACS_FITS_REC);
	if ((status = c_iraferr())) {
	    ref->pedigree[0] = '\0';
	    clear_cvoserr();
        status = 0;
	}

	c_tbhgtt (tp, "DESCRIP", ref->descrip, ACS_FITS_REC);
    trlmessage("DESCRIP read in as: ");
    sprintf(MsgText, "%s", ref->descrip);
    trlmessage(MsgText);
        if ((status = c_iraferr())) {
	    ref->descrip[0] = '\0';
	    clear_cvoserr();
        status = 0;
	}

	/* Is this a dummy reference file? */
	if (strncmp (ref->pedigree, "DUMMY", 5) == 0)
	    ref->goodPedigree = DUMMY_PEDIGREE;
	else
	    ref->goodPedigree = GOOD_PEDIGREE;

	/* Done with this table for the time being. */
	c_tbtclo (tp);

	return (status);
}
