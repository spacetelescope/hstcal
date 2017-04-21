# include <stdio.h>
# include <string.h>
# include "xtables.h"
# include "c_iraf.h"	/* clear_cvoserr */
# include "wf3.h"
# include "hstcalerr.h"

/* This routine opens the reference table header and gets pedigree, descrip,
   and filetype. If the input name is null (or the first character is a blank),
   or if the file cannot be opened, the file will be flagged as not
   existing, but the error return code will be zero (OK).
   It is not an error if pedigree and/or descrip is not found in the header,
   but it is an error if filetype is not found.
   If pedigree is present, if the first five letters are DUMMY, then
   goodPedigree will be set to zero; otherwise goodPedigree will be set
   to one.

   Note that for tables the pedigree, descrip, and filetype keywords are 
   read from the primary header of the FITS file.

   Phil Hodge, 1997 Nov 13:
	Don't print error message if the table can't be opened.

   Phil Hodge, 1998 Jan 14:
	Use GotFileName instead of explicitly checking ref->name[0];
	initialize ref->goodPedigree to GOOD_PEDIGREE.

   Warren Hack, 1998 May 26:
	Revised for ACS

   Howard Bushouse, 2009 Jan 08:
	Modified to read keywords from the primary header of the file and
        to also retrieve the FILETYPE keyword from the header.

   Howard Bushouse, 2010 Jun 25:
	Reset status to zero if table can't be opened, before returning.
*/

int TabPedigree (RefTab *ref) {

	extern int status;

	IRAFPointer tp;		/* for the reference table */
	int GotFileName (char *);
	char filename[SZ_LINE+1];

	ref->goodPedigree = GOOD_PEDIGREE;	/* initial value */

	if (!GotFileName (ref->name)) {
	    ref->exists = EXISTS_NO;
	    return (status);
	}

	/* Open the reference table. */
	strcpy (filename, ref->name);
	strcat (filename, "[0]");
	tp = c_tbtopn (filename, IRAF_READ_ONLY, 0);
	if ( (status = c_iraferr()) ) {
	    ref->exists = EXISTS_NO;
	    clear_cvoserr();
	    status = 0;
	    return (status);
	}
	ref->exists = EXISTS_YES;

	/* Get pedigree and descrip.  If either or both are missing,
	   that's not an error in this case.
	*/
	c_tbhgtt (tp, "PEDIGREE", ref->pedigree, SZ_FITS_REC);
	if ( (status = c_iraferr())) {
	    ref->pedigree[0] = '\0';
	    clear_cvoserr();
	    status = 0;
	}

	c_tbhgtt (tp, "DESCRIP", ref->descrip, SZ_FITS_REC);
	if ( (status = c_iraferr())) {
	    ref->descrip[0] = '\0';
	    clear_cvoserr();
	    status = 0;
	}

	/* Is this a dummy reference file? */
	if (strncmp (ref->pedigree, "DUMMY", 5) == 0)
	    ref->goodPedigree = DUMMY_PEDIGREE;
	else
	    ref->goodPedigree = GOOD_PEDIGREE;

	/* Get filetype. Clear the error if the keyword is not found, so
	   that processing of other files can continue, but print an
	   error message. 
	*/
	c_tbhgtt (tp, "FILETYPE", ref->type, SZ_FITS_REC);
	if ( (status = c_iraferr()) ) {
	    ref->type[0] = '\0';
	    clear_cvoserr();
	    status = 0;
	    /*sprintf (MsgText, "FILETYPE keyword not found in %s.", ref->name);
	    trlwarn (MsgText);*/
	}

	/* Done with this table for the time being. */
	c_tbtclo (tp);

	return (status);
}
