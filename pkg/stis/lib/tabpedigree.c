# include <stdio.h>
# include <stdlib.h>
# include <string.h>

# include "c_iraf.h"
# include "hstio.h"
# include "xtables.h"

# include "stis.h"
# include "hstcalerr.h"
# include "stisdef.h"

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
   
   Robert Jedrzejewski, 24 Feb 2017:
        Open the primary HDU (instead of the table HDU), in order to get
 	keywords PEDIGREE and DESCRIP.
*/

int TabPedigree (RefTab *ref) {

	IRAFPointer tp;		/* for the reference table */
        char *primary_hdu;      /* name including trailing "[0]" */

	ref->goodPedigree = GOOD_PEDIGREE;	/* initial value */

	if (!GotFileName (ref->name)) {
	    ref->exists = EXISTS_NO;
	    return (0);
	}

	/* Append "[0]" to the file name, for the primary HDU. */
	primary_hdu = calloc(strlen(ref->name) + 11, sizeof(char));
	if (primary_hdu == NULL) {
	    printf("ERROR    Out of memory.\n");
	    return(-1);
	}
	strcpy(primary_hdu, ref->name);
	strcat(primary_hdu, "[0]");

	/* Open the primary HDU of the reference table.  */
	tp = c_tbtopn(primary_hdu, IRAF_READ_ONLY, 0);
	if (c_iraferr()) {
	    ref->exists = EXISTS_NO;
	    clear_cvoserr();
	    free(primary_hdu);
	    return (0);
	}
	ref->exists = EXISTS_YES;

	/* Get pedigree and descrip.  If either or both are missing,
	   that's not an error in this case.
	*/
	c_tbhgtt (tp, "PEDIGREE", ref->pedigree, STIS_FITS_REC);
	if (c_iraferr()) {
	    ref->pedigree[0] = '\0';
	    clear_cvoserr();
	}

	c_tbhgtt (tp, "DESCRIP", ref->descrip, STIS_FITS_REC);
	if (c_iraferr()) {
	    ref->descrip[0] = '\0';
	    clear_cvoserr();
	}

	/* Is this a dummy reference file? */
	if (strncmp (ref->pedigree, "DUMMY", 5) == 0)
	    ref->goodPedigree = DUMMY_PEDIGREE;
	else
	    ref->goodPedigree = GOOD_PEDIGREE;

	/* Done with this table for the time being. */
	c_tbtclo (tp);

	free(primary_hdu);

	return (0);
}
