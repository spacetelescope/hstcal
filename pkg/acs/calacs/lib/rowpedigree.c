# include <stdio.h>
# include <string.h>
# include "xtables.h"
# include "acs.h"
# include "err.h"

/* This routine gets pedigree and descrip from the current table row.

   The pedigree and descrip columns need not be defined.  If not, this
   is flagged by their column pointers being zero, in which case no change
   will be made to the ref struct for this table.

   If the pedigree column is present, the value in the current row
   will replace any value previously gotten from the header.  As with
   the value from the header, if the first five letters of pedigree are
   DUMMY, then goodPedigree will be set to DUMMY_PEDIGREE; if not,
   goodPedigree will not be changed from its previous value (either as
   initialized or as gotten from the header).

   If the descrip column is present, the value is assigned to a second
   string, descrip2, rather than overwriting the one from the header.
   If the column is not found, descrip2 will be set to null.
*/

int RowPedigree (RefTab *ref, int row,
	IRAFPointer tp, IRAFPointer cp_pedigree, IRAFPointer cp_descrip) {

	extern int status;

	/* Get pedigree and descrip.  If either or both are missing,
	   that's not an error in this case.
	*/
	if (cp_pedigree > 0) {
	    c_tbegtt (tp, cp_pedigree, row, ref->pedigree, ACS_FITS_REC);
	    if (c_iraferr())
		return (status = TABLE_ERROR);
	    /* Is this row flagged as dummy? */
	    if (strncmp (ref->pedigree, "DUMMY", 5) == 0)
		ref->goodPedigree = DUMMY_PEDIGREE;
	    else
		ref->goodPedigree = GOOD_PEDIGREE;
	}

	if (cp_descrip > 0) {
	    c_tbegtt (tp, cp_descrip, row, ref->descrip2, ACS_FITS_REC);
	    if (c_iraferr())
		return (status = TABLE_ERROR);
	} else {
	    ref->descrip2[0] = '\0';
	}


	return (status);
}
