# include <fitsio.h>
# include "ctables.h"

void c_tbtclo (IRAFPointer tp) {

/* Close a table.
argument:
IRAFPointer tp          i: table descriptor
*/
	if (!tp)
		return;

	TableDescr * tbl_descr = (TableDescr *)tp;
	fitsfile * fptr = tbl_descr->fptr; //CFITSIO pointer
	int status = 0;

	/* fits_close_file = ffclos */
	fits_close_file (fptr, &status);

	free_tp (tp);
}

void c_tbtClose (IRAFPointer * tp)
{
	c_tbtclo(*tp);
	*tp = NULL;
}
