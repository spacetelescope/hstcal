# include <fitsio.h>
# include "ctables.h"

void c_tbtclo (IRAFPointer tp) {

	TableDescr *t_descr;
	fitsfile *fptr;		/* CFITSIO pointer */
	int status = 0;

	t_descr = getTableDescr (tp);
	fptr = t_descr->fptr;

	/* fits_close_file = ffclos */
	fits_close_file (fptr, &status);

	free_tp (tp);
}
