# include <fitsio.h>
# include "ctables.h"

void tbSaveInfo (IRAFPointer tp, char *tablename, char *fullname,
		fitsfile *fptr, int hdunum, int *status);

IRAFPointer c_tbtopn (char *tablename, int iomode, IRAFPointer template) {

	IRAFPointer tp;
	TableDescr *t_descr;
	char *fullname;		/* name with environment variable resolved */
	char *nameonly;		/* without trailing expression in brackets */
	char *extname;		/* FITS keyword for extension name */
	int extver;		/* FITS keyword for version number */
	int hdunum;		/* HDU number (1 for primary HDU) */
	int hdutype;		/* type of HDU */
	fitsfile *fptr;		/* CFITSIO pointer */
	int cf_iomode;		/* CFITSIO I/O mode */
	int status = 0;

	tp = init_tp();

	/* expand environment variable (if any) in tablename */
	fullname = expandfn (tablename);

	if (iomode == IRAF_READ_ONLY || iomode == IRAF_READ_WRITE) {
	    if (iomode == IRAF_READ_ONLY)
		cf_iomode = READONLY;
	    else
		cf_iomode = READWRITE;

	    /* fits_open_file = ffopen */
	    fits_open_file (&fptr, fullname, cf_iomode, &status);
	    if (status != 0) {
		setError (status, "");
		free_tp (tp);
		return -1;
	    }

	    /* fits_get_hdu_num = ffghdn */
	    fits_get_hdu_num (fptr, &hdunum);	/* primary HDU is hdunum = 1 */
	    if (hdunum == 1) {
		hdunum = 2;
		/* fits_movabs_hdu = ffmahd */
		fits_movabs_hdu (fptr, hdunum, &hdutype, &status);
		if (status != 0) {
		    setError (status, "");
		    free_tp (tp);
		    return -1;
		}
	    }
	    tbSaveInfo (tp, tablename, fullname, fptr, hdunum, &status);
	    if (status != 0) {
		setError (status, "");
		free_tp (tp);
		return -1;
	    }

	} else {
	    /* fits_create_file = ffinit */
	    fits_create_file (&fptr, fullname, &status);
	    if (status != 0) {
		setError (status, "");
		free_tp (tp);
		return -1;
	    }
	    /* xxx check for NEW_COPY */
	    t_descr = getTableDescr (tp);
	    /* the table hasn't actually been created yet */
	    t_descr->table_exists = 0;
	}

	return tp;
}
