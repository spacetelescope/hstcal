# include <stdlib.h>
# include <string.h>
# include "ctables.h"

# define EXTRA_COL_SPACE 100

void tb_save_info (IRAFPointer tp, char *tablename, char *fullname,
		fitsfile *fptr, int hdunum, int *status) {

	TableDescr *t_descr;
	ColumnDescr *c_descr;
	TblInfoPtr *tptr;
	IRAFPointer cp;
	char *colname;		/* column name */
	char cn_buffer[81];	/* column number as a string */
	int colnum;		/* column number */
	int hdutype, ncols;
	int typecode;
	long nrows;
	long repeat;		/* number of elements in column array */
	long width;
	int i;

	*status = 0;

	/* fits_get_hdu_type = ffghdt */
	fits_get_hdu_type (fptr, &hdutype, status);
	if (hdutype != ASCII_TBL && hdutype != BINARY_TBL) {
	    *status = ERR_NOT_A_TABLE;
	    return;
	}

	t_descr = getTableDescr (tp);
	strcpy (t_descr->tablename, tablename);
	strcpy (t_descr->filename, fullname);
	t_descr->brackets = strpbrk (fullname, "[");
	t_descr->fptr = fptr;
	t_descr->hdunum = hdunum;
	t_descr->hdutype = hdutype;

	/* fits_get_num_rows = ffgnrw */
	fits_get_num_rows (fptr, &nrows, status);
	t_descr->nrows = nrows;

	/* fits_get_num_cols = ffgncl*/
	fits_get_num_cols (fptr, &ncols, status);
	t_descr->ncols = ncols;
	if (*status != 0)
	    return;

	t_descr->alloc_cols = ncols + EXTRA_COL_SPACE;
	t_descr->columns = (IRAFPointer *)calloc (t_descr->alloc_cols,
			sizeof(IRAFPointer));

	colname = (char *)calloc (SZ_FITS_STR, sizeof(char));
	for (i = 0;  i < ncols;  i++) {

	    cp = init_cp();
	    c_descr = getColumnDescr (cp);

	    sprintf (cn_buffer, "%d", i+1);	/* one indexed */
	    /* fits_get_colname = ffgcnn */
	    fits_get_colname (fptr, CASEINSEN, cn_buffer, colname,
			&colnum, status);
	    strcpy (c_descr->name, colname);

	    /* fits_get_eqcoltype = ffeqty */
	    fits_get_eqcoltype (fptr, colnum, &typecode, &repeat,
			&width, status);
	    c_descr->dtype = typecode;
	    c_descr->repeat = repeat;
	    t_descr->columns[i] = cp;
	}

	return;
}
