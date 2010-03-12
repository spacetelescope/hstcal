# include <stdlib.h>
# include <string.h>
# include "ctables.h"

# define EXTRA_COL_SPACE 100

void tbSaveInfo (IRAFPointer tp, char *tablename, char *fullname,
		fitsfile *fptr, int hdunum, int *status) {

	TableDescr *t_descr;
	ColumnDescr *c_descr;
	IRAFPointer cp;
	char value[SZ_FITS_STR];	/* column name, units, print format */
	char comment[SZ_FITS_STR];	/* comment field in a FITS keyword */
	char keyword[SZ_FITS_STR];	/* keyword name, or column number */
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
	strcpy (t_descr->brackets, strpbrk (fullname, "["));
	t_descr->fptr = fptr;
	t_descr->hdunum = hdunum;
	t_descr->hdutype = hdutype;
	t_descr->table_exists = 1;	/* the table does exist */

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

	/* Get and save info for each column in the table. */
	for (i = 0;  i < ncols;  i++) {

	    cp = init_cp();
	    c_descr = getColumnDescr (cp);

	    /* get column name */
	    sprintf (keyword, "TTYPE%d", i+1);	/* one indexed column number */
	    /* fits_read_keyword = ffgkey */
	    fits_read_keyword (fptr, keyword, value, comment, status);
	    if (*status == 0) {
		strcpy (c_descr->name, value);
	    } else {
		c_descr->name[0] = '\0';
		status = 0;
	    }

	    /* get column units */
	    sprintf (keyword, "TUNIT%d", i+1);	/* keyword for units */
	    fits_read_keyword (fptr, keyword, value, comment, status);
	    if (*status == 0) {
		strcpy (c_descr->tunit, value);
	    } else {
		c_descr->tunit[0] = '\0';
		status = 0;
	    }

	    /* get display format */
	    sprintf (keyword, "TDISP%d", i+1);	/* keyword for print format */
	    fits_read_keyword (fptr, keyword, value, comment, status);
	    if (*status == 0) {
		strcpy (c_descr->tdisp, value);
	    } else {
		c_descr->tdisp[0] = '\0';
		status = 0;
	    }

	    /* fits_get_eqcoltype = ffeqty */
	    fits_get_eqcoltype (fptr, colnum, &typecode, &repeat,
			&width, status);
	    c_descr->colnum = i + 1;		/* one indexed */
	    c_descr->dtype = typecode;
	    c_descr->repeat = repeat;
	    t_descr->columns[i] = cp;
	}

	return;
}
