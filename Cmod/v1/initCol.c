# include <stdlib.h>
# include <string.h>
# include <fitsio.h>
# include "ctables.h"

void initCol (IRAFPointer tp, IRAFPointer *cp,
	char *colname, char *colunits, char *colfmt, int datatype, int nelem) {

	TableDescr *t_descr;
	ColumnDescr *c_descr;
	int ncols;		/* number of columns in the table */
	int colnum;		/* column number */
	int typecode;
	char *errmsg;

	/* Convert the IRAF data type to CFITSIO typecode for tables. */
	switch (datatype) {
	case IRAF_REAL:
	    typecode = TFLOAT;
	    break;
	case IRAF_DOUBLE:
	    typecode = TDOUBLE;
	    break;
	case IRAF_SHORT:
	    typecode = TSHORT;
	    break;
	case IRAF_INT:
	    typecode = TINT;
	    break;
	case IRAF_BOOL:
	    typecode = TLOGICAL;
	    break;
	case IRAF_CHAR:
	    typecode = TSTRING;
	    break;
	default:
	    errmsg = (char *)calloc (SZ_ERRMESS, sizeof(char));
	    sprintf (errmsg, "Data type code %d not recognized.", datatype);
	    setError (ERR_DATATYPE_UNKNOWN, errmsg);
	    free (errmsg);
	    *cp = -1;
	    return;
	}

	*cp = init_cp();
	c_descr = getColumnDescr (*cp);

	strcpy (c_descr->name, colname);
	strcpy (c_descr->tunit, colunits);
	strcpy (c_descr->tdisp, colfmt);

	c_descr->colnum = ncols + 1;		/* one indexed */
	c_descr->dtype = typecode;
	c_descr->repeat = nelem;

	t_descr = getTableDescr (tp);
	/* make sure we have space for one more column */
	columnSpace (t_descr, 1);
	ncols = t_descr->ncols;
	t_descr->columns[ncols] = *cp;
	t_descr->ncols = ncols + 1;
}
