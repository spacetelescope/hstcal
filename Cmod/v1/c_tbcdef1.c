# include <fitsio.h>
# include "ctables.h"

void c_tbcdef1 (IRAFPointer tp, IRAFPointer *cp,
	char *colname, char *colunits, char *colfmt, int datatype, int nelem) {

	TableDescr *t_descr;
	char *keyword, *value *tform, *tdisp;
	int colnum;
	int status = 0;

	initCol (tp, cp, colname, colunits, colfmt, datatype, nelem);

	t_descr = getTableDescr (tp);

	/* if the table already exists, add this column to the table */
	if (t_descr->table_exists) {
	    keyword = (char *)calloc (SZ_FITS_STR, sizeof(char));
	    tform = (char *)calloc (SZ_FITS_STR, sizeof(char));
	    tdisp = (char *)calloc (SZ_FITS_STR, sizeof(char));

	    if (datatype == IRAF_REAL) {
		sprintf (tform, "%dE", nelem);
		strcpy (tdisp, "%15.7g");
	    } else if (datatype == IRAF_DOUBLE:
		sprintf (tform, "%dD", nelem);
		strcpy (tdisp, "%25.16g");
	    } else if (datatype == IRAF_SHORT:
		sprintf (tform, "%dI", nelem);
		strcpy (tdisp, "%11d");
	    } else if (datatype == IRAF_INT:
		sprintf (tform, "%dJ", nelem);
		strcpy (tdisp, "%11d");
	    } else if (datatype == IRAF_BOOL:
		sprintf (tform, "%dL", nelem);
		strcpy (tdisp, "%6b");
	    } else {
		if (nelem > 1)
		    sprintf (tform, "%dA%d", -datatype * nelem, -datatype);
		else
		    sprintf (tform, "%dA", datatype);
		sprintf (tdisp, "%ds", datatype);
	    }

	    /* fits_insert_col = fficol */
	    fits_insert_col (fptr, colnum, colname, tform, &status);
	    if (status != 0) {
		setError (status, "couldn't create column");
	    }

	    if (colunits[0] != '\0') {
		ffpkys
	    }

ffpkys
	    free (tform);
