# include <fitsio.h>
# include "ctables.h"

# define MAXDIM 7

void c_tbciga (IRAFPointer tp, IRAFPointer cp, int *ndim, int *axlen,
                int maxdim) {

/* For a table column that contains arrays, this function gets the
   dimension of the array and the length of each axis.
arguments:
IRAFPointer tp          i: table descriptor
IRAFPointer cp          i: column descriptor
int *ndim               o: dimension of array (1 to maxdim)
int *axlen              o: array with the length of each axis
int maxdim              i: maximum size of the 'axlen' array (may be up to 7)
*/

        fitsfile *fptr;
        TableDescr *tbl_descr;
        ColumnDescr *col_descr;
        int colnum;
        int naxis;
        long naxes[MAXDIM];
        int i;
        int status = 0;

        if (maxdim > MAXDIM) {
            setError (ERR_ARRAY_TOO_LARGE, "c_tbciga:  maxdim is too large");
            return;
        }

        tbl_descr = (TableDescr *)tp;
        col_descr = (ColumnDescr *)cp;
        colnum = col_descr->colnum;

        fptr = tbl_descr->fptr;

        /* fits_read_tdim = ffgtdm */
        fits_read_tdim (fptr, colnum, maxdim, &naxis, naxes, &status);
        if (status != 0) {
            setError (status, "c_tbciga:  couldn't get dimensions");
            return;
        }

        *ndim = naxis;
        for (i = 0;  i < naxis;  i++)
             axlen[i] = naxes[i];
}
