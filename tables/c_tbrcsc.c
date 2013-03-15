# include <stdlib.h>
# include <string.h>
# include <fitsio.h>
# include "ctables.h"

void c_tbrcsc (IRAFPointer itp, IRAFPointer otp,
                IRAFPointer *icp, IRAFPointer *ocp,
                int irow, int orow, int ncols) {

/* Copy values from one or more columns in a row in one table to a row
   in another (or the same) table.
arguments:
IRAFPointer itp         i: descriptor for input table
IRAFPointer otp         i: descriptor for output table
IRAFPointer *icp        i: array of column descriptors in input table
IRAFPointer *ocp        i: array of column descriptors in output table
int irow                i: row number (one indexed) in input table
int orow                i: row number (one indexed) in output table
int ncols               i: number of columns to be copied (this is the length
                           of arrays 'icp' and 'ocp')
*/

        TableDescr *itbl_descr, *otbl_descr;
        ColumnDescr *icol_descr, *ocol_descr;

        /* buffers for copying data */
        double dx, *d_array;
        float rx, *r_array;
        int ix, *i_array;
        short sx, *s_array;
        Bool bx, *b_array;
        char *tx, **t_array;

        int j;
        int datatype, nelem, width;
        int status = 0;

        itbl_descr = (TableDescr *)itp;
        otbl_descr = (TableDescr *)otp;

        for (j = 0;  j < ncols;  j++) {
            icol_descr = (ColumnDescr *)icp[j];
            ocol_descr = (ColumnDescr *)ocp[j];
            datatype = c_tbcigi (icp[j], TBL_COL_DATATYPE);
            nelem = c_tbcigi (icp[j], TBL_COL_LENDATA);
            if (nelem > 1) {
                if (datatype == IRAF_REAL) {
                    r_array = (float *)malloc (nelem * sizeof(float));
                    if (r_array == NULL) {
                        setError (ERR_OUT_OF_MEMORY, "c_tbrcsc");
                        return;
                    }
                    c_tbagtr (itp, icp[j], irow, r_array, 1, nelem);
                    c_tbaptr (otp, ocp[j], orow, r_array, 1, nelem);
                    free (r_array);
                } else if (datatype == IRAF_DOUBLE) {
                    d_array = (double *)malloc (nelem * sizeof(double));
                    if (d_array == NULL) {
                        setError (ERR_OUT_OF_MEMORY, "c_tbrcsc");
                        return;
                    }
                    c_tbagtd (itp, icp[j], irow, d_array, 1, nelem);
                    c_tbaptd (otp, ocp[j], orow, d_array, 1, nelem);
                    free (d_array);
                } else if (datatype == IRAF_SHORT) {
                    s_array = (short *)malloc (nelem * sizeof(short));
                    if (s_array == NULL) {
                        setError (ERR_OUT_OF_MEMORY, "c_tbrcsc");
                        return;
                    }
                    c_tbagts (itp, icp[j], irow, s_array, 1, nelem);
                    c_tbapts (otp, ocp[j], orow, s_array, 1, nelem);
                    free (s_array);
                } else if (datatype == IRAF_INT) {
                    i_array = (int *)malloc (nelem * sizeof(int));
                    if (i_array == NULL) {
                        setError (ERR_OUT_OF_MEMORY, "c_tbrcsc");
                        return;
                    }
                    c_tbagti (itp, icp[j], irow, i_array, 1, nelem);
                    c_tbapti (otp, ocp[j], orow, i_array, 1, nelem);
                    free (i_array);
                } else if (datatype == IRAF_BOOL) {
                    b_array = (Bool *)malloc (nelem * sizeof(Bool));
                    if (b_array == NULL) {
                        setError (ERR_OUT_OF_MEMORY, "c_tbrcsc");
                        return;
                    }
                    c_tbagtb (itp, icp[j], irow, b_array, 1, nelem);
                    c_tbaptb (otp, ocp[j], orow, b_array, 1, nelem);
                    free (b_array);
                } else if (datatype < 0) {
                    int i;
                    width = icol_descr->width;
                    t_array = (char **)calloc (nelem, sizeof(char *));
                    if (t_array == NULL) {
                        setError (ERR_OUT_OF_MEMORY, "c_tbrcsc");
                        return;
                    }
                    for (i = 0;  i < nelem;  i++) {
                        t_array[i] = (char *)calloc (width+1, sizeof(char));
                        if (t_array[i] == NULL) {
                            setError (ERR_OUT_OF_MEMORY, "c_tbrcsc");
                            return;
                        }
                    }
                    c_tbagtt (itp, icp[j], irow, t_array, 1, nelem, width);
                    c_tbaptt (otp, ocp[j], orow, t_array, width, 1, nelem);
                    for (i = 0;  i < nelem;  i++)
                        free (t_array[i]);
                    free (t_array);
                } else {
                    setError (ERR_DATATYPE_UNKNOWN,
                        "c_tbrcsc:  data type not supported");
                    return;
                }
            } else {
                if (datatype == IRAF_REAL) {
                    c_tbegtr (itp, icp[j], irow, &rx);
                    c_tbeptr (otp, ocp[j], orow, rx);
                } else if (datatype == IRAF_DOUBLE) {
                    c_tbegtd (itp, icp[j], irow, &dx);
                    c_tbeptd (otp, ocp[j], orow, dx);
                } else if (datatype == IRAF_SHORT) {
                    c_tbegts (itp, icp[j], irow, &sx);
                    c_tbepts (otp, ocp[j], orow, sx);
                } else if (datatype == IRAF_INT) {
                    c_tbegti (itp, icp[j], irow, &ix);
                    c_tbepti (otp, ocp[j], orow, ix);
                } else if (datatype == IRAF_BOOL) {
                    c_tbegtb (itp, icp[j], irow, &bx);
                    c_tbeptb (otp, ocp[j], orow, bx);
                } else if (datatype < 0) {
                    width = icol_descr->width;
                    tx = (char *)calloc (width+1, sizeof(char));
                    c_tbegtt (itp, icp[j], irow, tx, width);
                    c_tbeptt (otp, ocp[j], orow, tx);
                    free (tx);
                } else {
                    setError (ERR_DATATYPE_UNKNOWN,
                        "c_tbrcsc:  data type not supported");
                    return;
                }
            }
            status = checkError();
            if (status) {
                setError (status, "c_tbrcsc:  couldn't copy a column");
                return;
            }
        }
}
