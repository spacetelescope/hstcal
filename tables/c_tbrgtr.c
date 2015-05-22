# include <fitsio.h>
# include "ctables.h"

void c_tbrgtr (IRAFPointer tp, IRAFPointer *cp, float *buffer, Bool *nullflag,
                int numcols, int row) {

/* Read one or more values of type float from a table column.
arguments:
IRAFPointer tp          i: table descriptor
IRAFPointer *cp         i: array of column descriptors
float *buffer           o: array of values read from table
Bool *nullflag          o: array of flags, True if the value is undefined
int numcols             i: the length of the arrays 'cp', 'buffer', 'nullflag'
int row                 i: row number (one indexed) from which to read values
*/

        int i;

        for (i = 0;  i < numcols;  i++) {
            c_tbegtr (tp, cp[i], row, &buffer[i]);
            if (buffer[i] >= 0.99999 * IRAF_INDEFR &&
                buffer[i] <= 1.00001 * IRAF_INDEFR)
                nullflag[i] = True;
            else
                nullflag[i] = False;
        }
}
