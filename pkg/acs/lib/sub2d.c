# include <math.h>
# include "hstio.h"
# include "hstcalerr.h"	/* SIZE_MISMATCH */

/* Subtract the second SingleGroup from the first, leaving the
   result in the first.  

   (*a) -= (*b)

   The science data arrays are subtracted; the error arrays are combined;
   the data quality arrays are combined.

   M.D. De La Pena: 05 June 2018
   Created this file/routine from sub1d.c.
   the data quality arrays are combined.  For all extensions, the entire 2D
   arrays are used.
*/
int sub2d (SingleGroup *a, SingleGroup *b) {

/* arguments:
SingleGroup *a   io: input data, output difference
SingleGroup *b   i: second input data
*/

    extern int status;

    float da, db;           /* errors for a and b */
    short dqa, dqb, dqab;   /* data quality for a, b, combined */

    if ((a->sci.data.nx != b->sci.data.nx) || (a->sci.data.ny != b->sci.data.ny))
        return (status = SIZE_MISMATCH);

    /* science, error, and DQ data */
    int dimx = a->sci.data.nx;
    int dimy = a->sci.data.ny;

    {unsigned int i, j;
    for (j = 0; j < dimy; j++) {
        for (i = 0;  i < dimx;  i++) {
            /* science array */
            Pix(a->sci.data, i, j) -= Pix(b->sci.data, i, j);

            /* error array */
            da = Pix (a->err.data, i, j);
            db = Pix (b->err.data, i, j);
            Pix (a->err.data, i, j) = sqrt (da * da + db * db);

            /* data quality */
            dqa = DQPix (a->dq.data, i, j);
            dqb = DQPix (b->dq.data, i, j);
            dqab = dqa | dqb;
            DQSetPix (a->dq.data, i, j, dqab);
       }
    }}

    return (status);
}
