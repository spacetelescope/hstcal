# include <fitsio.h>
# include "ctables.h"

void c_tbhcal (IRAFPointer itp, IRAFPointer otp) {

/* Copy most keywords from one table header to another.  Keywords describing
   the table itself (e.g. NAXISi, TTYPEi) will not be copied.
arguments:
IRAFPointer itp         i: descriptor for input table
IRAFPointer otp         i: descriptor for output table
*/

        TableDescr *itbl_descr, *otbl_descr;
        int status = 0;

        itbl_descr = (TableDescr *)itp;
        otbl_descr = (TableDescr *)otp;

        tbCopyHeader (itbl_descr->fptr, otbl_descr->fptr, &status);
        if (status != 0)
            setError (status, "c_tbhcal:  couldn't copy keywords");
}
