# include <string.h>
# include <fitsio.h>
#include "hstcal.h"
# include "ctables.h"

Bool c_tbhgtb (IRAFPointer tp, char *keyword) {

/* Get a boolean keyword from a header.
arguments:
IRAFPointer tp          i: table descriptor
char *keyword           i: keyword name

function value          o: value (True or False) for the keyword
*/

        TableDescr *tbl_descr;
        int value;
        int status = 0;

        tbl_descr = (TableDescr *)tp;

        /* fits_read_key = ffgky */
        fits_read_key (tbl_descr->fptr, TLOGICAL, keyword,
                        &value, NULL, &status);
        if (status != 0)
            setError (status, "c_tbhgtb:  error reading keyword");

        if (value)
            return True;
        else
            return False;
}

double c_tbhgtd (IRAFPointer tp, char *keyword) {

/* Get a keyword of type double from a header.
arguments:
IRAFPointer tp          i: table descriptor
char *keyword           i: keyword name

function value          o: value for the keyword
*/

        TableDescr *tbl_descr;
        double value;
        int status = 0;

        tbl_descr = (TableDescr *)tp;

        /* fits_read_key = ffgky */
        fits_read_key (tbl_descr->fptr, TDOUBLE, keyword,
                        &value, NULL, &status);
        if (status != 0)
            setError (status, "c_tbhgtd:  error reading keyword");

        return value;
}

float c_tbhgtr (IRAFPointer tp, char *keyword) {

/* Get a keyword of type float from a header.
arguments:
IRAFPointer tp          i: table descriptor
char *keyword           i: keyword name

function value          o: value for the keyword
*/

        TableDescr *tbl_descr;
        float value;
        int status = 0;

        tbl_descr = (TableDescr *)tp;

        /* fits_read_key = ffgky */
        fits_read_key (tbl_descr->fptr, TFLOAT, keyword,
                        &value, NULL, &status);
        if (status != 0)
            setError (status, "c_tbhgtr:  error reading keyword");

        return value;
}

int c_tbhgti (IRAFPointer tp, char *keyword) {

/* Get an integer keyword from a header.
arguments:
IRAFPointer tp          i: table descriptor
char *keyword           i: keyword name

function value          o: value for the keyword
*/

        TableDescr *tbl_descr;
        int value;
        int status = 0;

        tbl_descr = (TableDescr *)tp;

        /* fits_read_key = ffgky */
        fits_read_key (tbl_descr->fptr, TINT, keyword,
                        &value, NULL, &status);
        if (status != 0)
            setError (status, "c_tbhgti:  error reading keyword");

        return value;
}

void c_tbhgtt (IRAFPointer tp, char *keyword, char *text, int maxch) {

/* Get a string keyword from a header.
arguments:
IRAFPointer tp          i: table descriptor
char *keyword           i: keyword name
char *text              o: value for the keyword
int maxch               i: maximum length of the string (not including '\0')
*/

        TableDescr *tbl_descr;
        char value[CHAR_FNAME_LENGTH+1]; /* hopefully large enough for any keyword */
        int status = 0;

        tbl_descr = (TableDescr *)tp;

        /* fits_read_key = ffgky */
        fits_read_key (tbl_descr->fptr, TSTRING, keyword,
                        value, NULL, &status);
        if (status != 0)
            setError (status, "c_tbhgtd:  error reading keyword");
        copyString (text, value, maxch);
}
