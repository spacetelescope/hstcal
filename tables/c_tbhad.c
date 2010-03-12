# include <string.h>
# include <fitsio.h>
# include "ctables.h"

void c_tbhadb (IRAFPointer tp, char *keyword, Bool value) {

/* Update or add a boolean keyword to a header.
arguments:
IRAFPointer tp          i: table descriptor
char *keyword           i: keyword name
Bool value              i: value (True or False) for the keyword
*/

        TableDescr *tbl_descr;
        int status = 0;

        tbl_descr = (TableDescr *)tp;

        /* fits_update_key = ffuky */
        fits_update_key (tbl_descr->fptr, TLOGICAL, keyword, &value,
                                NULL, &status);
        if (status != 0)
            setError (status, "c_tbhadb:  error updating keyword");
}

void c_tbhadd (IRAFPointer tp, char *keyword, double value) {

/* Update or add a keyword of type double to a header.
arguments:
IRAFPointer tp          i: table descriptor
char *keyword           i: keyword name
double value            i: value for the keyword
*/

        TableDescr *tbl_descr;
        int status = 0;

        tbl_descr = (TableDescr *)tp;

        fits_update_key (tbl_descr->fptr, TDOUBLE, keyword, &value,
                                NULL, &status);
        if (status != 0)
            setError (status, "c_tbhadd:  error updating keyword");
}

void c_tbhadr (IRAFPointer tp, char *keyword, float value) {

/* Update or add a keyword of type float to a header.
arguments:
IRAFPointer tp          i: table descriptor
char *keyword           i: keyword name
float value             i: value for the keyword
*/

        TableDescr *tbl_descr;
        int status = 0;

        tbl_descr = (TableDescr *)tp;

        fits_update_key (tbl_descr->fptr, TFLOAT, keyword, &value,
                                NULL, &status);
        if (status != 0)
            setError (status, "c_tbhadr:  error updating keyword");
}

void c_tbhadi (IRAFPointer tp, char *keyword, int value) {

/* Update or add an integer keyword to a header.
arguments:
IRAFPointer tp          i: table descriptor
char *keyword           i: keyword name
int value               i: value for the keyword
*/

        TableDescr *tbl_descr;
        int status = 0;

        tbl_descr = (TableDescr *)tp;

        fits_update_key (tbl_descr->fptr, TINT, keyword, &value,
                                NULL, &status);
        if (status != 0)
            setError (status, "c_tbhadi:  error updating keyword");
}

void c_tbhadt (IRAFPointer tp, char *keyword, char *text) {

/* Update or add a string keyword to a header.
arguments:
IRAFPointer tp          i: table descriptor
char *keyword           i: keyword name
char *text              i: value for the keyword
*/

        TableDescr *tbl_descr;
        int status = 0;

        tbl_descr = (TableDescr *)tp;

        fits_update_key (tbl_descr->fptr, TSTRING, keyword, text,
                                NULL, &status);
        if (status != 0)
            setError (status, "c_tbhadt:  error updating keyword");
}
