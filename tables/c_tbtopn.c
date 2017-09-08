# include <stdlib.h>
# include <string.h>
# include <fitsio.h>
#include "hstcal.h"
# include "ctables.h"

IRAFPointer c_tbtopn (char *tablename, int iomode, IRAFPointer template) {

/* Open a table.
arguments:
char *tablename         i: name of table; may include extname or HDU number
                           in brackets (HDU = 0 is the primary HDU)
int iomode              i: IRAF I/O mode:  IRAF_READ_ONLY, IRAF_READ_WRITE,
                           IRAF_NEW_FILE, IRAF_NEW_COPY
IRAFPointer template    i: usually NULL (or 0), but if iomode is IRAF_NEW_COPY
                           this is the table descriptor for the template
function value          o: table descriptor
*/

        IRAFPointer tp = NULL;
        TableDescr *tbl_descr = NULL;
        char *fullname = NULL;         /* name with environment variable resolved */
        char *brackets = NULL;         /* expression in brackets, or NULL */
        int hdunum;                    /* HDU number (1 for primary HDU) */
        int hdutype;                   /* type of HDU */
        fitsfile *fptr = NULL;         /* CFITSIO pointer */
        int cf_iomode;                 /* CFITSIO I/O mode */
        int status = 0;

        clearError();

        /* expand environment variable (if any) in tablename */
        fullname = calloc (CHAR_FNAME_LENGTH+1, sizeof(*fullname));
        status = c_vfn2osfn (tablename, fullname);
        if (status != 0) {
            setError (status, "c_tbtopn:  error from c_vfn2osfn");
            return NULL;
        }

        tp = init_tp();
        tbl_descr = (TableDescr *)tp;

        strcpy (tbl_descr->tablename, tablename);
        strcpy (tbl_descr->fullname, fullname);
        brackets = strchr (fullname, '[');
        if (brackets == NULL) {
            tbl_descr->brackets[0] = '\0';
            strcpy (tbl_descr->filename, fullname);
        } else {
            int n;
            strcpy (tbl_descr->brackets, brackets);
            n = brackets - fullname;
            copyString (tbl_descr->filename, fullname, n);
        }
        tbl_descr->iomode = iomode;
        free (fullname);

        if (iomode == IRAF_READ_ONLY || iomode == IRAF_READ_WRITE) {
            if (iomode == IRAF_READ_ONLY)
                cf_iomode = READONLY;
            else
                cf_iomode = READWRITE;

            /* fits_open_file = ffopen */
            fits_open_file (&fptr, tbl_descr->fullname, cf_iomode, &status);
            if (status != 0) {
                setError (status, "c_tbtopn:  couldn't open file");
                free_tp (tp);
                return NULL;
            }
            tbl_descr->fptr = fptr;
            tbl_descr->template = NULL;         /* no template */
            tbl_descr->table_exists = 1;        /* the table does exist */

            /* fits_get_hdu_num = ffghdn */
            fits_get_hdu_num (fptr, &hdunum);   /* primary HDU is hdunum = 1 */
            if (tbl_descr->brackets != NULL &&
                strncmp (tbl_descr->brackets, "[0]", 3) == 0) {
                tbl_descr->hdunum = hdunum;
                tbl_descr->hdutype = IMAGE_HDU;
                tbl_descr->nrows = 0;
                tbl_descr->ncols = 0;
            } else {
                if (hdunum == 1) {
                    hdunum = 2;
                    /* fits_movabs_hdu = ffmahd */
                    fits_movabs_hdu (fptr, hdunum, &hdutype, &status);
                    if (status != 0) {
                        setError (status, "c_tbtopn:");
                        free_tp (tp);
                        return NULL;
                    }
                }
                tbl_descr->hdunum = hdunum;
                tbl_descr->hdutype = hdutype;

                tbSaveInfo (tp, &status);
                if (status != 0) {
                    setError (status, "c_tbtopn:");
                    free_tp (tp);
                    return NULL;
                }
            }

        } else {
            tbl_descr->template = template;
            /* the table hasn't actually been created yet */
            tbl_descr->table_exists = 0;
            tbl_descr->nrows = 0;
            if (iomode == IRAF_NEW_COPY)
                tbCopyTmpl (tp);        /* copy column info from template */
        }

        return tp;
}
