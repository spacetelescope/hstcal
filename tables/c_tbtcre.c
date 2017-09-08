# include <stdlib.h>
# include <string.h>
# include <fitsio.h>
#include "hstcal.h"
# include "ctables.h"

void c_tbtcre (IRAFPointer tp) {

/* In order to create a new table (either in a new file or appended to an
   existing file), the following procedure would be followed:
        call c_tbtopn to set up the table descriptor
        optionally call c_tbcdef1 one or more times to create columns
        call c_tbtcre to actually open the file and create the table
        optionally call c_tbhad[bdrit] to add header keywords
        call c_tbept[bdrist] or c_tbapt[bdrist] to write values to the table
        call c_tbtclo to close the table and file
argument:
IRAFPointer tp          i: table descriptor
*/

        fitsfile *fptr, *template_fptr;
	TableDescr *tbl_descr, *template_descr;
        IRAFPointer cp;
        ColumnDescr *col_descr;
	int file_exists;
	int i, ncols, nrows;
        int hdunum, hdutype;
        int bitpix=8, naxis=0;
        long naxes[2]={0,0};
	char **ttype, **tform, **tunit;
        char *filename;
        char keyword[SZ_FITS_STR+1];    /* keyword name */
        int status = 0;

	tbl_descr = (TableDescr *)tp;
	if (tbl_descr->table_exists)
	    return;
        if (tbl_descr->template != NULL) {
            template_descr = (TableDescr *)tbl_descr->template;
            template_fptr = template_descr->fptr;
        }

	file_exists = checkExists (tbl_descr->filename);

	if (file_exists) {
            /* fits_open_file = ffopen */
            fits_open_file (&fptr, tbl_descr->fullname, READWRITE, &status);
            if (status != 0) {
                setError (status, "c_tbtcre:  couldn't open file");
                return;
            }
	    tbl_descr->fptr = fptr;
	    /* find out how many extensions there are (hdunum) */
            /* fits_get_num_hdus = ffthdu */
            fits_get_num_hdus (fptr, &hdunum, &status);

            i = 1;
            /* fits_update_key = ffuky */
            fits_update_key (fptr, TLOGICAL, "EXTEND", &i,
                                "FITS file may contain extensions", &status);
            i = hdunum;
            fits_update_key (fptr, TINT, "NEXTEND", &i,
                                "Number of extensions", &status);

	    /* move to the last extension in the file */
            /* fits_movabs_hdu = ffmahd */
            fits_movabs_hdu (fptr, hdunum, &hdutype, &status);

	} else {

	    /* fits_create_file = ffinit */
	    fits_create_file (&fptr, tbl_descr->fullname, &status);
	    if (status != 0) {
		setError (status, "c_tbtcre:  couldn't create file");
		return;
	    }
	    tbl_descr->fptr = fptr;

            /* create primary header unit (with no data) */
            /* fits_write_imghdr = ffphps */
            fits_write_imghdr (fptr, bitpix, naxis, naxes, &status);
	    if (status != 0) {
		setError (status,
                        "c_tbtcre:  couldn't create primary header");
		return;
	    }
            /* copy keywords from template */
            if (tbl_descr->template != NULL)
                tbCopyPrimary (template_fptr, fptr, &status);

	    /* add or update FILENAME in the primary header */
            filename = (char *)calloc (CHAR_FNAME_LENGTH+1, sizeof(char));
            for (i = strlen (tbl_descr->filename) - 1;  i >= 0;  i--) {
                if (tbl_descr->filename[i] == '/')
                    break;
            }
            if (i > 0)
                strcpy (filename, tbl_descr->filename+i);
            else
                strcpy (filename, tbl_descr->filename);
            fits_update_key (fptr, TSTRING, "FILENAME", filename,
                                "File name", &status);
            free (filename);
            i = 1;
            fits_update_key (fptr, TLOGICAL, "EXTEND", &i,
                                "FITS file may contain extensions", &status);
            fits_update_key (fptr, TINT, "NEXTEND", &i,
                                "Number of extensions", &status);
        }

	ncols = tbl_descr->ncols;
	ttype = (char **)calloc (ncols, sizeof(char *));
	tform = (char **)calloc (ncols, sizeof(char *));
	tunit = (char **)calloc (ncols, sizeof(char *));
	for (i = 0;  i < ncols;  i++) {
            cp = tbl_descr->columns[i];
            col_descr = (ColumnDescr *)cp;
	    ttype[i] = (char *)calloc (FLEN_VALUE+1, sizeof(char));
	    tform[i] = (char *)calloc (FLEN_VALUE+1, sizeof(char));
	    tunit[i] = (char *)calloc (FLEN_VALUE+1, sizeof(char));
	    strcpy (ttype[i], col_descr->name);
	    strcpy (tform[i], col_descr->tform);
	    strcpy (tunit[i], col_descr->tunit);
	}

	nrows = 0;
	/* fits_create_tbl = ffcrtb */
	fits_create_tbl (fptr, BINARY_TBL, nrows, ncols,
		ttype, tform, tunit, NULL, &status);
	if (status != 0) {
	    setError (status, "c_tbtcre:  couldn't create table");
	    return;
	}
	for (i = 0;  i < ncols;  i++) {
	    free (ttype[i]);
	    free (tform[i]);
	    free (tunit[i]);
	}
	free (ttype);
	free (tform);
	free (tunit);

	for (i = 0;  i < ncols;  i++) {
            cp = tbl_descr->columns[i];
            col_descr = (ColumnDescr *)cp;
            if (col_descr->tdisp[0] != '\0') {
                sprintf (keyword, "TDISP%d", i+1);
                /* fits_update_key = ffuky */
                fits_update_key (fptr, TSTRING, keyword, col_descr->tdisp,
                        "display format for column", &status);
                if (status != 0) {
                    setError (status,
                        "c_tbtcre:  couldn't update TDISPi keyword");
                    return;
                }
            }
            if (col_descr->datatype == IRAF_SHORT ||
                col_descr->datatype == IRAF_INT) {
                sprintf (keyword, "TNULL%d", i+1);
                if (col_descr->datatype == IRAF_INT) {
                    int indef_int = IRAF_INDEFI;
                    fits_update_key (tbl_descr->fptr, TINT, keyword,
                        &indef_int, "undefined value for column", &status);
                } else {        /* short */
                    int indef_short = IRAF_INDEFS;
                    fits_update_key (tbl_descr->fptr, TINT, keyword,
                        &indef_short, "undefined value for column", &status);
                }
                if (status != 0) {
                    setError (status,
                        "c_tbtcre:  couldn't update TNULLi keyword");
                    return;
                }
            }
        }
        /* copy keywords from template */
        if (tbl_descr->template != NULL) {
            tbCopyHeader (template_fptr, fptr, &status);
            if (status != 0) {
                setError (status,
                        "c_tbtcre:  couldn't copy keywords from template");
                return;
            }
        }
}
