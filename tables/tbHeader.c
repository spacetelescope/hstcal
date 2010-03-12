# include <string.h>
# include <fitsio.h>
# include "ctables.h"

static void extractKeyword (char *card, char *keyword);
static int keywordExists (fitsfile *fptr, char *keyword, int *status);

void tbCopyPrimary (fitsfile *template_fptr, fitsfile *fptr, int *status) {

/* Copy the primary header from template_fptr to fptr.  
   The current HDU does not need to be the primary HDU, either in the
   template or in the output file.  This function will save the current
   HDUs, move to the primary for both files, copy the header, and then
   move back to the saved HDUs.
arguments:
fitsfile *template_fptr i: CFITSIO descriptor for template table
fitsfile *fptr          i: CFITSIO descriptor for table to receive keywords
int *status             o: CFITSIO status (0 is OK)
*/

        char card[SZ_FITS_STR+1];       /* FITS header record */
        char keyword[SZ_FITS_STR+1];    /* keyword from card */
        int nkeywords;                  /* number of keywords in header */
        int keynum;                     /* keyword number (one indexed) */
        int hdunum, hdutype;
        int save_hdunum, save_template_hdunum;

        /* save the current HDU numbers of the template and new table */
        /* fits_get_hdu_num = ffghdn */
        fits_get_hdu_num (template_fptr, &save_template_hdunum);
        fits_get_hdu_num (fptr, &save_hdunum);

        /* move to the primary header of each file */
        hdunum = 1;
        /* fits_movabs_hdu = ffmahd */
        fits_movabs_hdu (template_fptr, hdunum, &hdutype, status);
        fits_movabs_hdu (fptr, hdunum, &hdutype, status);

        /* fits_get_hdrspace = ffghsp */
        fits_get_hdrspace (template_fptr, &nkeywords, NULL, status);

        /* Copy keywords from the template to the current header, starting
           with keyword number 3 in order to skip SIMPLE and BITPIX.
        */
        for (keynum = 3;  keynum <= nkeywords;  keynum++) {
            /* fits_read_record = ffgrec */
            fits_read_record (template_fptr, keynum, card, status);
            if (strncmp (card, "NAXIS", 5) == 0)
                continue;
            if (strncmp (card, "ORIGIN", 6) == 0)
                continue;
            if (strncmp (card, "IRAF-TLM", 8) == 0)
                continue;
            if (strncmp (card,
                "COMMENT   FITS (Flexible Image Transport System)", 48) == 0)
                continue;
            if (strncmp (card,
                "COMMENT   and Astrophysics', volume 376, page 359", 49) == 0)
                continue;
            extractKeyword (card, keyword);
            if (keywordExists (fptr, keyword, status))
                continue;
            /* fits_write_record = ffprec */
            fits_write_record (fptr, card, status);
        }

        /* fits_write_date = ffpdat */
        fits_write_date (fptr, status);         /* write the DATE keyword */

        /* return to the saved HDU numbers */
        fits_movabs_hdu (template_fptr, save_template_hdunum,
                &hdutype, status);
        fits_movabs_hdu (fptr, save_hdunum, &hdutype, status);
}

void tbCopyHeader (fitsfile *template_fptr, fitsfile *fptr, int *status) {

/* Copy an extension header from template_fptr to fptr.  Keywords in the
   current template HDU will be copied to the current output HDU.
arguments:
fitsfile *template_fptr i: CFITSIO descriptor for template table
fitsfile *fptr          i: CFITSIO descriptor for table to receive keywords
int *status             o: CFITSIO status (0 is OK)
*/

        char card[SZ_FITS_STR+1];       /* FITS header record */
        char keyword[SZ_FITS_STR+1];    /* keyword from card */
        int nkeywords;                  /* number of keywords in header */
        int keynum;                     /* keyword number (one indexed) */

        /* fits_get_hdrspace = ffghsp */
        fits_get_hdrspace (template_fptr, &nkeywords, NULL, status);

        /* Copy keywords from the template to the current header, starting
           with keyword number 3 in order to skip XTENSION and BITPIX.
        */
        for (keynum = 3;  keynum <= nkeywords;  keynum++) {
            /* fits_read_record = ffgrec */
            fits_read_record (template_fptr, keynum, card, status);
            if (strncmp (card, "NAXIS", 5) == 0)
                continue;
            if (strncmp (card, "TTYPE", 5) == 0)
                continue;
            if (strncmp (card, "TFORM", 5) == 0)
                continue;
            if (strncmp (card, "TUNIT", 5) == 0)
                continue;
            if (strncmp (card, "TDISP", 5) == 0)
                continue;
            extractKeyword (card, keyword);
            if (strcmp (keyword, "PCOUNT") == 0)
                continue;
            if (strcmp (keyword, "GCOUNT") == 0)
                continue;
            if (strcmp (keyword, "TFIELDS") == 0)
                continue;
            if (strcmp (keyword, "ORIGIN") == 0)
                continue;
            if (keywordExists (fptr, keyword, status))
                continue;
            /* fits_write_record = ffprec */
            fits_write_record (fptr, card, status);
        }
}

static void extractKeyword (char *card, char *keyword) {

/* Copy the first eight characters of card to keyword, then trim
   trailing blanks.
arguments:
char *card              i: a header record
char *keyword           o: the keyword extracted from 'card'
*/

        int i;

        copyString (keyword, card, 8);
        for (i = strlen (keyword) - 1;  i > 0;  i--) {
            if (keyword[i] == ' ')
                keyword[i] = '\0';
        }
}

static int keywordExists (fitsfile *fptr, char *keyword, int *status) {

/* Does the specified keyword already exist in the current HDU?
   This is not a general-purpose function; it's used in this file only.
   A peculiarity of this function (which is appropriate in the context of
   its use in this file) is that for keywords HISTORY, COMMENT or blank,
   the function value will be 0, i.e. the keywords are considered to
   not exist, even if they do.  This is so that all HISTORY, all COMMENT,
   and all blank keywords will be copied.
arguments:
fitsfile *fptr          i: CFITSIO descriptor for table
char *keyword           i: a keyword name
int *status             o: CFITSIO status (0 is OK); it is not an error
                           if the keyword does not exist

function value          o: 1 if a non-commentary keyword exists, 0 otherwise
*/

        char value[SZ_FITS_STR+1];
        char comment[SZ_FITS_STR+1];
        int flag = 0;

        /* History, comment, and blank keywords are flagged as not existing,
           because we want all of them.
        */
        if (strcmp (keyword, "HISTORY") == 0 ||
            strcmp (keyword, "COMMENT") == 0 ||
            keyword[0] == ' ')
            return 0;

        /* fits_read_keyword = ffgkey */
        fits_read_keyword (fptr, keyword, value, comment, status);
        if (*status == 0) {
            flag = 1;
        } else if (*status == KEY_NO_EXIST) {
            flag = 0;
            *status = 0;
        } else {
            flag = 0;
        }

        return flag;
}
