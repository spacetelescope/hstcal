# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <fitsio.h>
# include "ctables.h"

# define SZ_TABLENAME 27    /* long enough for labq01hzq_lampflash.fits[1] */
# define MAXDIM 7

static int testError (char *message);

/* I typically use 'labq01hzq_lampflash.fits[1]' as the command-line argument.
   g.fits is opened (read-write).  junk.fits is created as a new copy of
   g.fits, and then array columns are added to junk.fits; values are written
   (not copied) to junk.fits using c_tbept and c_tbapt.
   junk3.fits is also created as a new copy of g.fits, and rows are copied
   from g.fits to junk3.fits using c_tbrcsc.
   j.fits is opened, and the character string column 'b' is gotten and
   printed; the maxch is only 4, which deliberately truncates "three".
   j4.fits is created, one column is defined and one element written; then
   c_tbhcal is used to copy most keywords from j.fits[1] to j4.fits[1].
   Then another table is appended to j4.fits (i.e. j4.fits[2]); one column
   is defined and one element is written.  Then specific keywords are copied
   one at a time from j.fits[1] to j4.fits[2] to test the keyword get and
   add functions, and the keyword comment get and put functions.
*/

int main (int argc, char **argv) {

        IRAFPointer imt;
        IRAFPointer tp, tp2, tp3, tp4;
        IRAFPointer cp;
        IRAFPointer cptr[6];
        IRAFPointer a_cp[6];
        IRAFPointer icp[50], ocp[50];   /* for c_tbrcsc to junk3.fits */
        TableDescr *tbl_descr;
        ColumnDescr *col_descr;
        int status;
        int copied;
        int lenfilename, done;
        int nrows, ncols, i, j;
        int ndim, axlen[MAXDIM];
        char tname[SZ_TABLENAME+1];
        int colnum, datatype, lendata, lenfmt;
        char colname[81], colunits[81], colfmt[81];
        char keyword[81], value[81], comment[81];
        int    a_value;
        short  b_value;
        double c_value;
        float  d_value;
        char   e_value[21];
        Bool   f_value;
        int    aa_value[10];
        short  bb_value[10];
        double cc_value[10];
        float  dd_value[10];
        char   **ee_value;
        Bool   ff_value[10];

        ee_value = (char **)calloc (10, sizeof(char *));
        for (i = 0;  i < 10;  i++) {
            ee_value[i] = (char *)calloc (7, sizeof(char));
        }
        sprintf (ee_value[0], "One");
        sprintf (ee_value[1], "Two");
        sprintf (ee_value[2], "Three");
        sprintf (ee_value[3], "Four");
        sprintf (ee_value[4], "Five");
        sprintf (ee_value[5], "Six");
        sprintf (ee_value[6], "Seven");
        sprintf (ee_value[7], "Eight");
        sprintf (ee_value[8], "Nine");
        sprintf (ee_value[9], "Ten");

        if (argc != 2) {
            printf ("syntax:  a.out filename\n");
            exit (2);
        }

        printf ("file name template functions:\n");
        imt = c_imtopen (argv[1]);
        printf ("  number of names in list = %d\n", c_imtlen (imt));
        done = 0;
        while (!done) {
            lenfilename = c_imtgetim (imt, comment, 80);
            done = (lenfilename == 0);
            if (!done)
                printf ("  file name = %s, len = %d\n", comment, lenfilename);
        }
        c_imtrew (imt);
        printf ("  get the list again, with a 5-character limit\n");
        done = 0;
        while (!done) {
            lenfilename = c_imtgetim (imt, comment, 5);
            if ((status = testError ("from c_imtgetim")))
                clear_cvoserr();
            done = (lenfilename == 0);
            if (!done)
                printf ("  file name = %s, len = %d\n", comment, lenfilename);
        }
        c_imtclose (imt);
        printf ("file name template closed\n");

        printf ("table exists?  %d\n", c_tbtacc (argv[1]));
        if ((status = testError ("from c_tbtacc")))
            exit (status);

        printf ("about to copy the primary header of %s\n", argv[1]);
        c_tbfpri (argv[1], "copy.fits", &copied);
        if (copied)
            printf ("primary header copied to copy.fits\n");
        else
            printf ("primary header was not copied to copy.fits\n");
        if ((status = testError ("from c_tbfpri")))
            exit (status);

        i = c_tbparse (argv[1], value, comment, 80, &j);
        printf ("from c_tbparse:  %d, '%s', '%s', %d\n", i, value, comment, j);

        tp = c_tbtopn (argv[1], IRAF_READ_ONLY, 0);
        if (tp != NULL)
            printf ("%s opened\n", argv[1]); fflush (stdout);
        if ((status = testError ("from c_tbtopn")))
            exit (status);

        c_tbtnam (tp, tname, SZ_TABLENAME);
        status = c_iraferr();
        printf ("table name is %s; status = %d\n", tname, status);
        if (status != 0)
            clear_cvoserr();

        c_tbcfnd1 (tp, "GROSS", &cp);
        if (cp == NULL) {
            printf ("GROSS column was not found\n");
        } else {
            c_tbciga (tp, cp, &ndim, axlen, MAXDIM);
            printf ("GROSS column was found (dimension = %d)\n", ndim);
            for (i = 0;  i < ndim;  i++)
                printf ("   axis %d, length = %d\n", i+1, axlen[i]);
        }

        /* begin section that relies on ctables.h */
        tbl_descr = (TableDescr *)tp;
        /* print info about the table */
        printf ("tablename = %s\n", tbl_descr->tablename);
        printf ("fullname = %s\n", tbl_descr->fullname);
        printf ("filename = %s\n", tbl_descr->filename);
        printf ("brackets = %s\n", tbl_descr->brackets);
        if (tbl_descr->template != NULL)
            printf ("template is not null\n");
        else
            printf ("template is null\n");
        printf ("iomode = %d\n", tbl_descr->iomode);
        printf ("hdunum = %d\n", tbl_descr->hdunum);
        printf ("hdutype = %d\n", tbl_descr->hdutype);
        printf ("nrows = %ld\n", tbl_descr->nrows);
        printf ("ncols = %d = %d\n", tbl_descr->ncols,
                                     c_tbpsta (tp, TBL_NCOLS));
        printf ("alloc_cols = %d\n", tbl_descr->alloc_cols);

        ncols = ((TableDescr *)tp)->ncols;
        for (i = 1;  i <= ncols;  i++) {
            cp = c_tbcnum (tp, i);
            if (cp == NULL) {
                printf ("%3d:  not found\n", i);
            } else {
                col_descr = (ColumnDescr *)cp;
                printf ("%3d:  '%s' '%s' '%s' '%s', %2d, %ld\n",
                        i, col_descr->name, col_descr->tform, col_descr->tunit,
                        col_descr->tdisp, col_descr->typecode,
                        col_descr->repeat);
                c_tbcinf (cp, &colnum, colname, colunits,
                        colfmt, &datatype, &lendata, &lenfmt);
                /* replace colname, colunits, colfmt */
                c_tbcigt (cp, TBL_COL_NAME, colname, 80);
                c_tbcigt (cp, TBL_COL_UNITS, colunits, 80);
                c_tbcigt (cp, TBL_COL_FMT, colfmt, 80);
                printf ("%3d:  '%s' '%s' '%s' '%s', %2d, %d %d\n",
                        colnum, colname, col_descr->tform, colunits,
                        colfmt, datatype, lendata, lenfmt);
            }
        }
        if ((status = testError ("from c_tbcinf, c_tbcigt")))
            exit (status);

        c_tbtclo (tp);
        if ((status = testError ("from c_tbtclo for <argv[1]>")))
            exit (status);
        else
            printf ("%s closed\n", argv[1]); fflush (stdout);

        printf ("about to open g.fits\n"); fflush (stdout);
        tp = c_tbtopn ("g.fits", IRAF_READ_WRITE, 0);
        printf ("g.fits opened\n"); fflush (stdout);
        if ((status = testError ("from c_tbtopn for g.fits")))
            exit (status);
        printf ("g.fits is type %d (expected %d)\n",
                c_tbpsta (tp, TBL_WHTYPE), TBL_TYPE_FITS);
        printf ("g.fits has %d header keywords\n", c_tbpsta (tp, TBL_NPAR));
        for (i = 0;  i < 50;  i++) {
            c_tbhgnp (tp, i+1, keyword, &datatype, value);
            if (c_iraferr() != 0) {
                clear_cvoserr();
                break;
            }
            c_tbhgcm (tp, keyword, comment, 80);
            printf ("%2d '%s' = '%s' [%d] '%s'\n", i+1, keyword, value,
                        datatype, comment);
        }
        c_tbcfnd1 (tp, "a", &cptr[0]);
        c_tbcfnd1 (tp, "b", &cptr[1]);
        c_tbcfnd1 (tp, "c", &cptr[2]);
        c_tbcfnd1 (tp, "d", &cptr[3]);
        c_tbcfnd1 (tp, "e", &cptr[4]);
        c_tbcfnd1 (tp, "f", &cptr[5]);
        nrows = c_tbpsta (tp, TBL_NROWS);
        for (i = 0;  i < nrows; i++) {
            c_tbegti (tp, cptr[0], i+1, &a_value);
            c_tbegts (tp, cptr[1], i+1, &b_value);
            c_tbegtd (tp, cptr[2], i+1, &c_value);
            c_tbegtr (tp, cptr[3], i+1, &d_value);
            c_tbrgtr (tp, cptr+3, dd_value, ff_value, 1, i+1);
            if (ff_value[0])
                strcpy (comment, "INDEF");
            else
                strcpy (comment, "OK");
            c_tbegtt (tp, cptr[4], i+1, e_value, 20);
            c_tbegtb (tp, cptr[5], i+1, &f_value);
            printf ("row %d:  %d %d %g %g (= %g, %s) %s %d\n",
                i+1, a_value, b_value, c_value,
                d_value, dd_value[0], comment,
                e_value, f_value);
        }

        /* create a new file and table */
        printf ("about to open junk.fits\n"); fflush (stdout);
        tp2 = c_tbtopn ("junk.fits", IRAF_NEW_COPY, tp);
        if ((status = testError ("from c_tbtopn for junk.fits")))
            exit (status);
        /* add columns */
        c_tbcdef1 (tp2, &cp, "aa", "", "%5d", IRAF_INT, 3);
        c_tbcdef1 (tp2, &cp, "bb", "", "%5d", IRAF_SHORT, 3);
        c_tbcdef1 (tp2, &cp, "cc", "", "%8.1f", IRAF_DOUBLE, 3);
        c_tbcdef1 (tp2, &cp, "dd", "", "%5.1f", IRAF_REAL, 3);
        c_tbcdef1 (tp2, &cp, "ee", "", "%-6s", -6, 3);
        c_tbcdef1 (tp2, &cp, "ff", "", "%5b", IRAF_BOOL, 3);
        status = c_iraferr();
        if (status == 0) {
            printf ("columns added, status = %d\n", status);
        } else {
            printf ("ERROR:  columns not added, %d, '%s'\n", status,
                        c_iraferrmsg());
            clear_cvoserr();
        }
        c_tbtcre (tp2);
        if ((status = testError ("from c_tbtcre for junk.fits")))
            exit (status);
        nrows = c_tbpsta (tp, TBL_NROWS);       /* number of rows in g.fits */
        c_tbcfnd1 (tp2, "a", &cptr[0]);         /* write to junk.fits */
        c_tbcfnd1 (tp2, "b", &cptr[1]);
        c_tbcfnd1 (tp2, "c", &cptr[2]);
        c_tbcfnd1 (tp2, "d", &cptr[3]);
        c_tbcfnd1 (tp2, "e", &cptr[4]);
        c_tbcfnd1 (tp2, "f", &cptr[5]);
        c_tbcfnd1 (tp2, "aa", &a_cp[0]);
        c_tbcfnd1 (tp2, "bb", &a_cp[1]);
        c_tbcfnd1 (tp2, "cc", &a_cp[2]);
        c_tbcfnd1 (tp2, "dd", &a_cp[3]);
        c_tbcfnd1 (tp2, "ee", &a_cp[4]);
        c_tbcfnd1 (tp2, "ff", &a_cp[5]);
        for (i = 0;  i < 6;  i++) {
            if (cptr[i] == NULL || a_cp[i] == NULL) {
                printf ("ERROR:  column not found\n");
                exit (4);
            }
        }
        for (i = 0;  i < 3;  i++) {
            aa_value[i] = i;
            bb_value[i] = -i;
            cc_value[i] = i;
            dd_value[i] = -i;
        }
        ff_value[0] = True;
        ff_value[1] = False;
        ff_value[2] = True;
        for (i = 0;  i < nrows; i++) {
            a_value = i+1 - 4;
            b_value = i+1 - 4;
            c_value = i+1. - 4.;
            d_value = i+1. - 4.;
            sprintf (e_value, "%d", i+1 - 4);
            if (i / 2 * 2 == i)
                f_value = False;
            else
                f_value = True;
            c_tbepti (tp2, cptr[0], i+1, a_value);
            c_tbepts (tp2, cptr[1], i+1, b_value);
            c_tbeptd (tp2, cptr[2], i+1, c_value);
            c_tbeptr (tp2, cptr[3], i+1, d_value);
            c_tbeptt (tp2, cptr[4], i+1, e_value);
            c_tbeptb (tp2, cptr[5], i+1, f_value);
            c_tbapti (tp2, a_cp[0], i+1, aa_value, 1, 3);
            c_tbapts (tp2, a_cp[1], i+1, bb_value, 1, 3);
            c_tbaptd (tp2, a_cp[2], i+1, cc_value, 1, 3);
            c_tbaptr (tp2, a_cp[3], i+1, dd_value, 1, 3);
            c_tbaptt (tp2, a_cp[4], i+1, ee_value, 6, 1, 3);
            c_tbaptb (tp2, a_cp[5], i+1, ff_value, 1, 3);
        }
        for (i = 0;  i < 3;  i++) {
            aa_value[i] = IRAF_INDEFI;
            bb_value[i] = IRAF_INDEFS;
            cc_value[i] = IRAF_INDEFD;
            dd_value[i] = IRAF_INDEFR;
            strcpy (ee_value[i], "INDEF");
            ff_value[i] = False;
        }
        for (i = 0;  i < nrows; i++) {
            j = c_tbagti (tp2, a_cp[0], i+1, aa_value, 1, 3);
            j = c_tbagts (tp2, a_cp[1], i+1, bb_value, 1, 3);
            j = c_tbagtd (tp2, a_cp[2], i+1, cc_value, 1, 3);
            j = c_tbagtr (tp2, a_cp[3], i+1, dd_value, 1, 3);
            j = c_tbagtt (tp2, a_cp[4], i+1, ee_value, 1, 3, 6);
            j = c_tbagtb (tp2, a_cp[5], i+1, ff_value, 1, 3);
            if ((status = testError ("getting array values from junk.fits")))
                exit (status);
            printf ("row %d:  %d  %d  %d\n",
                        i, aa_value[0], aa_value[1], aa_value[2]);
            printf ("row %d:  %d  %d  %d\n",
                        i, bb_value[0], bb_value[1], bb_value[2]);
            printf ("row %d:  %g  %g  %g\n",
                        i, cc_value[0], cc_value[1], cc_value[2]);
            printf ("row %d:  %g  %g  %g\n",
                        i, dd_value[0], dd_value[1], dd_value[2]);
            printf ("row %d:  %s  %s  %s\n",
                        i, ee_value[0], ee_value[1], ee_value[2]);
            printf ("row %d:  %d  %d  %d\n",
                        i, ff_value[0], ff_value[1], ff_value[2]);
        }
        /* set row 3 (a, b, c, d, e, f) to INDEF */
        c_tbrudf (tp2, cptr, 6, 3);
        if ((status = testError ("from c_tbrudf for junk.fits")))
            exit (status);

        c_tbtclo (tp2);
        status = c_iraferr();
        printf ("junk.fits closed, %d\n", status); fflush (stdout);

        printf ("about to open junk3.fits\n"); fflush (stdout);
        tp3 = c_tbtopn ("junk3.fits", IRAF_NEW_COPY, tp);
        if ((status = testError ("from c_tbtopn for junk3.fits")))
            exit (status);
        printf ("junk3.fits opened\n"); fflush (stdout);
        c_tbtcre (tp3);
        if ((status = testError ("from c_tbtcre for junk3.fits")))
            exit (status);
        printf ("junk3.fits created\n"); fflush (stdout);
        ncols = ((TableDescr *)tp)->ncols;      /* columns in g.fits */
        printf ("ncols = %d\n", ncols); fflush (stdout);
        for (j = 0;  j < ncols; j++) {
            icp[j] = c_tbcnum (tp, j+1);
            ocp[j] = c_tbcnum (tp3, j+1);
        }
        nrows = c_tbpsta (tp, TBL_NROWS);       /* number of rows in g.fits */
        printf ("nrows = %d\n", nrows); fflush (stdout);
        for (i = 1;  i <= nrows; i++) {
            c_tbrcsc (tp, tp3, icp, ocp, i, i, ncols);
            if ((status = testError ("from c_tbrcsc to junk3.fits")))
                exit (status);
        }
        c_tbtclo (tp3);
        if ((status = testError ("from c_tbtclo for junk3.fits")))
            exit (status);
        printf ("junk3.fits closed\n"); fflush (stdout);

        c_tbtclo (tp);
        if ((status = testError ("from c_tbtclo for g.fits")))
            exit (status);
        printf ("g.fits closed\n"); fflush (stdout);

        /* get an array of character strings */
        printf ("open j.fits\n"); fflush (stdout);
        tp = c_tbtopn ("j.fits", IRAF_READ_ONLY, 0);
        if ((status = testError ("from c_tbtopn for j.fits")))
            exit (status);
        nrows = c_tbpsta (tp, TBL_NROWS);
        c_tbcfnd1 (tp, "b", &cp);
        if (cp == NULL) {
            printf ("'b' column was not found\n");
            exit (2);
        }
        for (i = 0;  i < nrows; i++) {
            j = c_tbagtt (tp, cp, i+1, ee_value, 1, 3, 4); /* 4 to truncate */
            printf ("row %d:  '%s'  '%s'  '%s'\n",
                        i+1, ee_value[0], ee_value[1], ee_value[2]);
        }
        if ((status = testError ("from c_tbagtt for j.fits")))
            exit (status);

        printf ("create j4.fits[1] with one column\n"); fflush (stdout);
        tp4 = c_tbtopn ("j4.fits", IRAF_NEW_FILE, 0);
        if ((status = testError ("from c_tbtopn for j4.fits")))
            exit (status);
        printf ("j4.fits opened\n"); fflush (stdout);
        c_tbcdef1 (tp4, &cp, "a", "", "%5d", IRAF_INT, 1);
        if (cp == NULL) {
            printf ("column 'a' not created\n"); fflush (stdout);
        } else {
            printf ("column 'a' was created OK\n"); fflush (stdout);
        }
        c_tbtcre (tp4);
        if ((status = testError ("from c_tbtcre for j4.fits[1]")))
            exit (status);
        c_tbepti (tp4, cp, 2, 11);
        if ((status = testError ("from c_tbepti for j4.fits[1] row 2")))
            exit (status);
        printf ("copy keywords from j.fits to j4.fits[1] using c_tbhcal\n");
        c_tbhcal (tp, tp4);
        if ((status = testError ("from c_tbhcal for j4.fits[1]")))
            exit (status);
        printf ("header keywords copied using c_tbhcal\n");
        c_tbtclo (tp4);
        if ((status = testError ("from c_tbtclo for j4.fits[1]")))
            exit (status);
        printf ("j4.fits[1] closed\n");
        fflush (stdout);

        printf ("create j4.fits[2] with one column\n"); fflush (stdout);
        tp4 = c_tbtopn ("j4.fits", IRAF_NEW_FILE, 0);
        if ((status = testError ("from c_tbtopn for j4.fits[2]")))
            exit (status);
        c_tbcdef1 (tp4, &cp, "a", "", "%5d", IRAF_INT, 1);
        c_tbtcre (tp4);
        c_tbepti (tp4, cp, 2, 22);
        if ((status = testError ("from c_tbepti for j4.fits[2] row 2")))
            exit (status);

        printf ("copy keywords from j.fits to j4.fits[2]\n");
        a_value = c_tbhgti (tp, "one");
        c_tbhadi (tp4, "one", a_value);
        c_tbhgcm (tp, "one", comment, 80);
        c_tbhpcm (tp4, "one", comment);

        c_value = c_tbhgtd (tp, "two");
        c_tbhadd (tp4, "two", c_value);
        c_tbhgcm (tp, "two", comment, 80);
        c_tbhpcm (tp4, "two", comment);

        d_value = c_tbhgtr (tp, "two_r");
        c_tbhadr (tp4, "two_r", d_value);
        c_tbhgcm (tp, "two_r", comment, 80);
        c_tbhpcm (tp4, "two_r", comment);

        c_tbhgtt (tp, "three", value, 80);
        c_tbhadt (tp4, "three", value);
        /* c_tbhptt (tp4, "three", value);  commented out after testing once */
        c_tbhgcm (tp, "three", comment, 80);
        c_tbhpcm (tp4, "three", comment);

        f_value = c_tbhgtb (tp, "four");
        c_tbhadb (tp4, "four", f_value);
        c_tbhgcm (tp, "four", comment, 80);
        c_tbhpcm (tp4, "four", comment);

        f_value = c_tbhgtb (tp, "five");
        c_tbhadb (tp4, "five", f_value);
        c_tbhgcm (tp, "five", comment, 80);
        c_tbhpcm (tp4, "five", comment);

        if ((status = testError ("from c_tbh* for j4.fits[2]")))
            exit (status);

        c_tbtclo (tp4);
        if ((status = testError ("from c_tbtclo for j4.fits[2]")))
            exit (status);
        printf ("j4.fits[2] closed\n");
        fflush (stdout);

        c_tbtclo (tp);
        if ((status = testError ("from c_tbtclo for j.fits")))
            exit (status);
        printf ("j.fits closed\n");
        fflush (stdout);

        exit (0);
}

static int testError (char *message) {

        int status;

        status = c_iraferr();
        if (status) {
            printf ("ERROR %d:  %s\n", status, message);
            fflush (stdout);
            printf ("    %s\n", c_iraferrmsg());
            fflush (stdout);
        }

        return status;
}
